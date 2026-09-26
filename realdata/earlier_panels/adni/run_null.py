#!/usr/bin/env python
"""ADNI 12-event population: comparisons of two groups drawn from e3/e3 carriers (Fig. 4e).

For each replicate (draw) two disjoint groups are drawn without replacement, within diagnosis,
from the e3/e3 participants of the frame built by run_observed.build_frame (core population,
12 events, 180-day window). Group A has 69/5/12 CN/MCIc/AD participants (close to the e2 group)
and group B 56/58/86 (the e4 proportions scaled to 200). Both groups are e3/e3 carriers, so any
ordering difference detected is a false positive. Each draw is tested once, as a single
two-group comparison: P = (1 + E)/(B + 1) with B = 599 and ties counted, rejection if E <= 29.
Exact-decision early stopping ends the permutations of a draw once its decision is fixed.

Two fits per draw (output keys):
  standard_U  separately fitted DEBM, unmodified pyebm search, unrestricted permutation
  concord_D   CONCORD (common proportions = smallest proportion of each diagnosis over the two
              drawn groups; continued search), within-diagnosis permutation

  --stage n2       the only stage (kept as a name for the output folder)
  --panel P12      12 events (default); P8 and P5 use subsets of the same adjusted columns and the
                   same draws (Fig. 4e, ADNI 8- and 5-event settings)
  --replicates a:b half-open range of draw indices; one JSON per draw (rep_XXXX.json)

Draw i uses numpy.random.default_rng([seed, 2, i]) and permutation seed seed*1000 + i
(seed 20260914). Requires the concord-ebm package (version 0.1.1, commit 40fc046) and
pyebm 2.0.3. The rates of Fig. 4e are aggregated by ../aggregate_earlier_panels.py.

  python run_null.py --table <derived>/adni_p12_deid.csv --out $CONCORD_WORK_DIR/earlier/adni/null \\
      --stage n2 --replicates 0:200
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from concord import compare

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_observed import GROUP_ORDER, PANELS, PLAN_VERSION, build_frame  # noqa: E402

# Design of the draws: group A close to the e2 group, group B the e4 proportions scaled to n = 200.
N2_DRAW = {"pA_e2comp": {"CN": 69, "MCIc": 5, "AD": 12},
           "pB_e4comp": {"CN": 56, "MCIc": 58, "AD": 86}}
# (label, estimator, search, permutation scheme) as passed to concord.compare
ARMS_FROZEN = [("standard_U", "standard", "original", "unrestricted"),
               ("concord_D", "invariant_min", "repaired", "diagnosis")]


def draw_within(frame, group, counts, rng, exclude=None):
    sub = frame[frame.APOE == group]
    if exclude is not None:
        sub = sub[~sub.index.isin(exclude)]
    picked = []
    for label, n in counts.items():
        idx = sub.index[sub.Diagnosis == label].to_numpy()
        if len(idx) < n:
            raise ValueError(f"{group}/{label}: need {n}, have {len(idx)}")
        picked.extend(rng.choice(idx, size=n, replace=False).tolist())
    return frame.loc[picked]


def result_summary(res):
    out = {"status": res.status, "distances": res.distances, "tests": {}}
    for scheme, test in res.tests.items():
        out["tests"][scheme] = {"pairs": {p: {"exceedances": c["exceedances"], "p": c["p"], "reject": c["reject"],
                                             "p_bounds": c.get("p_bounds")} for p, c in test["pairs"].items()},
                                "reject_any_pair": any(c["reject"] is True for c in test["pairs"].values())}
    out["fits"] = res.fits
    return out


def rng_for(seed, stage, rep):
    code = {"n2": 2}[stage]
    return np.random.default_rng([seed, code, rep])


def n2_pseudo(frame, rep, args):
    """Two disjoint groups from e3/e3; the same replicate index gives the same draw."""
    rng = rng_for(args.seed, "n2", rep)
    (name_a, counts_a), (name_b, counts_b) = N2_DRAW.items()
    part_a = draw_within(frame, "e33", counts_a, rng)
    part_b = draw_within(frame, "e33", counts_b, rng, exclude=part_a.index)
    pseudo = pd.concat([part_a.assign(APOE=name_a), part_b.assign(APOE=name_b)]).reset_index(drop=True)
    return pseudo, (name_a, name_b)


def n2_replicate(frame, labels, events, rep, args, arms=ARMS_FROZEN, stage="n2"):
    pseudo, names = n2_pseudo(frame, rep, args)
    out = {"stage": stage, "replicate": rep, "n": int(len(pseudo)),
           "composition": pd.crosstab(pseudo.APOE, pseudo.Diagnosis).to_dict()}
    for arm, estimator, consensus, scheme in arms:
        res = compare(pseudo, group_column="APOE", labels=labels, biomarkers=events,
                      group_order=names, estimator=estimator, consensus=consensus,
                      schemes=(scheme,), B=args.B, alpha=0.05, seed=args.seed * 1000 + rep,
                      workers=args.workers, stop_when_decided=True, stability_resamples=0, verbose=False)
        out[arm] = result_summary(res)          # one pair: reject iff E <= 29 at B=599
    return out


STAGES = {"n2": n2_replicate}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--stage", choices=tuple(STAGES), default="n2")
    ap.add_argument("--replicates", default="0:1", help="half-open range a:b of replicate indices")
    ap.add_argument("--B", type=int, default=599)
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--panel", choices=("P12", "P5", "P8"), default="P12",
                    help="event subset of the 12-event frame; same rows, same draws")
    args = ap.parse_args(argv)
    a, b = (int(x) for x in args.replicates.split(":"))
    out = args.out / args.stage
    out.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.table)
    frame, labels, events, _ = build_frame(table, "core", "P12", 180, "DXSUM")   # core population only
    if args.panel != "P12":
        events = [e for e in events if e in PANELS[args.panel]]   # same rows, same draws; per-event adjustment unchanged
    print(f"{args.stage}: frame n={len(frame)} labels={labels} panel={args.panel} events={len(events)} "
          f"replicates {a}:{b} definition {PLAN_VERSION}")
    if args.dry_run:
        pseudo, _ = n2_pseudo(frame, a, args)
        print(f"dry run: groups drawn, {pseudo.APOE.value_counts().to_dict()} disjoint rows; no fit performed")
        return 0
    for rep in range(a, b):
        path = out / f"rep_{rep:04d}.json"
        if path.exists():
            print(f"skip {path.name} (exists)")
            continue
        result = STAGES[args.stage](frame, labels, events, rep, args)
        result["B"] = args.B
        path.write_text(json.dumps(result, indent=1, default=str))
        print(f"wrote {path.name}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
