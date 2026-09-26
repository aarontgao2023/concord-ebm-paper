#!/usr/bin/env python
"""NACC amyloid PET population: analysis frame and comparisons of two groups drawn from e3/e3
carriers (Fig. 4e, NACC amyloid PET-cognition setting).

build_frame() adjusts the four cognitive scores for age, sex and education and the Centiloid
value for age and sex, by OLS in the CN participants, once on the full frame (the same rule as
the ADNI 12-event population). run_nacc_stability.py uses the same frame.

  --stage n2   one JSON per draw (--replicates a:b): two disjoint groups drawn without
               replacement, within diagnosis, from the e3/e3 participants. Group A has the
               composition of the e2 group and group B the e4 proportions, scaled to the largest
               multiple of 10 that the remaining e3/e3 participants allow (n2_sizes). Each draw is
               tested once as a two-group comparison: P = (1 + E)/(B + 1), B = 599, rejection if
               E <= 29; permutations stop once the decision is fixed.
  --stage a0   descriptive counts only (no fits); not reported in the paper

Four fits per draw (output keys); Fig. 4e uses standard_U and concord_D:
  standard_U     separately fitted DEBM, unmodified pyebm search, unrestricted permutation
  standardrep_D  separately fitted DEBM, continued search, within-diagnosis permutation
  concord_U      CONCORD, unrestricted permutation
  concord_D      CONCORD (common proportions = smallest proportion of each diagnosis over the two
                 drawn groups; continued search), within-diagnosis permutation

Draw i uses numpy.random.default_rng([seed, 2, i]) and permutation seed seed*1000 + i
(seed 20260917). Requires the concord-ebm package (version 0.1.1, commit 40fc046) and pyebm 2.0.3.

  python run_nacc.py --table <derived>/nacc_p5_deid.csv --out $CONCORD_WORK_DIR/earlier/nacc \\
      --stage n2 --replicates 0:200
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from concord import compare, prepare_data
from concord.invariant import composition_weights

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "adni"))
from run_null import result_summary, draw_within  # noqa: E402  (same draw and summary code as ADNI)

EVENTS = ["CENTILOID", "MOCA", "MEM", "EXF", "LAN"]
COGNITIVE = ["MOCA", "MEM", "EXF", "LAN"]
GROUP_ORDER = ("e2", "e33", "e4")
LABELS = ("CN", "MCI", "AD")
DIAG_CODE = {"CN": 1, "MCI": 2, "AD": 3}
# (label, estimator, search, permutation scheme) as passed to concord.compare
ARMS = [("standard_U", "standard", "original", "unrestricted"),
        ("standardrep_D", "standard", "repaired", "diagnosis"),
        ("concord_U", "invariant_min", "repaired", "unrestricted"),
        ("concord_D", "invariant_min", "repaired", "diagnosis")]


def build_frame(table: pd.DataFrame):
    df = table[table.APOE.isin(GROUP_ORDER)].copy()
    df["SexMale"] = (df.Sex == "Male").astype(float)
    control = (df.Diagnosis == "CN").to_numpy()
    log = {}

    def adjust(cols, factors, tag):
        X_all = df[factors].to_numpy(dtype=float)
        dev = X_all - np.nanmean(X_all, axis=0)
        dev[np.isnan(dev)] = 0.0
        for c in cols:
            y = df[c].to_numpy(dtype=float)
            rows = control & ~np.isnan(y) & ~np.isnan(X_all).any(axis=1)
            X = np.column_stack([np.ones(rows.sum()), X_all[rows]])
            beta, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
            df[c] = y - dev @ beta[1:]
            log[f"{tag}:{c}"] = {"n_control": int(rows.sum()), "slopes": dict(zip(factors, map(float, beta[1:])))}

    adjust(COGNITIVE, ["Age", "SexMale", "Education"], "cognition")
    adjust(["CENTILOID"], ["Age", "SexMale"], "amyloid")
    return df[["PTID", "Diagnosis", "APOE"] + EVENTS].reset_index(drop=True), log


def n2_sizes(frame):
    tab = pd.crosstab(frame.APOE, frame.Diagnosis).reindex(index=list(GROUP_ORDER), columns=list(LABELS))
    e2 = {c: int(tab.loc["e2", c]) for c in LABELS}
    e4p = tab.loc["e4"] / tab.loc["e4"].sum()
    remaining = tab.loc["e33"] - pd.Series(e2)
    n4 = (int(min(remaining[c] / e4p[c] for c in LABELS)) // 10) * 10
    return e2, {c: int(round(n4 * e4p[c])) for c in LABELS}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, type=Path, help="nacc_p5_deid.csv from build_cohort_nacc.py")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--stage", choices=("a0", "n2"), required=True)
    ap.add_argument("--replicates", default="0:1")
    ap.add_argument("--B", type=int, default=599)
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    frame, log = build_frame(pd.read_csv(args.table))
    prepare_data(frame, group_column="APOE", labels=LABELS, biomarkers=EVENTS, group_order=GROUP_ORDER)
    print(f"NACC frame n={len(frame)} events={EVENTS}")
    (args.out / "meta_nacc.json").write_text(json.dumps({"n": int(len(frame)), "events": EVENTS, "labels": LABELS,
                                                         "B": args.B, "seed": args.seed, "correction": log}, indent=2))
    if args.stage == "a0":
        tab = pd.crosstab(frame.APOE, frame.Diagnosis).reindex(index=list(GROUP_ORDER), columns=list(LABELS))
        chi2, p, dof, _ = stats.chi2_contingency(tab.to_numpy())
        diag = frame.Diagnosis.map(DIAG_CODE).to_numpy()
        _, _, ess_min = composition_weights(diag, frame.APOE.to_numpy(), np.ones(len(frame), bool), "min")
        _, _, ess_pool = composition_weights(diag, frame.APOE.to_numpy(), np.ones(len(frame), bool), "pooled")
        e2, e4 = n2_sizes(frame)
        desc = {"table": tab.to_dict(), "chi2": chi2, "dof": dof, "p": p, "ess_min": ess_min, "ess_pooled": ess_pool,
                "availability": {g: {e: int(frame.loc[frame.APOE == g, e].notna().sum()) for e in EVENTS} for g in GROUP_ORDER},
                "n2_draws": {"pseudo_e2": e2, "pseudo_e4": e4}}
        (args.out / "a0_nacc.json").write_text(json.dumps(desc, indent=2, default=float))
        print(tab.to_string()); print(f"chi2 {chi2:.1f} p {p:.2e}; ESS min {ess_min}; draw sizes {desc['n2_draws']}")
        return 0
    if args.dry_run:
        print("dry run: inputs validated; no fit performed")
        return 0
    a, b = (int(x) for x in args.replicates.split(":"))
    outdir = args.out / args.stage
    outdir.mkdir(exist_ok=True)
    e2_draw, e4_draw = n2_sizes(frame)
    for rep in range(a, b):
        path = outdir / f"rep_{rep:04d}.json"
        if path.exists():
            continue
        rng = np.random.default_rng([args.seed, {"n2": 2}[args.stage], rep])
        pa = draw_within(frame, "e33", e2_draw, rng)
        pb = draw_within(frame, "e33", e4_draw, rng, exclude=pa.index)
        data = pd.concat([pa.assign(APOE="pA_e2comp"), pb.assign(APOE="pB_e4comp")]).reset_index(drop=True)
        order, alpha = ("pA_e2comp", "pB_e4comp"), 0.05
        out = {"stage": args.stage, "replicate": rep, "n": int(len(data)),
               "composition": pd.crosstab(data.APOE, data.Diagnosis).to_dict()}
        for arm, estimator, consensus, scheme in ARMS:
            res = compare(data, group_column="APOE", labels=LABELS, biomarkers=EVENTS, group_order=order,
                          estimator=estimator, consensus=consensus, schemes=(scheme,), B=args.B, alpha=alpha,
                          seed=args.seed * 1000 + rep, workers=args.workers, stop_when_decided=True,
                          stability_resamples=0, verbose=False)
            out[arm] = result_summary(res)
        path.write_text(json.dumps(out, indent=1, default=str))
        print(f"wrote {path.name}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
