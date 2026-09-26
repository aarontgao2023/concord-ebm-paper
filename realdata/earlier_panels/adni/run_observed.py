#!/usr/bin/env python
"""ADNI 12-event population: analysis frame and bootstrap refits (Supplementary Fig. 3).

Input is the de-identified wide table written by build_cohort.py. build_frame() selects the rows,
labels, events and measurement window and applies the covariate adjustment once; run_null.py uses
the same frame for the ADNI rows of Fig. 4e.

  --stage a2s      bootstrap refits: R resamples within genotype-by-diagnosis cells of the frame,
                   each refitted from scratch (pooled mixture, weights and ordering; the covariate
                   adjustment stays fixed). Two fits per resample: concord_min (CONCORD, smallest-
                   proportion common proportions of the three genotype groups, continued search)
                   and standard_original (separately fitted DEBM, unmodified pyebm search). The
                   paper uses concord_min with --R 500 (Supplementary Fig. 3).
  --stage a0       descriptive counts only (no fits): genotype-by-diagnosis table, chi-square,
                   event availability and effective sample sizes; not reported in the paper
  --dry-run        validate the inputs up to the first model fit, then stop

Options not used in the paper are kept for completeness: --cohort S1 (all MCI), the panels P11,
P12-FS7 and P12-ComBat, --window 365, --diagnosis-source PHC and --correction paper.

Outputs are JSON aggregates under --out (a2s_<tag>.json holds the resampled orderings, meta_<tag>.json
the frame size and the adjustment coefficients); no per-participant data are written.
Requires the concord-ebm package (version 0.1.1, commit 40fc046) and pyebm 2.0.3.

  python run_observed.py --table <derived>/adni_p12_deid.csv --out $CONCORD_WORK_DIR/earlier/adni/main \\
      --stage a2s --R 500
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from concord import prepare_data
from concord.core import CompareSpec, fit_once, resample_within_strata
from concord.invariant import composition_weights

PLAN_VERSION = "v2 2026-09-14"   # version of the population definition
EVENTS = ["ABETA", "TAU", "PTAU", "MMSE", "MEM", "Hippocampus", "Entorhinal", "Fusiform",
          "MiddleTemporal", "Precuneus", "WholeBrain", "Ventricles"]
VOLUMES = EVENTS[5:]
COGNITIVE = ["MMSE", "MEM"]
PANELS = {"P12": EVENTS, "P11": [e for e in EVENTS if e != "MEM"],
          "P12-FS7": EVENTS, "P12-ComBat": EVENTS,
          # reduced panels of the 12-event frame (Fig. 4e, ADNI 5- and 8-event settings)
          "P5": ["ABETA", "MMSE", "MEM", "Hippocampus", "Entorhinal"],
          "P8": ["ABETA", "TAU", "PTAU", "MMSE", "MEM", "Hippocampus", "Entorhinal", "WholeBrain"]}
GROUP_ORDER = ("e2", "e33", "e4")
DIAG_CODE = {"CN": 1, "MCI": 2, "MCIc": 2, "AD": 3}


# ----------------------------------------------------------------------------- frame construction

def select_cohort(table: pd.DataFrame, cohort: str, panel: str, window: int, diagnosis_source: str):
    """Rows, labels, panel columns and the measurement window."""
    df = table.copy()
    if cohort == "core":
        df = df[df.cohort == "core"]
        labels = ("CN", "MCIc", "AD")
    elif cohort == "S1":                                   # all MCI, three-label pyebm default
        df = df.copy()
        df["Diagnosis"] = df.Diagnosis.replace({"MCIc": "MCI"})
        labels = ("CN", "MCI", "AD")
    else:
        raise ValueError(cohort)
    if diagnosis_source == "PHC":                          # PHC harmonised diagnosis at baseline
        code_to_label = {1: "CN", 2: "MCI", 3: "AD"}
        df = df[df.PHC_Diagnosis_bl.isin([1, 2, 3])].copy()
        phc = df.PHC_Diagnosis_bl.map(code_to_label)
        if cohort == "core":                               # keep the converter definition for MCI
            keep = (phc == df.Diagnosis.replace({"MCIc": "MCI"}))
            df = df[keep]
        else:
            df["Diagnosis"] = phc
    df = df[df.APOE.isin(GROUP_ORDER)].copy()

    # panel: volume source
    if panel == "P12-FS7":
        for e in VOLUMES:
            df[e] = df[f"{e}_FS7"]
        df["ICV"] = df["ICV_FS7"]
        off = df["FS7_offset_days"].abs()
        df.loc[~off.le(window), VOLUMES] = np.nan
    elif panel == "P12-ComBat":
        for e in VOLUMES:
            df[e] = df[f"{e}_ComBat"]
        df["ICV"] = df["ICV_ComBat"]
        off = df["ComBat_offset_days"].abs()
        df.loc[~off.le(window), VOLUMES] = np.nan
    # FS6 (P12 / P11): BL_ columns, no date matching

    # window for dated modalities
    for offset, cols in [("CSF_offset_days", ["ABETA", "TAU", "PTAU"]),
                         ("MMSE_offset_days", ["MMSE"]), ("MEM_offset_days", ["MEM"])]:
        df.loc[~df[offset].abs().le(window), cols] = np.nan
    events = PANELS[panel]
    df = df[df[events].notna().any(axis=1)].copy()        # at least one event observed
    return df.reset_index(drop=True), labels, events


def confounder_correction(df: pd.DataFrame, events: list, labels: tuple, mode: str = "authors"):
    """Covariate adjustment, applied once on the full analysis frame, outside every resampling loop.

    pyebm's CorrectConfounders rule: OLS of each biomarker on the factors in controls (first label),
    then subtract (factor - mean over all rows) . slope, with a missing factor deviation set to 0.
    Because controls are the same set under every relabelling of the genotype groups, applying the
    rule once equals pyebm applying it inside every permutation fit. For the resampling analyses
    (the draws of run_null.py and the bootstrap refits) the adjustment is therefore fixed from the
    full frame; it is not re-estimated inside each subsample.

    mode='authors' (used in the paper): Factors=['Age','Sex','ICV'] on all events, after the two
    cognitive events are first adjusted for education by the same rule (sequential; not a joint
    regression).
    mode='paper' (not used in the paper): modality-specific joint regressions, volumes ~ age + sex +
    ICV; CSF ~ age + sex; cognition ~ age + sex + education.
    """
    out = df.copy()
    out["SexMale"] = (out.Sex == "Male").astype(float)
    control = out.Diagnosis == labels[0]
    log = {"mode": mode}

    def adjust(cols, factors, tag):
        X_all = out[factors].to_numpy(dtype=float)
        means = np.nanmean(X_all, axis=0)
        dev = X_all - means
        dev[np.isnan(dev)] = 0.0
        for c in cols:
            y = out[c].to_numpy(dtype=float)
            rows = control.to_numpy() & ~np.isnan(y) & ~np.isnan(X_all).any(axis=1)
            X = np.column_stack([np.ones(rows.sum()), X_all[rows]])
            beta, *_ = np.linalg.lstsq(X, y[rows], rcond=None)
            out[c] = y - dev @ beta[1:]
            log[f"{tag}:{c}"] = {"n_control": int(rows.sum()), "factors": list(factors),
                                 "slopes": dict(zip(factors, map(float, beta[1:])))}

    if mode == "authors":
        adjust([c for c in COGNITIVE if c in events], ["Education"], "education")
        adjust(events, ["Age", "SexMale", "ICV"], "age_sex_icv")
    elif mode == "paper":
        adjust([c for c in VOLUMES if c in events], ["Age", "SexMale", "ICV"], "volumes")
        adjust([c for c in ("ABETA", "TAU", "PTAU") if c in events], ["Age", "SexMale"], "csf")
        adjust([c for c in COGNITIVE if c in events], ["Age", "SexMale", "Education"], "cognition")
    else:
        raise ValueError(mode)
    return out, log


def build_frame(table, cohort, panel, window, diagnosis_source, correction="authors"):
    df, labels, events = select_cohort(table, cohort, panel, window, diagnosis_source)
    df, correction_log = confounder_correction(df, events, labels, correction)
    frame = df[["PTID", "Diagnosis", "APOE"] + events].copy()
    return frame, labels, events, correction_log


def missingness_diagnostics(df: pd.DataFrame, events: list) -> dict:
    """Genotype x diagnosis x event availability, observed-count distribution and the most
    frequent missingness patterns, per group. Counts only."""
    out = {"availability_by_group_dx": {}, "observed_count_by_group_dx": {}, "patterns_by_group": {},
           "strata_diagnosis_count": {}}
    obs = df[events].notna()
    count = obs.sum(axis=1)
    for g in GROUP_ORDER:
        sel = df.APOE == g
        out["availability_by_group_dx"][g] = {d: {e: int(obs.loc[sel & (df.Diagnosis == d), e].sum())
                                                 for e in events} for d in sorted(df.Diagnosis.unique())}
        out["observed_count_by_group_dx"][g] = {d: count[sel & (df.Diagnosis == d)].value_counts().sort_index().to_dict()
                                                for d in sorted(df.Diagnosis.unique())}
        pat = obs[sel].apply(lambda r: "".join("1" if v else "0" for v in r), axis=1).value_counts().head(8)
        out["patterns_by_group"][g] = {k: int(v) for k, v in pat.items()}
    strata = df.Diagnosis.astype(str) + ":" + count.astype(str)
    tab = pd.crosstab(strata, df.APOE)
    out["strata_diagnosis_count"] = {s: {g: int(tab.loc[s, g]) if g in tab.columns else 0 for g in GROUP_ORDER}
                                     for s in tab.index}
    out["n_strata"] = int(len(tab)); out["n_strata_single_group"] = int(((tab > 0).sum(axis=1) == 1).sum())
    out["event_names"] = events
    return out


# ----------------------------------------------------------------------------- a0 (no fits)

def a0_descriptives(frame, labels, events):
    tab = pd.crosstab(frame.APOE, frame.Diagnosis).reindex(index=list(GROUP_ORDER), columns=list(labels)).fillna(0).astype(int)
    chi2, p, dof, _ = stats.chi2_contingency(tab.to_numpy())
    avail = {g: {e: int(frame.loc[frame.APOE == g, e].notna().sum()) for e in events} for g in GROUP_ORDER}
    diag = frame.Diagnosis.map(DIAG_CODE).to_numpy()
    gv = frame.APOE.to_numpy()
    mask = np.ones(len(frame), dtype=bool)
    ess = {}
    weights_summary = {}
    for ref in ("min", "pooled"):
        gvals, w, e = composition_weights(diag, gv, mask, ref)
        ess[ref] = e
        weights_summary[ref] = {str(g): {"min": float(w[g].min()), "max": float(w[g].max()),
                                         "unique": sorted(set(np.round(w[g], 4).tolist()))} for g in gvals}
    comp = tab.div(tab.sum(axis=1), axis=0).round(3)
    return {"n": int(len(frame)), "table": tab.to_dict(), "proportions": comp.to_dict(),
            "chi2": float(chi2), "chi2_dof": int(dof), "chi2_p": float(p),
            "availability": avail, "ess": ess, "weights": weights_summary}


# ----------------------------------------------------------------------------- a2s (bootstrap refits)

def a2s_stability(frame, labels, events, R, seed, workers_unused):
    """Resampling within genotype-by-diagnosis cells; both models refitted in every resample."""
    prepared, group_values, names = prepare_data(frame, group_column="APOE", labels=labels,
                                                  biomarkers=events, group_order=GROUP_ORDER)
    out = {"R": R, "events": list(names), "groups": list(GROUP_ORDER), "arms": {}}
    for arm, estimator, consensus in (("standard_original", "standard", "original"),
                                      ("concord_min", "invariant_min", "repaired")):
        spec = CompareSpec(group_column="APOE", group_values=group_values, labels=labels, biomarkers=names,
                           estimator=estimator, consensus=consensus, seed=seed)
        orderings, failed = [], 0
        for r in range(R):
            sub = resample_within_strata(prepared, spec, r)
            fit = fit_once(sub, spec)
            if not fit.ok:
                failed += 1
                continue
            orderings.append([[names[i] for i in o] for o in fit.orderings])
            if (r + 1) % 25 == 0:
                print(f"A2s {arm} {r + 1}/{R}", flush=True)
        out["arms"][arm] = {"completed": len(orderings), "failed": failed, "orderings": orderings,
                            "note": "each resample refits the whole pipeline (pooled mixture, weights, consensus); "
                                    "correction parameters frozen from the full frame"}
    return out


# ----------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, type=Path, help="adni_p12_deid.csv from build_cohort.py")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--stage", choices=("a0", "a2s"), required=True)
    ap.add_argument("--correction", choices=("authors", "paper"), default="authors")
    ap.add_argument("--cohort", choices=("core", "S1"), default="core")
    ap.add_argument("--panel", choices=tuple(PANELS), default="P12")
    ap.add_argument("--window", type=int, choices=(180, 365), default=180)
    ap.add_argument("--diagnosis-source", choices=("DXSUM", "PHC"), default="DXSUM")
    ap.add_argument("--B", type=int, default=599, help="recorded in the meta file only")
    ap.add_argument("--R", type=int, default=200, help="number of bootstrap refits (the paper used 500)")
    ap.add_argument("--seed", type=int, default=20260914)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    table = pd.read_csv(args.table)
    frame, labels, events, correction_log = build_frame(table, args.cohort, args.panel, args.window,
                                                        args.diagnosis_source, args.correction)
    tag = f"{args.cohort}_{args.panel}_w{args.window}_{args.diagnosis_source}"
    if args.correction != "authors":
        tag += f"_corr{args.correction}"
    meta = {"plan": PLAN_VERSION, "cohort": args.cohort, "panel": args.panel, "window": args.window,
            "diagnosis_source": args.diagnosis_source, "correction": args.correction, "labels": labels,
            "events": events, "n": int(len(frame)), "B": args.B, "seed": args.seed,
            "correction_log": correction_log}
    (args.out / f"meta_{tag}.json").write_text(json.dumps(meta, indent=2))
    print(f"frame: n={len(frame)} labels={labels} events={len(events)} tag={tag}")

    # input validation identical to what compare() does, without fitting
    prepare_data(frame, group_column="APOE", labels=labels, biomarkers=events, group_order=GROUP_ORDER)

    if args.stage == "a0":
        desc = a0_descriptives(frame, labels, events)
        desc["missingness"] = missingness_diagnostics(frame, events)
        (args.out / f"a0_{tag}.json").write_text(json.dumps(desc, indent=2))
        miss = desc["missingness"]
        print(f"diagnosis x observed-count strata: {miss['n_strata']} (single-group strata: {miss['n_strata_single_group']})")
        for g in GROUP_ORDER:
            print(f"  {g} top missingness patterns (1=observed, order {events}):", miss["patterns_by_group"][g])
        tab =pd.DataFrame(desc["table"]).reindex(index=list(GROUP_ORDER))
        print("composition (counts):\n" + tab.to_string())
        print(f"chi-square = {desc['chi2']:.1f}, dof {desc['chi2_dof']}, p = {desc['chi2_p']:.2e}")
        print("ESS min-rule proportions:", {k: round(v, 1) for k, v in desc["ess"]["min"].items()})
        print("ESS pooled proportions:", {k: round(v, 1) for k, v in desc["ess"]["pooled"].items()})
        print("availability by group:\n" + pd.DataFrame(desc["availability"]).to_string())
        return 0

    if args.dry_run:
        print("dry run: inputs validated; no fit performed")
        return 0

    if args.stage == "a2s":
        res = a2s_stability(frame, labels, events, args.R, args.seed, args.workers)
        (args.out / f"a2s_{tag}.json").write_text(json.dumps(res))
        print({arm: {k: v for k, v in a.items() if k != "orderings"} for arm, a in res["arms"].items()})
        return 0


if __name__ == "__main__":
    sys.exit(main())
