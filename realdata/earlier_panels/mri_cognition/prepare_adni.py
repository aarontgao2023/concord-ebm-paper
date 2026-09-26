#!/usr/bin/env python3
"""ADNI MRI-cognition population (1,025 participants; Supplementary Methods 2): write the input.

Earlier analysis population behind the ADNI MRI-cognition rows of Fig. 4e. Starts from the core
ADNI 12-event population of ../adni/build_cohort.py (same baseline, MCI = converters within 3
years, APOE groups) and keeps participants with all eight events (MMSE, ADSP-PHC memory and six
ADSP-PHC ComBat MRI volumes, each within 180 days of baseline) plus age, sex, education and ICV.
No regression or model fit is done here; run_shared_panel.correct() adjusts the events later.

  python prepare_adni.py --data-dir $CONCORD_ADNI_DIR --out $CONCORD_WORK_DIR/mri_cognition/adni

--derived-dir (default <data-dir>/derived) holds adni_p12_deid.csv and adni_cohort_full.csv from
../adni/build_cohort.py. --out must lie outside the repository. Only fresh pseudonyms go into
raw_input.csv; the identifying linkage stays in <out>/local_only. The script asserts the
population sizes of Supplementary Methods 2 (1,070 and 1,025); cell counts can be checked with
--expected-counts, a local JSON file {"core": [[..3 x 3..]], "final": [[..3 x 3..]]} (rows
e2, e3/e3, e4; columns CN, MCIc, AD). An existing output can be regenerated only when its bytes
are unchanged. The next step is freeze_campaign.py.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
OLD_BUILDER = HERE.parent / "adni/build_cohort.py"
OLD_RUNNER = HERE.parent / "adni/run_observed.py"
EVENTS = ["MMSE", "MEM", "Hippocampus", "Entorhinal", "Fusiform",
          "MiddleTemporal", "Precuneus", "Ventricles"]
MRI = EVENTS[2:]
MAP = {
    "Hippocampus": ("Left.Hippocampus_combat", "Right.Hippocampus_combat"),
    "Entorhinal": ("lh_entorhinal_volume_combat", "rh_entorhinal_volume_combat"),
    "Fusiform": ("lh_fusiform_volume_combat", "rh_fusiform_volume_combat"),
    "MiddleTemporal": ("lh_middletemporal_volume_combat", "rh_middletemporal_volume_combat"),
    "Precuneus": ("lh_precuneus_volume_combat", "rh_precuneus_volume_combat"),
    "Ventricles": ("Left.Lateral.Ventricle_combat", "Right.Lateral.Ventricle_combat"),
}
ICV_FIELD = "EstimatedTotalIntraCranialVol_combat"
GROUPS = ["e2", "e33", "e4"]
DIAGNOSES = ["CN", "MCIc", "AD"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def frozen_write(path, content, private=False):
    payload = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != payload:
            raise RuntimeError(f"Refusing to overwrite changed frozen file: {path.name}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    path.chmod(0o600 if private else 0o444)


def csv_text(frame):
    return frame.to_csv(index=False, lineterminator="\n", float_format="%.17g")


def counts(frame):
    return pd.crosstab(frame.APOE, frame.Diagnosis).reindex(
        index=GROUPS, columns=DIAGNOSES, fill_value=0).astype(int)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=os.environ.get("CONCORD_ADNI_DIR"),
                    help="folder with the ADNI files (default: $CONCORD_ADNI_DIR)")
    ap.add_argument("--derived-dir", type=Path, help="outputs of build_cohort.py (default: <data-dir>/derived)")
    ap.add_argument("--out", type=Path, help="output folder outside the repository "
                    "(default: $CONCORD_WORK_DIR/mri_cognition/adni)")
    ap.add_argument("--expected-counts", type=Path, help="optional local JSON with expected cell counts")
    args = ap.parse_args(argv)
    if args.out is None and os.environ.get("CONCORD_WORK_DIR"):
        args.out = Path(os.environ["CONCORD_WORK_DIR"]) / "mri_cognition/adni"
    if args.data_dir is None or args.out is None:
        ap.error("--data-dir and --out are required (or set CONCORD_ADNI_DIR and CONCORD_WORK_DIR)")
    DATA = args.data_dir.resolve()
    derived = (args.derived_dir or DATA / "derived").resolve()
    OUT = args.out.resolve()
    if OUT == REPO or REPO in OUT.parents:
        raise SystemExit(f"The output folder must lie outside the repository: {OUT}")
    PRIVATE = OUT / "local_only"
    DEID = derived / "adni_p12_deid.csv"
    FULL = derived / "adni_cohort_full.csv"
    COMBAT = DATA / "ADSP_PHC_T1_FS_18Oct2025.csv"
    expected = json.loads(args.expected_counts.read_text()) if args.expected_counts else {}
    deid = pd.read_csv(DEID)
    full = pd.read_csv(FULL).sort_values("RID").reset_index(drop=True)
    assert len(deid) == len(full)
    assert full.RID.is_unique and deid.PTID.is_unique
    # Verify the original builder's deterministic pseudonym-to-full-row mapping.
    assert deid.PTID.tolist() == [f"a{i:05d}" for i in range(len(full))]
    for col, raw_col in [("Diagnosis", "Diagnosis"), ("APOE", "G"), ("cohort", "cohort"), ("Sex", "Sex")]:
        assert deid[col].fillna("").equals(full[raw_col].fillna("")), col
    assert np.allclose(deid.Age, full.Age.round(2), equal_nan=True)
    deid["_RID"] = full.RID.to_numpy()
    core = deid[(deid.cohort == "core") & deid.APOE.isin(GROUPS)].copy()
    assert len(core) == 1070
    if "core" in expected:
        assert counts(core).values.tolist() == expected["core"]
    assert set(core.Diagnosis) == set(DIAGNOSES)

    # Reconstruct ONLY the old ComBat measurement choice to verify provenance,
    # including the old first-record tie convention. Do not choose a newer scan.
    source_cols = sorted({v for pair in MAP.values() for v in pair} | {ICV_FIELD, "BrainSegVolNotVentSurf_combat"})
    scans = pd.read_csv(COMBAT, usecols=["RID", "PHC_SCANDATE"] + source_cols)
    scans["_source_csv_row"] = np.arange(len(scans)) + 2
    scans["PHC_SCANDATE"] = pd.to_datetime(scans.PHC_SCANDATE, errors="coerce")
    for col in source_cols:
        scans[col] = pd.to_numeric(scans[col], errors="coerce")
    scans = scans.dropna(subset=["PHC_SCANDATE", MAP["Hippocampus"][0]])
    base = full[["RID", "BLDATE"]].copy()
    base["BLDATE"] = pd.to_datetime(base.BLDATE)
    scans = scans.merge(base, on="RID", validate="many_to_one")
    scans["delta"] = (scans.PHC_SCANDATE - scans.BLDATE).dt.days
    scans = scans[scans.delta.abs().le(365)].copy()
    scans["absd"] = scans.delta.abs()
    selected = scans.sort_values(["RID", "absd"]).drop_duplicates("RID").set_index("RID")
    selected_core = selected.reindex(core._RID)
    assert np.allclose(selected_core.delta, core.ComBat_offset_days, equal_nan=True)
    for event, pair in MAP.items():
        # Match the stored old mean exactly, then independently audit both sides.
        reconstructed = selected_core[list(pair)].mean(axis=1)
        assert np.allclose(reconstructed, core[event + "_ComBat"], atol=1e-8, rtol=1e-12, equal_nan=True), event
    assert np.allclose(selected_core[ICV_FIELD], core.ICV_ComBat, atol=1e-8, rtol=1e-12, equal_nan=True)

    frame = core.copy()
    for event in MRI:
        frame[event] = frame[event + "_ComBat"]
    frame["ICV"] = frame.ICV_ComBat
    for offset, events in [("ComBat_offset_days", MRI), ("MMSE_offset_days", ["MMSE"]), ("MEM_offset_days", ["MEM"])]:
        frame.loc[~frame[offset].abs().le(180), events] = np.nan
    features_complete = np.isfinite(frame[EVENTS]).all(axis=1)
    numeric_complete = np.isfinite(frame[EVENTS + ["Age", "Education", "ICV"]]).all(axis=1)
    complete = numeric_complete & frame.Sex.isin(["Male", "Female"])
    # The existing mean(skipna=True) is safe here only if both hemispheres exist.
    bilateral_ok = np.isfinite(selected_core[[v for p in MAP.values() for v in p]]).all(axis=1).to_numpy()
    assert bilateral_ok[complete.to_numpy()].all(), "Unilateral source value in a complete candidate"
    final = frame.loc[complete].copy().sort_values("PTID").reset_index(drop=True)
    assert len(final) == 1025
    if "final" in expected:
        assert counts(final).values.tolist() == expected["final"]
    assert final.MMSE.between(0, 30).all()
    assert (final.ICV > 0).all() and (final[MRI] > 0).all().all()
    assert (final.Age > 0).all() and (final.Education >= 0).all()

    new_ids = [f"adni_shared_{i:06d}" for i in range(1, len(final) + 1)]
    raw = final[["PTID", "APOE", "Diagnosis", "Age", "Sex", "Education", "ICV"] + EVENTS].copy()
    raw["PTID"] = new_ids
    # No regression, direction flip, imputation, resampling or model fit here.
    raw_path = OUT / "raw_input.csv"
    frozen_write(raw_path, csv_text(raw))
    final_full = full.set_index("RID").loc[final._RID]
    link = pd.DataFrame({"PTID": new_ids, "source_pseudonym": final.PTID.to_numpy(),
                         "source_RID": final._RID.to_numpy(), "source_PTID": final_full.PTID.to_numpy(),
                         "baseline_date": final_full.BLDATE.to_numpy(),
                         "selected_MRI_date": selected.loc[final._RID, "PHC_SCANDATE"].dt.strftime("%Y-%m-%d").to_numpy(),
                         "selected_MRI_source_csv_row": selected.loc[final._RID, "_source_csv_row"].astype(int).to_numpy(),
                         "MRI_offset_days": final.ComBat_offset_days.to_numpy(),
                         "MMSE_offset_days": final.MMSE_offset_days.to_numpy(),
                         "MEM_offset_days": final.MEM_offset_days.to_numpy()})
    frozen_write(PRIVATE / "subject_linkage.csv", csv_text(link), private=True)
    cell_table = counts(final)
    cells = cell_table.rename_axis("APOE").reset_index()
    frozen_write(OUT / "group_diagnosis_counts.csv", csv_text(cells))
    old_a, old_b = [69, 5, 12], [56, 58, 86]
    capacities = cell_table.loc["e33"].tolist()
    n2 = {d: {"available_e33": capacity, "legacy_arm_a": a, "legacy_arm_b": b,
              "legacy_required_total": a + b, "legacy_feasible": bool(a + b <= capacity)}
          for d, capacity, a, b in zip(DIAGNOSES, capacities, old_a, old_b)}
    source_paths = [DEID, FULL, COMBAT, OLD_BUILDER, OLD_RUNNER,
                    DATA / "MMSE_18Jan2026.csv", DATA / "ADSP_PHC_COGN_14Sep2026.csv",
                    DATA / "DXSUM_27Aug2025.csv", DATA / "APOERES_14Sep2026.csv",
                    DATA / "PTDEMOG_27Aug2025.csv", Path(__file__)]
    hashes = {str(p.relative_to(REPO)) if p.is_relative_to(REPO) else p.name: sha(p) for p in source_paths}
    manifest = {
        "cohort": "ADNI", "release": "nc_shared_panel_20260920", "n": len(raw),
        "input_csv": "raw_input.csv", "input_sha256": sha(raw_path), "input_columns": raw.columns.tolist(),
        "panels": {"K4": EVENTS[:4], "K8": EVENTS}, "same_people_for_both_panels": True,
        "diagnosis_labels": DIAGNOSES, "group_order": GROUPS,
        "selection": {"base": "Existing core ADNI1/GO/2, original DXSUM baseline and MCIc converter definition, e2/e33/e4 only; unchanged",
                      "MRI": "Original nearest PHC_SCANDATE ComBat record; absolute baseline lag <=180 days; both hemispheres finite",
                      "cognition": "Original nearest MMSE and PHC_MEM records, each absolute baseline lag <=180 days",
                      "completeness": "All 8 events plus Age, Sex, Education and ComBat ICV; no outcome-dependent exclusions",
                      "MRI_ICV_source_field": ICV_FIELD, "MRI_fields": {k: list(v) for k, v in MAP.items()},
                      "covariate_age": "Existing frozen deidentified Age rounded to two decimals; no recomputation",
                      "MRI_scale": "Existing PHC ComBat volumes, bilateral arithmetic means; no extra scaling"},
        "flow": {"core": 1070, "MRI_record_within180": int(core.ComBat_offset_days.abs().le(180).sum()),
                 "complete8": int(features_complete.sum()), "complete8_and_covariates": len(raw),
                 "excluded_from_core": int(1070 - len(raw)), "unilateral_candidates_removed": 0},
        "cells": {g: {d: int(cell_table.loc[g, d]) for d in DIAGNOSES} for g in GROUPS},
        "preprocessing_state": "RAW selected PHC/ComBat values; no regression or sign flip applied",
        "agreed_runner_preprocessing": {"fit_population": "Pooled CN in this complete K8 cohort, once before resampling/permutations",
                 "cognition_factors": ["Age", "SexMale", "Education"], "MRI_factors": ["Age", "SexMale", "ICV"],
                 "rule": "Joint CN OLS with intercept; y_adjusted = y - (X - mean_full_complete_cohort) @ beta_slopes",
                 "K4": "Select K4 columns from the already-adjusted K8 table; do not re-estimate corrections",
                 "direction": "No explicit sign flip; CN/AD measurement initialization determines direction"},
        "conventional_abnormal_direction": {e: "higher" if e == "Ventricles" else "lower" for e in EVENTS},
        "legacy_N2_capacity": n2,
        "legacy_N2_joint_draw_feasible": all(x["legacy_feasible"] for x in n2.values()),
        "legacy_N2_note": "Legacy N2 AD counts require 98 disjoint e33 records, but only 93 remain. Do not silently sample with replacement or reuse original draw counts.",
        "privacy": {"upload_allowlist": ["raw_input.csv", "manifest.json", "group_diagnosis_counts.csv", "validation.json", "README.md"],
                    "local_only": "local_only/subject_linkage.csv; never upload or display individual rows",
                    "identifier": "Fresh sequential analysis pseudonym; original RID/PTID and dates only in local linkage"},
        "source_sha256": hashes,
    }
    validation = {"source_core_counts_match": True, "source_ComBat_offsets_match": True,
                  "source_ComBat_bilateral_means_match": True, "both_hemispheres_finite": True,
                  "all_8_events_and_covariates_complete": True, "unique_fresh_pseudonyms": bool(raw.PTID.is_unique),
                  "MRI_and_ICV_positive": True, "MMSE_in_0_30": True, "new_model_fits": 0,
                  "new_simulations": 0, "new_permutations": 0, "n": len(raw)}
    frozen_write(OUT / "manifest.json", json.dumps(manifest, indent=2) + "\n")
    frozen_write(OUT / "validation.json", json.dumps(validation, indent=2) + "\n")
    cell_text = ", ".join(f"{g}:" + "/".join(str(int(cell_table.loc[g, d])) for d in DIAGNOSES) for g in GROUPS)
    readme = f"""# ADNI common MRI/cognition frozen input

`raw_input.csv` contains 1,025 complete participants from the unchanged original
ADNI core. K4 is the first four event columns of K8. MRI values and ICV use the
existing PHC ComBat source and its original visit selection, restricted to
absolute baseline lag <=180 days. Cognition retains its original ±180-day rule.
MCIc remains the original three-year converter definition.

No CN regression or sign reversal has been applied. The shared runner must use
the joint pooled-CN correction specified in `manifest.json`, once on this full
complete cohort; K4 then subsets the adjusted K8 columns.

Genotype by CN/MCIc/AD counts are {cell_text}.
The draw counts of the 12-event population cannot be copied unchanged: their two
AD cells require 98 disjoint e33 profiles, but this pool contains 93. The draws
are therefore scaled by 93/98 in freeze_campaign.py.

Only the manifest's upload allowlist may be copied to HPC. The `local_only/`
directory contains identifying linkage and acquisition dates; never upload it.
Outputs are protected against changed-byte overwrites; SHA256 values in the
manifest identify the exact raw input and sources. Running the preparation
again verifies and retains byte-identical files.
"""
    frozen_write(OUT / "README.md", readme)
    print(json.dumps({"n": len(raw), "cells": manifest["cells"], "input_sha256": manifest["input_sha256"],
                      "legacy_N2_joint_draw_feasible": False, "validation_passed": True}))


if __name__ == "__main__":
    main()
