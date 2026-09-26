#!/usr/bin/env python
"""Build the ADNI 12-event population (1,070 participants; Supplementary Methods 2).

Earlier analysis population behind the ADNI rows of Fig. 4e and the ADNI panel of Supplementary
Fig. 3. Baseline is each participant's first DXSUM diagnosis in ADNI1, ADNIGO or ADNI2; MCI
participants are those who progressed to dementia within 3 years (label MCIc); APOE e2 = e2/e2 or
e2/e3, e3/e3, e4 = e3/e4 or e4/e4 (e2/e4 excluded). Events: CSF Abeta42, total tau and p-tau181
(ADSP-PHC), MMSE, ADSP-PHC memory composite and seven FreeSurfer 6.0 volumes; incomplete profiles
are allowed. The definition was fixed on 14 September 2026.

Reads the ADNI files listed in FILES from --data-dir (default $CONCORD_ADNI_DIR) and writes to
--out-dir (default <data-dir>/derived; keep it outside the repository):

  adni_cohort_full.csv      one row per participant of the core population or of the MCI
                            non-converter extension (not used in the paper); RID/PTID/dates
                            kept; stays on this machine
  adni_p12_deid.csv         de-identified wide table for the compute cluster: sequential id,
                            Diagnosis, APOE, Age, Sex, ICV, Education, 12 events + sensitivity
                            volume sources and measurement offsets; no RID/PTID/dates
  fs6_baseline_extract.csv  cached columns from the FreeSurfer 6.0 workbook (RID-keyed, local)
  build_cohort_report.txt   counts and availability only (no identifiers)

The script exits non-zero if the population does not have 1,070 participants. More detailed
expected counts can be checked with --expected-counts, a local JSON file (keep it out of the
repository) of the form {"table": {"e2": {"CN": n, "MCIc": n, "AD": n}, "e33": {...},
"e4": {...}}, "avail_180": {"MMSE": n, "FS6": n, "CSF": n}}.

Nothing here fits an ordering or computes a test statistic.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PLAN_VERSION = "v2 2026-09-14"   # version of the population definition
THREE_YEARS_DAYS = 1095.75
WINDOW_DAYS_MAX = 365          # store the closest measurement within this; run scripts apply 180/365

FILES = {
    "dxsum": "DXSUM_27Aug2025.csv",
    "demodx": "ADSP_PHC_DEMODX_14Sep2026.csv",
    "apoe": "APOERES_14Sep2026.csv",
    "ptdemog": "PTDEMOG_27Aug2025.csv",
    "csf": "ADSP_PHC_CSF_14Sep2026.csv",
    "mmse": "MMSE_18Jan2026.csv",
    "cogn": "ADSP_PHC_COGN_14Sep2026.csv",
    "fs6": "ADNI1GO23_FS6_alltimepoints_04282021.xlsx",
    "fs7": "UCSFFSX7_27Feb2026.csv",
    "combat": "ADSP_PHC_T1_FS_18Oct2025.csv",
}

# Expected size of the core population (Supplementary Methods 2). Cell counts are not stored
# here; see --expected-counts.
EXPECTED = {"n": 1070}

# FS6 workbook columns; bilateral totals are halved to give the left-right mean.
FS6_COLS = ["SubjID", "RID", "BL_ScanType", "BL_ICV", "BL_LHippVol", "BL_RHippVol",
            "BL_entorhinal_volume", "BL_fusiform_volume", "BL_middletemporal_volume",
            "BL_precuneus_volume", "BL_BrainSegVolNotVent", "BL_LLatVent", "BL_RLatVent"]

# UCSF FS7 ST codes (sensitivity volumes): (left, right) or single.
FS7_MAP = {
    "Hippocampus": ("ST29SV", "ST88SV"), "Entorhinal": ("ST24CV", "ST83CV"),
    "Fusiform": ("ST26CV", "ST85CV"), "MiddleTemporal": ("ST40CV", "ST99CV"),
    "Precuneus": ("ST52CV", "ST111CV"), "Ventricles": ("ST37SV", "ST96SV"),
    "WholeBrain": ("ST128SV",), "ICV": ("ST10CV",),
}
COMBAT_MAP = {
    "Hippocampus": ("Left.Hippocampus_combat", "Right.Hippocampus_combat"),
    "Entorhinal": ("lh_entorhinal_volume_combat", "rh_entorhinal_volume_combat"),
    "Fusiform": ("lh_fusiform_volume_combat", "rh_fusiform_volume_combat"),
    "MiddleTemporal": ("lh_middletemporal_volume_combat", "rh_middletemporal_volume_combat"),
    "Precuneus": ("lh_precuneus_volume_combat", "rh_precuneus_volume_combat"),
    "Ventricles": ("Left.Lateral.Ventricle_combat", "Right.Lateral.Ventricle_combat"),
    "WholeBrain": ("BrainSegVolNotVentSurf_combat",), "ICV": ("EstimatedTotalIntraCranialVol_combat",),
}
VOLUME_EVENTS = ["Hippocampus", "Entorhinal", "Fusiform", "MiddleTemporal", "Precuneus",
                 "WholeBrain", "Ventricles"]
EVENTS_P12 = ["ABETA", "TAU", "PTAU", "MMSE", "MEM"] + VOLUME_EVENTS


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def genotype_group(g):
    if pd.isna(g):
        return "NA"
    alleles = set(str(g).split("/"))
    if alleles == {"2", "4"}:
        return "e2e4"
    if "2" in alleles:
        return "e2"
    if "4" in alleles:
        return "e4"
    return "e33"


# ----------------------------------------------------------------------------- diagnosis / cohort

def baseline_and_class(dxsum: pd.DataFrame) -> pd.DataFrame:
    """First DXSUM record with EXAMDATE and DIAGNOSIS is the baseline; MCIc = AD within 3 years."""
    d = dxsum.copy()
    d["EXAMDATE"] = pd.to_datetime(d["EXAMDATE"], errors="coerce")
    d = d.dropna(subset=["EXAMDATE", "DIAGNOSIS"]).sort_values(["RID", "EXAMDATE"])
    first = d.groupby("RID", sort=False).first().reset_index()
    first = first.rename(columns={"EXAMDATE": "BLDATE", "DIAGNOSIS": "BL_DIAGNOSIS",
                                  "PHASE": "BL_PHASE", "VISCODE2": "BL_VISCODE2"})
    later = d.merge(first[["RID", "BLDATE"]], on="RID")
    later = later[(later.EXAMDATE > later.BLDATE)
                  & (later.EXAMDATE <= later.BLDATE + pd.Timedelta(days=THREE_YEARS_DAYS))]
    conv = later.groupby("RID")["DIAGNOSIS"].apply(lambda s: bool((s == 3).any()))
    first["conv3"] = first.RID.map(conv).fillna(False).astype(bool)
    fu = d.groupby("RID").EXAMDATE.max().reset_index().rename(columns={"EXAMDATE": "LASTDATE"})
    first = first.merge(fu, on="RID")
    first["followup_days"] = (first.LASTDATE - first.BLDATE).dt.days
    cls = np.select(
        [first.BL_DIAGNOSIS == 1, (first.BL_DIAGNOSIS == 2) & first.conv3,
         (first.BL_DIAGNOSIS == 2) & ~first.conv3, first.BL_DIAGNOSIS == 3],
        ["CN", "MCIc", "MCInc", "AD"], "other")
    first["cls"] = cls
    keep = ["RID", "PTID", "BL_PHASE", "BL_VISCODE2", "BLDATE", "BL_DIAGNOSIS", "conv3",
            "followup_days", "cls"]
    return first[keep]


def phc_baseline_diagnosis(demodx: pd.DataFrame, base: pd.DataFrame) -> pd.Series:
    """Cross-check: PHC_Diagnosis at the DEMODX visit closest to BLDATE (≤180 d)."""
    d = demodx[["RID", "EXAMDATE", "PHC_Diagnosis"]].copy()
    d["EXAMDATE"] = pd.to_datetime(d["EXAMDATE"], errors="coerce")
    d = d.dropna(subset=["EXAMDATE", "PHC_Diagnosis"]).merge(base[["RID", "BLDATE"]], on="RID")
    d["delta"] = (d.EXAMDATE - d.BLDATE).dt.days.abs()
    d = d[d.delta <= 180].sort_values(["RID", "delta"]).drop_duplicates("RID")
    return base.RID.map(d.set_index("RID").PHC_Diagnosis)


def demographics(ptdemog: pd.DataFrame) -> pd.DataFrame:
    """First PTDEMOG record per participant with usable fields; PTDOB is MM/YYYY (day = 15)."""
    p = ptdemog.copy()
    p["VISDATE"] = pd.to_datetime(p["VISDATE"], errors="coerce")
    p = p.sort_values(["RID", "VISDATE"])
    p["PTEDUCAT"] = pd.to_numeric(p["PTEDUCAT"], errors="coerce")
    p.loc[p.PTEDUCAT < 0, "PTEDUCAT"] = np.nan
    p["PTGENDER"] = pd.to_numeric(p["PTGENDER"], errors="coerce")
    dob = pd.to_datetime(p["PTDOB"].astype(str).str.strip() + "/15", format="%m/%Y/%d", errors="coerce")
    p["DOB"] = dob
    p.loc[p.DOB.isna() & p.PTDOBYY.notna(), "DOB"] = pd.to_datetime(
        p.loc[p.DOB.isna() & p.PTDOBYY.notna(), "PTDOBYY"].astype(int).astype(str) + "-07-01", errors="coerce")
    out = {}
    for col in ["PTGENDER", "DOB", "PTEDUCAT"]:
        out[col] = p.dropna(subset=[col]).drop_duplicates("RID").set_index("RID")[col]
    res = pd.DataFrame(out).reset_index()
    res["Sex"] = res.PTGENDER.map({1: "Male", 2: "Female"})   # pyebm expects these strings
    return res[["RID", "Sex", "DOB", "PTEDUCAT"]].rename(columns={"PTEDUCAT": "Education"})


# ----------------------------------------------------------------------------- measurements

def closest_within(df: pd.DataFrame, datecol: str, base: pd.DataFrame, cols: list, prefix: str):
    """Closest record to BLDATE within WINDOW_DAYS_MAX with any of `cols` observed; keeps the offset."""
    d = df[["RID", datecol] + cols].copy()
    d[datecol] = pd.to_datetime(d[datecol], errors="coerce")
    d = d.dropna(subset=[datecol]).dropna(subset=cols, how="all")
    d = d.merge(base[["RID", "BLDATE"]], on="RID")
    d["delta"] = (d[datecol] - d.BLDATE).dt.days
    d = d[d.delta.abs() <= WINDOW_DAYS_MAX].assign(absd=lambda x: x.delta.abs())
    d = d.sort_values(["RID", "absd"]).drop_duplicates("RID")
    d = d[["RID", "delta"] + cols].rename(columns={"delta": f"{prefix}_offset_days"})
    return d


def load_fs6(path: Path, cache: Path) -> pd.DataFrame:
    """FS 6.0 workbook, sheet 'allscans'; one row per RID (scan type, then scan letter)."""
    if cache.exists():
        fs = pd.read_csv(cache)
    else:
        import openpyxl  # noqa: WPS433  (only needed once; cached afterwards)
        wb = openpyxl.load_workbook(path, read_only=True)
        ws = wb["allscans"]
        rows = ws.iter_rows(values_only=True)
        header = list(next(rows))
        idx = [header.index(c) for c in FS6_COLS]
        fs = pd.DataFrame([[r[i] for i in idx] for r in rows], columns=FS6_COLS)
        fs.to_csv(cache, index=False)
    fs = fs.dropna(subset=["BL_LHippVol"]).copy()
    fs["scan_letter"] = fs.SubjID.astype(str).str[-1]
    fs["type_rank"] = fs.BL_ScanType.fillna(9)
    fs = fs.sort_values(["RID", "type_rank", "scan_letter"]).drop_duplicates("RID")
    out = pd.DataFrame({
        "RID": fs.RID.astype(int),
        "FS6_ScanType": fs.BL_ScanType,
        "ICV": fs.BL_ICV,
        "Hippocampus": (fs.BL_LHippVol + fs.BL_RHippVol) / 2,
        "Entorhinal": fs.BL_entorhinal_volume / 2,
        "Fusiform": fs.BL_fusiform_volume / 2,
        "MiddleTemporal": fs.BL_middletemporal_volume / 2,
        "Precuneus": fs.BL_precuneus_volume / 2,
        "WholeBrain": fs.BL_BrainSegVolNotVent,
        "Ventricles": (fs.BL_LLatVent + fs.BL_RLatVent) / 2,
    })
    return out


def volumes_from_map(df: pd.DataFrame, datecol: str, base: pd.DataFrame, mapping: dict, suffix: str):
    cols = sorted({c for v in mapping.values() for c in v})
    d = df[["RID", datecol] + cols].copy()
    for c in cols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    key = mapping["Hippocampus"][0]
    d = closest_within(d.dropna(subset=[key]), datecol, base, cols, suffix)
    out = pd.DataFrame({"RID": d.RID, f"{suffix}_offset_days": d[f"{suffix}_offset_days"]})
    for name, src in mapping.items():
        out[f"{name}_{suffix}"] = d[list(src)].mean(axis=1) if len(src) == 2 else d[src[0]]
    return out


# ----------------------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=os.environ.get("CONCORD_ADNI_DIR"),
                    help="folder with the ADNI files (default: $CONCORD_ADNI_DIR)")
    ap.add_argument("--out-dir", type=Path, default=None, help="default: <data-dir>/derived")
    ap.add_argument("--expected-counts", type=Path, default=None,
                    help="optional local JSON with expected cell counts and availability")
    ap.add_argument("--no-hash", action="store_true", help="skip SHA256 of inputs (slow for the DTI-sized files)")
    args = ap.parse_args(argv)
    if args.data_dir is None:
        ap.error("--data-dir is required (or set CONCORD_ADNI_DIR)")
    data = args.data_dir
    out = args.out_dir or (data / "derived")
    out.mkdir(parents=True, exist_ok=True)
    report = []

    def say(*parts):
        line = " ".join(str(p) for p in parts)
        report.append(line)
        print(line)

    say(f"build_cohort  definition {PLAN_VERSION}")
    if not args.no_hash:
        for key, name in FILES.items():
            say(f"  input {name}: {sha256(data / name)}")

    dxsum = pd.read_csv(data / FILES["dxsum"], low_memory=False)
    base = baseline_and_class(dxsum)
    apoe = pd.read_csv(data / FILES["apoe"], low_memory=False).drop_duplicates("RID")[["RID", "GENOTYPE"]]
    base = base.merge(apoe, on="RID", how="left")
    base["G"] = base.GENOTYPE.map(genotype_group)
    base["PHC_Diagnosis_bl"] = phc_baseline_diagnosis(pd.read_csv(data / FILES["demodx"], low_memory=False), base)

    old = base.BL_PHASE.isin(["ADNI1", "ADNIGO", "ADNI2"])
    grp_ok = base.G.isin(["e2", "e33", "e4"])
    core_mask = old & grp_ok & base.cls.isin(["CN", "MCIc", "AD"])
    s1_mask = old & grp_ok & base.cls.isin(["MCInc"])          # S1 adds MCI non-converters
    cohort = base[core_mask | s1_mask].copy()
    cohort["cohort"] = np.where(core_mask[core_mask | s1_mask], "core", "S1_only")
    cohort["Diagnosis"] = cohort.cls.replace({"MCInc": "MCI", "MCIc": "MCIc"})

    # covariates
    demo = demographics(pd.read_csv(data / FILES["ptdemog"], low_memory=False))
    cohort = cohort.merge(demo, on="RID", how="left")
    cohort["Age"] = (cohort.BLDATE - cohort.DOB).dt.days / 365.25

    # events
    csf = closest_within(pd.read_csv(data / FILES["csf"], low_memory=False), "DRAWDATE", cohort,
                         ["PHC_AB42", "PHC_Tau", "PHC_pTau181"], "CSF")
    csf = csf.rename(columns={"PHC_AB42": "ABETA", "PHC_Tau": "TAU", "PHC_pTau181": "PTAU"})
    mmse = pd.read_csv(data / FILES["mmse"], low_memory=False)
    mmse["MMSCORE"] = pd.to_numeric(mmse["MMSCORE"], errors="coerce")
    mmse = closest_within(mmse, "VISDATE", cohort, ["MMSCORE"], "MMSE").rename(columns={"MMSCORE": "MMSE"})
    cogn = closest_within(pd.read_csv(data / FILES["cogn"], low_memory=False), "EXAMDATE", cohort,
                          ["PHC_MEM"], "MEM").rename(columns={"PHC_MEM": "MEM"})
    fs6 = load_fs6(data / FILES["fs6"], out / "fs6_baseline_extract.csv")
    fs7 = volumes_from_map(pd.read_csv(data / FILES["fs7"], low_memory=False), "EXAMDATE", cohort, FS7_MAP, "FS7")
    combat = volumes_from_map(pd.read_csv(data / FILES["combat"], low_memory=False), "PHC_SCANDATE", cohort,
                              COMBAT_MAP, "ComBat")
    for part in (csf, mmse, cogn, fs6, fs7, combat):
        cohort = cohort.merge(part, on="RID", how="left")

    # ---------------------------------------------------------------- report (counts only)
    core = cohort[cohort.cohort == "core"]
    say(f"\ncore cohort n = {len(core)}   (expected {EXPECTED['n']})")
    tab = pd.crosstab(core.G, core.cls).reindex(index=["e2", "e33", "e4"], columns=["CN", "MCIc", "AD"])
    say(tab.to_string())
    ok = len(core) == EXPECTED["n"]
    detail = json.loads(args.expected_counts.read_text()) if args.expected_counts else {}
    for g, cells in detail.get("table", {}).items():
        for c, n in cells.items():
            ok &= int(tab.loc[g, c]) == int(n)
    say(f"class totals {core.cls.value_counts().to_dict()}  group totals {core.G.value_counts().to_dict()}")
    say(f"baseline phase {core.BL_PHASE.value_counts().to_dict()}  visit code {core.BL_VISCODE2.value_counts().to_dict()}")

    def n_within(col_offset, value_col, days):
        m = core[col_offset].abs().le(days) & core[value_col].notna()
        return int(m.sum())
    avail = {"MMSE": n_within("MMSE_offset_days", "MMSE", 180),
             "CSF": n_within("CSF_offset_days", "ABETA", 180),
             "MEM": n_within("MEM_offset_days", "MEM", 180),
             "FS6": int(core.Hippocampus.notna().sum()),
             "FS7": n_within("FS7_offset_days", "Hippocampus_FS7", 180),
             "ComBat": n_within("ComBat_offset_days", "Hippocampus_ComBat", 180)}
    say(f"availability at baseline (±180 d; FS6 = BL_ columns): {avail}")
    for k, v in detail.get("avail_180", {}).items():
        ok &= avail[k] == int(v)
    say("per-group availability (±180 d):")
    for g in ["e2", "e33", "e4"]:
        sub = core[core.G == g]
        say(f"  {g:3s} n={len(sub):3d} CSF={int((sub.CSF_offset_days.abs().le(180) & sub.ABETA.notna()).sum()):3d}"
            f" MMSE={int((sub.MMSE_offset_days.abs().le(180) & sub.MMSE.notna()).sum()):3d}"
            f" MEM={int((sub.MEM_offset_days.abs().le(180) & sub.MEM.notna()).sum()):3d}"
            f" FS6={int(sub.Hippocampus.notna().sum()):3d}")
    say(f"covariates missing: Age {int(core.Age.isna().sum())}, Sex {int(core.Sex.isna().sum())}, "
        f"Education {int(core.Education.isna().sum())}, ICV(FS6) {int(core.ICV.isna().sum())}")
    say(f"Age (years) range {core.Age.min():.1f}–{core.Age.max():.1f}; Sex {core.Sex.value_counts().to_dict()}")
    dis = core.PHC_Diagnosis_bl.notna() & (core.PHC_Diagnosis_bl != core.BL_DIAGNOSIS)
    say(f"PHC_Diagnosis cross-check: matched {int(core.PHC_Diagnosis_bl.notna().sum())}, disagree {int(dis.sum())}")
    say(f"S1 extension (MCI non-converters): {int((cohort.cohort == 'S1_only').sum())}")
    say(f"\nEXPECTED COUNTS {'MATCH' if ok else 'DO NOT MATCH'}"
        + ("" if detail else " (population size only; no --expected-counts given)"))

    # ---------------------------------------------------------------- outputs
    cohort = cohort.sort_values("RID").reset_index(drop=True)
    cohort.to_csv(out / "adni_cohort_full.csv", index=False)
    deid_cols = (["cohort", "Diagnosis", "G", "Age", "Sex", "Education", "ICV", "FS6_ScanType", "BL_PHASE"]
                 + EVENTS_P12 + ["CSF_offset_days", "MMSE_offset_days", "MEM_offset_days",
                                 "FS7_offset_days", "ComBat_offset_days", "PHC_Diagnosis_bl"]
                 + [f"{e}_FS7" for e in VOLUME_EVENTS] + ["ICV_FS7"]
                 + [f"{e}_ComBat" for e in VOLUME_EVENTS] + ["ICV_ComBat"])
    deid = cohort[deid_cols].copy()
    deid.insert(0, "PTID", [f"a{i:05d}" for i in range(len(deid))])   # sequential; mapping stays in the full table
    deid = deid.rename(columns={"G": "APOE"})
    deid["Age"] = deid.Age.round(2)
    deid.to_csv(out / "adni_p12_deid.csv", index=False)
    (out / "build_cohort_report.txt").write_text("\n".join(report) + "\n")
    say(f"\nwrote {out / 'adni_cohort_full.csv'} (identifiers; local only)")
    say(f"wrote {out / 'adni_p12_deid.csv'} (de-identified; {len(deid)} rows, {deid.shape[1]} columns)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
