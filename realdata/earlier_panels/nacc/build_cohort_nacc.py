#!/usr/bin/env python
"""Build the NACC amyloid PET population (1,200 participants; Supplementary Methods 2). Fits nothing.

Earlier analysis population behind the NACC amyloid PET row of Fig. 4e and the NACC panel of
Supplementary Fig. 3. Each participant's first ADSP-PHC amyloid PET scan with a Centiloid value is
linked to the UDS visit closest in age to the ADSP-PHC cognition age (at most 1 year apart);
diagnoses CN, prevalent MCI and AD dementia with AD aetiology; APOE e2 (codes 3, 6), e3/e3 (1) and
e4 (2, 4), e2/e4 excluded. Events: Centiloid, MoCA and the ADSP-PHC memory, executive function
and language composites; incomplete cognitive profiles are allowed. The definition was fixed on
17 September 2026.

Inputs (--data-dir is the NACC ADSP-PHC release folder, default $CONCORD_NACC_DIR):
  --pet-file   default <data-dir>/Imaging_PET/NACC_ADSP_PHC_Amyloid_Simple_2024.csv
  --cog-file   default <data-dir>/Cognition/NACC_ADSP_PHC_Cognition_2024.csv
  --uds-file   CSV export of the NACC UDS file (default $CONCORD_NACC_UDS)
Writes to --out-dir (required; keep it outside the repository):
<out>/nacc_cohort_full.csv (identifiers; local only), <out>/nacc_p5_deid.csv (cluster input) and
<out>/build_cohort_report.txt (counts only).
"""
from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

EVENTS = ["CENTILOID", "MOCA", "MEM", "EXF", "LAN"]
APOE_GROUP = {1: "e33", 2: "e4", 3: "e2", 4: "e4", 5: "e2e4", 6: "e2"}


def sha256(p: Path):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=os.environ.get("CONCORD_NACC_DIR"),
                    help="NACC ADSP-PHC release folder (default: $CONCORD_NACC_DIR)")
    ap.add_argument("--pet-file", type=Path, help="default: <data-dir>/Imaging_PET/NACC_ADSP_PHC_Amyloid_Simple_2024.csv")
    ap.add_argument("--cog-file", type=Path, help="default: <data-dir>/Cognition/NACC_ADSP_PHC_Cognition_2024.csv")
    ap.add_argument("--uds-file", type=Path, default=os.environ.get("CONCORD_NACC_UDS"),
                    help="CSV export of the NACC UDS file (default: $CONCORD_NACC_UDS)")
    ap.add_argument("--out-dir", type=Path, required=True, help="output folder outside the repository")
    args = ap.parse_args(argv)
    d = args.data_dir
    if args.uds_file is None or (d is None and None in (args.pet_file, args.cog_file)):
        ap.error("give --uds-file and either --data-dir or both --pet-file and --cog-file")
    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    rep = []

    def say(*p):
        line = " ".join(str(x) for x in p)
        rep.append(line)
        print(line)

    files = {"pet": args.pet_file or d / "Imaging_PET/NACC_ADSP_PHC_Amyloid_Simple_2024.csv",
             "cog": args.cog_file or d / "Cognition/NACC_ADSP_PHC_Cognition_2024.csv",
             "uds": args.uds_file}
    for k, f in files.items():
        say(f"input {Path(f).name}: {sha256(Path(f))}")
    pet = pd.read_csv(files["pet"], low_memory=False)
    pet = pet[pet.PHC_CENTILOIDS.notna()].sort_values(["NACCID", "PHC_Age_PET"]).drop_duplicates("NACCID")
    uds = pd.read_csv(files["uds"], low_memory=False)
    cog = pd.read_csv(files["cog"], low_memory=False)

    m = pet[["NACCID", "PHC_Age_Cognition", "PHC_Age_PET", "PHC_CENTILOIDS", "PHC_TRACER"]].merge(uds, on="NACCID", how="inner")
    m["age_diff"] = (m.NACCAGE - m.PHC_Age_Cognition).abs()
    m = m.sort_values(["NACCID", "age_diff"]).drop_duplicates("NACCID")
    m = m[m.age_diff <= 1.0]
    m["APOE"] = m.NACCAPOE.map(APOE_GROUP).fillna("NA")
    cls = np.select([m.NACCUDSD == 1, m.NACCUDSD == 3, (m.NACCUDSD == 4) & (m.NACCALZD == 1)],
                    ["CN", "MCI", "AD"], "excl")
    m["Diagnosis"] = cls
    say(f"PET subjects with a matched visit (|age diff| <= 1 y): {len(m)}")
    say("excluded: impaired-not-MCI", int((m.NACCUDSD == 2).sum()), "| dementia not coded AD",
        int(((m.NACCUDSD == 4) & (m.NACCALZD != 1)).sum()), "| APOE 2/4 or unknown", int((~m.APOE.isin(["e2", "e33", "e4"])).sum()))
    core = m[m.APOE.isin(["e2", "e33", "e4"]) & m.Diagnosis.isin(["CN", "MCI", "AD"])].copy()
    core = core.merge(cog, on=["NACCID", "NACCVNUM"], how="left")
    core["CENTILOID"] = core.PHC_CENTILOIDS
    core["MOCA"] = pd.to_numeric(core.NACCMOCA, errors="coerce").where(lambda s: s.between(0, 30))
    core["MEM"], core["EXF"], core["LAN"] = core.PHC_MEM, core.PHC_EXF, core.PHC_LAN
    core["Age"] = core.NACCAGE
    core["Sex"] = core.SEX.map({1: "Male", 2: "Female"})
    core["Education"] = pd.to_numeric(core.EDUC, errors="coerce").where(lambda s: s < 99)
    core = core[core[EVENTS].notna().any(axis=1)].copy()

    tab = pd.crosstab(core.APOE, core.Diagnosis).reindex(index=["e2", "e33", "e4"], columns=["CN", "MCI", "AD"]).fillna(0).astype(int)
    say(f"\ncore cohort n = {len(core)}")
    say(tab.to_string())
    if (tab == 0).any().any():
        say("EMPTY group x class cell: stop")
        return 1
    say("availability:", {e: int(core[e].notna().sum()) for e in EVENTS})
    say("covariates missing: Age", int(core.Age.isna().sum()), "Sex", int(core.Sex.isna().sum()), "Education", int(core.Education.isna().sum()))
    say("tracers:", core.PHC_TRACER.value_counts().to_dict())
    say("composition (CN/MCI/AD):", {g: tuple((tab.loc[g] / tab.loc[g].sum()).round(2)) for g in tab.index})
    # draw sizes for the Fig. 4e comparison (recomputed by run_nacc.n2_sizes)
    e2 = tab.loc["e2"].to_dict()
    e4p = (tab.loc["e4"] / tab.loc["e4"].sum())
    remaining = tab.loc["e33"] - pd.Series(e2)
    n_max = int(min(remaining[c] / e4p[c] for c in ["CN", "MCI", "AD"]))
    n4 = (n_max // 10) * 10
    draw4 = {c: int(round(n4 * e4p[c])) for c in ["CN", "MCI", "AD"]}
    say(f"draw sizes for Fig. 4e: group A = {e2} (n={sum(e2.values())}); group B = {draw4} (n={sum(draw4.values())}; max feasible {n_max})")

    core = core.sort_values("NACCID").reset_index(drop=True)
    core.to_csv(out / "nacc_cohort_full.csv", index=False)
    deid = core[["Diagnosis", "APOE", "Age", "Sex", "Education", "age_diff"] + EVENTS].copy()
    deid.insert(0, "PTID", [f"n{i:05d}" for i in range(len(deid))])
    deid["Age"] = deid.Age.round(2)
    deid.to_csv(out / "nacc_p5_deid.csv", index=False)
    (out / "build_cohort_report.txt").write_text("\n".join(rep) + "\n")
    say(f"\nwrote {out / 'nacc_p5_deid.csv'} ({len(deid)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
