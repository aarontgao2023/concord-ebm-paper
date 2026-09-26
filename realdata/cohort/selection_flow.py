#!/usr/bin/env python3
"""Six-event panel: selection-flow counts of Supplementary Table 3 (no model fits, no tests).

Reads the two aggregate audit files written by the cohort audits and writes
S5_selection_flow.csv (columns cohort, stage, participants; 8 rows), the file that
figures/tables/make_si_cohort_tables.py reads from figures/inputs/realdata/.

  python realdata/cohort/selection_flow.py \\
      --adni-audit $CONCORD_WORK_DIR/audit/adni/coverage_audit.json \\
      --nacc-audit $CONCORD_WORK_DIR/audit/nacc/coverage_summary.json \\
      --out $CONCORD_WORK_DIR/audit/S5_selection_flow.csv

Both inputs contain aggregate counts only. The rows and the keys they come from:

  stage                                ADNI (coverage_audit.json,          NACC (coverage_summary.json,
                                       adni_audit.py)                      nacc_audit.py)
  Supplied PHC CSF participants        source_CSF_people                   source.csf_participants
  Eligible CSF anchor diagnosis APOE   anchor_people                       CSF.anchor_n
  Valid age sex education              panels.CSFcog6.covariate_complete_n CSF.core_covariate_n
  Complete six events and covariates   panels.CSFcog6.                     CSF.CSFcog6.complete_n
                                       complete_events_and_covariates_n

Supplied: participants with at least one row in the ADSP-PHC CSF table. Eligible: participants
with an anchor visit (CSF measurement with an eligible APOE genotype and a CN, MCI or AD
diagnosis). Valid: anchors with age, recorded sex and years of education. Complete: anchors with
these covariates and all six events (ABETA, TAU, PTAU, MEM, EXF, LAN); these are the 1,203 ADNI
and 1,584 NACC participants of the analysis.
"""
import argparse
import csv
import json
import os
from pathlib import Path

STAGES = [
    ("Supplied PHC CSF participants",
     ("source_CSF_people",),
     ("source", "csf_participants")),
    ("Eligible CSF anchor diagnosis APOE",
     ("anchor_people",),
     ("CSF", "anchor_n")),
    ("Valid age sex education",
     ("panels", "CSFcog6", "covariate_complete_n"),
     ("CSF", "core_covariate_n")),
    ("Complete six events and covariates",
     ("panels", "CSFcog6", "complete_events_and_covariates_n"),
     ("CSF", "CSFcog6", "complete_n")),
]


def lookup(data, keys, source):
    value = data
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise SystemExit(f"{source}: key {'.'.join(keys)} not found")
        value = value[key]
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise SystemExit(f"{source}: {'.'.join(keys)} is not a count: {value!r}")
    return value


def selection_rows(adni, nacc, adni_name="coverage_audit.json", nacc_name="coverage_summary.json"):
    rows = []
    for cohort, data, name, column in (("ADNI", adni, adni_name, 1), ("NACC", nacc, nacc_name, 2)):
        counts = [lookup(data, stage[column], name) for stage in STAGES]
        if any(later > earlier for earlier, later in zip(counts, counts[1:])):
            raise SystemExit(f"{cohort}: counts increase along the selection steps: {counts}")
        rows.extend((cohort, stage[0], n) for stage, n in zip(STAGES, counts))
    # The ADNI anchor count is recorded twice; both must agree.
    if adni["panels"]["CSFcog6"].get("anchor_n", adni["anchor_people"]) != adni["anchor_people"]:
        raise SystemExit(f"{adni_name}: panels.CSFcog6.anchor_n differs from anchor_people")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    work = os.environ.get("CONCORD_WORK_DIR")
    default = (lambda *p: Path(work, *p)) if work else (lambda *p: None)
    ap.add_argument("--adni-audit", type=Path, default=default("audit", "adni", "coverage_audit.json"),
                    help="coverage_audit.json of adni_audit.py (default: $CONCORD_WORK_DIR/audit/adni/)")
    ap.add_argument("--nacc-audit", type=Path, default=default("audit", "nacc", "coverage_summary.json"),
                    help="coverage_summary.json of nacc_audit.py (default: $CONCORD_WORK_DIR/audit/nacc/)")
    ap.add_argument("--out", type=Path, default=default("audit", "S5_selection_flow.csv"),
                    help="output CSV (default: $CONCORD_WORK_DIR/audit/S5_selection_flow.csv)")
    args = ap.parse_args(argv)
    if args.adni_audit is None or args.nacc_audit is None or args.out is None:
        ap.error("--adni-audit, --nacc-audit and --out are required (or set CONCORD_WORK_DIR)")
    adni = json.loads(args.adni_audit.read_text())
    nacc = json.loads(args.nacc_audit.read_text())
    rows = selection_rows(adni, nacc, args.adni_audit.name, args.nacc_audit.name)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["cohort", "stage", "participants"])
        writer.writerows(rows)
    for cohort, stage, n in rows:
        print(f"{cohort:5s} {stage:36s} {n:6,d}")
    print(args.out)


if __name__ == "__main__":
    main()
