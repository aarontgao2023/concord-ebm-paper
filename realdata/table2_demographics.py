#!/usr/bin/env python3
"""Table 2: participants of the six-event panel by cohort and APOE group (aggregate output only).

Reads the prepared, unadjusted analysis inputs written by prepare_campaign.py
(runtime/inputs/<cohort>/raw_deid.csv; columns PTID, APOE, Diagnosis, Age, Sex, Education and the
six events). These files hold participant-level data and are available only under the ADNI and
NACC data use agreements. This script prints and writes group-level aggregates only.

For each cohort and APOE group (e2, e3/e3, e4 and all participants) it computes the number of
participants, the CN/MCI/AD counts, age and years of education (mean and SD, n - 1 denominator)
and the number and percentage of women, and, per cohort, the Pearson chi-square test of the
3 x 3 genotype-by-diagnosis table (4 d.f.). Means, SDs and percentages are rounded to one decimal
as in Table 2.

  python table2_demographics.py --inputs-dir realdata/runtime/inputs --out $CONCORD_WORK_DIR/table2

or give the two files with --adni and --nacc. --out must lie outside the repository.
Outputs: Table2_demographics_aggregate.csv and Table2_chi_square.csv.
"""
import argparse
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
GROUPS = ['e2', 'e33', 'e4']
DX = ['CN', 'MCI', 'AD']
COHORTS = {'adni': 'ADNI', 'nacc': 'NACC'}
COLUMNS = ['cohort', 'genotype', 'n', 'CN', 'MCI', 'AD', 'age_mean', 'age_sd', 'female_n', 'female_pct',
           'education_mean', 'education_sd']


def one_decimal(x):
    return f'{x:.1f}'


def group_row(cohort, label, frame):
    n = len(frame)
    dx = frame.Diagnosis.value_counts().reindex(DX, fill_value=0)
    female = int(frame.Sex.eq('Female').sum())
    return {'cohort': cohort, 'genotype': label, 'n': n, 'CN': int(dx['CN']), 'MCI': int(dx['MCI']),
            'AD': int(dx['AD']), 'age_mean': one_decimal(frame.Age.mean()), 'age_sd': one_decimal(frame.Age.std(ddof=1)),
            'female_n': female, 'female_pct': one_decimal(100 * female / n),
            'education_mean': one_decimal(frame.Education.mean()),
            'education_sd': one_decimal(frame.Education.std(ddof=1))}


def chi_square(frame):
    """Pearson chi-square of the genotype x diagnosis table; survival function for 4 d.f."""
    table = pd.crosstab(frame.APOE, frame.Diagnosis).reindex(index=GROUPS, columns=DX, fill_value=0).to_numpy(float)
    expected = table.sum(axis=1, keepdims=True) * table.sum(axis=0, keepdims=True) / table.sum()
    statistic = float(((table - expected) ** 2 / expected).sum())
    df = (table.shape[0] - 1) * (table.shape[1] - 1)
    assert df == 4
    return statistic, df, math.exp(-statistic / 2) * (1 + statistic / 2)


def read_input(path):
    frame = pd.read_csv(path, float_precision='round_trip')
    required = {'APOE', 'Diagnosis', 'Age', 'Sex', 'Education'}
    missing = required - set(frame.columns)
    if missing:
        raise SystemExit(f'{path.name}: missing columns {sorted(missing)}')
    if not set(frame.APOE) <= set(GROUPS) or not set(frame.Diagnosis) <= set(DX):
        raise SystemExit(f'{path.name}: unexpected APOE or diagnosis labels')
    if not set(frame.Sex) <= {'Male', 'Female'}:
        raise SystemExit(f'{path.name}: Sex must be Male or Female')
    if not np.isfinite(frame[['Age', 'Education']].to_numpy(float)).all():
        raise SystemExit(f'{path.name}: age and education must be complete')
    return frame


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--inputs-dir', type=Path, help='folder with adni/raw_deid.csv and nacc/raw_deid.csv')
    ap.add_argument('--adni', type=Path, help='ADNI raw_deid.csv (overrides --inputs-dir)')
    ap.add_argument('--nacc', type=Path, help='NACC raw_deid.csv (overrides --inputs-dir)')
    ap.add_argument('--out', type=Path, required=True, help='output folder outside the repository')
    args = ap.parse_args(argv)
    paths = {'adni': args.adni, 'nacc': args.nacc}
    for cohort in paths:
        if paths[cohort] is None:
            if args.inputs_dir is None:
                ap.error('give --inputs-dir or both --adni and --nacc')
            paths[cohort] = args.inputs_dir / cohort / 'raw_deid.csv'
    out = args.out.resolve()
    if out == REPO or REPO in out.parents:
        raise SystemExit(f'The output folder must lie outside the repository: {out}')
    out.mkdir(parents=True, exist_ok=True)

    rows, tests = [], []
    for cohort, name in COHORTS.items():
        frame = read_input(paths[cohort])
        for g in GROUPS:
            rows.append(group_row(name, g, frame[frame.APOE.eq(g)]))
        rows.append(group_row(name, 'all', frame))
        statistic, df, p = chi_square(frame)
        tests.append({'cohort': name, 'chi2': round(statistic, 4), 'df': df, 'p': f'{p:.3g}'})

    with open(out / 'Table2_demographics_aggregate.csv', 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    with open(out / 'Table2_chi_square.csv', 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['cohort', 'chi2', 'df', 'p'])
        writer.writeheader()
        writer.writerows(tests)
    for row in rows:
        print(','.join(str(row[c]) for c in COLUMNS))
    for test in tests:
        print(f"{test['cohort']}: chi-square {test['chi2']:.1f}, {test['df']} d.f., P = {test['p']}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
