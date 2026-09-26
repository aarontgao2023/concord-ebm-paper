#!/usr/bin/env python3
"""Earlier analysis populations: aggregate the stored results for Fig. 4e and Supplementary Fig. 3.

Standard library only. No model is fitted. Reads the per-draw JSON files of adni/run_null.py and
nacc/run_nacc.py and the bootstrap file of adni/run_observed.py (--stage a2s), and writes:

  fig3_constructed_null_summary.csv  rejections out of 200 draws, ADNI 12-event and NACC amyloid
                                     PET settings (Fig. 4e); Wilson 95% intervals with
                                     z = NormalDist().inv_cdf(0.975)
  precision_adni_rates.csv           rejections out of 200 draws, ADNI 12-, 8- and 5-event settings
                                     (Fig. 4e uses the 8- and 5-event rows); Wilson 95% intervals
                                     with z = 1.96; mean and SD of the between-group distance
  precision_adni_provenance.json     checks that the three ADNI settings used the same draws
  fig4_ADNI_bootstrap_orderings.csv  the 500 bootstrap orderings per genotype of the ADNI 12-event
                                     panel (input of ../apoe_refit_summaries.py; Supplementary Fig. 3)

The four MRI-cognition rows of Fig. 4e come from mri_cognition/aggregate_completed_results.py
(n2_rejection_rates.csv). In every file, standard_U is the separately fitted DEBM (unmodified
pyebm search) with unrestricted permutation and concord_D is CONCORD with within-diagnosis
permutation; a draw counts as a rejection when E <= 29 of B = 599 permutations (reject_any_pair).

Expected layout (as written by the sbatch files in adni/hpc and nacc/hpc):
  <adni-results>/null/n2/rep_0000.json ... rep_0199.json          12 events
  <adni-results>/null_P8/n2/rep_*.json, <adni-results>/null_P5/n2/rep_*.json
  <adni-results>/main/a2s_core_P12_w180_DXSUM.json               500 bootstrap refits
  <nacc-results>/n2/rep_0000.json ... rep_0199.json

  python aggregate_earlier_panels.py --adni-results $CONCORD_WORK_DIR/earlier/adni \\
      --nacc-results $CONCORD_WORK_DIR/earlier/nacc --out $CONCORD_WORK_DIR/earlier

This replaces the corresponding parts of the scripts used for the paper (an R figure script and a
precision-table script); the numbers are the same, the CSV formatting of floats may differ.
"""
import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from statistics import NormalDist

ARMS = ('standard_U', 'concord_D')
LABELS = {'standard_U': ('separately fitted DEBM', 'unrestricted'),
          'concord_D': ('CONCORD', 'within-diagnosis')}
R = 200


def wilson(k, n, z):
    p = k / n
    scale = 1 + z * z / n
    center = (p + z * z / (2 * n)) / scale
    radius = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / scale
    return center - radius, center + radius


def scheme(arm):
    return 'unrestricted' if arm.endswith('_U') else 'diagnosis'


def read_draws(folder):
    files = sorted(Path(folder).glob('rep_*.json'))
    if len(files) != R:
        raise SystemExit(f'{folder}: expected {R} draw files, found {len(files)}')
    records = [json.loads(f.read_text()) for f in files]
    if sorted(d['replicate'] for d in records) != list(range(R)):
        raise SystemExit(f'{folder}: replicate indices are not 0..{R - 1}')
    return records


def rejections(records, arm):
    for d in records:
        if arm not in d:
            raise SystemExit(f'Missing result {arm} in replicate {d["replicate"]}')
    return sum(bool(d[arm]['tests'][scheme(arm)]['reject_any_pair']) for d in records)


def write_csv(path, rows, lineterminator='\r\n'):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator=lineterminator)
        writer.writeheader()
        writer.writerows(rows)


def constructed_null_summary(adni_results, nacc_results):
    """Port of the N2 part of the constructed-null summary (12-event ADNI, NACC amyloid PET)."""
    z = NormalDist().inv_cdf(0.975)
    rows = []
    for cohort, folder in (('ADNI', adni_results / 'null/n2'), ('NACC', nacc_results / 'n2')):
        records = read_draws(folder)
        for arm in sorted(ARMS):
            k = rejections(records, arm)
            low, high = wilson(k, len(records), z)
            method, permutation = LABELS[arm]
            rows.append({'cohort': cohort, 'null': 'N2', 'arm': arm, 'reject': k, 'R': len(records),
                         'low': low, 'high': high, 'rate': k / len(records), 'method': method,
                         'permutation': permutation})
    return rows


def adni_panel_rates(adni_results):
    """Port of the ADNI reduced-panel rate table (12, 8 and 5 events; same draws)."""
    rows, matched, sources = [], {}, []
    for panel in ('P5', 'P8', 'P12'):
        base = adni_results / ('null' if panel == 'P12' else 'null_' + panel)
        records = read_draws(base / 'n2')
        for d in records:
            assert d['B'] == 599 and d['n'] == 286, (panel, d['replicate'])
            identity = (d['n'], d['composition']); r = d['replicate']
            if r in matched:
                assert matched[r] == identity, (panel, r)
            else:
                matched[r] = identity
        for arm in ARMS:
            assert all(d[arm]['status'] == 'ok' and len(d[arm]['tests'][scheme(arm)]['pairs']) == 1 for d in records)
            k = rejections(records, arm)
            distances = [d[arm]['distances']['max'] for d in records]
            lo, hi = wilson(k, R, 1.96)
            rows.append(dict(panel=panel, stage='n2', arm=arm, K=int(panel[1:]), source_frame_n=1070, pseudo_draw_n=286,
                             planned=R, reject=k, rate=k / R, wilson95_low=lo, wilson95_high=hi,
                             distance_mean=statistics.mean(distances), distance_sd=statistics.stdev(distances)))
        sources.append(dict(panel=panel, stage='n2', folder=str(Path(base.name) / 'n2'), replicates=R, draw_n=286))
    provenance = {'sources': sources,
                  'draw_identity_check': 'All 200 replicate-indexed compositions and n agree across the three panels; '
                                         'the participant frame is fixed before the events are subset.',
                  'composition': matched[0][1]}
    return rows, provenance


def adni_bootstrap_orderings(adni_results):
    """Port of the extraction of the 500 ADNI 12-event bootstrap orderings (1-based positions)."""
    boot = json.loads((adni_results / 'main/a2s_core_P12_w180_DXSUM.json').read_text())
    groups = boot['groups']
    rows = []
    for armname, a in boot['arms'].items():
        assert a['completed'] == 500 and a['failed'] == 0, armname
        for r, ordering in enumerate(a['orderings'], start=1):
            for g, group in enumerate(groups):
                for position, event in enumerate(ordering[g], start=1):
                    rows.append({'arm': armname, 'replicate': r, 'group': group, 'position': position, 'event': event})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    work = Path(os.environ['CONCORD_WORK_DIR']) if os.environ.get('CONCORD_WORK_DIR') else None
    ap.add_argument('--adni-results', type=Path, default=work / 'earlier/adni' if work else None,
                    help='folder with null/, null_P8/, null_P5/ and main/ (default: $CONCORD_WORK_DIR/earlier/adni)')
    ap.add_argument('--nacc-results', type=Path, default=work / 'earlier/nacc' if work else None,
                    help='folder with n2/ (default: $CONCORD_WORK_DIR/earlier/nacc)')
    ap.add_argument('--out', type=Path, default=work / 'earlier' if work else None,
                    help='output folder (default: $CONCORD_WORK_DIR/earlier)')
    args = ap.parse_args(argv)
    if None in (args.adni_results, args.nacc_results, args.out):
        ap.error('--adni-results, --nacc-results and --out are required (or set CONCORD_WORK_DIR)')
    args.out.mkdir(parents=True, exist_ok=True)

    summary = constructed_null_summary(args.adni_results, args.nacc_results)
    write_csv(args.out / 'fig3_constructed_null_summary.csv', summary, lineterminator='\n')
    rates, provenance = adni_panel_rates(args.adni_results)
    write_csv(args.out / 'precision_adni_rates.csv', rates)
    (args.out / 'precision_adni_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    orderings = adni_bootstrap_orderings(args.adni_results)
    write_csv(args.out / 'fig4_ADNI_bootstrap_orderings.csv', orderings, lineterminator='\n')
    for row in summary:
        print(f"{row['cohort']} {row['arm']}: {row['reject']}/{row['R']} ({100 * row['rate']:.1f}%)")
    for row in rates:
        print(f"ADNI {row['K']} events {row['arm']}: {row['reject']}/{row['planned']} ({100 * row['rate']:.1f}%)")
    print(f"wrote {len(orderings)} bootstrap ordering rows")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
