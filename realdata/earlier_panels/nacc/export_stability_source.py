"""Collect the NACC bootstrap refits of run_nacc_stability.py into aggregate tables.

Checks that both runs (common proportions min and pooled) completed all 500 refits without
failure, then writes to --out:
  r2_nacc_position_probabilities.csv, r2_nacc_bootstrap_orderings.csv, r2_nacc_pair_precedence.csv
  (both runs stacked) and r2_nacc_stability_summary.csv (Centiloid first/last and rank summaries).
r2_nacc_bootstrap_orderings.csv is the input of ../../apoe_refit_summaries.py (Supplementary Fig. 3).
No participant-level records are read or written. A completion report goes to
<results-dir>/completion.json.

  python export_stability_source.py --results-dir $CONCORD_WORK_DIR/earlier/nacc_stability \\
      --out $CONCORD_WORK_DIR/earlier
"""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--results-dir', type=Path, required=True,
                    help='folder with the min/ and pooled/ outputs of run_nacc_stability.py')
    ap.add_argument('--out', type=Path, required=True, help='output folder for the aggregate tables')
    args = ap.parse_args(argv)
    RES, DATA = args.results_dir, args.out
    DATA.mkdir(parents=True, exist_ok=True)
    rows = []
    for reference in ('min', 'pooled'):
        d = RES / reference
        validation = json.loads((d / 'validation.json').read_text())
        assert validation['passed'] and validation['complete'] == 500 and validation['failed'] == 0
        pos = pd.read_csv(d / 'position_probabilities.csv'); assert len(pos) == 75
        assert (pos.groupby(['group', 'event'])['count'].sum() == 500).all()
        summary = json.loads((d / 'summary.json').read_text())
        for g in ('e2', 'e33', 'e4'):
            rows.append({'reference': reference, 'group': g, 'R': 500, **summary['centiloid'][g],
                         **{k: v for k, v in summary['groups'][g].items() if k != 'position_sd'}})
    for filename in ('position_probabilities.csv', 'bootstrap_orderings.csv', 'pair_precedence.csv'):
        pd.concat([pd.read_csv(RES / r / filename) for r in ('min', 'pooled')],
                  ignore_index=True).to_csv(DATA / ('r2_nacc_' + filename), index=False)
    pd.DataFrame(rows).to_csv(DATA / 'r2_nacc_stability_summary.csv', index=False)
    report = {'references': 2, 'planned_resamples_per_reference': 500, 'completed_per_reference': 500, 'failed': 0,
              'source_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('r2_nacc_*.csv')}}
    (RES / 'completion.json').write_text(json.dumps(report, indent=2) + '\n')
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
