"""Summarize the fixed-sequence simulation runs (Fig. 2a,c,d): one row per run, cell and model.

Reads the per-dataset records of run_v2.py (RUNS_DIR/<run>/cell_*/chunk_*/rows.jsonl) and of
scripts/dev/saebm_mechanism_dev.py (RUNS_DIR/method_saebm_a/rows_chunk*.jsonl). Only the fitted
orderings of the three simulated groups are used; nothing is refitted.

Groups: g1 = CN-heavy (code e2), g2 = intermediate (e33), g3 = AD-heavy (e4). Positions are ranks in a
group's fitted ordering. Output columns:
  R                                   datasets
  mean_error_g*                       error against the true sequence (normalized Kendall distance), mean over datasets
  mean_distance_g*_g*                 normalized Kendall distance between two groups' orderings, mean over datasets
  mean_position_gap_g1_g3             mean-position gap: mean over events of |mean position in g1 - mean position in g3|
  mean_per_cohort_position_gap_g1_g3  the same absolute differences averaged within each dataset
  majority_wrong_pairs_g*             systematically misordered event pairs: pairs that the group's ordering
                                      reverses relative to the true sequence in more than half of the datasets
  majority_pairwise_error_g*          the same count divided by the number of event pairs (91)

The default runs are those of the published mechanism_all_cells.csv: the five runs used in the paper and
two diagnostic runs with known event models (method_oracle_a, method_oracle_estprior_a) that no figure
uses and whose configurations are not included here. Runs without records are skipped.

Usage:
  python simulation/summaries/build_mechanism_source_data.py --runs-dir RUNS_DIR --out mechanism_all_cells.csv
"""
import argparse, csv, itertools, json, os
from pathlib import Path

PAPER_RUNS = ['method_mechanism_a', 'method_mechanism_16n_a', 'method_saebm_reference_a',
              'method_saebm_a', 'method_kde_mechanism_a']
RUNS = [
    'method_mechanism_a', 'method_mechanism_16n_a', 'method_oracle_a',
    'method_oracle_estprior_a', 'method_saebm_reference_a',
    'method_saebm_a', 'method_kde_mechanism_a',
]
FIXED_ORDER = [5,9,1,12,11,10,3,8,13,0,6,7,2,4]


def summarize(runs_dir, runs):
    acc = {}
    for run in runs:
        if not (runs_dir / run).is_dir():
            print(f'{run}: no records in {runs_dir}; skipped')
            continue
        for path in (runs_dir / run).rglob('*rows*.jsonl'):
            for line in path.open():
                row = json.loads(line)
                assert row['status'] == 'ok', (run, row['cell'], row['seed'])
                obs = row.get('observed', row)
                orderings = obs['orderings']
                errors = obs['distance_to_truth']
                engine = row.get('engine', row.get('variant'))
                key = run, row['cell'], engine
                k = len(orderings[0])
                if 'truth' in row:
                    n = sum(row['truth']['group_sizes'].values())
                    assert row['truth']['config']['fixed_base_order'] is not None
                    missing = row['truth']['config']['missing']
                    assert row['truth']['base_order'] == FIXED_ORDER
                else:
                    n = 3884 if 'N4' in row['cell'] else 971
                    missing = False
                pairs=list(itertools.combinations(range(k),2))
                truth_positions=[FIXED_ORDER.index(i) for i in range(k)]
                value = acc.setdefault(key, dict(R=0, n=n, K=k, missing=missing,
                     pair_errors=[[0]*len(pairs) for _ in range(3)],
                     errors=[0.0]*3, pair_distances=[0.0]*3, positions=[[0.0]*k for _ in range(3)],
                     per_cohort_position_gap=0.0, seeds=set()))
                assert row['seed'] not in value['seeds'], key
                value['seeds'].add(row['seed']); value['R'] += 1
                positions=[]
                for g, order in enumerate(orderings):
                    pos=[0]*k
                    for rank, event in enumerate(order):pos[int(event)]=rank
                    positions.append(pos)
                    value['positions'][g]=[a+b for a,b in zip(value['positions'][g],pos)]
                    value['errors'][g]+=errors[g]
                    pair_errors=[int((pos[i]-pos[j])*(truth_positions[i]-truth_positions[j])<0) for i,j in pairs]
                    assert abs(sum(pair_errors)/len(pairs)-errors[g])<1e-12,(run,key,row['seed'],g)
                    value['pair_errors'][g]=[a+b for a,b in zip(value['pair_errors'][g],pair_errors)]
                value['per_cohort_position_gap']+=sum(abs(a-b) for a,b in zip(positions[0],positions[2]))/k
                for j,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
                    value['pair_distances'][j]+=sum((positions[a][i]-positions[a][l])*(positions[b][i]-positions[b][l])<0 for i in range(k) for l in range(i+1,k))/(k*(k-1)/2)
    rows=[]
    for (run,cell,engine),v in sorted(acc.items()):
        r=v['R']
        rows.append(dict(run_id=run,cell=cell,engine=engine,R=r,n_subjects=v['n'],K_events=v['K'],missingness_enabled=v['missing'],
            true_order_design='fixed across cohorts',
            mean_error_g1=v['errors'][0]/r,mean_error_g2=v['errors'][1]/r,mean_error_g3=v['errors'][2]/r,
            mean_distance_g1_g2=v['pair_distances'][0]/r,mean_distance_g1_g3=v['pair_distances'][1]/r,mean_distance_g2_g3=v['pair_distances'][2]/r,
            mean_position_gap_g1_g3=sum(abs(a-b) for a,b in zip(v['positions'][0],v['positions'][2]))/r/v['K'],
            mean_per_cohort_position_gap_g1_g3=v['per_cohort_position_gap']/r,
            majority_wrong_pairs_g1=sum(x>r/2 for x in v['pair_errors'][0]),
            majority_wrong_pairs_g2=sum(x>r/2 for x in v['pair_errors'][1]),
            majority_wrong_pairs_g3=sum(x>r/2 for x in v['pair_errors'][2]),
            majority_pairwise_error_g1=sum(x>r/2 for x in v['pair_errors'][0])/len(v['pair_errors'][0]),
            majority_pairwise_error_g2=sum(x>r/2 for x in v['pair_errors'][1])/len(v['pair_errors'][1]),
            majority_pairwise_error_g3=sum(x>r/2 for x in v['pair_errors'][2])/len(v['pair_errors'][2])))
    return rows


def main():
    ap = argparse.ArgumentParser(description='Summarize the fixed-sequence simulation runs (Fig. 2a,c,d).')
    ap.add_argument('--runs-dir', type=Path, default=os.environ.get('CONCORD_RUNS_DIR'),
                    help='directory with one subdirectory of records per run (default: $CONCORD_RUNS_DIR)')
    ap.add_argument('--out', type=Path, required=True, help='output CSV (mechanism_all_cells.csv)')
    ap.add_argument('--runs', nargs='+', default=RUNS,
                    help='runs to summarize (default: the runs of the published file; the paper uses '
                         + ', '.join(PAPER_RUNS) + ')')
    args = ap.parse_args()
    if args.runs_dir is None:
        ap.error('give --runs-dir or set CONCORD_RUNS_DIR')
    rows = summarize(Path(args.runs_dir), args.runs)
    if not rows:
        raise SystemExit(f'No records found in {args.runs_dir}')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(f'{len(rows)} run/cell/model summaries -> {args.out}')


if __name__ == '__main__':
    main()
