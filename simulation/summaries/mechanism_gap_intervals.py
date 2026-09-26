"""Bootstrap intervals for the mean-position gap of the fixed-sequence simulation runs (Fig. 2a,d).

For every run, cell and model listed in the plan (the output of build_mechanism_source_data.py), the
script reads each dataset's fitted orderings from the run records and computes the event-position shift
of each event (its position in the CN-heavy group g1 minus its position in the AD-heavy group g3). The
mean-position gap is the mean over events of the absolute mean shift. Its 95% interval is given by the
2.5th and 97.5th percentiles of the gap over 2,000 bootstrap resamples of datasets (multinomial weights)
drawn from one random stream (NumPy default_rng, seed 20260919) that runs through the groups in sorted
order of (run, cell, model). The gap can lie outside its interval because |mean| is biased upward.

Because the random stream is shared, a group's interval depends on the groups sorted before it. The
published intervals were computed over the 59 groups of the published mechanism_all_cells.csv, which
include two diagnostic runs that no figure uses (method_oracle_a, method_oracle_estprior_a); pass that
file as --plan to reproduce them. A plan group without records is left out of the output, but its
resamples are still drawn, so that the intervals of all other groups are unchanged.

Usage:
  python simulation/summaries/mechanism_gap_intervals.py --runs-dir RUNS_DIR \
      --plan mechanism_all_cells.csv --out v2_mechanism_gap_intervals.csv [--fit-metrics fits.csv]
"""
import argparse, json, os, re
from pathlib import Path
import numpy as np
import pandas as pd

SEED = 20260919
RESAMPLES = 2000
K = 14


def load_fits(runs_dir, plan, cache=None):
    """Per-dataset positions and errors for the plan groups, read from the run records.

    With cache set, an existing file is read instead of the records, and a new one is written.
    """
    if cache is not None and cache.exists(): return pd.read_csv(cache)
    allowed=set(zip(plan.run_id,plan.cell,plan.engine)); rows=[]; dec=json.JSONDecoder()
    for run in sorted(plan.run_id.unique()):
        if not (runs_dir/run).is_dir():
            print('No records:',run,flush=True)
            continue
        for f in sorted((runs_dir/run).rglob('*rows*.jsonl')):
            with f.open() as stream:
                for line in stream:
                    if run=='method_saebm_a':
                        ob=json.loads(line);cell=ob['cell'];engine=ob['variant'];seed=ob['seed']
                        assert ob['status']=='ok'
                    else:
                        head=line[:line.find('"truth"')]
                        if '"status": "ok"' not in head: continue
                        cell=re.search(r'"cell"\s*:\s*"([^"]+)"',head).group(1)
                        engine=re.search(r'"engine"\s*:\s*"([^"]+)"',head).group(1)
                        seed=int(re.search(r'"seed"\s*:\s*([0-9]+)',head).group(1))
                        match=re.search(r'"observed"\s*:\s*(\{\s*"orderings")',line)
                        if match is None: raise ValueError((run,cell,engine,seed,'missing saved orderings'))
                        ob,_=dec.raw_decode(line,match.start(1))
                    if (run,cell,engine) not in allowed: continue
                    orders=np.asarray(ob['orderings'],int)
                    assert orders.shape==(3,K) and all(np.array_equal(np.sort(o),np.arange(K)) for o in orders)
                    pos=np.argsort(orders,axis=1)+1
                    row=dict(run_id=run,cell=cell,engine=engine,seed=seed)
                    row.update({f'error_g{g+1}':float(x) for g,x in enumerate(ob['distance_to_truth'])})
                    ii,jj=np.triu_indices(K,1)
                    row.update({f'distance_{a+1}{b+1}':float(np.mean((pos[a,ii]-pos[a,jj])*(pos[b,ii]-pos[b,jj])<0)) for a,b in [(0,1),(0,2),(1,2)]})
                    row.update({f'delta_{e}':int(pos[0,e]-pos[2,e]) for e in range(K)})
                    row.update({f'position_g1_{e}':int(pos[0,e]) for e in range(K)})
                    rows.append(row)
        print('Read records:',run,flush=True)
    columns=['run_id','cell','engine','seed']
    fits=pd.DataFrame(rows,columns=None if rows else columns)
    assert not fits.duplicated(columns).any()
    if cache is not None:
        cache.parent.mkdir(parents=True,exist_ok=True)
        tmp=cache.with_suffix('.tmp');fits.to_csv(tmp,index=False);tmp.replace(cache)
    return fits


def gap_summary(fits, plan):
    """Mean-position gap and bootstrap interval for each plan group, in sorted order."""
    counts=fits.groupby(['run_id','cell','engine']).size() if len(fits) else pd.Series(dtype=int)
    groups=dict(tuple(fits.groupby(['run_id','cell','engine'],sort=True))) if len(fits) else {}
    expected={(r.run_id,r.cell,r.engine):int(r.R) for r in plan.itertuples()}
    rows=[];rng=np.random.default_rng(SEED)
    for key in sorted(expected):
        n_plan=expected[key]
        if key not in groups:
            # no records: draw (and discard) this group's resamples so later groups keep their intervals
            rng.multinomial(n_plan,np.full(n_plan,1/n_plan),size=RESAMPLES)
            print('No records, left out:',*key,flush=True)
            continue
        assert counts[key]==n_plan,(key,int(counts[key]),n_plan)
        g=groups[key]
        d=g[[f'delta_{i}' for i in range(K)]].to_numpy(float);n=len(d)
        weights=rng.multinomial(n,np.full(n,1/n),size=RESAMPLES)
        samples=np.abs(weights@d/n).mean(axis=1)
        rows.append(dict(zip(['run_id','cell','engine'],key))|dict(R=n,gap=np.abs(d.mean(axis=0)).mean(),low=np.quantile(samples,.025),high=np.quantile(samples,.975)))
    return pd.DataFrame(rows,columns=['run_id','cell','engine','R','gap','low','high'])


def main():
    ap=argparse.ArgumentParser(description='Bootstrap intervals for the mean-position gap (Fig. 2a,d).')
    ap.add_argument('--runs-dir',type=Path,default=os.environ.get('CONCORD_RUNS_DIR'),
                    help='directory with one subdirectory of records per run (default: $CONCORD_RUNS_DIR)')
    ap.add_argument('--plan',type=Path,required=True,
                    help='mechanism_all_cells.csv listing the run/cell/model groups and their numbers of datasets')
    ap.add_argument('--out',type=Path,required=True,help='output CSV (v2_mechanism_gap_intervals.csv)')
    ap.add_argument('--fit-metrics',type=Path,default=None,
                    help='optional per-dataset table; read instead of the records if it exists, otherwise written')
    args=ap.parse_args()
    if args.runs_dir is None and (args.fit_metrics is None or not args.fit_metrics.exists()):
        ap.error('give --runs-dir or set CONCORD_RUNS_DIR')
    plan=pd.read_csv(args.plan)
    fits=load_fits(Path(args.runs_dir) if args.runs_dir is not None else None,plan,args.fit_metrics)
    gaps=gap_summary(fits,plan)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    gaps.to_csv(args.out,index=False)
    print(f'{len(gaps)} run/cell/model intervals -> {args.out}')


if __name__=='__main__':main()
