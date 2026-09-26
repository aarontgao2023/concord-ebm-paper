"""SA-EBM (pysaebm; Hao et al. 2025) on the complete-data fixed-sequence cells (run method_saebm_a; Fig. 2a,c).

Each dataset is simulated with design_v2 with one fixed true sequence and missing=False, because pysaebm
needs a complete participant-by-biomarker matrix. Each simulated group is fitted separately with
run_ebm(algorithm='conjugate_priors'); CN participants are non-diseased and MCI and AD participants are
diseased. The reported ordering is pysaebm's order_with_highest_ll. The pysaebm seed is the dataset seed
modulo 100000. Each dataset gives one JSON line with the three group orderings, their normalized Kendall
distances to the true sequence and the fit times (rows_chunkNNN.jsonl in --output), which
simulation/summaries/ reads.

The run in the paper used cells IID_S, REF_S, BAL_S and REF_S_N4, seeds 53400000-53400999, 10,000
iterations and 2,500 burn-in iterations (hpc/saebm_mechanism_freeze.sbatch); the argument defaults below
are those of an earlier development run. Run it in an environment with pysaebm (the paper used pysaebm
7.7.7). design_v2 is imported from the directory in the environment variable EBM_SNAPSHOT if it is set,
otherwise from simulation/scripts.
"""
import argparse, json, os, sys, time, tempfile, shutil, warnings, logging
from multiprocessing import get_context
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
SNAP = os.environ.get('EBM_SNAPSHOT', str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, SNAP)
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]
MECH_CELLS = {'IID_S': ('IID_H0', {}), 'REF_S': ('REF_H0', {}), 'BAL_S': ('BALANCED_H0', {}),
              'IID_S_N4': ('IID_H0', {'sample_scale': 4}), 'REF_S_N4': ('REF_H0', {'sample_scale': 4}),
              'BAL_S_N4': ('BALANCED_H0', {'sample_scale': 4})}

def kendall_norm(a, b):
    ra = np.argsort(np.asarray(a)); rb = np.argsort(np.asarray(b)); n = len(ra); d = 0
    for i in range(n):
        for j in range(i + 1, n): d += (ra[i] - ra[j]) * (rb[i] - rb[j]) < 0
    return d / (n * (n - 1) / 2)

def fit_group(sub, names, n_iter, burn_in, seed):
    from pysaebm import run_ebm
    sub = sub.reset_index(drop=True)
    rows = [dict(participant=i, biomarker=b, measurement=float(r[b]), diseased=bool(r['Diagnosis'] != 'CN'))
            for i, r in sub.iterrows() for b in names]
    tmp = tempfile.mkdtemp(prefix='saebm_')
    try:
        f = os.path.join(tmp, 'g.csv'); pd.DataFrame(rows).to_csv(f, index=False)
        t0 = time.monotonic()
        res = run_ebm(algorithm='conjugate_priors', data_file=f, output_dir=os.path.join(tmp, 'out'), n_iter=n_iter,
                      burn_in=burn_in, skip_heatmap=True, skip_traceplot=True, save_results=False, save_details=True, seed=seed)
        secs = time.monotonic() - t0
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    pos = res['order_with_highest_ll']                      # biomarker name -> position
    order = [names.index(b) for b, _ in sorted(pos.items(), key=lambda kv: kv[1])]
    return order, secs, float(res.get('healthy_ratio', float('nan')))

def job(args):
    cell, seed, n_iter, burn_in = args
    from design_v2 import resolve, simulate
    rec = dict(mode='mechanism', cell=cell, seed=seed, variant='saebm_conjugate')
    try:
        name, ov = MECH_CELLS[cell]
        df, truth = simulate(resolve(name, fixed_base_order=tuple(FIXED_S), missing=False, **ov), seed)
        names = [c for c in df.columns if c not in ('PTID', 'Diagnosis', 'EXAMDATE', 'APOE')]
        orders, secs = [], []
        for g in sorted(df.APOE.unique()):
            o, s, hr = fit_group(df[df.APOE == g], names, n_iter, burn_in, seed % 100000); orders.append(o); secs.append(round(s, 1))
        S = truth['base_order']
        rec.update(status='ok', orderings=orders, distance_to_truth=[kendall_norm(o, S) for o in orders], fit_seconds=secs)
    except Exception as e:
        rec.update(status='error', error=f'{type(e).__name__}: {e}')
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', default='IID_S,REF_S,BAL_S')
    ap.add_argument('--seed-start', type=int, default=31800000)
    ap.add_argument('--datasets', type=int, default=200)
    ap.add_argument('--n-iter', type=int, default=2000); ap.add_argument('--burn-in', type=int, default=500)
    ap.add_argument('--chunk', type=int, default=0); ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--workers', type=int, default=16)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    seeds = list(range(a.seed_start, a.seed_start + a.datasets))[a.chunk::a.nchunks]
    a.output.mkdir(parents=True, exist_ok=True)
    out = a.output / f'rows_chunk{a.chunk:03d}.jsonl'
    done = set()
    if out.exists():
        for line in open(out):
            r = json.loads(line)
            if r.get('status') == 'ok': done.add((r['cell'], r['seed']))
    jobs = [(c, s, a.n_iter, a.burn_in) for c in a.cells.split(',') for s in seeds if (c, s) not in done]
    print(f'saebm chunk {a.chunk}/{a.nchunks}: {len(jobs)} jobs ({len(done)} done)', flush=True)
    t0 = time.monotonic(); n = 0
    with open(out, 'a', buffering=1) as fh, get_context('spawn').Pool(a.workers) as pool:
        for rec in pool.imap_unordered(job, jobs, chunksize=1):
            fh.write(json.dumps(rec) + '\n'); n += 1
            if n % 10 == 0: print(f'  {n}/{len(jobs)}  {(time.monotonic()-t0)/60:.1f} min', flush=True)
    print(f'done {n} jobs in {(time.monotonic()-t0)/3600:.2f} h', flush=True)

if __name__ == '__main__':
    main()
