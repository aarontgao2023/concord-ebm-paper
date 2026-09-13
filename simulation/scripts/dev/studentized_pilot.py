"""Development pilot: a studentised ordering statistic under unrestricted permutation.

Question: does dividing the Kendall distance by each group's own estimation scale restore
calibration WITHOUT stratifying the permutation (the Chung-Romano / Behrens-Fisher route)?

Per dataset (REF_H0, development seeds): one observed fit, B unrestricted permutations,
B diagnosis-stratified permutations. Every fit uses the repaired engine and records
  taus[3]   plain normalised Kendall distances between fitted group orderings
  vhat[3]   per-group scale: mean d_K between the fitted ordering and NBOOT orderings
            refitted at the consensus step from bootstrap resamples of that group's
            subjects (mixture parameters fixed, same repaired consensus)
so plain and studentised statistics come from the same permutations. Decisions are taken
offline. Full budget, no early stopping. Appends JSON lines; safe to resume.
This script does not touch any frozen run and uses only development seeds."""
import argparse, json, os, signal, sys, time, hashlib
from multiprocessing import get_context
from pathlib import Path
import numpy as np

SNAP = os.environ.get('EBM_SNAPSHOT', '/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2')
sys.path.insert(0, SNAP)

def kendall_norm(a, b):
    ra = np.argsort(np.asarray(a)); rb = np.argsort(np.asarray(b)); n = len(ra)
    disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            disc += (ra[i] - ra[j]) * (rb[i] - rb[j]) < 0
    return disc / (n * (n - 1) / 2)

def permute(df, seed, scheme, pid):
    tag = int.from_bytes(hashlib.sha256(scheme.encode()).digest()[:4], 'little')
    rng = np.random.default_rng(np.random.SeedSequence([seed, tag, pid, 62000000]))
    labels = df.APOE.to_numpy().copy()
    strata = np.zeros(len(df), dtype=int) if scheme == 'unrestricted' else df.Diagnosis.to_numpy()
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s); labels[idx] = rng.permutation(labels[idx])
    out = df.copy(); out['APOE'] = labels; return out

class _Timeout(Exception): pass
def _alarm(*_): raise _Timeout()

def fit_with_scale(df, nboot, boot_seed):
    """Repaired-engine fit; capture per-group fitMallows inputs/outputs; bootstrap scale."""
    from engine_v2 import EngineConfig, _temporary_engine
    from pyebm import debm
    from pyebm.central_ordering import generalized_mallows as gm
    captured = []
    orig = gm.weighted_mallows.__dict__['fitMallows']          # descriptor, for restoration
    orig_fn = gm.weighted_mallows.fitMallows                    # bound classmethod: (p_yes, mixing)
    def rec(p_yes, mixing):
        pi0, ecs, het = orig_fn(p_yes, mixing)
        captured.append((np.array(p_yes, copy=True), np.array(mixing, copy=True), list(pi0)))
        return pi0, ecs, het
    gm.weighted_mallows.fitMallows = staticmethod(rec)
    diag = {}
    try:
        with _temporary_engine(EngineConfig(mode='repaired', expected_events=14), diag):
            t0 = time.monotonic()
            model, _, _ = debm.fit(df, Factors=[], Labels=['CN', 'MCI', 'AD'], Groups=['APOE'])
            t_fit = time.monotonic() - t0
            if len(captured) != 3: raise RuntimeError(f'captured {len(captured)} groups')
            pis = [c[2] for c in captured]
            taus = [kendall_norm(pis[0], pis[1]), kendall_norm(pis[0], pis[2]), kendall_norm(pis[1], pis[2])]
            rng = np.random.default_rng(boot_seed); vhat = []; t1 = time.monotonic()
            for p_yes, mixing, pi0 in captured:
                n = p_yes.shape[0]; d = []
                for _ in range(nboot):
                    idx = rng.integers(0, n, n)
                    pib, _, _ = orig_fn(p_yes[idx], mixing)   # patched repaired consensus applies
                    d.append(kendall_norm(pib, pi0))
                vhat.append(float(np.mean(d)))
            t_boot = time.monotonic() - t1
    finally:
        gm.weighted_mallows.fitMallows = orig
    return dict(status='ok', taus=[float(t) for t in taus], vhat=vhat, n_groups=[int(c[0].shape[0]) for c in captured],
                fit_seconds=round(t_fit, 1), boot_seconds=round(t_boot, 1),
                optimizer_failures=int(sum(not i.get('success', True) for i in diag.get('optimizer', []))))

def job(args):
    seed, scheme, pid, nboot, timeout = args
    from design_v2 import resolve, simulate
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(timeout)
    rec = dict(seed=seed, scheme=scheme, pid=pid)
    try:
        df, truth = simulate(resolve('REF_H0'), seed)
        if scheme is not None: df = permute(df, seed, scheme, pid)
        rec.update(fit_with_scale(df, nboot, boot_seed=(seed * 1000003 + (pid or 0) * 7919 + (0 if scheme is None else 1 if scheme == 'unrestricted' else 2))))
    except _Timeout:
        rec.update(status='timeout')
    except Exception as e:
        rec.update(status='error', error=f'{type(e).__name__}: {e}')
    finally:
        signal.alarm(0)
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed-start', type=int, default=31500000)
    ap.add_argument('--datasets', type=int, default=200)
    ap.add_argument('--chunk', type=int, default=0); ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--bperm', type=int, default=199); ap.add_argument('--nboot', type=int, default=10)
    ap.add_argument('--workers', type=int, default=16); ap.add_argument('--timeout', type=int, default=900)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    seeds = list(range(a.seed_start, a.seed_start + a.datasets))[a.chunk::a.nchunks]
    a.output.mkdir(parents=True, exist_ok=True)
    out = a.output / f'rows_chunk{a.chunk:03d}.jsonl'
    done = set()
    if out.exists():
        for line in open(out):
            r = json.loads(line)
            if r.get('status') == 'ok': done.add((r['seed'], r['scheme'], r['pid']))
    jobs = []
    for s in seeds:
        if (s, None, None) not in done: jobs.append((s, None, None, a.nboot, a.timeout))
        for sc in ('unrestricted', 'diagnosis'):
            for b in range(a.bperm):
                if (s, sc, b) not in done: jobs.append((s, sc, b, a.nboot, a.timeout))
    print(f'chunk {a.chunk}/{a.nchunks}: {len(seeds)} seeds, {len(jobs)} jobs to run ({len(done)} done)', flush=True)
    t0 = time.monotonic(); n = 0
    with open(out, 'a', buffering=1) as fh, get_context('spawn').Pool(a.workers) as pool:
        for rec in pool.imap_unordered(job, jobs, chunksize=1):
            fh.write(json.dumps(rec) + '\n'); n += 1
            if n % 50 == 0: print(f'  {n}/{len(jobs)}  {(time.monotonic()-t0)/60:.1f} min', flush=True)
    print(f'done {n} jobs in {(time.monotonic()-t0)/3600:.2f} h', flush=True)

if __name__ == '__main__':
    main()
