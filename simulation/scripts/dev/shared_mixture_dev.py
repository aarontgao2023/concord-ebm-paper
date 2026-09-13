"""Development: shared-mixture DEBM.

One two-component mixture per biomarker is fitted on the POOLED sample (no group labels), the
resulting parameters are handed to every group, and only the per-group Mallows consensus uses
the labels. This removes the composition dependence of the mixture step by construction; the
question is whether the fitted group orderings then share one target under a common true
ordering (mechanism mode) and whether unrestricted / stratified permutation calibrate
(calibration mode). Because the pooled mixture does not depend on labels, it is cached from the
observed fit and every permutation refits only the consensus.

Implemented as a patch on pyebm.core_utilities.do_mixturemodel inside the repaired-engine
context of engine_v2; nothing else in the pipeline changes. Development seeds only."""
import argparse, json, os, signal, sys, time, hashlib
from multiprocessing import get_context
from pathlib import Path
import numpy as np

SNAP = os.environ.get('EBM_SNAPSHOT', '/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2')
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

def permute(df, seed, scheme, pid):
    tag = int.from_bytes(hashlib.sha256(scheme.encode()).digest()[:4], 'little')
    rng = np.random.default_rng(np.random.SeedSequence([seed, tag, pid, 63000000]))
    labels = df.APOE.to_numpy().copy()
    strata = np.zeros(len(df), dtype=int) if scheme == 'unrestricted' else df.Diagnosis.to_numpy()
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s); labels[idx] = rng.permutation(labels[idx])
    out = df.copy(); out['APOE'] = labels; return out

class _Timeout(Exception): pass
def _alarm(*_): raise _Timeout()

def fit_shared(df, cache):
    """cache: dict; if empty, the pooled mixture is fitted and stored; else reused."""
    from engine_v2 import EngineConfig, _temporary_engine, map_orderings
    from pyebm import debm
    import pyebm.core_utilities as cu
    orig = cu.do_mixturemodel
    def shared(DMO, data_AD_raw, data_CN_raw, Data_all, Groups, GroupValues, GroupValues_cn, GroupValues_ad,
               HyperParams=1, flag_init_together=1, only_init=0):
        if 'params' not in cache:
            BP, p_yes, p_no, lpost, lpre = orig(DMO, data_AD_raw, data_CN_raw, Data_all, [], [], [], [],
                                                HyperParams=HyperParams, flag_init_together=1, only_init=only_init)
            cache['params'] = (np.array(BP.Control, copy=True), np.array(BP.Disease, copy=True), np.array(BP.Mixing, copy=True))
            cache['post'] = (np.array(p_yes, copy=True), np.array(p_no, copy=True), np.array(lpost, copy=True), np.array(lpre, copy=True))
            cache['n'] = Data_all.shape[0]
        if cache['n'] != Data_all.shape[0]: raise RuntimeError('cached posterior does not match this dataset')
        C, D, M = cache['params']; p_yes, p_no, lpost, lpre = cache['post']
        k = len(np.unique(GroupValues[0])) if len(Groups) else 1
        from collections import namedtuple
        BPn = namedtuple('BiomarkerParams', 'Control Disease Mixing')
        BPn.Control = [C.copy() for _ in range(k)]; BPn.Disease = [D.copy() for _ in range(k)]; BPn.Mixing = [M.copy() for _ in range(k)]
        return BPn, p_yes.copy(), p_no.copy(), lpost.copy(), lpre.copy()
    diag = {}
    cu.do_mixturemodel = shared
    try:
        with _temporary_engine(EngineConfig(mode='repaired', expected_events=14), diag):
            t0 = time.monotonic()
            model, _, _ = debm.fit(df, Factors=[], Labels=['CN', 'MCI', 'AD'], Groups=['APOE'])
            secs = time.monotonic() - t0
            names = [c for c in df.columns if c not in ('PTID', 'Diagnosis', 'EXAMDATE', 'APOE')]
            orders = [list(map(int, o)) for o in map_orderings(model.MeanCentralOrdering, list(model.BiomarkerList), names, 3)]
    finally:
        cu.do_mixturemodel = orig
    return orders, secs

def job(args):
    mode, cell, seed, scheme, pid, timeout = args
    from design_v2 import resolve, simulate
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(timeout)
    rec = dict(mode=mode, cell=cell, seed=seed, scheme=scheme, pid=pid)
    try:
        if mode == 'mechanism':
            name, ov = MECH_CELLS[cell]
            df, truth = simulate(resolve(name, fixed_base_order=tuple(FIXED_S), **ov), seed)
            orders, secs = fit_shared(df, {})
            S = truth['base_order']
            rec.update(status='ok', orderings=orders, distance_to_truth=[kendall_norm(o, S) for o in orders], fit_seconds=round(secs, 1))
        else:  # calibration on REF: observed + all permutations of one dataset in ONE job (cache reuse)
            df, truth = simulate(resolve('REF_H0'), seed)
            cache = {}
            obs, secs = fit_shared(df, cache)
            taus_obs = [kendall_norm(obs[0], obs[1]), kendall_norm(obs[0], obs[2]), kendall_norm(obs[1], obs[2])]
            B = pid  # number of permutations per scheme
            out = {'observed': taus_obs, 'mixture_seconds': round(secs, 1), 'perm': {}}
            t1 = time.monotonic()
            for sc in ('unrestricted', 'diagnosis'):
                T = []
                for b in range(B):
                    o, _ = fit_shared(permute(df, seed, sc, b), cache)
                    T.append([kendall_norm(o[0], o[1]), kendall_norm(o[0], o[2]), kendall_norm(o[1], o[2])])
                out['perm'][sc] = T
            out['perm_seconds'] = round(time.monotonic() - t1, 1)
            # consistency check: cached path with original labels must reproduce the observed fit
            chk, _ = fit_shared(df, cache)
            out['cache_consistent'] = (chk == obs)
            rec.update(status='ok', **out)
    except _Timeout:
        rec.update(status='timeout')
    except Exception as e:
        rec.update(status='error', error=f'{type(e).__name__}: {e}')
    finally:
        signal.alarm(0)
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['mechanism', 'calibration'], required=True)
    ap.add_argument('--seed-start', type=int, default=31800000)
    ap.add_argument('--datasets', type=int, default=200)
    ap.add_argument('--bperm', type=int, default=199)
    ap.add_argument('--chunk', type=int, default=0); ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--workers', type=int, default=16); ap.add_argument('--timeout', type=int, default=7000)
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
    jobs = []
    if a.mode == 'mechanism':
        for cell in MECH_CELLS:
            for s in seeds:
                if (cell, s) not in done: jobs.append(('mechanism', cell, s, None, None, a.timeout))
    else:
        for s in seeds:
            if ('REF_H0', s) not in done: jobs.append(('calibration', 'REF_H0', s, None, a.bperm, a.timeout))
    print(f'{a.mode} chunk {a.chunk}/{a.nchunks}: {len(jobs)} jobs ({len(done)} done)', flush=True)
    t0 = time.monotonic(); n = 0
    with open(out, 'a', buffering=1) as fh, get_context('spawn').Pool(a.workers) as pool:
        for rec in pool.imap_unordered(job, jobs, chunksize=1):
            fh.write(json.dumps(rec) + '\n'); n += 1
            if n % 20 == 0: print(f'  {n}/{len(jobs)}  {(time.monotonic()-t0)/60:.1f} min', flush=True)
    print(f'done {n} jobs in {(time.monotonic()-t0)/3600:.2f} h', flush=True)

if __name__ == '__main__':
    main()
