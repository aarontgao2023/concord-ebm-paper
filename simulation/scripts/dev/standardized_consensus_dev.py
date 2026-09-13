"""Development: composition-standardised consensus (direct standardisation inside weighted Mallows).

Follow-up to shared_mixture_dev.py. That run showed the pooled mixture removes ~half of the
composition-dependent target under REF; the residual lives in the consensus step. Here each
subject's contribution to the consensus objective (the mean over subjects of the weighted
Kendall score) is reweighted so that every group has a common diagnostic composition:
    w_i = pi_ref[d_i] / pi_g[d_i],   d_i in {CN, MCI, AD}
with pi_ref either the pooled composition of the whole sample ('pooled') or the renormalised
element-wise minimum over groups ('min'; the sparsest common composition, lowest weight
variance). Weights are recomputed from the labels at every (pseudo-)group fit, so under
permutation they are a function of the permuted labels as required.

Variants: pergroup_std_pooled | shared_std_pooled | shared_std_min (mechanism);
          shared_std_pooled | shared_std_min (calibration, REF_H0, U vs D, B permutations).
Implementation: patches on cu.parse_inputs (capture Diagnosis + Data_all), cu.find_central_ordering
(build per-group weight queue; asserts row-order agreement), gm.weighted_mallows.fitMallows
(pop the group's weights), gm.weighted_mallows.totalconsensus (weighted mean). All inside the
repaired-engine context of engine_v2. Development seeds only."""
import argparse, json, os, signal, sys, time, hashlib
from multiprocessing import get_context
from pathlib import Path
import numpy as np
EXTRA = {}

SNAP = os.environ.get('EBM_SNAPSHOT', '/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2')
sys.path.insert(0, SNAP)
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]
MECH_CELLS = {'IID_S': ('IID_H0', {}), 'REF_S': ('REF_H0', {}), 'BAL_S': ('BALANCED_H0', {}),
              'IID_S_N4': ('IID_H0', {'sample_scale': 4}), 'REF_S_N4': ('REF_H0', {'sample_scale': 4}),
              'BAL_S_N4': ('BALANCED_H0', {'sample_scale': 4})}
VARIANTS = {'plain': (False, None), 'pergroup_std_uniform': (False, 'uniform'),
            'pergroup_std_pooled': (False, 'pooled'), 'shared_std_pooled': (True, 'pooled'),
            'shared_std_min': (True, 'min'), 'shared': (True, None)}

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

def composition_weights(diag, gv, mask, ref):
    """diag: int codes 1/2/3 per row; gv: group value per row; mask: bool rows used.
    Returns {group value: weight vector over that group's masked rows (in row order)} and ESS."""
    gvals = np.unique(gv); gvals = gvals[~(gvals != gvals)] if gvals.dtype.kind == 'f' else gvals
    codes = np.array([1, 2, 3])
    comp = {}
    for g in gvals:
        rows = (gv == g) & mask
        comp[g] = np.array([(diag[rows] == c).mean() for c in codes])
    if ref == 'uniform': pi_ref = None
    elif ref == 'pooled': pi_ref = np.array([(diag[mask] == c).mean() for c in codes])
    elif ref == 'min':
        m = np.min(np.stack([comp[g] for g in gvals]), axis=0); pi_ref = m / m.sum()
    else: raise ValueError(ref)
    W = {}; ess = {}
    for g in gvals:
        rows = (gv == g) & mask; d = diag[rows]
        if pi_ref is None: w = np.ones(rows.sum())
        else:
            w = np.zeros(rows.sum())
            for k, c in enumerate(codes):
                sel = d == c
                if sel.any(): w[sel] = pi_ref[k] / comp[g][k]
        W[g] = w; ess[str(g)] = float(w.sum() ** 2 / (w ** 2).sum())
    return gvals, W, ess

def fit_variant(df, cache, shared, ref):
    from engine_v2 import EngineConfig, _temporary_engine, map_orderings
    from pyebm import debm
    import pyebm.core_utilities as cu
    from pyebm.central_ordering import generalized_mallows as gm
    from collections import namedtuple
    orig_mm = cu.do_mixturemodel; orig_parse = cu.parse_inputs; orig_fco = cu.find_central_ordering
    fm_desc = gm.weighted_mallows.__dict__['fitMallows']; fm_orig = gm.weighted_mallows.fitMallows
    tc_desc = gm.weighted_mallows.__dict__['totalconsensus']; tc_orig = gm.weighted_mallows.totalconsensus
    holder = {}; queue = []; cur = {}; info = {}
    def shared_mm(DMO, data_AD_raw, data_CN_raw, Data_all, Groups, GroupValues, GroupValues_cn, GroupValues_ad,
                  HyperParams=1, flag_init_together=1, only_init=0):
        if 'params' not in cache:
            BP, p_yes, p_no, lpost, lpre = orig_mm(DMO, data_AD_raw, data_CN_raw, Data_all, [], [], [], [],
                                                   HyperParams=HyperParams, flag_init_together=1, only_init=only_init)
            cache['params'] = (np.array(BP.Control, copy=True), np.array(BP.Disease, copy=True), np.array(BP.Mixing, copy=True))
            cache['post'] = (np.array(p_yes, copy=True), np.array(p_no, copy=True), np.array(lpost, copy=True), np.array(lpre, copy=True))
            cache['n'] = Data_all.shape[0]
        if cache['n'] != Data_all.shape[0]: raise RuntimeError('cached posterior does not match this dataset')
        C, D, M = cache['params']; p_yes, p_no, lpost, lpre = cache['post']
        k = len(np.unique(GroupValues[0])) if len(Groups) else 1
        BPn = namedtuple('BiomarkerParams', 'Control Disease Mixing')
        BPn.Control = [C.copy() for _ in range(k)]; BPn.Disease = [D.copy() for _ in range(k)]; BPn.Mixing = [M.copy() for _ in range(k)]
        return BPn, p_yes.copy(), p_no.copy(), lpost.copy(), lpre.copy()
    def parse_cap(*a, **k):
        out = orig_parse(*a, **k)
        holder['diag'] = out[7]['Diagnosis'].to_numpy().astype(int); holder['Data_all'] = out[6]
        return out
    def fco(Data_all, p_yes, BP, Groups, GroupValues, DMO, algo_type, maskidx):
        A = holder['Data_all']
        if A.shape != Data_all.shape or not np.allclose(np.nan_to_num(A[:, :, 0], nan=-9e9), np.nan_to_num(Data_all[:, :, 0], nan=-9e9)):
            raise RuntimeError('row order of Data_all differs from parse_inputs output')
        gv = np.asarray(GroupValues[0]); mask = np.asarray(maskidx, bool) if len(maskidx) else np.ones(len(gv), bool)
        gvals, W, ess = composition_weights(holder['diag'], gv, mask, ref)
        queue[:] = [W[g] for g in gvals]; info['ess'] = ess
        return orig_fco(Data_all, p_yes, BP, Groups, GroupValues, DMO, algo_type, maskidx)
    def fm(p_yes, mixing):
        w = queue.pop(0)
        if len(w) != p_yes.shape[0]: raise RuntimeError(f'weight length {len(w)} != group rows {p_yes.shape[0]}')
        cur['w'] = w
        try: return fm_orig(p_yes, mixing)
        finally: cur.pop('w', None)
    def tc(pi0, D, prob):
        tscore, score, score_indv = tc_orig(pi0, D, prob)
        w = cur.get('w')
        if w is None: raise RuntimeError('totalconsensus called outside a weighted fitMallows')
        if len(w) != len(score): raise RuntimeError('weight/score length mismatch')
        return np.float64(np.sum(w * np.asarray(score, float)) / np.sum(w)), score, score_indv
    diag = {}
    if shared: cu.do_mixturemodel = shared_mm
    if ref is not None:
        cu.parse_inputs = parse_cap; cu.find_central_ordering = fco
        gm.weighted_mallows.fitMallows = staticmethod(fm); gm.weighted_mallows.totalconsensus = staticmethod(tc)
    try:
        with _temporary_engine(EngineConfig(mode='repaired', expected_events=14), diag):
            t0 = time.monotonic()
            model, _, _ = debm.fit(df, Factors=[], Labels=['CN', 'MCI', 'AD'], Groups=['APOE'])
            secs = time.monotonic() - t0
            names = [c for c in df.columns if c not in ('PTID', 'Diagnosis', 'EXAMDATE', 'APOE')]
            orders = [list(map(int, o)) for o in map_orderings(model.MeanCentralOrdering, list(model.BiomarkerList), names, 3)]
    finally:
        cu.do_mixturemodel = orig_mm; cu.parse_inputs = orig_parse; cu.find_central_ordering = orig_fco
        gm.weighted_mallows.fitMallows = fm_desc; gm.weighted_mallows.totalconsensus = tc_desc
    if ref is not None and queue: raise RuntimeError('unconsumed group weights')
    return orders, secs, info.get('ess')

def job(args):
    mode, cell, seed, variant, pid, timeout, extra = args
    global EXTRA; EXTRA = extra
    from design_v2 import resolve, simulate
    signal.signal(signal.SIGALRM, _alarm); signal.alarm(timeout)
    shared, ref = VARIANTS[variant]
    rec = dict(mode=mode, cell=cell, seed=seed, variant=variant, pid=pid)
    try:
        if mode == 'mechanism':
            name, ov = MECH_CELLS[cell]
            df, truth = simulate(resolve(name, fixed_base_order=tuple(FIXED_S), **ov, **EXTRA), seed)
            orders, secs, ess = fit_variant(df, {}, shared, ref)
            S = truth['base_order']
            rec.update(status='ok', orderings=orders, distance_to_truth=[kendall_norm(o, S) for o in orders], ess=ess, fit_seconds=round(secs, 1))
        else:
            df, truth = simulate(resolve('REF_H0'), seed)
            cache = {}
            obs, secs, ess = fit_variant(df, cache, shared, ref)
            taus_obs = [kendall_norm(obs[0], obs[1]), kendall_norm(obs[0], obs[2]), kendall_norm(obs[1], obs[2])]
            out = {'observed': taus_obs, 'ess': ess, 'mixture_seconds': round(secs, 1), 'perm': {}}
            t1 = time.monotonic()
            for sc in ('unrestricted', 'diagnosis'):
                T = []
                for b in range(pid):
                    o, _, _ = fit_variant(permute(df, seed, sc, b), cache, shared, ref)
                    T.append([kendall_norm(o[0], o[1]), kendall_norm(o[0], o[2]), kendall_norm(o[1], o[2])])
                out['perm'][sc] = T
            out['perm_seconds'] = round(time.monotonic() - t1, 1)
            chk, _, _ = fit_variant(df, cache, shared, ref)
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
    ap.add_argument('--mode', choices=['mechanism', 'calibration', 'smoke'], required=True)
    ap.add_argument('--variants', default=None, help='comma list; default per mode')
    ap.add_argument('--seed-start', type=int, default=31800000)
    ap.add_argument('--datasets', type=int, default=200)
    ap.add_argument('--bperm', type=int, default=199)
    ap.add_argument('--chunk', type=int, default=0); ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--workers', type=int, default=16); ap.add_argument('--timeout', type=int, default=7000)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--no-missing', action='store_true', help='simulate complete data (missing=False); for the SA-EBM comparison')
    a = ap.parse_args()
    global EXTRA; EXTRA = {'missing': False} if a.no_missing else {}
    if a.mode == 'smoke':  # plumbing check: uniform weights must reproduce the plain fit exactly
        from design_v2 import resolve, simulate
        df, truth = simulate(resolve('REF_H0', fixed_base_order=tuple(FIXED_S)), a.seed_start)
        o1, s1, _ = fit_variant(df, {}, False, None); o2, s2, e2 = fit_variant(df, {}, False, 'uniform')
        print('plain  ', o1, round(s1, 1)); print('uniform', o2, round(s2, 1), e2)
        print('IDENTICAL' if o1 == o2 else 'MISMATCH'); 
        for v in ('pergroup_std_pooled', 'shared_std_pooled', 'shared_std_min'):
            o, s, e = fit_variant(df, {}, *VARIANTS[v]); print(v, [round(kendall_norm(x, truth['base_order']), 3) for x in o], round(s, 1), e)
        return
    seeds = list(range(a.seed_start, a.seed_start + a.datasets))[a.chunk::a.nchunks]
    variants = a.variants.split(',') if a.variants else (['pergroup_std_pooled', 'shared_std_pooled', 'shared_std_min'] if a.mode == 'mechanism' else ['shared_std_pooled', 'shared_std_min'])
    a.output.mkdir(parents=True, exist_ok=True)
    out = a.output / f'rows_chunk{a.chunk:03d}.jsonl'
    done = set()
    if out.exists():
        for line in open(out):
            r = json.loads(line)
            if r.get('status') == 'ok': done.add((r['cell'], r['seed'], r['variant']))
    jobs = []
    if a.mode == 'mechanism':
        for v in variants:
            for cell in MECH_CELLS:
                for s in seeds:
                    if (cell, s, v) not in done: jobs.append(('mechanism', cell, s, v, None, a.timeout, EXTRA))
    else:
        for v in variants:
            for s in seeds:
                if ('REF_H0', s, v) not in done: jobs.append(('calibration', 'REF_H0', s, v, a.bperm, a.timeout, EXTRA))
    print(f'{a.mode} chunk {a.chunk}/{a.nchunks}: {len(jobs)} jobs ({len(done)} done)', flush=True)
    t0 = time.monotonic(); n = 0
    with open(out, 'a', buffering=1) as fh, get_context('spawn').Pool(a.workers) as pool:
        for rec in pool.imap_unordered(job, jobs, chunksize=1):
            fh.write(json.dumps(rec) + '\n'); n += 1
            if n % 20 == 0: print(f'  {n}/{len(jobs)}  {(time.monotonic()-t0)/60:.1f} min', flush=True)
    print(f'done {n} jobs in {(time.monotonic()-t0)/3600:.2f} h', flush=True)

if __name__ == '__main__':
    main()
