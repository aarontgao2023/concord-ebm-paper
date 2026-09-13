"""Composition-standardised consensus dev analysis (runs/v2/dev/stdcons/).
Mechanism: per variant, the §3-style table (dist / maj-d / pairs>0.5 / max per group) and the
mean |e2-e4| target-position gap, alongside the per-group and shared-mixture baselines.
Calibration: per variant, U vs D rejection rates under the B=199 dev rule, p-shift, ESS."""
import json, glob, sys, itertools, collections
import numpy as np
G = ['e2', 'e33', 'e4']; CELLS = ['IID_S', 'REF_S', 'BAL_S', 'IID_S_N4', 'REF_S_N4', 'BAL_S_N4']
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]; PAIRS = ['e2-e33', 'e2-e4', 'e33-e4']
def rank(order): rk = np.empty(len(order), int); rk[np.asarray(order)] = np.arange(len(order)); return rk
def load(pattern, key):
    rows = collections.defaultdict(list); bad = collections.Counter()
    for f in glob.glob(pattern):
        for line in open(f):
            r = json.loads(line)
            if r['status'] == 'ok': rows[key(r)].append(r)
            else: bad[r['status']] += 1
    return rows, bad
def mech_table(rows, label, get):
    S = rank(FIXED_S); P = list(itertools.combinations(range(14), 2))
    print(f"\n[{label}]  cell  R   " + '  '.join(f"{g}: dist/maj-d/n>0.5/max" for g in G) + "   |e2-e4| gap   ESS e2/e33/e4")
    for cell in CELLS:
        if cell not in rows: continue
        data = rows[cell]; out = []; meanpos = []
        for gi in range(3):
            rv = np.zeros(len(P)); d = []; pos = []
            for r in data:
                ords, dist = get(r); rk = rank(ords[gi]); d.append(dist[gi]); pos.append(rk)
                for k, (i, j) in enumerate(P): rv[k] += ((rk[i]-rk[j])*(S[i]-S[j])) < 0
            rv /= len(data); meanpos.append(np.mean(pos, axis=0))
            out.append(f"{np.mean(d):.3f}/{(rv>0.5).sum()/len(P):.3f}/{int((rv>0.5).sum()):2d}/{rv.max():.2f}")
        gap = np.abs(meanpos[0]-meanpos[2]).mean()
        ess = [r.get('ess') for r in data if r.get('ess')]
        esss = '/'.join(f"{np.mean([e[k] for e in ess]):.0f}" for k in ('0', '1', '2')) if ess else '—'
        print(f"  {cell:9s} {len(data):3d}  " + '   '.join(f"{o:>22s}" for o in out) + f"   {gap:.2f}   {esss}")
pg, _ = load('runs/v2/hpc_results/dev_mechanism_fixedS_a/cell_*/chunk_*/rows.jsonl', lambda r: r['cell'])
mech_table(pg, 'baseline: per-group mixture', lambda r: (r['observed']['orderings'], r['observed']['distance_to_truth']))
sh, _ = load('runs/v2/dev/shared_mixture/shared_mixture_mechanism/rows_chunk*.jsonl', lambda r: r['cell'])
mech_table(sh, 'baseline: shared mixture', lambda r: (r['orderings'], r['distance_to_truth']))
st, bad = load('runs/v2/dev/stdcons/stdcons_mechanism/rows_chunk*.jsonl', lambda r: (r['variant'], r['cell']))
variants = sorted({v for v, _ in st})
for v in variants:
    mech_table({c: st[(v, c)] for (vv, c) in st if vv == v}, f'{v}  non-ok={dict(bad)}', lambda r: (r['orderings'], r['distance_to_truth']))
# calibration
cal, bad = load('runs/v2/dev/stdcons/stdcons_calibration/rows_chunk*.jsonl', lambda r: r['variant'])
def ci(k, n):
    if n == 0: return '—'
    p = k/n; se = np.sqrt(p*(1-p)/n); return f"{p:.3f} [{max(0,p-1.96*se):.3f},{min(1,p+1.96*se):.3f}]"
for v, recs in sorted(cal.items()):
    n = len(recs); print(f"\n[calibration REF_H0, {v}]  datasets ok={n}  non-ok={dict(bad)}  cache_consistent all={all(r['cache_consistent'] for r in recs)}")
    res = {}
    for sc in ('unrestricted', 'diagnosis'):
        pv = []; rej_any = 0; rej_max = 0; rej_pair = np.zeros(3, int)
        for r in recs:
            obs = np.array(r['observed']); T = np.array(r['perm'][sc]); B = len(T)
            cnt = (T >= obs[None, :] - 1e-12).sum(0); pv.append((cnt+1)/(B+1))
            rej_pair += cnt <= 2; rej_any += bool((cnt <= 2).any())
            rej_max += (T.max(1) >= obs.max() - 1e-12).sum() <= 9
        pv = np.array(pv); res[sc] = pv
        print(f"  {sc:12s} B={B}  mean p: " + '  '.join(f"{PAIRS[i]} {pv[:,i].mean():.3f}" for i in range(3)))
        print(f"  {'':12s} pair rej (p<=0.015): " + '  '.join(f"{PAIRS[i]} {ci(rej_pair[i], n)}" for i in range(3)))
        print(f"  {'':12s} family any-pair {ci(rej_any, n)}   max-stat (<=9) {ci(rej_max, n)}   ceilings 0.015 / 0.045 / 0.050")
    d = res['diagnosis'] - res['unrestricted']
    print("  p-shift (D - U): " + '  '.join(f"{PAIRS[i]} {d[:,i].mean():+.3f}" for i in range(3)) + "   mean observed d_K: " + '  '.join(f"{np.mean([r['observed'][i] for r in recs]):.3f}" for i in range(3)))
    print(f"  ESS e2/e33/e4: " + '/'.join(f"{np.mean([r['ess'][k] for r in recs]):.0f}" for k in ('0','1','2')) + f"   perms {np.mean([r['perm_seconds'] for r in recs])/60:.0f} min/dataset")
