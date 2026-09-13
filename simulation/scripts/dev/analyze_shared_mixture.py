"""Shared-mixture DEBM dev analysis.
Mechanism: same table as analyze_mechanism_fixedS.py (dist / maj-d / pairs>0.5 / max) plus the
mean |e2-e4| target-position gap over the 14 biomarkers, computed identically for the per-group
mixture run (dev_mechanism_fixedS_a) and the shared-mixture run so the two are comparable.
Calibration (REF_H0, B=199, U vs D): per-pair p, dev rule (pair rejects iff #{T_b>=T_obs}<=2;
family = any pair; max iff <=9), p-shift, and cache consistency."""
import json, glob, sys, itertools, collections
import numpy as np
G = ['e2', 'e33', 'e4']; CELLS = ['IID_S', 'REF_S', 'BAL_S', 'IID_S_N4', 'REF_S_N4', 'BAL_S_N4']
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]
def rank(order): rk = np.empty(len(order), int); rk[np.asarray(order)] = np.arange(len(order)); return rk
def load_pergroup(root):
    rows = collections.defaultdict(list)
    for f in glob.glob(f'{root}/cell_*/chunk_*/rows.jsonl'):
        for line in open(f):
            r = json.loads(line)
            if r['status'] == 'ok': rows[r['cell']].append((r['observed']['orderings'], r['observed']['distance_to_truth']))
    return rows
def load_shared(root):
    rows = collections.defaultdict(list); bad = collections.Counter()
    for f in glob.glob(f'{root}/rows_chunk*.jsonl'):
        for line in open(f):
            r = json.loads(line)
            if r['status'] == 'ok': rows[r['cell']].append((r['orderings'], r['distance_to_truth']))
            else: bad[r['status']] += 1
    return rows, bad
def mech_table(rows, label):
    S = rank(FIXED_S); P = list(itertools.combinations(range(14), 2))
    print(f"\n[{label}]  cell  R   " + '  '.join(f"{g}: dist/maj-d/n>0.5/max" for g in G) + "   |e2-e4| gap")
    for cell in CELLS:
        if cell not in rows: continue
        data = rows[cell]; out = []; meanpos = []
        for gi in range(3):
            rv = np.zeros(len(P)); d = []; pos = []
            for ords, dist in data:
                rk = rank(ords[gi]); d.append(dist[gi]); pos.append(rk)
                for k, (i, j) in enumerate(P): rv[k] += ((rk[i]-rk[j])*(S[i]-S[j])) < 0
            rv /= len(data); meanpos.append(np.mean(pos, axis=0))
            out.append(f"{np.mean(d):.3f}/{(rv>0.5).sum()/len(P):.3f}/{int((rv>0.5).sum()):2d}/{rv.max():.2f}")
        gap = np.abs(meanpos[0]-meanpos[2]).mean()
        print(f"  {cell:9s} {len(data):3d}  " + '   '.join(f"{o:>22s}" for o in out) + f"   {gap:.2f}")
pg = load_pergroup('runs/v2/hpc_results/dev_mechanism_fixedS_a')
sh, bad = load_shared('runs/v2/dev/shared_mixture/shared_mixture_mechanism')
mech_table(pg, 'per-group mixture (dev_mechanism_fixedS_a)'); mech_table(sh, f'SHARED mixture; non-ok={dict(bad)}')
# ---- calibration
recs = []; bad = collections.Counter()
for f in glob.glob('runs/v2/dev/shared_mixture/shared_mixture_calibration/rows_chunk*.jsonl'):
    for line in open(f):
        r = json.loads(line)
        if r['status'] == 'ok': recs.append(r)
        else: bad[r['status']] += 1
n = len(recs); print(f"\n[calibration REF_H0, shared mixture]  datasets ok={n}  non-ok={dict(bad)}  cache_consistent all={all(r['cache_consistent'] for r in recs)}")
PAIRS = ['e2-e33', 'e2-e4', 'e33-e4']
def ci(k, n):
    if n == 0: return '—'
    p = k/n; se = np.sqrt(p*(1-p)/n); return f"{p:.3f} [{max(0,p-1.96*se):.3f},{min(1,p+1.96*se):.3f}]"
res = {}
for sc in ('unrestricted', 'diagnosis'):
    pv = []; rej_any = 0; rej_max = 0; rej_pair = np.zeros(3, int)
    for r in recs:
        obs = np.array(r['observed']); T = np.array(r['perm'][sc]); B = len(T)
        cnt = (T >= obs[None, :] - 1e-12).sum(0); pv.append((cnt+1)/(B+1))
        rej_pair += cnt <= 2; rej_any += bool((cnt <= 2).any())
        cmax = (T.max(1) >= obs.max() - 1e-12).sum(); rej_max += cmax <= 9
    pv = np.array(pv); res[sc] = pv
    print(f"  {sc:12s} B={B}  mean p: " + '  '.join(f"{PAIRS[i]} {pv[:,i].mean():.3f}" for i in range(3)))
    print(f"  {'':12s} pair rej (p<=0.015): " + '  '.join(f"{PAIRS[i]} {ci(rej_pair[i], n)}" for i in range(3)))
    print(f"  {'':12s} family any-pair {ci(rej_any, n)}   max-stat (<=9) {ci(rej_max, n)}   ceilings: pair 0.015, family 0.045, max 0.050")
if n:
    d = res['diagnosis'] - res['unrestricted']
    print("  p-shift (D - U): " + '  '.join(f"{PAIRS[i]} {d[:,i].mean():+.3f} [{d[:,i].mean()-1.96*d[:,i].std(ddof=1)/np.sqrt(n):+.3f},{d[:,i].mean()+1.96*d[:,i].std(ddof=1)/np.sqrt(n):+.3f}]" for i in range(3)))
    print(f"  mean observed d_K: " + '  '.join(f"{PAIRS[i]} {np.mean([r['observed'][i] for r in recs]):.3f}" for i in range(3)) + f"   mixture {np.mean([r['mixture_seconds'] for r in recs]):.0f}s  perms {np.mean([r['perm_seconds'] for r in recs])/60:.0f} min/dataset")
