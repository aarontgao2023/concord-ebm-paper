"""Offline analysis of the studentised pilot (development).

For each dataset and scheme, p-values for the plain statistic T and the studentised
statistic T* = T / (vhat_g + vhat_h), from the same permutations. Inclusive rule, +1
convention, Bonferroni alpha/3 per pair. Reports family-wise rejection with Wilson
intervals and the dataset-level mean p-shift, for the 2x2 of scheme x statistic."""
import json, glob, sys, math, collections
import numpy as np
root = sys.argv[1] if len(sys.argv) > 1 else 'runs/v2/dev/studentized_pilot'
ALPHA = 0.05; PAIRS = ['e2-e33', 'e2-e4', 'e33-e4']; PIDX = [(0, 1), (0, 2), (1, 2)]
rows = collections.defaultdict(dict)   # (seed, scheme) -> pid -> rec
obs = {}
for f in glob.glob(f'{root}/rows_chunk*.jsonl'):
    for line in open(f):
        r = json.loads(line)
        if r.get('status') != 'ok': continue
        if r['scheme'] is None: obs[r['seed']] = r
        else: rows[(r['seed'], r['scheme'])][r['pid']] = r
def stats(r):
    t = np.array(r['taus']); v = np.array(r['vhat'])
    ts = np.array([t[k] / (v[i] + v[j]) for k, (i, j) in enumerate(PIDX)])
    return t, ts
def wilson(k, n, z=1.96):
    if n == 0: return (float('nan'),) * 2
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)
out = collections.defaultdict(lambda: dict(n=0, rej=0, p=[], B=[]))
pshift = collections.defaultdict(list)
for seed, o in obs.items():
    to, tso = stats(o)
    pv = {}
    for scheme in ('unrestricted', 'diagnosis'):
        perms = rows.get((seed, scheme), {})
        if not perms: continue
        B = len(perms)
        T = np.array([stats(p)[0] for p in perms.values()]); TS = np.array([stats(p)[1] for p in perms.values()])
        for name, mat, o_ in (('plain', T, to), ('studentised', TS, tso)):
            p = (1 + (mat >= o_ - 1e-12).sum(0)) / (B + 1)      # ties count as exceedances
            rej = bool((p <= ALPHA / 3).any())
            k = (scheme, name); out[k]['n'] += 1; out[k]['rej'] += rej; out[k]['p'].append(p); out[k]['B'].append(B)
            pv[k] = p
    for name in ('plain', 'studentised'):
        if ('unrestricted', name) in pv and ('diagnosis', name) in pv:
            pshift[name].append(float(np.mean(pv[('diagnosis', name)] - pv[('unrestricted', name)])))
    if ('unrestricted', 'plain') in pv and ('unrestricted', 'studentised') in pv:
        pshift['studentised - plain, unrestricted'].append(float(np.mean(pv[('unrestricted', 'studentised')] - pv[('unrestricted', 'plain')])))
print(f"datasets with observed fit: {len(obs)}")
print(f"{'scheme':14s}{'statistic':13s}{'n':>5s}{'FWER':>8s}{'95% Wilson':>18s}{'mean p':>9s}{'B range':>10s}")
for (scheme, name), d in sorted(out.items()):
    lo, hi = wilson(d['rej'], d['n']); mp = float(np.mean(np.concatenate(d['p']))) if d['p'] else float('nan')
    print(f"{scheme:14s}{name:13s}{d['n']:5d}{d['rej']/d['n']:8.3f}   [{lo:.3f}, {hi:.3f}]{mp:9.3f}   {min(d['B'])}-{max(d['B'])}")
print("\nmean paired p-shift (dataset-level), with bootstrap 95% CI:")
rng = np.random.default_rng(0)
for name, vals in pshift.items():
    v = np.array(vals)
    if len(v) < 3: print(f"  {name}: n={len(v)} {v}"); continue
    bs = [rng.choice(v, len(v)).mean() for _ in range(2000)]
    print(f"  {name:38s} n={len(v):4d}  {v.mean():+.4f}  [{np.percentile(bs,2.5):+.4f}, {np.percentile(bs,97.5):+.4f}]")
# observed scale by group
if obs:
    V = np.array([o['vhat'] for o in obs.values()])
    print(f"\nobserved vhat by group (mean over datasets): e2 {V[:,0].mean():.3f}  e33 {V[:,1].mean():.3f}  e4 {V[:,2].mean():.3f}")
