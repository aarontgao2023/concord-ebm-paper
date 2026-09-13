"""Complete-data (missing=False) estimator comparison: per-group DEBM (plain), shared mixture +
composition-standardised consensus (ours), and SA-EBM conjugate_priors (Hao et al. 2025).
Same fixed S, same seeds, same simulated cohorts. Table per estimator: dist / maj-d / pairs>0.5 /
max per group, and the mean |e2-e4| target-position gap."""
import json, glob, itertools, collections
import numpy as np
G = ['e2', 'e33', 'e4']; CELLS = ['IID_S', 'REF_S', 'BAL_S', 'IID_S_N4', 'REF_S_N4', 'BAL_S_N4']
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]
def rank(order): rk = np.empty(len(order), int); rk[np.asarray(order)] = np.arange(len(order)); return rk
def load(pattern):
    rows = collections.defaultdict(list); bad = collections.Counter()
    for f in glob.glob(pattern):
        for line in open(f):
            r = json.loads(line)
            if r['status'] == 'ok': rows[(r['variant'], r['cell'])].append(r)
            else: bad[(r['variant'], r['status'])] += 1
    return rows, bad
def table(rows, variant, label):
    S = rank(FIXED_S); P = list(itertools.combinations(range(14), 2))
    print(f"\n[{label}]  cell  R   " + '  '.join(f"{g}: dist/maj-d/n>0.5/max" for g in G) + "   |e2-e4| gap   fit s")
    for cell in CELLS:
        data = rows.get((variant, cell), [])
        if not data: continue
        out = []; meanpos = []
        for gi in range(3):
            rv = np.zeros(len(P)); d = []; pos = []
            for r in data:
                rk = rank(r['orderings'][gi]); d.append(r['distance_to_truth'][gi]); pos.append(rk)
                for k, (i, j) in enumerate(P): rv[k] += ((rk[i]-rk[j])*(S[i]-S[j])) < 0
            rv /= len(data); meanpos.append(np.mean(pos, axis=0))
            out.append(f"{np.mean(d):.3f}/{(rv>0.5).sum()/len(P):.3f}/{int((rv>0.5).sum()):2d}/{rv.max():.2f}")
        gap = np.abs(meanpos[0]-meanpos[2]).mean()
        secs = np.mean([np.mean(r['fit_seconds']) if isinstance(r['fit_seconds'], list) else r['fit_seconds'] for r in data])
        print(f"  {cell:9s} {len(data):3d}  " + '   '.join(f"{o:>22s}" for o in out) + f"   {gap:.2f}   {secs:.0f}")
a, ba = load('runs/v2/dev/stdcons_nomiss/stdcons_nomiss_mechanism/rows_chunk*.jsonl')
b, bb = load('runs/v2/dev/saebm/saebm_mechanism/rows_chunk*.jsonl')
rows = {**a, **b}
print('non-ok:', dict(ba), dict(bb))
table(rows, 'plain', 'complete data: per-group DEBM (current practice)')
table(rows, 'shared_std_pooled', 'complete data: shared mixture + standardised consensus (ours)')
table(rows, 'saebm_conjugate', 'complete data: SA-EBM conjugate_priors (Hao et al. 2025), per group')
