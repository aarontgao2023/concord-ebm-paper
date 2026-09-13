"""Mechanism: precision or target?  Fixed true ordering S; compositions IID / REF / BALANCED;
sample scale 1 and 4. For each cell and group: distance to S, per-pair reversal probability
P(pair (i,j) fitted in the opposite order to S), the pairwise-majority (Condorcet) ordering
and its distance to S, and how many of the 91 pairs are reversed in the majority.
If reversal probabilities shrink toward 0 with N (and majority ordering -> S), the effect is
precision; pairs whose reversal probability stays high or exceeds 0.5 at 4N indicate a
composition-dependent estimation target."""
import json, glob, sys, itertools, collections
import numpy as np
root = sys.argv[1] if len(sys.argv) > 1 else 'runs/v2/hpc_results/dev_mechanism_fixedS_a'
G = ['e2', 'e33', 'e4']
rows = collections.defaultdict(list)   # cell -> list of (orderings, truth S)
for f in glob.glob(f'{root}/cell_*/chunk_*/rows.jsonl'):
    for line in open(f):
        r = json.loads(line)
        if r['status'] != 'ok': continue
        rows[r['cell']].append((r['observed']['orderings'], r['truth']['base_order'], r['observed']['distance_to_truth']))
def rank(order):  # order: list of biomarker ids by position -> rank array
    rk = np.empty(len(order), int); rk[np.asarray(order)] = np.arange(len(order)); return rk
def pairs(n): return list(itertools.combinations(range(n), 2))
order = ['IID_S', 'REF_S', 'BAL_S', 'IID_S_N4', 'REF_S_N4', 'BAL_S_N4']
print(f"{'cell':10s}{'R':>4s}  " + '  '.join(f"{g:>22s}" for g in G))
print(f"{'':16s}" + '  '.join(f"{'dist':>6s}{'maj-d':>6s}{'n>0.5':>5s}{'max':>5s}" for g in G))
rev = {}
for cell in order:
    if cell not in rows: continue
    data = rows[cell]; n = len(data[0][1]); P = pairs(n)
    out = []
    for gi, g in enumerate(G):
        S = rank(data[0][1]); rv = np.zeros(len(P)); d = []
        for ords, base, dist in data:
            rk = rank(ords[gi]); d.append(dist[gi])
            for k, (i, j) in enumerate(P):
                rv[k] += ((rk[i] - rk[j]) * (S[i] - S[j])) < 0
        rv /= len(data); rev[(cell, g)] = rv
        maj_reversed = int((rv > 0.5).sum())
        out.append(f"{np.mean(d):6.3f}{maj_reversed/len(P):6.3f}{maj_reversed:5d}{rv.max():5.2f}")
    print(f"{cell:10s}{len(data):4d}  " + '  '.join(out))
print("\ndist = mean normalised Kendall distance of fitted ordering to S; maj-d = distance of pairwise-majority ordering to S;")
print("n>0.5 = pairs reversed in the majority (a composition-dependent target if >0 and stable with N); max = largest per-pair reversal probability.")
# does the reversal profile differ by composition, beyond scale?  correlation of per-pair reversal vectors
print("\nper-pair reversal-probability profiles: correlation across compositions (same group, same N)")
for N in ('', '_N4'):
    for g in G:
        a, b, c = (rev.get((f'IID_S{N}', g)), rev.get((f'REF_S{N}', g)), rev.get((f'BAL_S{N}', g)))
        if a is None or b is None or c is None: continue
        print(f"  N{'4' if N else '1'} {g:4s}  corr(IID,REF)={np.corrcoef(a,b)[0,1]:.2f}  corr(IID,BAL)={np.corrcoef(a,c)[0,1]:.2f}  corr(REF,BAL)={np.corrcoef(b,c)[0,1]:.2f}   mean|REF-IID|={np.abs(b-a).mean():.3f}")
# persistence: pairs with reversal prob >= 0.4 at N4, by cell/group
print("\npairs with reversal probability >= 0.40 at 4N (persistent reversals):")
for cell in ('IID_S_N4', 'REF_S_N4', 'BAL_S_N4'):
    for g in G:
        rv = rev.get((cell, g))
        if rv is None: continue
        n = len(rows[cell][0][1]); P = pairs(n); base = rows[cell][0][1]; S = rank(base)
        hits = [(i, j, round(rv[k], 2)) for k, (i, j) in enumerate(P) if rv[k] >= 0.40]
        print(f"  {cell:10s} {g:4s} {len(hits):2d} pairs: {hits[:8]}{' ...' if len(hits) > 8 else ''}")
