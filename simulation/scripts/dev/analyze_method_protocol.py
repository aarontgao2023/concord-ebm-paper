"""Second confirmation protocol, mechanism part (frozen rows, seeds 53.xM).

method_mechanism_a      fixed-S IID/REF/BAL at n and 4n; engines repaired / shared / invariant_min / invariant_pooled
method_saebm_reference_a complete-data cells, engines repaired / shared / invariant_min
method_saebm_a          SA-EBM (pysaebm conjugate priors), same seeds and cells, dev-script rows
Per engine x cell: mean distance to truth, majority-ordering distance, number of majority-wrong
pairs, max pair error, per group; the mean |e2-e4| target-position gap; ESS (invariant engines)."""
import json, glob, itertools, collections, sys
import numpy as np
G = ['e2', 'e33', 'e4']
FIXED_S = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]
ROOT = 'runs/v2/hpc_results'
def rank(order):
    rk = np.empty(len(order), int); rk[np.asarray(order)] = np.arange(len(order)); return rk
def table(rows, label):
    """rows: {cell: [(orderings, distances, ess)]}"""
    S = rank(FIXED_S); P = list(itertools.combinations(range(14), 2))
    print(f"\n[{label}]")
    print(f"{'cell':12s}{'R':>5s}  " + '  '.join(f"{g}: dist/maj-d/n>0.5/max" for g in G) + "   |e2-e4| gap   ESS e2/e33/e4")
    for cell in sorted(rows, key=lambda c: (('N4' in c), c)):
        data = rows[cell]; out = []; meanpos = []
        for gi in range(3):
            rv = np.zeros(len(P)); d = []; pos = []
            for ords, dist, _ in data:
                rk = rank(ords[gi]); d.append(dist[gi]); pos.append(rk)
                for k, (i, j) in enumerate(P): rv[k] += ((rk[i]-rk[j])*(S[i]-S[j])) < 0
            rv /= len(data); meanpos.append(np.mean(pos, axis=0))
            out.append(f"{np.mean(d):.3f}/{(rv>0.5).sum()/len(P):.3f}/{int((rv>0.5).sum()):2d}/{rv.max():.2f}")
        gap = np.abs(meanpos[0]-meanpos[2]).mean()
        ess = [e for _, _, e in data if e]
        esss = '/'.join(f"{np.mean([float(e[k]) for e in ess]):.0f}" for k in ('0', '1', '2')) if ess else '-'
        print(f"{cell:12s}{len(data):5d}  " + '   '.join(f"{o:>22s}" for o in out) + f"   {gap:.2f}   {esss}")
def load_run(run):
    rows = collections.defaultdict(lambda: collections.defaultdict(list)); bad = collections.Counter()
    for f in glob.glob(f'{ROOT}/{run}/cell_*/chunk_*/rows.jsonl'):
        for line in open(f):
            r = json.loads(line)
            if r['status'] != 'ok': bad[r['status']] += 1; continue
            assert r['truth']['base_order'] == FIXED_S
            rows[r['engine']][r['cell']].append((r['observed']['orderings'], r['observed']['distance_to_truth'],
                                                 r['diagnostics'].get('effective_sample_size')))
    return rows, bad
for run in ('method_mechanism_a', 'method_saebm_reference_a'):
    rows, bad = load_run(run)
    print(f"\n===== {run}  non-ok fits: {dict(bad) or 0}")
    for eng in ('repaired', 'shared', 'invariant_pooled', 'invariant_min'):
        if eng in rows: table(rows[eng], f'{run} / {eng}')
sa = collections.defaultdict(list); bad = collections.Counter()
for f in glob.glob(f'{ROOT}/method_saebm_a/rows_chunk*.jsonl'):
    for line in open(f):
        r = json.loads(line)
        if r['status'] != 'ok': bad[r['status']] += 1; continue
        sa[r['cell']].append((r['orderings'], r['distance_to_truth'], None))
print(f"\n===== method_saebm_a  non-ok: {dict(bad) or 0}")
table(sa, 'SA-EBM conjugate priors, complete data, same seeds as method_saebm_reference_a')
