"""Review response, part B: analyses that need no new simulation.

B2  per-cohort gap between the e2 and e4 orderings (mean over cohorts of the per-cohort mean absolute
    position difference) next to the population gap (difference of mean positions) reported so far.
B3  D-all under partial nulls: how often the pair that does NOT involve the displaced group is rejected,
    and decomposition of family rejections into 'a truly different pair rejected' vs 'only the innocent pair'.
B4  per-pair D-all rejection at the complete null against the exact per-pair bound 1/60.
Exact decisions under the frozen rule (B=599, le)."""
import json, glob, collections, math
import numpy as np
ROOT = 'runs/v2/hpc_results'; B, PT = 599, 10
PAIR = ['e2-e33', 'e2-e4', 'e33-e4']
def rows(run):
    for f in glob.glob(f'{ROOT}/{run}/cell_*/chunk_*/rows.jsonl'):
        for line in open(f):
            r = json.loads(line)
            if r.get('status') == 'ok': yield r
def dec(count, nperm, T=PT):
    rem = B - nperm
    return True if count + rem <= T - 1 else False if count >= T else None
def wilson(k, n, z=1.96):
    if not n: return float('nan'), float('nan'), float('nan')
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return p, c-h, c+h
def rank(o):
    rk = np.empty(len(o), int); rk[np.asarray(o)] = np.arange(len(o)); return rk

print("== B2  gap definitions, method_mechanism_a (R=1000 per cell)")
print(f"{'engine':16s}{'cell':10s}  population gap  per-cohort gap  per-cohort Kendall d(e2,e4)")
acc = collections.defaultdict(list)
for r in rows('method_mechanism_a'):
    o = r['observed']['orderings']; acc[(r['engine'], r['cell'])].append((rank(o[0]), rank(o[2]), r['observed']['taus'][1]))
for (eng, cell) in sorted(acc):
    if eng not in ('repaired', 'invariant_min', 'shared', 'invariant_pooled'): continue
    R2 = np.array([a for a, _, _ in acc[(eng, cell)]]); R4 = np.array([b for _, b, _ in acc[(eng, cell)]])
    pop = np.abs(R2.mean(0) - R4.mean(0)).mean(); per = np.abs(R2 - R4).mean(1).mean(); kd = np.mean([t for _, _, t in acc[(eng, cell)]])
    print(f"{eng:16s}{cell:10s}  {pop:14.2f}  {per:14.2f}  {kd:12.3f}")

print("\n== B3  D-all under partial nulls: innocent pair and decomposition (exact decisions)")
innocent = {'PWR_E2_K14': 2, 'PWR_E2_K27': 2, 'PWR_E2_K55': 2, 'PWR_E4_K14': 0, 'PWR_E4_K27': 0, 'PWR_E4_K55': 0}
print(f"{'run':32s}{'cell':12s}{'engine':14s}{'det':>5s}  innocent-pair rej [95%]   family rej   true-pair rej   only-innocent rej")
for run in ('confirm_power_global_a', 'confirm_power_global_k14_k55_a', 'method_power_a'):
    T = collections.defaultdict(collections.Counter)
    for r in rows(run):
        v = r['schemes'].get('diagnosis'); cell = r['cell']
        if v is None or cell not in innocent: continue
        d = [dec(v['gt'][i] + v['eq'][i], v['nperm']) for i in range(3)]
        if any(x is None for x in d): continue
        t = T[(cell, r['engine'])]; t['det'] += 1; inn = innocent[cell]; true_pairs = [i for i in range(3) if i != inn]
        t['inn'] += d[inn]; t['fam'] += any(d); t['true'] += any(d[i] for i in true_pairs); t['only_inn'] += d[inn] and not any(d[i] for i in true_pairs)
    for (cell, eng), t in sorted(T.items()):
        p, lo, hi = wilson(t['inn'], t['det'])
        print(f"{run:32s}{cell:12s}{eng:14s}{t['det']:5d}  {p:.3f} [{lo:.3f}, {hi:.3f}]   {t['fam']/t['det']:.3f}        {t['true']/t['det']:.3f}          {t['only_inn']/t['det']:.3f}")

print("\n== B4  per-pair D-all rejection at the complete null (REF), bound 1/60 = 0.0167")
for run, engines in (('confirm_core_a', ('original', 'repaired')), ('method_calibration_a', ('shared', 'invariant_min')), ('method_calibration_pooled_a', ('invariant_pooled',))):
    T = collections.defaultdict(lambda: np.zeros(3, int)); N = collections.Counter()
    for r in rows(run):
        if r['cell'] != 'REF_H0': continue
        v = r['schemes']['diagnosis']; d = [dec(v['gt'][i] + v['eq'][i], v['nperm']) for i in range(3)]
        if any(x is None for x in d): continue
        T[r['engine']] += np.array(d, int); N[r['engine']] += 1
    for eng in engines:
        print(f"{run:28s}{eng:16s} n={N[eng]:5d}  " + "  ".join(f"{PAIR[i]} {T[eng][i]/N[eng]:.4f}" for i in range(3)))
