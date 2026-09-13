"""Completeness of the second protocol under the frozen rule (same logic as tally.py)."""
import json, glob, collections
root = 'runs/v2/hpc_results'
plan = {'method_mechanism_a': 24000, 'method_saebm_reference_a': 12000, 'method_calibration_a': 2000,
        'method_power_a': 3000, 'method_boundary_a': 1500}
B, PAIR_T, MAX_T = 599, 10, 30
def determined(v):
    """Exact decision under the frozen rule over every completion of the planned budget:
    reject determined if G+E+(B-nperm) <= T-1; not-reject determined if G+E >= T."""
    d = v['definition']; rem = B - v['nperm']
    ge = [g+e for g,e in zip(v['gt'], v['eq'])]
    idx = d.get('tested_pair_index'); need = [idx] if idx is not None else range(len(ge))
    pairs_ok = all(ge[i] >= PAIR_T or ge[i] + rem <= PAIR_T-1 for i in need)
    m = v['maxgt']+v['maxeq']
    max_ok = (not d.get('require_max')) or (m >= MAX_T or m + rem <= MAX_T-1)
    return pairs_ok and max_ok
def rejects(v):
    d = v['definition']; rem = B - v['nperm']; ge=[g+e for g,e in zip(v['gt'], v['eq'])]; idx=d.get('tested_pair_index')
    need = [idx] if idx is not None else range(len(ge))
    return any(ge[i] + rem <= PAIR_T-1 for i in need)
print(f"{'run':28s}{'plan':>6s}{'rows':>6s}{'ok':>6s}{'determ':>7s}{'%':>6s}  per scheme determined (full-B, rejects) | by engine determined")
tp = td = 0
for run, P in plan.items():
    T = collections.Counter(); E = collections.Counter(); En = collections.Counter(); S = set()
    for f in glob.glob(f'{root}/{run}/cell_*/chunk_*/rows.jsonl'):
        for line in open(f):
            r = json.loads(line); T['rows'] += 1
            if r.get('status') != 'ok': T['fail'] += 1; continue
            T['ok'] += 1; En[r['engine']] += 1
            alld = True
            for s, v in r['schemes'].items():
                S.add(s); dd = determined(v); T[s+':det'] += dd; alld &= dd
                T[s+':full'] += v['nperm'] >= B; T[s+':rej'] += rejects(v)
            T['det'] += alld; E[r['engine']] += alld
    tp += P; td += T['det']
    per = ' | '.join(f"{s}={T[s+':det']} ({T[s+':full']}f,{T[s+':rej']}r)" for s in sorted(S))
    eng = ', '.join(f"{e}={E[e]}/{En[e]}" for e in sorted(En))
    print(f"{run:28s}{P:6d}{T['rows']:6d}{T['ok']:6d}{T['det']:7d}{100*T['det']/P:6.1f}  {per} | {eng}")
print(f"{'TOTAL':28s}{tp:6d}{'':12s}{td:7d}{100*td/tp:6.1f}")
n = sum(1 for f in glob.glob(f'{root}/method_saebm_a/rows_chunk*.jsonl') for l in open(f) if '"status": "ok"' in l)
print(f"method_saebm_a (dev script): {n}/4000 ok rows")
