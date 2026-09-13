"""Completeness by *determined decision* under the frozen rule (B=599, le).
pair reject iff G+E<=9  -> determined-not-reject once G+E>=10
max  reject iff G+E<=29 -> determined-not-reject once G+E>=30
A scheme is determined if nperm==599, or every required endpoint is determined."""
import json, glob, sys, collections
root = sys.argv[1] if len(sys.argv) > 1 else 'runs/v2/hpc_results'
plan = {'confirm_core_a': 4000, 'confirm_stage_a': 1000, 'confirm_power_global_a': 1000,
 'confirm_pair_complete_ref_a': 1000, 'confirm_pair_power_k27_a': 1000,
 'confirm_pair_partial_e2_topup_a': 500, 'confirm_pair_partial_e4_topup_a': 500,
 'confirm_pair_partial_e33_a': 1000, 'confirm_power_global_k14_k55_a': 2000}
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
T = collections.defaultdict(collections.Counter); S = collections.defaultdict(set)
for f in glob.glob(f'{root}/confirm_*/cell_*/chunk_*/rows.jsonl'):
    run = f.split('/')[-4]
    if run not in plan: continue
    for line in open(f):
        r = json.loads(line); t = T[run]; t['rows'] += 1
        if r.get('status') != 'ok': t['rowfail'] += 1; continue
        alld = True
        for s, v in r['schemes'].items():
            S[run].add(s); dd = determined(v); t[s+':det'] += dd; alld &= dd
            t[s+':rej'] += rejects(v); t[s+':fail'] += v['failed']; t[s+':full'] += v['nperm'] >= B
        t['det'] += alld
print(f"{'run':34s}{'plan':>6s}{'rows':>6s}{'determ':>7s}{'%plan':>6s}  per-scheme determined/rows (full-B, rejects, failed fits)")
tp = td = 0
for run, P in plan.items():
    t = T[run]; tp += P; td += t['det']
    per = ' | '.join(f"{s.replace('diagnosis_pair_','pair_')}={t[s+':det']}/{t['rows']-t['rowfail']} ({t[s+':full']}f,{t[s+':rej']}r,{t[s+':fail']}x)" for s in sorted(S[run]))
    print(f"{run:34s}{P:6d}{t['rows']:6d}{t['det']:7d}{100*t['det']/P:6.1f}  {per}")
print(f"{'TOTAL':34s}{tp:6d}{'':6s}{td:7d}{100*td/tp:6.1f}")
