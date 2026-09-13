"""Family-wise and max rejection under the frozen rule with EXACT decisions (early-determined
rejects included), vs the 'full-budget-only' convention, per cell x engine x scheme."""
import json, glob, math, sys, collections
root = 'runs/v2/hpc_results'; run = sys.argv[1] if len(sys.argv) > 1 else 'confirm_core_a'
B, PT, MT = 599, 10, 30
def wilson(k, n, z=1.96):
    if n == 0: return (float('nan'),)*3
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return p, c-h, c+h
def dec(count, nperm, T):
    rem = B-nperm
    if count + rem <= T-1: return True
    if count >= T: return False
    return None
agg = collections.defaultdict(lambda: collections.Counter())
for f in glob.glob(f'{root}/{run}/cell_*/chunk_*/rows.jsonl'):
    for line in open(f):
        r = json.loads(line)
        if r.get('status') != 'ok': continue
        for s, v in r['schemes'].items():
            key = (r['cell'], r['engine'], s); a = agg[key]; a['rows'] += 1
            d = v['definition']; idx = d.get('tested_pair_index')
            need = [idx] if idx is not None else range(3)
            pd = [dec(v['gt'][i]+v['eq'][i], v['nperm'], PT) for i in need]
            md = dec(v['maxgt']+v['maxeq'], v['nperm'], MT) if d.get('require_max') else False
            fam_det = all(x is not None for x in pd); max_det = md is not None
            if fam_det:
                a['fam_det'] += 1; a['fam_rej'] += any(pd)
                # old convention: reject only counted at full budget
                a['fam_rej_old'] += any(pd) and v['nperm'] >= B
                a['fam_det_old'] += (v['nperm'] >= B) or all(x is False for x in pd)
            if max_det:
                a['max_det'] += 1; a['max_rej'] += bool(md)
                a['max_rej_old'] += bool(md) and v['nperm'] >= B
                a['max_det_old'] += (v['nperm'] >= B) or (md is False)
print(f"{'cell':10s}{'engine':16s}{'scheme':14s}{'rows':>5s} | exact: det  fam rate [95%]        max rate | old: det  fam   max")
for key in sorted(agg):
    a = agg[key]
    p, lo, hi = wilson(a['fam_rej'], a['fam_det']); pm = wilson(a['max_rej'], a['max_det'])[0]
    po = a['fam_rej_old']/a['fam_det_old'] if a['fam_det_old'] else float('nan')
    pmo = a['max_rej_old']/a['max_det_old'] if a['max_det_old'] else float('nan')
    print(f"{key[0]:10s}{key[1]:16s}{key[2]:14s}{a['rows']:5d} | {a['fam_det']:5d}  {p:.3f} [{lo:.3f}, {hi:.3f}]  {pm:.3f}     | {a['fam_det_old']:5d}  {po:.3f} {pmo:.3f}")
