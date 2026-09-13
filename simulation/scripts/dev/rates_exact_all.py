"""Exact frozen-rule decisions for every run: per cell x engine x scheme, determined count,
family/pair rejection rate with Wilson 95% CI, and max-statistic rate where required."""
import json, glob, math, sys, collections
root = 'runs/v2/hpc_results'
runs = sys.argv[1:] or sorted({p.split('/')[-4] for p in glob.glob(f'{root}/*/cell_*/chunk_*/rows.jsonl')})
B, PT, MT = 599, 10, 30
def wilson(k, n, z=1.96):
    if n == 0: return float('nan'), float('nan'), float('nan')
    p = k/n; d = 1+z*z/n; c = (p+z*z/(2*n))/d; h = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return p, c-h, c+h
def dec(count, nperm, T):
    rem = B-nperm
    if count + rem <= T-1: return True
    if count >= T: return False
    return None
for run in runs:
    agg = collections.defaultdict(collections.Counter); truth = {}
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
                if all(x is not None for x in pd):
                    a['det'] += 1; a['rej'] += any(pd)
                if d.get('require_max') and md is not None:
                    a['mdet'] += 1; a['mrej'] += bool(md)
                a['full'] += v['nperm'] >= B
    print(f"== {run}")
    print(f"{'cell':14s}{'engine':16s}{'scheme':22s}{'rows':>5s}{'det':>5s}{'rej':>5s}  rate [95% CI]         | max: det rate")
    for key in sorted(agg):
        a = agg[key]; p, lo, hi = wilson(a['rej'], a['det'])
        m = f"{a['mdet']:5d} {a['mrej']/a['mdet']:.3f}" if a['mdet'] else "    -     -"
        print(f"{key[0]:14s}{key[1]:16s}{key[2]:22s}{a['rows']:5d}{a['det']:5d}{a['rej']:5d}  {p:.3f} [{lo:.3f}, {hi:.3f}] | {m}")
