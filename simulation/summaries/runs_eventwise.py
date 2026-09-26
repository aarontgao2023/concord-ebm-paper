"""Event-position shifts in the random-sequence simulation runs (Fig. 2b); no fitting.

For each run, cell and model, reads the fitted orderings stored in the run_v2.py records
(RUNS_DIR/<run>/cell_*/chunk_*/rows.jsonl) and writes, for each simulated biomarker, the event-position
shift (position in the CN-heavy group g1 minus position in the AD-heavy group g3, mean over datasets) with
its standard error sd/sqrt(R). Positive values mean that the event is placed later in the CN-heavy group.
Also written: the mean normalized Kendall distance between the g1 and g3 orderings and the mean errors
against the true sequence in g1 and g3. The biomarker column holds the code names of design_v2.BIOMARKERS,
in input order = biomarkers 1-14 of the paper (1-5 modeled on CSF markers, 6-7 on cognitive scores,
8-14 on imaging measures). A summary line per run, cell and model is printed, including the Pearson
correlation of the 14 shifts with the component AUCs.

Fig. 2b uses confirm_core_a (REF_H0, repaired), method_calibration_a (REF_H0, shared and invariant_min)
and method_kde_calibration_a (REF_H0_CD, kde_gmm). The default runs are those of the published
runs_random_truth_eventwise.csv, which also include two runs that the paper does not use
(method_calibration_pooled_a, method_kde_precision_a; configurations not included). Runs without
records are skipped.

Usage:
  python simulation/summaries/runs_eventwise.py --runs-dir RUNS_DIR --out runs_random_truth_eventwise.csv \
      [--datahash benchmark_gaussian_null_datahash.json]
"""
import argparse, json, os, math, csv, glob
from collections import defaultdict

PAPER_RUNS=['confirm_core_a','method_calibration_a','method_kde_calibration_a']
RUNS=['confirm_core_a','method_calibration_a','method_calibration_pooled_a','method_kde_calibration_a','method_kde_precision_a']
AUC={'ABETA':0.88,'PTAU':0.85,'TAU':0.80,'NG':0.68,'NFL':0.72,'ADAS13':0.92,'MMSE':0.86,'Hippocampus':0.82,'Entorhinal':0.80,'MidTemp':0.78,'Fusiform':0.76,'WholeBrain':0.74,'Ventricles':0.72,'Precuneus':0.70}
def pearson(x,y):
    mx=sum(x)/len(x); my=sum(y)/len(y)
    sxy=sum((a-mx)*(b-my) for a,b in zip(x,y)); sx=math.sqrt(sum((a-mx)**2 for a in x)); sy=math.sqrt(sum((b-my)**2 for b in y))
    return sxy/(sx*sy)

def main():
    ap=argparse.ArgumentParser(description='Event-position shifts in the random-sequence simulation runs (Fig. 2b).')
    ap.add_argument('--runs-dir',default=os.environ.get('CONCORD_RUNS_DIR'),
                    help='directory with one subdirectory of records per run (default: $CONCORD_RUNS_DIR)')
    ap.add_argument('--out',required=True,help='output CSV (runs_random_truth_eventwise.csv)')
    ap.add_argument('--runs',nargs='+',default=RUNS,
                    help='runs to summarize (default: the runs of the published file; Fig. 2b uses '+', '.join(PAPER_RUNS)+')')
    ap.add_argument('--datahash',default=None,
                    help='optional JSON mapping "GAUSSIAN_H0|<seed>" to the data SHA-256 of the benchmark datasets; '
                         'if given, checks that method_calibration_a simulated the same datasets')
    args=ap.parse_args()
    if args.runs_dir is None:
        ap.error('give --runs-dir or set CONCORD_RUNS_DIR')
    BASE=args.runs_dir
    datahash=json.load(open(args.datahash)) if args.datahash else None
    os.makedirs(os.path.dirname(os.path.abspath(args.out)),exist_ok=True)
    out=open(args.out,'w')
    w=csv.writer(out)
    w.writerow(['run_id','cell','engine','R_ok','R_rows','biomarker','mean_pos_g1_minus_g3','se','mean_tau_g1g3','err_g1','err_g3'])
    for run in args.runs:
        acc=defaultdict(lambda: defaultdict(list)); cnt=defaultdict(int); names=None; seen=set(); hashmatch=[0,0]
        files=sorted(glob.glob(os.path.join(BASE,run,'cell_*','chunk_*','rows.jsonl')))
        if not files:
            print('%s: no records in %s; skipped'%(run,BASE)); continue
        for f in files:
            for line in open(f):
                r=json.loads(line); k=(r['cell'],r['engine']); cnt[k]+=1
                if (k,r['seed']) in seen: continue
                seen.add((k,r['seed']))
                if r['status']!='ok': continue
                t=r['truth']; names=t['biomarker_names']; tp={b:i+1 for i,b in enumerate(t['base_order'])}
                if run=='method_calibration_a' and datahash is not None:
                    h=datahash.get('GAUSSIAN_H0|%d'%r['seed'])
                    if h is not None:
                        hashmatch[1]+=1; hashmatch[0]+= (h==t['manifest'].get('data_sha256'))
                o=r['observed']; pos=[{b:i+1 for i,b in enumerate(x)} for x in o['orderings']]
                a=acc[k]
                a['tau'].append(o['taus'][1]); a['e1'].append(o['distance_to_truth'][0]); a['e3'].append(o['distance_to_truth'][2])
                for b in range(len(names)): a['c_%d'%b].append(pos[0][b]-pos[2][b])
        if run=='method_calibration_a' and datahash is not None: print('method_calibration_a data_sha256 matches the benchmark GAUSSIAN_H0 datasets: %d/%d'%tuple(hashmatch))
        for k,a in sorted(acc.items()):
            R=len(a['tau']); ms=[]
            for b in range(len(names)):
                x=a['c_%d'%b]; m=sum(x)/R; se=math.sqrt(sum((v-m)**2 for v in x)/(R-1))/math.sqrt(R); ms.append((m,se))
                w.writerow([run,k[0],k[1],R,cnt[k],names[b],round(m,4),round(se,4),round(sum(a['tau'])/R,4),round(sum(a['e1'])/R,4),round(sum(a['e3'])/R,4)])
            gap=sum(abs(m) for m,_ in ms)/len(ms); deb=sum(math.sqrt(max(m*m-s*s,0)) for m,s in ms)/len(ms)
            print('%-28s %-12s %-16s R=%4d (rows %d) d13=%.3f err1=%.3f err3=%.3f eventwise gap=%.3f debiased=%.3f r(AUC)=%+.2f'%(run,k[0],k[1],R,cnt[k],sum(a['tau'])/R,sum(a['e1'])/R,sum(a['e3'])/R,gap,deb,pearson([m for m,_ in ms],[AUC[n] for n in names])))
            print('     ',' '.join('%s:%+.2f'%(n[:6],m) for n,(m,_) in zip(names,ms)))
    out.close()

if __name__=='__main__':
    main()
