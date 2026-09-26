"""No fits: verify every paired historical dataset and frozen permutation stream."""
from pathlib import Path
import argparse,ast,hashlib,importlib.util,json,sys,time,os
from concurrent.futures import ProcessPoolExecutor
import numpy as np
HERE=Path(__file__).resolve().parent;SIM_SCRIPTS=HERE.parents[1]/'simulation/scripts'
# Working directory for outputs; CONCORD_HISTORICAL_ROOT holds <config>/cell_*/chunk_*/rows.jsonl written by
# simulation/scripts/run_v2.py for confirm_power_global_a, method_alt_geometry_a and method_calibration_a.
WORK=Path(os.environ.get('CONCORD_WORK_DIR',str(HERE.parent/'results')))
OUT=Path(os.environ.get('CONCORD_REUSE_OUT',str(WORK/'reuse_audit')))
CELLS=[('confirm_power_global_a','PWR_E2_K27',['repaired']),('confirm_power_global_a','PWR_E4_K27',['repaired']),('method_alt_geometry_a','K13_SINGLE_E4',['repaired','invariant_min']),('method_alt_geometry_a','K6_DISJOINT_E2',['repaired','invariant_min']),('method_calibration_a','REF_H0',['shared','invariant_min'])]

def sha(x):return hashlib.sha256(x).hexdigest()
def csha(x):return sha(json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode())
def load_design(tag):
 path=HERE/'vendor'/tag/'design_v2.py';name='frozen_'+tag
 if name in sys.modules:return sys.modules[name]
 sp=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(sp);sys.modules[name]=mod;sp.loader.exec_module(mod);return mod

def perm_source(path):
 tree=ast.parse(path.read_text());f=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='permute')
 class DiagnosisPath(ast.NodeTransformer):
  def visit_If(self,node):
   self.generic_visit(node)
   # The only historical extension is a branch never selected by within-diagnosis permutation ('diagnosis').
   if ast.unparse(node.test)=="stratify == 'diagnosis_observed_count'":
    return node.orelse
   return node
 f=DiagnosisPath().visit(f)
 return ast.dump(f,include_attributes=False)

def permutation_stream(df,seed,budget=599):
 # Identical old permutation algorithm; avoids DataFrame copies for hashing.
 labels=df.APOE.to_numpy();strata=df.Diagnosis.to_numpy();ix=[np.flatnonzero(strata==d) for d in np.unique(strata)]
 tag=int.from_bytes(hashlib.sha256(b'diagnosis').digest()[:4],'little');h=hashlib.sha256()
 for pid in range(budget):
  rng=np.random.default_rng(np.random.SeedSequence([seed,tag,pid,61000000]));p=labels.copy()
  for indexes in ix:p[indexes]=rng.permutation(p[indexes])
  h.update(np.asarray(p,dtype='<i8').tobytes())
 return h.hexdigest()

def work(task):
 run,cell,seed,records=task;d=load_design(run);saved=records[0]['truth'];cfg=d.DesignConfig.from_dict(saved['config']);df,truth=d.simulate(cfg,seed)
 # Historical runner added these deterministic fields after generation.
 distances=[truth['pair_truth'][p]['normalized_distance'] for p in ('e2_vs_e33','e2_vs_e4','e33_vs_e4')]
 truth.update(pair_null=[x==0 for x in distances],pair_distances=distances,h1_distance=max(distances))
 checks={'data_sha':truth['manifest']['data_sha256']==saved['manifest']['data_sha256'],'config_sha':truth['manifest']['config_sha256']==saved['manifest']['config_sha256'],'generator_sha':truth['manifest']['source_sha256']==saved['manifest']['source_sha256'],'truth_sha':csha(truth)==csha(saved),'all_saved_arm_truth_equal':all(csha(x['truth'])==csha(saved) for x in records)}
 # File path copied under vendor preserves __file__ bytes, hence generator SHA.
 if not all(checks.values()):raise AssertionError((run,cell,seed,checks))
 return {'run':run,'cell':cell,'seed':seed,'data_sha256':truth['manifest']['data_sha256'],'config_sha256':truth['manifest']['config_sha256'],'truth_sha256':csha(truth),'permutation_label_stream_sha256':permutation_stream(df,seed),'checks':checks,'saved_engines':{r['engine']:{'status':r['status'],'computed_permutations':r['schemes']['diagnosis']['nperm'],'failed':r['schemes']['diagnosis']['failed'],'pair_reject':r['schemes']['diagnosis']['exact_decision']['le']['pair_reject']} for r in records}}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=4);ap.add_argument('--limit',type=int);ap.add_argument('--only-run');args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 tasks=[];source_checks=[];historical=[]
 for run,cell,engines in CELLS:
  if args.only_run and run!=args.only_run:continue
  p=Path(os.environ.get('CONCORD_HISTORICAL_ROOT',str(WORK/'runs')))/run;by={};expected_seeds=range(42200000,42200500) if run.startswith('confirm') else (range(53100000,53101000) if run=='method_calibration_a' else range(53800000,53800500))
  files=list(p.glob('cell_*/chunk_*/rows.jsonl'))
  for f in files:
   manifest=json.loads(f.with_name('manifest.json').read_text());v=HERE/'vendor'/run
   for name in ('design_v2.py','run_v2.py'):
    assert manifest['source_sha256'][name]==sha((v/name).read_bytes()),(f,name)
   assert manifest['config']['bperm']==599 and manifest['config']['rule']=='le'
   for line in f.open():
    r=json.loads(line)
    if r['cell']!=cell or r['engine']not in engines:continue
    k=(r['seed'],r['engine']);assert k not in by,k;by[k]=r
  for seed in expected_seeds:
   records=[by[(seed,e)] for e in engines]
   for r in records:
    s=r['schemes']['diagnosis'];assert r['status']=='ok' and all(x is not None for x in s['exact_decision']['le']['pair_reject']);assert s['requested_nperm']==599
   historical.extend(records);tasks.append((run,cell,seed,records))
  source_checks.append({'run':run,'cell':cell,'rows':len(by),'permutation_function_same_as_current':perm_source(v/'run_v2.py')==perm_source(SIM_SCRIPTS/'run_v2.py')})
 assert all(x['permutation_function_same_as_current'] for x in source_checks)
 if args.limit is not None:tasks=tasks[:args.limit]
 started=time.monotonic();dest=OUT/('reuse_audit_pretest.jsonl' if args.limit else 'reuse_audit.jsonl')
 with dest.open('w') as out,ProcessPoolExecutor(max_workers=args.workers) as pool:
  for i,r in enumerate(pool.map(work,tasks,chunksize=4),1):
   out.write(json.dumps(r,separators=(',',':'))+'\n')
   if i%100==0:out.flush();print(i,'datasets verified',round(time.monotonic()-started,1),flush=True)
 summary={'status':'PASS','formal_reuse_gate':not bool(args.limit or args.only_run),'datasets':len(tasks),'planned_datasets':3000,'source_checks':source_checks,'elapsed_seconds':time.monotonic()-started,'method':'Regenerate each full dataset/truth from its exact frozen source. Hash all599 deterministic within-diagnosis label arrays; frozen permutation function AST equal after excluding the unused diagnosis_observed_count extension. Historical labels were not separately archived, so the label stream hash is a reconstructed provenance identity, not a comparison to a preexisting historical hash.','output_sha256':sha(dest.read_bytes())}
 (OUT/('reuse_pretest.json' if args.limit else 'reuse_summary.json')).write_text(json.dumps(summary,indent=2)+'\n')
 if not args.limit:
  with (OUT/'historical_reused_rows.jsonl').open('w') as f:
   for r in historical:f.write(json.dumps(r,separators=(',',':'))+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
