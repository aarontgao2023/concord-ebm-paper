"""Append-only single-dataset runner; formal campaigns require external gate."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,hashlib,json,sys,time,signal,fcntl,traceback
import numpy as np
HERE=Path(__file__).resolve().parent
# Engine and runner modules (engine_v2, invariant_engine_v2, run_v2, ...) are loaded from
# simulation/scripts of this repository; the benchmark ran with byte-identical copies (see README).
sys.path.insert(0,str(HERE.parents[1]/'simulation/scripts'))
from engine_v2 import EngineConfig,fit_orderings
from invariant_engine_v2 import fit_invariant_orderings,_CACHE,_matrix_key
from fast_likelihood_v2 import fast_likelihood_context
from dependency_v2 import verify_pyebm
from run_v2 import pair_taus,kendall,permute,atomic_json,FitTimeout,timeout_handler,decision_bounds,compact_diagnostics
from robust_design import simulate_robust,canonical_sha
from audit_reuse import load_design,permutation_stream
from allgroup_engine import allgroup_stopping
_CACHE_PROVENANCE={}

def historical_index():
 p=HERE/'inputs/historical_reused_rows.jsonl'
 return {(r['run_id'],r['cell'],r['seed'],r['engine']):r for r in map(json.loads,p.open())}

def generate(campaign,spec,seed):
 if campaign=='robustness':return simulate_robust(spec['mechanism'],spec['cell'],seed),None
 rows=historical_index();r=rows[(spec['historical_run'],spec['historical_cell'],seed,spec.get('historical_anchor_engine','repaired'))];d=load_design(spec['historical_run']);df,t=d.simulate(d.DesignConfig.from_dict(r['truth']['config']),seed)
 assert t['manifest']['data_sha256']==r['truth']['manifest']['data_sha256']
 return (df,t),{a:rows.get((spec['historical_run'],spec['historical_cell'],seed,a)) for a in ('repaired','shared','invariant_min')}

def save_cache(df,names,path):
 c=_CACHE.get(_matrix_key(df,names))
 if c is None:return
 payload={f'params_{i}':x for i,x in enumerate(c['params'])}|{f'post_{i}':x for i,x in enumerate(c['post'])}
 tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,**payload);tmp.replace(path)

def arm_state(records,arm,budget):
 observed=records.get((arm,None))
 if observed is None or observed.get('status')!='ok':return None
 good=[r for (e,pid),r in records.items() if e==arm and pid is not None and r.get('status')=='ok']
 o=np.asarray(observed['taus']);gt=np.zeros(3,dtype=int);eq=np.zeros(3,dtype=int)
 for r in good:
  v=np.asarray(r['taus']);gt+=v>o+1e-9;eq+=np.abs(v-o)<=1e-9
 decisions=decision_bounds(gt,eq,len(good),budget,.05)['le']['pair_reject'];e=gt+eq
 return {'computed_permutations':len(good),'exceedances':e.tolist(),'pair_reject':decisions,'all_decisions_resolved':all(v is not None for v in decisions),'adjusted_p_lower':np.minimum(1,3*(1+e)/(budget+1)).tolist(),'adjusted_p_upper':np.minimum(1,3*(1+e+budget-len(good))/(budget+1)).tolist()}

def run(campaign,index,seed_index,out,benchmark_b=None,all_arms=False):
 cfg=json.loads((HERE/'campaign.json').read_text());spec=cfg[campaign][index];seed=spec['base_seed']+seed_index
 if benchmark_b is None:
  gate=json.loads((HERE/'EXECUTION_GATE.json').read_text())
  assert gate.get('formal_launch_approved') is True and gate.get('engine')=='allgroup_original_tolerance'
  assert gate['campaign_sha256']==hashlib.sha256((HERE/'campaign.json').read_bytes()).hexdigest()
  assert gate['inputs_sha256']==hashlib.sha256((HERE/'inputs/historical_reused_rows.jsonl').read_bytes()).hexdigest()
  assert 0<=seed_index<spec['datasets']
  budget=cfg['bperm'];arms=spec['new_arms']
 else:
  assert 0<=benchmark_b<=8;budget=benchmark_b;arms=['repaired','shared','invariant_min'] if all_arms else spec['new_arms']
  if campaign=='robustness':seed=92310000+index*1000+seed_index
 safe_stopping=benchmark_b is None
 out=Path(out)/campaign/spec['id']/f'seed_{seed}';out.mkdir(parents=True,exist_ok=True)
 with (out/'lock').open('w') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  environment=verify_pyebm();(df,truth),old=generate(campaign,spec,seed)
  td=pair_taus([truth['group_orderings'][g] for g in ('e2','e33','e4')]);truth.update(pair_null=[x==0 for x in td],pair_distances=td,h1_distance=max(td))
  source={str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(HERE.rglob('*.py')) if '__pycache__'not in str(p)}
  identity={'campaign':campaign,'cell':spec,'seed':seed,'budget':budget,'arms':arms,'benchmark':benchmark_b is not None,'exact_fixed_budget_stopping':safe_stopping,'campaign_sha256':hashlib.sha256((HERE/'campaign.json').read_bytes()).hexdigest(),'source_sha256':source,'environment':environment,'data_sha256':truth['manifest']['data_sha256'],'truth_sha256':canonical_sha(truth),'full599_permutation_label_stream_sha256':permutation_stream(df,seed)}
  if benchmark_b is None:
   assert gate['source_sha256']==source,'Production source differs from frozen gate'
   assert gate['environment_versions']==environment['versions'],'Production numerical environment differs from gate'
  ident=canonical_sha(identity)
  if (out/'manifest.json').exists():assert json.loads((out/'manifest.json').read_text())['identity_sha256']==ident,'Changed environment/config/source: use a versioned output'
  else:atomic_json(out/'manifest.json',identity|{'identity_sha256':ident});atomic_json(out/'truth.json',truth)
  events=out/'fit_records.jsonl';records={}
  if events.exists():
   raw=events.read_bytes();lines=raw.splitlines(keepends=True);offset=0
   for li,line in enumerate(lines):
    try:r=json.loads(line)
    except (json.JSONDecodeError,UnicodeDecodeError):
     if li!=len(lines)-1 or line.endswith(b'\n'):raise
     (out/('torn_tail_'+hashlib.sha256(line).hexdigest()+'.bin')).write_bytes(line)
     with events.open('r+b') as recovery:recovery.truncate(offset);recovery.flush();os.fsync(recovery.fileno())
     break
    key=(r['engine'],r['perm_id']);assert key not in records;records[key]=r;offset+=len(line)
    if li==len(lines)-1 and not line.endswith(b'\n'):
     with events.open('ab') as recovery:recovery.write(b'\n');recovery.flush();os.fsync(recovery.fileno())
  ec=EngineConfig(mode='repaired',expected_events=14,biomarker_names=tuple(truth['biomarker_names']),audit_original_neighbors=True)
  signal.signal(signal.SIGALRM,timeout_handler)
  begin=time.monotonic()
  with events.open('a',buffering=1) as stream:
   for pid in [None]+list(range(budget)):
    x=df if pid is None else permute(df,truth,seed,'diagnosis',pid,{'stratify':'diagnosis','require_max':True})
    labels_sha=hashlib.sha256(np.asarray(x.APOE,dtype='<i8').tobytes()).hexdigest()
    for arm in arms:
     if (arm,pid) in records:continue
     state=arm_state(records,arm,budget)
     if pid is not None and safe_stopping and state is not None and state['all_decisions_resolved']:continue
     started=time.monotonic();record={'seed':seed,'engine':arm,'perm_id':pid,'kind':'observed' if pid is None else 'permutation','permuted_labels_sha256':labels_sha}
     signal.alarm(cfg['fit_timeout_s']);measurement=None
     try:
      with allgroup_stopping(detailed=pid is None) as measurement, fast_likelihood_context(True):
       result=fit_orderings(x,ec) if arm=='repaired' else fit_invariant_orderings(x,ec,variant=arm)
      record['measurement']=measurement.summary()
      if not result.ok:_CACHE.clear();_CACHE_PROVENANCE.clear()
      if arm!='repaired' and result.ok:
       cache_key=_matrix_key(x,list(truth['biomarker_names']))
       if measurement.termination=='all_groups_at_upstream_tolerance':
        _CACHE_PROVENANCE[cache_key]={'certified_allgroup01':True,'origin_seed':seed,'origin_perm_id':pid,'outer_iterations':len(measurement.iterations),'termination':measurement.termination}
       if cache_key not in _CACHE_PROVENANCE:raise RuntimeError('Uncertified pooled score cache; stopping provenance missing')
       record['pooled_score_provenance']=_CACHE_PROVENANCE[cache_key]
      record.update(result.to_dict());record['taus']=pair_taus(result.orderings) if result.ok else None
      if pid is None and result.ok:
       record['distance_to_truth']=[float(kendall(o,truth['group_orderings'][g])) for o,g in zip(result.orderings,('e2','e33','e4'))]
       if old is not None and old.get(arm) is not None:
        record['historical_observed_order_match']=record['orderings']==old[arm]['observed']['orderings']
        if arm!='repaired' and not record['historical_observed_order_match']:raise ValueError('Unchanged pooled arm failed historical observed bridge')
       if arm!='repaired':save_cache(x,list(truth['biomarker_names']),out/'pooled_score_cache.npz')
      if pid is not None and result.ok:
       original_diagnostics=record['diagnostics'];record['diagnostics']=compact_diagnostics(original_diagnostics)
       optimizer=original_diagnostics.get('optimizer',[])
       record['diagnostics']['optimizer_iterations_total']=sum(max(0,v.get('nit',0)) for v in optimizer)
       record['diagnostics']['optimizer_iterations_max']=max([v.get('nit',0) for v in optimizer],default=0)
     except Exception as exc:
      # An invariant fit may populate _CACHE before the outer cap check raises.
      # Never allow the next pooled arm to consume those uncertified scores.
      _CACHE.clear();_CACHE_PROVENANCE.clear()
      record.update(status='error',error=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
      if measurement is not None:
       try:record['failed_measurement']=measurement.summary(include_trajectory=True)
       except TypeError:record['failed_measurement']=measurement.summary()
     finally:signal.alarm(0)
     record['fit_seconds']=time.monotonic()-started;stream.write(json.dumps(record,allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno());records[(arm,pid)]=record
    if safe_stopping and pid is not None and all(arm_state(records,arm,budget) and arm_state(records,arm,budget)['all_decisions_resolved'] for arm in arms):break
    if pid is None or (pid+1)%20==0:
     atomic_json(out/'progress.json',{'records':len(records),'expected':len(arms)*(budget+1),'elapsed_seconds':time.monotonic()-begin,'last_perm_id':pid})
  summary={'seed':seed,'cell':spec['id'],'budget':budget,'benchmark':benchmark_b is not None,'execution_policy':'exact_fixed_budget_decision_bounds' if safe_stopping else 'all_permutations','maximum_planned_records':len(arms)*(budget+1),'records':len(records),'failed':sum(r['status']!='ok' for r in records.values()),'arms':{}}
  for arm in arms:
   observed=records[(arm,None)];perms=[r for (e,pid),r in records.items() if e==arm and pid is not None];good=[r for r in perms if r['status']=='ok'];state=arm_state(records,arm,budget)
   a={'observed':observed,'computed_permutations':len(good),'attempted_permutations':len(perms),'failed_permutations':len(perms)-len(good),'computed_perm_ids':sorted(r['perm_id'] for r in good),'fit_seconds_sum':sum(r['fit_seconds'] for r in [observed]+perms),'mean_permutation_fit_seconds':float(np.mean([r['fit_seconds'] for r in perms])) if perms else None,'full_permutation_budget_complete':len(good)==budget}
   if state is not None:
    a.update(state)
    if len(good)==budget:a['adjusted_p']=state['adjusted_p_lower']
    else:a['adjusted_p']=None
    if state['all_decisions_resolved']:
     dec=np.array(state['pair_reject']);a.update(any_reject=bool(np.any(dec)),true_pair_detection=bool(np.any(dec&(~np.array(truth['pair_null'])))))
   summary['arms'][arm]=a
  summary['complete']=all(a['observed']['status']=='ok' and (a.get('all_decisions_resolved',False) if safe_stopping else a['full_permutation_budget_complete']) for a in summary['arms'].values())
  summary['unresolved_arms']=[arm for arm,a in summary['arms'].items() if a['observed']['status']!='ok' or not a.get('all_decisions_resolved',False)]
  atomic_json(out/'result.json',summary)
  print(json.dumps({k:summary[k] for k in ('cell','seed','complete','records','failed')},sort_keys=True),flush=True)
  return summary

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--campaign',choices=['power','robustness'],required=True);ap.add_argument('--cell-index',type=int,required=True);ap.add_argument('--seed-index',type=int,default=0);ap.add_argument('--output',required=True);ap.add_argument('--benchmark-b',type=int);ap.add_argument('--all-arms',action='store_true');a=ap.parse_args()
 run(a.campaign,a.cell_index,a.seed_index,a.output,a.benchmark_b,a.all_arms)
if __name__=='__main__':main()
