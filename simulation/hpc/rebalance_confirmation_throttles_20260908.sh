#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import sys,json,re,fcntl,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
root=Path('/N/slate/tg11/ebmcal_v2')
sys.path.insert(0,str(root/'toolkit/hpc/v2'))
from dispatch_campaign import durable_json
targets={'10263134':'confirm_stage_a','10263135':'confirm_power_global_a','10263136':'confirm_pair_complete_ref_a','10263137':'confirm_pair_power_k27_a'}
ledger=root/'control/submission_ledger.json'
path=root/'control/campaign_versions/confirmation_20260908_expand1024/throttle_rebalance.json'
def command(argv):return subprocess.run(argv,check=True,capture_output=True,text=True,timeout=30).stdout
with (root/'campaign.dispatch.lock').open('a+') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 assert Path(json.loads((root/'campaign.dispatch.binding.json').read_text())['ledger_path']).resolve()==ledger.resolve()
 assert not path.exists(),'One-shot; inspect saved evidence rather than repeat.'
 raw=ledger.read_bytes();data=json.loads(raw);entries={};before={}
 for job,run in targets.items():
  matches=[e for e in data['submissions'] if e.get('job_id')==job]
  assert len(matches)==1
  entry=matches[0];assert entry['run_id']==run and entry['state']=='submitted' and entry['max_concurrent']==8
  detail=command(['scontrol','show','job','-o',job])
  assert set(re.findall(r'\bArrayTaskThrottle=(\d+)',detail))=={'8'}
  assert set(re.findall(r'\bComment=(\S+)',detail))=={'ebmcal_v2:'+entry['plan_id']}
  entries[job]=entry;before[job]=detail
 evidence={'created_utc':datetime.now(timezone.utc).isoformat(),'state':'intent','original_ledger_sha256':hashlib.sha256(raw).hexdigest(),'purpose':'Redistribute project concurrency across all ten frozen required runs; prioritize first coverage for remaining e4 topup/e33 partial-null/K14-K55. No running task is canceled; current allocations drain naturally.','old_array_throttle':8,'new_array_throttle':4,'scientific_changes':False,'before':before,'outcomes':[]}
 durable_json(path,evidence)
 for job,run in targets.items():
  item={'job_id':job,'run_id':run,'state':'update_unknown'};evidence['outcomes'].append(item);durable_json(path,evidence)
  result=subprocess.run(['scontrol','update','JobId='+job,'ArrayTaskThrottle=4'],capture_output=True,text=True,timeout=30,check=False)
  item.update(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr)
  detail=command(['scontrol','show','job','-o',job]);item['after']=detail
  item['state']='verified' if result.returncode==0 and set(re.findall(r'\bArrayTaskThrottle=(\d+)',detail))=={'4'} else 'needs_reconciliation'
  durable_json(path,evidence)
  entries[job].setdefault('execution_adjustments',[]).append({'kind':'array_throttle_reduction','old_max_concurrent':8,'new_max_concurrent':4,'state':item['state'],'evidence_file':str(path),'utc':datetime.now(timezone.utc).isoformat()})
  durable_json(ledger,data)
  assert item['state']=='verified','Inspect uncertain update before further actions'
 evidence['state']='verified';evidence['finished_utc']=datetime.now(timezone.utc).isoformat();durable_json(path,evidence)
 print(json.dumps({'state':'verified','evidence':str(path),'new_throttles':{j:4 for j in targets}}))
PY
