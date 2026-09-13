#!/bin/bash
# One-shot execution amendment. Existing scientific snapshots and receipt identities stay immutable.
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import sys,json,fcntl,hashlib,subprocess
from pathlib import Path
from datetime import datetime,timezone
root=Path('/N/slate/tg11/ebmcal_v2')
sys.path.insert(0,str(root/'toolkit/hpc/v2'))
from dispatch_campaign import durable_json
targets={'10263134':'confirm_stage_a','10263135':'confirm_power_global_a','10263136':'confirm_pair_complete_ref_a','10263137':'confirm_pair_power_k27_a','10263138':'confirm_pair_partial_e2_topup_a'}
ledger=root/'control/submission_ledger.json'
evidence_path=root/'control/campaign_versions/confirmation_20260908_expand1024/pending_walltime_amendment.json'
def query():
 return subprocess.run(['squeue','-r','-u','tg11','-h','-o','%i|%T|%l|%k'],check=True,capture_output=True,text=True,timeout=30).stdout
with (root/'campaign.dispatch.lock').open('a+') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 assert Path(json.loads((root/'campaign.dispatch.binding.json').read_text())['ledger_path']).resolve()==ledger.resolve()
 assert not evidence_path.exists(), 'One-shot evidence already exists; inspect rather than rerun.'
 raw=ledger.read_bytes(); data=json.loads(raw)
 records={}
 for job_id,run in targets.items():
  matches=[e for e in data['submissions'] if e.get('job_id')==job_id]
  assert len(matches)==1
  e=matches[0]; assert e['run_id']==run and e['state']=='submitted' and e['walltime_s']==7200 and e['cpus_per_task']==16
  snapshot=root/'snapshots'/run
  assert hashlib.sha256((snapshot/'SHA256SUMS').read_bytes()).hexdigest()==e['snapshot_sha256']
  for line in (snapshot/'SHA256SUMS').read_text().splitlines():
   digest,name=line.split(None,1); p=snapshot/name.lstrip('*')
   assert p.resolve().is_relative_to(snapshot.resolve()) and hashlib.sha256(p.read_bytes()).hexdigest()==digest
  cfg=json.loads((snapshot/'configs/v2'/f'{run}.json').read_text())
  assert cfg['engines']==['repaired'] and cfg['fit_timeout_s']==300 and not cfg.get('paired_standard_engines')
  script=(snapshot/'hpc/v2/run_array.sbatch').read_text()
  assert '#SBATCH --signal=B:USR1@1500' in script and 'exec "$HOME/ebmcal-env/bin/python" scripts/v2/run_v2.py' in script
  records[job_id]=e
 before=query(); selected={j:[] for j in targets}; skipped=[]
 for line in before.splitlines():
  job,state,limit,comment=line.strip().split('|',3)
  if '_' not in job or job.split('_')[0] not in targets:continue
  parent,task=job.split('_'); e=records[parent]
  assert task.isdigit() and int(task) in e['array_task_ids'] and comment=='ebmcal_v2:'+e['plan_id']
  if state=='PENDING':
   assert limit=='2:00:00', (job,limit)
   selected[parent].append(job)
  else:skipped.append({'job_id':job,'state':state,'time_limit':limit})
 assert any(selected.values())
 stamp=datetime.now(timezone.utc).isoformat()
 evidence={'created_utc':stamp,'status':'intent','original_ledger_sha256':hashlib.sha256(raw).hexdigest(),'scope':'Only individually enumerated members observed PENDING of the five exact accepted receipts. A scheduler start race during this bounded call can only shorten a newly started member with nearly its full hour remaining. Previously running members are excluded.','old_walltime_s':7200,'new_walltime_s':3600,'scientific_changes':False,'signal':'B:USR1@1500','nominal_max_wave_s':600,'nominal_extra_drain_margin_s':900,'limitation':'Nominal timeout-wave bound excludes arbitrary native/system blocking; checkpoint and fit status audits still required.','basis':['https://slurm.schedmd.com/sbatch.html','https://slurm.schedmd.com/scontrol.html'],'selected':selected,'skipped':skipped,'queue_before':before,'outcomes':[]}
 durable_json(evidence_path,evidence)
 for parent,ids in selected.items():
  if not ids:continue
  argv=['scontrol','update','JobId='+','.join(ids),'TimeLimit=01:00:00']
  outcome={'parent_job_id':parent,'job_ids':ids,'argv':argv,'started_utc':datetime.now(timezone.utc).isoformat(),'state':'update_unknown'}
  evidence['outcomes'].append(outcome);durable_json(evidence_path,evidence)
  try:
   result=subprocess.run(argv,capture_output=True,text=True,timeout=30,check=False)
   outcome.update(returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,state='command_returned')
  except subprocess.TimeoutExpired as exc:
   outcome['error']=str(exc);durable_json(evidence_path,evidence);raise
  after=query(); outcome['queue_after']=after
  observed={line.strip().split('|',3)[0]:line.strip().split('|',3) for line in after.splitlines()}
  verified=all(j in observed and observed[j][2]=='1:00:00' and observed[j][3]=='ebmcal_v2:'+records[parent]['plan_id'] for j in ids)
  outcome['state']='verified' if result.returncode==0 and verified else 'needs_reconciliation'
  outcome['finished_utc']=datetime.now(timezone.utc).isoformat();durable_json(evidence_path,evidence)
  records[parent].setdefault('execution_adjustments',[]).append({'kind':'pending_member_walltime_reduction','evidence_file':str(evidence_path),'job_ids':ids,'old_walltime_s':7200,'new_walltime_s':3600,'state':outcome['state'],'utc':outcome['finished_utc']})
  durable_json(ledger,data)
  assert outcome['state']=='verified','Inspect evidence; do not automatically retry uncertain scheduler update.'
 evidence['status']='verified';evidence['finished_utc']=datetime.now(timezone.utc).isoformat();durable_json(evidence_path,evidence)
 print(json.dumps({'status':'verified','evidence':str(evidence_path),'changed':{k:len(v) for k,v in selected.items()},'skipped':skipped}))
PY
