#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import sys, json, fcntl, hashlib, subprocess
from pathlib import Path
from datetime import datetime,timezone
root=Path('/N/slate/tg11/ebmcal_v2')
sys.path.insert(0,str(root/'toolkit/hpc/v2'))
from dispatch_campaign import durable_json
plan_id='0603409cc8bed311a74875527680406fcfc556ad64adc7b5373306dbf2b93b79'
ledger=root/'control/submission_ledger.json'
receipt=root/'control/campaign_versions/confirmation_20260908_1455/dispatch_apply.json'
with (root/'campaign.dispatch.lock').open('a+') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX)
 assert Path(json.loads((root/'campaign.dispatch.binding.json').read_text())['ledger_path']).resolve()==ledger.resolve()
 original=ledger.read_bytes(); data=json.loads(original)
 candidates=[e for e in data['submissions'] if e['plan_id']==plan_id]
 assert len(candidates)==1
 entry=candidates[0]
 assert entry['state']=='submission_unknown' and entry['run_id']=='confirm_core_a' and not entry.get('job_id')
 outcome=[e for e in json.loads(receipt.read_text())['outcomes'] if e['plan_id']==plan_id][0]
 expected='sbatch: error: QOSMaxSubmitJobPerUserLimit\nsbatch: error: Batch job submission failed: Job violates accounting/QOS policy (job submit limit, user\'s size and/or time limits)\n'
 assert outcome['returncode']==1 and outcome['stdout']=='' and outcome['stderr']==expected and not outcome.get('job_id')
 assert entry['returncode']==outcome['returncode'] and entry['stderr']==outcome['stderr'] and entry['stdout']==''
 queue=subprocess.run(['squeue','-r','-u','tg11','-h','-o','%i|%T|%k'],check=True,capture_output=True,text=True,timeout=30).stdout
 assert plan_id not in queue
 accounting=subprocess.run(['sacct','-X','-n','-P','-u','tg11','-S','2026-09-08T10:54:00','-o','JobID%80,State%40,Comment%100,Submit,Start,End'],check=True,capture_output=True,text=True,timeout=30).stdout
 assert plan_id not in accounting
 qos=subprocess.run(['sacctmgr','-n','-P','show','qos','where','name=allocated','format=Name,MaxJobsPU,MaxSubmitPU'],check=True,capture_output=True,text=True,timeout=30).stdout
 assert 'allocated|1000|1000' in qos
 stamp=datetime.now(timezone.utc).isoformat()
 evidence={'kind':'explicit_scheduler_submission_rejection','verified_utc':stamp,'plan_id':plan_id,'basis':'Synchronous sbatch return code 1, empty stdout, no job ID, and exact QOSMaxSubmitJobPerUserLimit + Batch job submission failed response. Queue/accounting absence is corroboration, not the basis for treating an unknown transport response as unsubmitted.','original_outcome':outcome,'original_ledger_sha256':hashlib.sha256(original).hexdigest(),'queue':queue,'accounting':accounting,'qos':qos,'query_scope':'all current tg11 jobs and allocations since 10:54 cluster local time','subprocess_invoked':True,'scheduler_rejected':True}
 evidence_path=root/'control/campaign_versions/confirmation_20260908_1455/core_qos_reconciliation.json'
 assert not evidence_path.exists()
 durable_json(evidence_path,evidence)
 entry['reconciliation_history']=entry.get('reconciliation_history',[])+[{'previous_state':entry['state'],'verified_utc':stamp,'evidence_file':str(evidence_path)}]
 entry['state']='not_submitted_verified';entry['updated_utc']=stamp;entry['reconciliation_evidence']=evidence
 durable_json(ledger,data)
 print(json.dumps({'plan_id':plan_id,'state':entry['state'],'evidence':str(evidence_path)}))
PY
