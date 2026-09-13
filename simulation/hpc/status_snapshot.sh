#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import json,subprocess
from datetime import datetime,timezone
from pathlib import Path
from collections import Counter,defaultdict
root=Path('/N/slate/tg11/ebmcal_v2')
start=datetime.now(timezone.utc)
registered={str(e['job_id']) for e in json.loads((root/'control/submission_ledger.json').read_text())['submissions'] if e.get('job_id')}
registered.update(str(e['job_id']) for e in json.loads((root/'control/job_registry.json').read_text()).values() if e.get('job_id'))
q=subprocess.run(['squeue','-r','-u','tg11','-h','-o','%i|%T|%C|%l|%k'],check=True,capture_output=True,text=True,timeout=30)
jobs=[];groups=defaultdict(Counter)
for line in q.stdout.splitlines():
 job,state,cpus,limit,comment=line.strip().split('|',4)
 entry={'job_id':job,'state':state,'cpus':int(cpus),'time_limit':limit,'comment':comment,'own_project':job.split('_')[0] in registered or comment.startswith('ebmcal_v2:')};jobs.append(entry)
 if entry['own_project']:
  groups[job.split('_')[0]][state]+=1;groups[job.split('_')[0]]['limit_'+limit]+=1
runs=[]
for run in sorted((root/'runs').iterdir()):
 ps=list(run.glob('cell_*/chunk_*/progress.json'))
 if not ps:continue
 ds=[(p,json.loads(p.read_text())) for p in ps]
 cells={}
 for cell in sorted({p.parent.parent.name for p,d in ds}):
  vals=[d for p,d in ds if p.parent.parent.name==cell]
  cells[cell]={'started_chunks':len(vals),'done_chunks':sum(d['done'] for d in vals),'fit_records':sum(d['fit_records'] for d in vals),'rows':sum(d['rows'] for d in vals)}
 runs.append({'run_id':run.name,'started_chunks':len(ds),'done_chunks':sum(d['done'] for p,d in ds),'fit_records':sum(d['fit_records'] for p,d in ds),'rows':sum(d['rows'] for p,d in ds),'statuses':dict(sum((Counter(d['statuses']) for p,d in ds),Counter())),'latest_utc':max(d['updated_utc'] for p,d in ds),'cells':cells})
result={'captured_start_utc':start.isoformat(),'captured_end_utc':datetime.now(timezone.utc).isoformat(),'capture_policy':'Independent live progress files; descriptive execution state, not a complete cohort analysis. Own-project resource totals include ledger/registry legacy IDs as well as plan comments.','queue_argv':q.args,'queue_stderr':q.stderr,'jobs':jobs,'own_jobs_by_parent':{k:dict(v) for k,v in groups.items()},'all_user_active_members':len(jobs),'own_running_allocated_cpus':sum(j['cpus'] for j in jobs if j['own_project'] and j['state']=='RUNNING'),'runs':runs}
directory=root/'diagnostics/status_snapshots';directory.mkdir(parents=True,exist_ok=True)
path=directory/(start.strftime('%Y%m%dT%H%M%S%f')+'.json')
with path.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'evidence_file':str(path),'captured_utc':result['captured_end_utc'],'all_user_active_members':result['all_user_active_members'],'own_running_allocated_cpus':result['own_running_allocated_cpus'],'jobs':result['own_jobs_by_parent'],'runs':runs}))
PY
