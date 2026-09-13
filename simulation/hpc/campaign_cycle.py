#!/usr/bin/env python3
"""One bounded HPC campaign cycle; no scheduler loop, daemon or automatic retry.

Explicit --apply dispatches a fresh plan using the existing bound ledger.
Each invocation stores a new receipt directory; old proposals/receipts are kept.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import uuid

ROOT = Path('/N/slate/tg11/ebmcal_v2')


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--campaign', type=Path, required=True)
    ap.add_argument('--weekly-used-percent', type=float, required=True)
    ap.add_argument('--usage-observed-utc', required=True)
    ap.add_argument('--apply', action='store_true')
    a = ap.parse_args(argv)
    now = datetime.now(timezone.utc)
    observed = datetime.fromisoformat(a.usage_observed_utc.replace('Z', '+00:00'))
    if observed.tzinfo is None or not 0 <= (now-observed).total_seconds() <= 1800:
        raise ValueError('Use an explicit account usage observation from the last30 minutes')
    if not 0 <= a.weekly_used_percent <= 100:
        raise ValueError('Usage must be a percentage')
    campaign = a.campaign.resolve()
    proposals = ROOT/'control/campaign_versions'
    if not campaign.is_relative_to(proposals) or not campaign.is_file():
        raise ValueError('Use a saved versioned proposal in this HPC project')
    ledger = ROOT/'control/submission_ledger.json'
    if not ledger.is_file():
        raise ValueError('Existing authoritative ledger is required; never create an empty replacement')
    code = ROOT/'toolkit/hpc/v2'
    cycle = ROOT/'control/cycles'/(now.strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8])
    cycle.mkdir(parents=True, exist_ok=False)
    evidence = {'created_utc':now.isoformat(),'campaign':str(campaign),
                'campaign_sha256':hashlib.sha256(campaign.read_bytes()).hexdigest(),
                'ledger':str(ledger),'weekly_used_percent':a.weekly_used_percent,
                'usage_observed_utc':a.usage_observed_utc,'apply_requested':a.apply,
                'tools_sha256':{name:hashlib.sha256((code/name).read_bytes()).hexdigest()
                    for name in ('campaign_cycle.py','reconcile_ledger.py','accounted_members.py','collect_queue.py','manage_campaign.py','dispatch_campaign.py','inventory_run.py')}}
    (cycle/'input.json').write_text(json.dumps(evidence,indent=2)+'\n')
    steps=[]
    def run(name, argv):
        result=subprocess.run(argv, capture_output=True, text=True, timeout=180, check=False)
        (cycle/(name+'.stdout')).write_text(result.stdout)
        (cycle/(name+'.stderr')).write_text(result.stderr)
        steps.append({'step':name,'returncode':result.returncode})
        if result.returncode:
            raise RuntimeError(f'{name} returned{result.returncode}; inspect saved evidence, do not blindly retry')
    status='incomplete'
    error=None
    try:
        run('reconcile',[sys.executable,str(code/'reconcile_ledger.py'),'--ledger',str(ledger)])
        run('queue',[sys.executable,str(code/'collect_queue.py'),'--registry',str(ROOT/'control/job_registry.json'),'--ledger',str(ledger),'--output',str(cycle/'queue.json')])
        run('plan',[sys.executable,str(code/'manage_campaign.py'),'--campaign',str(campaign),'--queue',str(cycle/'queue.json'),'--ledger',str(ledger),'--weekly-used-percent',str(a.weekly_used_percent),'--output',str(cycle/'plan.json')])
        # Keep the original planning evidence and timestamps intact. Expensive
        # inventory planning does not authorize dispatch against its old queue.
        run('queue_dispatch',[sys.executable,str(code/'collect_queue.py'),'--registry',str(ROOT/'control/job_registry.json'),'--ledger',str(ledger),'--output',str(cycle/'queue_dispatch.json')])
        run('dispatch',[sys.executable,str(code/'dispatch_campaign.py'),'--project-root',str(ROOT),'--plan',str(cycle/'plan.json'),'--queue',str(cycle/'queue_dispatch.json'),'--ledger',str(ledger),'--weekly-used-percent',str(a.weekly_used_percent),'--plan-max-age-s','300','--output',str(cycle/'dispatch.json')]+(['--apply'] if a.apply else []))
        status='cycle_finished'
    except (RuntimeError,subprocess.TimeoutExpired,OSError) as exc:
        error=str(exc)
    result={'cycle':str(cycle),'status':status,'steps':steps,'error':error}
    if (cycle/'plan.json').exists():
        plan=json.loads((cycle/'plan.json').read_text())
        result.update(cpu_budget=plan.get('cpu_budget'),user_job_budget=plan.get('user_job_budget'),
                      planned=[{'run_id':x['run_id'],'members':len(x['array_task_ids']),'parallel':x['max_concurrent'],'generation':x['dispatch_generation']} for x in plan.get('submission_plans',[])],
                      runs=[{'run_id':x['run_id'],'state':x['state'],'counts':x.get('state_counts'),'issues':x['issues']} for x in plan.get('runs',[])],
                      plan_issues=plan.get('issues',[]))
    if (cycle/'dispatch.json').exists():
        dispatch=json.loads((cycle/'dispatch.json').read_text())
        result.update(dispatch_status=dispatch.get('status'),outcomes=dispatch.get('outcomes'))
    (cycle/'cycle.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    return 0 if status=='cycle_finished' else 3


if __name__=='__main__':
    sys.exit(main())
