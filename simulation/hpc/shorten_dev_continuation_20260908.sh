#!/bin/bash
# Exact one-shot execution amendment; never retry an existing evidence record.
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import fcntl, hashlib, json, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path
root = Path('/N/slate/tg11/ebmcal_v2').resolve()
sys.path.insert(0, str(root/'toolkit/hpc/v2'))
import manage_campaign as M
from dispatch_campaign import durable_json, submission_identity
ledger = root/'control/submission_ledger.json'
evidence_path = root/'control/campaign_versions/confirmation_20260908_dev_continue_1h/pending_dev_walltime_amendment.json'
parent, member = '10264821', '10264821_13'
plan_id = '6d9d64a72f26dad6f28d1ae32deb257a2771b96969384d6405762deb774ba495'
def utc(): return datetime.now(timezone.utc).isoformat()
def require(value, message):
    if not value:
        raise ValueError(message)
def query():
    p = subprocess.run(['squeue','-r','-u','tg11','-h','-o','%i|%T|%l|%k'],
                       capture_output=True,text=True,timeout=30,check=True)
    return p.stdout
with (root/'campaign.dispatch.lock').open('a+') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
    binding = json.loads((root/'campaign.dispatch.binding.json').read_text())
    require(binding.get('project') == 'ebmcal_v2' and
            Path(binding['ledger_path']).resolve() == ledger.resolve() and ledger.resolve().is_relative_to(root),
            'Project root/ledger binding differs')
    require(evidence_path.resolve().is_relative_to(root/'control/campaign_versions'), 'Evidence path escapes project versions')
    require(not evidence_path.exists(), 'One-shot evidence exists; inspect, do not rerun.')
    original = ledger.read_bytes(); data = json.loads(original)
    require(data.get('schema_version') == 'campaign_ledger_v1' and data.get('project') == 'ebmcal_v2', 'Wrong authoritative project ledger')
    matches = [e for e in data['submissions'] if e.get('job_id') == parent]
    require(len(matches) == 1, 'Exact accepted parent receipt must be unique')
    receipt = matches[0]
    require(receipt['plan_id'] == plan_id and M.I.digest(submission_identity(receipt)) == plan_id
            and receipt['run_id'] == 'dev_calibration_a', 'Original receipt identity hash differs')
    require(receipt['state'] == 'submitted' and receipt['array_task_ids'] == [13], 'Exact accepted member differs')
    require(receipt['walltime_s'] == 14400 and receipt['cpus_per_task'] == 16 and receipt['nchunks'] == 8,
            'Original execution dimensions differ')
    adjustments = receipt.get('execution_adjustments', [])
    require(isinstance(adjustments, list) and all(isinstance(a, dict) for a in adjustments), 'Malformed adjustment history')
    require(not any(a.get('kind') == 'pending_member_walltime_reduction' and member in a.get('job_ids', [])
                    for a in adjustments), 'Prior member amendment exists; inspect, do not rerun')
    snapshot = root/'snapshots/dev_calibration_a'
    require(snapshot.resolve().is_relative_to(root/'snapshots'), 'Snapshot escapes bound project root')
    frozen = M.snapshot_info(snapshot,'dev_calibration_a',receipt['snapshot_sha256'])
    cfg = json.loads(frozen['config_path'].read_text())
    proof = M.validate_execution_walltime(frozen,cfg,16,3600,M.DEV_CALIBRATION_1H_POLICY,nchunks=8)
    query_started = time.monotonic()
    before = query()
    rows = [line.strip().split('|',3) for line in before.splitlines() if line.split('|',1)[0] == member]
    require(len(rows) <= 1 and all(len(row) == 4 for row in rows), 'Malformed/duplicate exact member queue rows')
    evidence = {'created_utc':utc(),'original_ledger_sha256':hashlib.sha256(original).hexdigest(),
                'manager_sha256':hashlib.sha256(Path(M.__file__).read_bytes()).hexdigest(),
                'job_id':member,'plan_id':plan_id,'old_walltime_s':14400,'new_walltime_s':3600,
                'scientific_changes':False,'execution_policy':M.DEV_CALIBRATION_1H_POLICY,
                'validated_evidence':proof,'queue_before':before,
                'scope':'Only this exact member observed pending. Already running or absent members are skipped. A start race can only shorten a newly started member with nearly its full hour remaining.'}
    if not rows or rows[0][1] != 'PENDING':
        evidence.update(status='skipped_no_change',observed_member=rows,finished_utc=utc())
        durable_json(evidence_path,evidence)
        print(json.dumps({'status':evidence['status'],'evidence':str(evidence_path),'observed':rows}))
        sys.exit(0)
    require(rows[0][2:] == ['4:00:00','ebmcal_v2:'+plan_id], 'Pending member walltime/comment differs')
    argv = ['scontrol','update','JobId='+member,'TimeLimit=01:00:00']
    evidence.update(status='update_unknown',argv=argv,subprocess_start_utc=utc())
    durable_json(evidence_path,evidence)
    query_age = time.monotonic()-query_started
    if not 0 <= query_age <= 60:
        evidence.update(status='skipped_stale_queue_no_change',query_age_s=query_age,
                        subprocess_invoked=False,finished_utc=utc())
        durable_json(evidence_path,evidence)
        print(json.dumps({'status':evidence['status'],'evidence':str(evidence_path)}))
        sys.exit(0)
    evidence.update(query_age_s=query_age,subprocess_invoked=True)
    try:
        response = subprocess.run(argv,capture_output=True,text=True,timeout=30,check=False)
        evidence.update(returncode=response.returncode,stdout=response.stdout,stderr=response.stderr)
        after = query(); evidence['queue_after'] = after
        current = [line.strip().split('|',3) for line in after.splitlines() if line.split('|',1)[0] == member]
        verified = response.returncode == 0 and len(current) == 1 and current[0][2:] == ['1:00:00','ebmcal_v2:'+plan_id]
        evidence.update(status='verified' if verified else 'needs_reconciliation',finished_utc=utc())
    except (subprocess.TimeoutExpired,subprocess.CalledProcessError,OSError) as exc:
        evidence.update(error=str(exc),exception_type=type(exc).__name__,finished_utc=utc())
    durable_json(evidence_path,evidence)
    receipt.setdefault('execution_adjustments',[]).append({'kind':'pending_member_walltime_reduction',
        'evidence_file':str(evidence_path),'job_ids':[member],'old_walltime_s':14400,'new_walltime_s':3600,
        'execution_walltime_policy':M.DEV_CALIBRATION_1H_POLICY,'state':evidence['status'],'utc':evidence['finished_utc']})
    durable_json(ledger,data)
    require(evidence['status'] == 'verified', 'Inspect evidence; no automatic retry of uncertain scheduler update.')
    print(json.dumps({'status':'verified','evidence':str(evidence_path),'job_id':member}))
PY
