#!/bin/bash
# Exact one-shot throttle3->2 request; no sbatch, cancellation, time/ID edits.
# Preparation only: do not deploy/execute before root's independent final queue-fill review.
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - "$0" <<'PY'
import fcntl, hashlib, json, os, re, subprocess, sys, tempfile, time
from datetime import datetime, timezone
from pathlib import Path
root = Path('/N/slate/tg11/ebmcal_v2').resolve()
helper = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(root/'toolkit/hpc/v2'))
import manage_campaign as M
import dispatch_campaign as D
ledger = root/'control/submission_ledger.json'
evidence_path = root/'control/campaign_versions/confirmation_20260908_final_queue_fill/reserve_final_queue_slot.json'
parent = '10264903'
plan_id = 'c997a2e98a610c02aabd45acfb81f8396cdcd7507ee61ff4b2219925b1e380bc'
members = list(range(414, 430))
snapshot_hash = 'e4701a864230fbc6faf64709defb4dec5fa5a3d07ac3625b95cc36cb06847464'
config_hash = '11b716e46bb5055f8b24a3ee5b756683508f66165fe3c361b359b477edf30a52'
kind = 'exact_final_queue_fill_array_throttle_reduction'
def utc(): return datetime.now(timezone.utc).isoformat()
def require(value, message):
    if not value:
        raise ValueError(message)
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def durable(path, value):
    # Random exclusive temp creation avoids following a pre-existing temp alias.
    with tempfile.NamedTemporaryFile(mode='w',dir=path.parent,prefix='.final-queue-',delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(value,stream,indent=2,allow_nan=False); stream.write('\n')
            stream.flush(); os.fsync(stream.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.replace(temporary,path)
        fd = os.open(path.parent,os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)
    finally: temporary.unlink(missing_ok=True)
def query(argv):
    p = subprocess.run(argv,capture_output=True,text=True,timeout=30,check=True)
    return p.stdout
queue_argv = ['squeue','-r','-j',parent,'-h','-o','%i|%T|%l|%k']
control_argv = ['scontrol','show','job',parent]
def parse_queue(raw, pending=True):
    rows = {}
    for line in raw.splitlines():
        if not line.strip(): continue
        values = line.strip().split('|')
        require(len(values) == 4, 'Malformed target queue row')
        job, state, wall, comment = values
        match = re.fullmatch(parent+r'_([0-9]+)',job)
        require(match is not None and str(int(match[1])) == match[1], 'Queue must contain exact expanded target member IDs')
        index = int(match[1])
        require(index in members and index not in rows, 'Unplanned or duplicate target queue member')
        require(wall in ('2:00:00','02:00:00') and comment == 'ebmcal_v2:'+plan_id, 'Target queue walltime/comment differs')
        require(state == 'PENDING' if pending else state in ('PENDING','RUNNING'), 'Target member no longer pending or unrecognized final state')
        rows[index] = state
    require(set(rows) == set(members), 'Not every original target member is present')
    return rows
def parse_control(raw, throttle, pending=True):
    values = {}
    for key, value in re.findall(r'(?:^|\s)([A-Za-z][A-Za-z0-9]*)=([^\s]*)',raw):
        if key in ('JobId','ArrayJobId','ArrayTaskId','ArrayTaskThrottle','JobState','Comment','TimeLimit','NumCPUs'):
            require(key not in values, 'Duplicate target controller identity field')
            values[key] = value
    require(values.get('JobId') == parent and values.get('ArrayJobId') == parent, 'Controller exact parent identity differs')
    require(values.get('Comment') == 'ebmcal_v2:'+plan_id, 'Controller comment differs')
    require(values.get('ArrayTaskThrottle') == str(throttle), 'Controller throttle differs')
    require(values.get('TimeLimit') in ('2:00:00','02:00:00') and values.get('NumCPUs') == '16', 'Controller walltime/CPU differs')
    require(values.get('JobState') == 'PENDING' if pending else values.get('JobState') in ('PENDING','RUNNING'), 'Controller target is no longer pending')
    expression = values.get('ArrayTaskId','')
    if '%' in expression:
        expression, suffix = expression.split('%',1)
        require(suffix == str(throttle), 'Controller member-list throttle differs')
    indexes = []
    for part in expression.split(','):
        match = re.fullmatch(r'([0-9]+)(?:-([0-9]+))?',part)
        require(match is not None, 'Unsupported controller member expression')
        first, last = int(match[1]), int(match[2] or match[1])
        require(0 <= first <= last < 1000, 'Controller member range outside frozen layout')
        indexes.extend(range(first,last+1))
    require(len(indexes) == len(set(indexes)) and sorted(indexes) == members, 'Controller array member allocation differs')
    return values
with (root/'campaign.dispatch.lock').open('a+') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX | fcntl.LOCK_NB)
    binding = json.loads((root/'campaign.dispatch.binding.json').read_text())
    require(binding.get('project') == 'ebmcal_v2' and Path(binding['ledger_path']).resolve() == ledger.resolve()
            and ledger.resolve().is_relative_to(root) and not ledger.is_symlink(), 'Project root/ledger binding differs')
    require(evidence_path.resolve().is_relative_to(root/'control/campaign_versions'), 'Evidence path escapes bound project versions')
    require(not evidence_path.exists() and not evidence_path.is_symlink(), 'One-shot evidence exists; inspect, do not rerun')
    original = ledger.read_bytes(); data = json.loads(original)
    require(data.get('schema_version') == 'campaign_ledger_v1' and data.get('project') == 'ebmcal_v2', 'Wrong authoritative project ledger')
    matches = [e for e in data['submissions'] if e.get('job_id') == parent]
    require(len(matches) == 1, 'Exact accepted parent receipt must be unique')
    receipt = matches[0]
    require(receipt.get('state') == 'submitted' and receipt.get('run_id') == 'confirm_core_a'
            and receipt.get('plan_id') == plan_id and M.I.digest(D.submission_identity(receipt)) == plan_id,
            'Original accepted receipt identity hash differs')
    identity_before = D.submission_identity(receipt)
    require(receipt.get('array_task_ids') == members and receipt.get('snapshot_sha256') == snapshot_hash
            and receipt.get('cpus_per_task') == 16 and receipt.get('max_concurrent') == 3
            and receipt.get('walltime_s') == 7200 and receipt.get('nchunks') == 500
            and receipt.get('mapping') == 'chunk_major', 'Original exact execution dimensions differ')
    adjustments = receipt.get('execution_adjustments',[])
    require(isinstance(adjustments,list) and all(isinstance(a,dict) for a in adjustments), 'Malformed adjustment history')
    require(not any(a.get('kind') == kind for a in adjustments), 'Prior throttle amendment exists; inspect, do not rerun')
    snapshot = root/'snapshots/confirm_core_a'
    require(snapshot.resolve().is_relative_to(root/'snapshots'), 'Snapshot escapes bound project root')
    frozen = M.snapshot_info(snapshot,'confirm_core_a',snapshot_hash)
    require(sha(frozen['config_path']) == config_hash, 'Frozen core configuration hash differs')
    cfg = json.loads(frozen['config_path'].read_text())
    M.validate_execution_walltime(frozen,cfg,16,7200,nchunks=500)
    source_hashes = {str(path):sha(path) for path in (helper,Path(M.__file__),Path(D.__file__))}
    frozen_hashes = {str(snapshot/'SHA256SUMS'):snapshot_hash,str(frozen['config_path']):config_hash}
    evidence = {'created_utc':utc(),'job_id':parent,'plan_id':plan_id,'array_task_ids':members,
        'original_ledger_sha256':hashlib.sha256(original).hexdigest(),'source_sha256':source_hashes,
        'frozen_input_sha256':frozen_hashes,
        'snapshot_sha256':snapshot_hash,'config_sha256':config_hash,'old_array_task_throttle':3,'new_array_task_throttle':2,
        'scientific_changes':False,'walltime_changes':False,'receipt_identity_changes':False,
        'queue_command':queue_argv,'control_read_command':control_argv,
        'scope':'Only this exact accepted array, all16 members observed PENDING. The sole update reduces throttle3 to2. A start race never terminates a running task.',
        'reservation_note':'Potential16CPU reservation reduction only. Fresh dispatcher queue accounting must use max(running,throttle); this helper does not claim capacity is available or submit a final queue-fill array.',
        'execution_scope':'Separate from the core1h pilot. A final stage-only proposal uses the previously reviewed1h policy; this throttle adjustment does not submit work or verify signal delivery.'}
    evidence_path.parent.mkdir(parents=True,exist_ok=True)
    query_started = time.monotonic()
    try:
        evidence['queue_before'] = query(queue_argv)
        evidence['parsed_queue_before'] = parse_queue(evidence['queue_before'])
        evidence['control_before'] = query(control_argv)
        evidence['parsed_control_before'] = parse_control(evidence['control_before'],3)
        # Recheck after controller query; don't proceed if a member has started.
        evidence['queue_preupdate'] = query(queue_argv)
        evidence['parsed_queue_preupdate'] = parse_queue(evidence['queue_preupdate'])
    except (ValueError,subprocess.TimeoutExpired,subprocess.CalledProcessError,OSError) as exc:
        evidence.update(status='skipped_no_change',error=str(exc),exception_type=type(exc).__name__,
                        subprocess_invoked=False,finished_utc=utc())
        durable(evidence_path,evidence)
        print(json.dumps({'status':evidence['status'],'evidence':str(evidence_path)}))
        sys.exit(0)
    argv = ['scontrol','update','JobId=10264903','ArrayTaskThrottle=2']
    evidence.update(status='update_unknown',argv=argv,subprocess_start_utc=utc(),subprocess_invoked=None)
    durable(evidence_path,evidence)  # Persist unknown before the sole mutation.
    unchanged = ledger.read_bytes() == original and all(sha(p) == digest for p,digest in (source_hashes | frozen_hashes).items())
    query_age = time.monotonic()-query_started
    if not 0 <= query_age <= 60 or not unchanged:
        evidence.update(status='skipped_stale_or_changed_no_change',query_age_s=query_age,sources_and_ledger_unchanged=unchanged,
                        subprocess_invoked=False,finished_utc=utc())
        durable(evidence_path,evidence)
        print(json.dumps({'status':evidence['status'],'evidence':str(evidence_path)}))
        sys.exit(0)
    evidence.update(query_age_s=query_age,subprocess_invoked=True)
    try:
        response = subprocess.run(argv,capture_output=True,text=True,timeout=30,check=False)
        evidence.update(returncode=response.returncode,stdout=response.stdout,stderr=response.stderr,status='needs_reconciliation')
        evidence['control_after'] = query(control_argv)
        evidence['parsed_control_after'] = parse_control(evidence['control_after'],2,pending=False)
        evidence['queue_after'] = query(queue_argv)
        evidence['parsed_queue_after'] = parse_queue(evidence['queue_after'],pending=False)
        evidence.update(status='verified' if response.returncode == 0 else 'needs_reconciliation',finished_utc=utc())
    except (ValueError,subprocess.TimeoutExpired,subprocess.CalledProcessError,OSError) as exc:
        evidence.update(error=str(exc),exception_type=type(exc).__name__,finished_utc=utc())
    durable(evidence_path,evidence)
    require(ledger.read_bytes() == original, 'Ledger changed despite bound lock; preserve evidence and reconcile, do not overwrite')
    receipt.setdefault('execution_adjustments',[]).append({'kind':kind,'evidence_file':str(evidence_path),
        'job_id':parent,'array_task_ids':members,'old_array_task_throttle':3,'new_array_task_throttle':2,
        'state':evidence['status'],'utc':evidence['finished_utc'],'receipt_identity_unchanged':True})
    require(D.submission_identity(receipt) == identity_before and M.I.digest(identity_before) == plan_id,
            'Adjustment unexpectedly changes original receipt identity')
    durable(ledger,data)
    require(evidence['status'] == 'verified', 'Inspect evidence; no automatic retry of uncertain scheduler update')
    print(json.dumps({'status':'verified','evidence':str(evidence_path),'job_id':parent,'requested_throttle':2}))
PY
