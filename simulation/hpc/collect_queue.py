"""Read this user's complete Slurm queue and normalize it for campaign planning.

No scheduler mutations. Unregistered ebmv2 jobs are errors, not free capacity.
Use on Quartz with an explicit local copy of the project job registry/ledger.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess

ROOT = Path('/N/slate/tg11/ebmcal_v2')


def command(argv):
    return subprocess.run(argv, check=True, text=True, capture_output=True,
                          timeout=30).stdout


def registrations(registry, ledger):
    registered = {}
    for name, entry in registry.items():
        if not entry.get('job_id'):
            continue
        snapshot = ROOT / 'snapshots' / Path(entry['snapshot']).name
        digest = json.loads((snapshot/'snapshot.json').read_text())['sha256_manifest']
        registered[str(entry['job_id'])] = {
            'run_id': name, 'snapshot_sha256': digest,
            'max_concurrent': entry.get('current_throttle'),
        }
    for entry in ledger.get('submissions', []):
        if not entry.get('job_id'):
            continue
        candidate = {'run_id': entry['run_id'], 'snapshot_sha256': entry['snapshot_sha256'],
                     'max_concurrent': entry['max_concurrent']}
        previous = registered.get(str(entry['job_id']))
        if previous and any(previous[k] != candidate[k] for k in ('run_id', 'snapshot_sha256')):
            raise ValueError('Registry and ledger disagree about a submitted job')
        registered[str(entry['job_id'])] = candidate
    return registered


def collect(registry, ledger):
    known = registrations(registry, ledger)
    raw = command(['squeue', '-r', '-u', os.environ['USER'], '-h',
                   '-o', '%i|%j|%T|%C|%Z|%k'])
    grouped = defaultdict(list)
    for line in raw.splitlines():
        parts = line.split('|')
        if len(parts) != 6:
            raise ValueError('Unrecognized squeue record')
        jid, name, state, cpus, workdir, comment = (x.strip() for x in parts)
        match = re.fullmatch(r'(\d+)(?:_(\d+))?', jid)
        if not match:
            raise ValueError(f'Unexpanded or unsupported Slurm job identity: {jid}')
        parent, task = match.groups()
        own_path = workdir == str(ROOT) or workdir.startswith(str(ROOT)+'/')
        own = parent in known or own_path or name.startswith('ebmv2')
        if own and parent not in known:
            raise ValueError(f'Unregistered project job {parent}; reconcile its submission first')
        grouped[parent].append({'job_id': jid, 'task': None if task is None else int(task),
                                'state': state, 'cpus': int(cpus), 'own': own,
                                'name': name, 'workdir': workdir, 'comment': comment})
    jobs = []
    for parent, rows in sorted(grouped.items()):
        if len({r['cpus'] for r in rows}) != 1 or len({r['own'] for r in rows}) != 1:
            raise ValueError('Unexpected mixed CPU/project array')
        tasks = [r['task'] for r in rows]
        is_array = all(t is not None for t in tasks)
        if not is_array and (len(rows) != 1 or tasks != [None]):
            raise ValueError('Mixed parent/task rows must not be double counted')
        states = {r['state'] for r in rows}
        entry = {'job_id': parent, 'state': next(iter(states)) if len(states) == 1 else 'MIXED',
                 'project': 'ebmcal_v2' if rows[0]['own'] else 'other_user_project',
                 'cpus_per_task': rows[0]['cpus'],
                 'array_task_ids': sorted(tasks) if is_array else None}
        if is_array:
            detail = command(['scontrol', 'show', 'job', rows[0]['job_id']])
            match = re.search(r'\bArrayTaskThrottle=(\d+)\b', detail)
            if not match:
                raise ValueError('Cannot determine the current array throttle')
            entry.update(max_concurrent=int(match[1]) or len(tasks),
                         running_tasks=sum(r['state'] != 'PENDING' for r in rows))
            if rows[0]['own']:
                entry.update(run_id=known[parent]['run_id'],
                             snapshot_sha256=known[parent]['snapshot_sha256'])
        jobs.append(entry)
    return {'schema_version': 'campaign_queue_v1', 'captured_utc': datetime.now(timezone.utc).isoformat(),
            'scope': 'all_user_jobs', 'complete': True, 'jobs': jobs,
            'source': 'squeue -r (all user jobs), plus scontrol current array throttle'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--registry', type=Path, required=True)
    ap.add_argument('--ledger', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    result = collect(json.loads(args.registry.read_text()), json.loads(args.ledger.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix('.tmp')
    with tmp.open('w') as f:
        json.dump(result, f, indent=2)
        f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, args.output)
    print(json.dumps({'output': str(args.output), 'jobs': len(result['jobs'])}))


if __name__ == '__main__':
    main()
