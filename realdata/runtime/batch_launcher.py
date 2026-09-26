"""Run a frozen empirical batch without exceeding its 16 allocated CPU slots."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mapping', type=Path, required=True)
    ap.add_argument('--index', type=int, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--gate', type=Path, required=True)
    args = ap.parse_args()
    gate = json.loads(args.gate.read_text())
    mapping_hash = hashlib.sha256(args.mapping.read_bytes()).hexdigest()
    if mapping_hash != gate['batch_mapping_sha256']:
        raise ValueError('Batch mapping differs from the approved frozen gate')
    batch = json.loads(args.mapping.read_text())[args.index]
    if not batch or len(batch)>16:
        raise ValueError('Batch must contain one to sixteen work units')
    if any(x['mode']=='apoe' for x in batch) and len(batch)!=1:
        raise ValueError('A parallel APOE fit family occupies its own batch')
    allocated = int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))
    if allocated < 16:
        raise ValueError('This frozen batch requires sixteen allocated CPUs')
    def run(item):
        cmd = [sys.executable, '-B', str(HERE/'run_common_panel.py'), '--config', str(HERE/'campaign.json'),
               '--mode', item['mode'], '--cohort', item['cohort'], '--out', str(args.out), '--gate', str(args.gate),
               '--workers', '16' if item['mode']=='apoe' else '1']
        if item['mode']=='apoe':
            cmd += ['--panel', item['panel']]
        else:
            cmd += ['--start', str(item['replicate']), '--stop', str(item['replicate']+1)]
        result = subprocess.run(cmd, check=False)
        return {'unit':item, 'returncode':result.returncode}
    outcomes=[]
    with ThreadPoolExecutor(max_workers=len(batch)) as executor:
        for f in as_completed([executor.submit(run,item) for item in batch]):
            outcome=f.result(); outcomes.append(outcome)
            print(json.dumps({'batch_index':args.index, **outcome}), flush=True)
    if any(x['returncode'] for x in outcomes):
        raise SystemExit('At least one frozen work unit failed; its records remain retained')

if __name__=='__main__':
    main()
