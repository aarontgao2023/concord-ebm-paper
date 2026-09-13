"""Freeze runtime and configuration before any remote execution."""
import argparse, hashlib, json, shutil
from pathlib import Path
from datetime import datetime, timezone
p=argparse.ArgumentParser(); p.add_argument('name'); args=p.parse_args()
root=Path(__file__).resolve().parents[2]
out=root/'runs'/'v2'/'snapshots'/args.name
if out.exists(): raise SystemExit('Snapshot already exists; choose a new name')
files=[*sorted((root/'scripts'/'v2').glob('*.py')), root/'configs'/'v2'/f'{args.name}.json', *sorted((root/'hpc'/'v2').glob('*.sbatch')), root/'hpc'/'v2'/'inventory_run.py']
for f in files:
    target=out/f.relative_to(root); target.parent.mkdir(parents=True,exist_ok=True); shutil.copy2(f,target)
checks=[]
for f in files:
    rel=f.relative_to(root); checks.append(f'{hashlib.sha256((out/rel).read_bytes()).hexdigest()}  {rel}')
(out/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
(out/'snapshot.json').write_text(json.dumps({'name':args.name,'created_utc':datetime.now(timezone.utc).isoformat(),'sha256_manifest':hashlib.sha256((out/'SHA256SUMS').read_bytes()).hexdigest()},indent=2)+'\n')
print(out)
