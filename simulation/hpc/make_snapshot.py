"""Copy the simulation code and one run configuration into a new, checksummed snapshot directory.

A run is executed from a snapshot so that the code and configuration cannot change while it runs.
The snapshot holds scripts/ (with scripts/dev/), configs/NAME.json and hpc/*.sbatch in the layout of
this directory, plus SHA256SUMS, which run_array.sbatch and saebm_mechanism_freeze.sbatch verify
before running, and snapshot.json. An existing snapshot is never overwritten.

Usage: python simulation/hpc/make_snapshot.py NAME [--snapshot-root DIR]
NAME is a configuration in simulation/configs. The snapshot is written to DIR/NAME; the default DIR is
$CONCORD_WORK_DIR/snapshots, or runs/snapshots in the repository if CONCORD_WORK_DIR is not set.
"""
import argparse, hashlib, json, os, shutil
from pathlib import Path
from datetime import datetime, timezone
p = argparse.ArgumentParser(description='Copy the simulation code and one configuration into a new snapshot.')
p.add_argument('name', help='configuration name (simulation/configs/NAME.json); also the snapshot name')
p.add_argument('--snapshot-root', type=Path, default=None,
               help='directory for snapshots (default: $CONCORD_WORK_DIR/snapshots, else runs/snapshots in the repository)')
args = p.parse_args()
sim = Path(__file__).resolve().parents[1]
if args.snapshot_root is not None:
    snapshot_root = args.snapshot_root
elif os.environ.get('CONCORD_WORK_DIR'):
    snapshot_root = Path(os.environ['CONCORD_WORK_DIR']) / 'snapshots'
else:
    snapshot_root = sim.parent / 'runs' / 'snapshots'
out = snapshot_root / args.name
config = sim / 'configs' / f'{args.name}.json'
if not config.exists(): raise SystemExit(f'No configuration {config.name} in {config.parent}')
if out.exists(): raise SystemExit('Snapshot already exists; choose a new name')
files = [*sorted((sim / 'scripts').rglob('*.py')), config, *sorted((sim / 'hpc').glob('*.sbatch'))]
for f in files:
    dest = out / f.relative_to(sim); dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(f, dest)
checks = []
for f in files:
    rel = f.relative_to(sim); checks.append(f'{hashlib.sha256((out / rel).read_bytes()).hexdigest()}  {rel.as_posix()}')
(out / 'SHA256SUMS').write_text('\n'.join(checks) + '\n')
(out / 'snapshot.json').write_text(json.dumps({'name': args.name, 'created_utc': datetime.now(timezone.utc).isoformat(),
                                               'sha256_manifest': hashlib.sha256((out / 'SHA256SUMS').read_bytes()).hexdigest()}, indent=2) + '\n')
print(out)
