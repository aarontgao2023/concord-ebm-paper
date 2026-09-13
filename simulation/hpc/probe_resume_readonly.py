"""Read-only identity/environment probe for an existing ebmcal_v2 HPC root.

Run through the already authenticated SSH master with Python stdin. It never
creates a ledger, changes a snapshot, fits a model or contacts the scheduler.
"""
import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

root = Path('/N/slate/tg11/ebmcal_v2')
ledger = root / 'control/submission_ledger.json'
protocol = root / 'control/campaign_versions/confirmation_20260908_1455/confirmation_protocol.json'
names = ('campaign_cycle.py', 'reconcile_ledger.py', 'accounted_members.py',
         'collect_queue.py', 'manage_campaign.py', 'dispatch_campaign.py', 'inventory_run.py')
def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

value = json.loads(ledger.read_text())
print(json.dumps({
    'captured_utc': datetime.now(timezone.utc).isoformat(),
    'root': str(root), 'stop_file_exists': (root / 'STOP').exists(),
    'ledger_sha256': digest(ledger), 'receipt_count': len(value['submissions']),
    'protocol_sha256': digest(protocol),
    'management_source_sha256': {name: digest(root / 'toolkit/hpc/v2' / name) for name in names},
    'snapshot_manifests_sha256': {p.parent.name: digest(p) for p in sorted((root / 'snapshots').glob('*/SHA256SUMS'))},
    'python_executable': sys.executable, 'python_version': sys.version,
    'distributions': {name: importlib.metadata.version(name) for name in ('numpy', 'scipy', 'pyebm')},
    'read_only': True,
}, indent=2))
