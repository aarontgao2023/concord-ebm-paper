#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
from pathlib import Path
import subprocess,sys
root=Path('/N/slate/tg11/ebmcal_v2')
alias=root/'runs/v2/snapshots'
alias.parent.mkdir(parents=True,exist_ok=True)
if not alias.exists():alias.symlink_to(root/'snapshots',target_is_directory=True)
assert alias.resolve()==root/'snapshots'
out=root/'diagnostics/stage_first_seed_20260908b'
assert not out.exists(), 'Use a new diagnostic output path rather than overwrite'
out.mkdir(parents=True)
subprocess.run([sys.executable,str(root/'toolkit/scripts/v2/analyze_stage_exchangeability.py'),'--protocol',str(root/'control/campaign_versions/confirmation_20260908_1455/confirmation_protocol.json'),'--workspace',str(root),'--chunk-dir',str(root/'runs/confirm_stage_a/cell_0/chunk_0'),'--seed','42100000','--output-json',str(out/'stage_blocks_seed42100000.json'),'--output-md',str(out/'stage_blocks_seed42100000.md')],check=True)
print((out/'stage_blocks_seed42100000.md').read_text())
PY
