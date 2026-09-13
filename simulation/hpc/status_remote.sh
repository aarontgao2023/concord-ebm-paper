#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
squeue -u "$USER" -h -o '%.18i %.16j %.9T %.10M %.6C %R'
"$HOME/ebmcal-env/bin/python" - <<'PY'
import json
from pathlib import Path
root=Path('/N/slate/tg11/ebmcal_v2')
for p in sorted((root/'runs').glob('*/cell_*/chunk_*/progress.json')):
 d=json.loads(p.read_text()); print(p.relative_to(root),json.dumps(d))
for p in sorted((root/'logs').glob('*.err')):
 if p.stat().st_size: print(str(p),p.read_text()[-1200:])
PY
