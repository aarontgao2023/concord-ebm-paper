#!/bin/bash
set -euo pipefail
squeue -u "$USER" -h -o '%.18i %.16j %.9T %.10M %.6C %R'
"$HOME/ebmcal-env/bin/python" - <<'PY'
import json
from pathlib import Path
from collections import Counter
root=Path('/N/slate/tg11/ebmcal_v2')
for run in sorted((root/'runs').iterdir()):
 ps=list(run.glob('cell_*/chunk_*/progress.json'))
 if not ps: continue
 ds=[json.loads(p.read_text()) for p in ps]
 print(json.dumps({'run':run.name,'started_chunks':len(ds),'done_chunks':sum(d['done'] for d in ds),'fit_records':sum(d['fit_records'] for d in ds),'rows':sum(d['rows'] for d in ds),'statuses':dict(sum((Counter(d['statuses']) for d in ds),Counter())),'latest_utc':max(d['updated_utc'] for d in ds)}))
for p in sorted((root/'logs').glob('*.err')):
 if p.stat().st_size: print(str(p),p.read_text()[-800:])
PY
