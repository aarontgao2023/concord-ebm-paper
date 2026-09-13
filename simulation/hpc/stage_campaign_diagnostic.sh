#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
from pathlib import Path
from datetime import datetime,timezone
import hashlib,subprocess,sys,json
root=Path('/N/slate/tg11/ebmcal_v2')
checks={'analyze_stage_campaign.py':'abc81c60efd4726d0e1a6e2ec6c94c04458df1dad6df162022e8a4b63de7ef82','analyze_stage_exchangeability.py':'e02acc932af7e15e606730f81c41824786aa8e5bd1a168f55b66bbecd95c1ac7'}
for name,digest in checks.items():assert hashlib.sha256((root/'toolkit/scripts/v2'/name).read_bytes()).hexdigest()==digest
assert (root/'runs/v2/snapshots').resolve()==root/'snapshots'
out=root/'diagnostics/stage_campaign'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
out.mkdir(parents=True,exist_ok=False)
subprocess.run([sys.executable,str(root/'toolkit/scripts/v2/analyze_stage_campaign.py'),'--protocol',str(root/'control/campaign_versions/confirmation_20260908_1455/confirmation_protocol.json'),'--workspace',str(root),'--raw-root',str(root/'runs'),'--output-json',str(out/'stage_campaign.json'),'--output-md',str(out/'stage_campaign.md')],check=True)
d=json.loads((out/'stage_campaign.json').read_text())
print(json.dumps({'output':str(out),**{k:d[k] for k in ('status','planned_datasets','verified_datasets','pending_datasets','unavailable_datasets','error_datasets','errors')},'cohort_statistics_available':d['cohort_summary']['available']}))
PY
