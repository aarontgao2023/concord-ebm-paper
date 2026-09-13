#!/bin/bash
set -euo pipefail
hostname
date -u '+%Y-%m-%dT%H:%M:%SZ'
df -h /N/slate /N/u/tg11/Quartz
squeue -u "$USER" -h -o '%.18i %.16j %.9T %.10M %.6C %R'
sacctmgr -n show assoc user="$USER" format=Account,Partition,GrpCPUs,GrpTRESMins
"$HOME/ebmcal-env/bin/python" - <<'PY'
import sys, importlib.metadata, inspect, hashlib, json
import pyebm.central_ordering.generalized_mallows as gm
import numpy, pandas, scipy
p = inspect.getsourcefile(gm)
print(json.dumps({'python':sys.version, 'pyebm':importlib.metadata.version('pyebm'),
 'numpy':numpy.__version__, 'pandas':pandas.__version__, 'scipy':scipy.__version__,
 'consensus_file':p, 'consensus_source_sha256':hashlib.sha256(open(p,'rb').read()).hexdigest()},indent=2))
print(inspect.getsource(gm.weighted_mallows.consensus))
PY
