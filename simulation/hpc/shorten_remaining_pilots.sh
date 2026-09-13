#!/bin/bash
set -euo pipefail
"$HOME/ebmcal-env/bin/python" - <<'PY'
import re,subprocess
for parent in ('10255861','10258286'):
 result=subprocess.run(['squeue','-r','-j',parent,'-h','-o','%i %T'],check=True,text=True,capture_output=True)
 for line in result.stdout.splitlines():
  job,state=line.split()
  if state=='PENDING' and re.fullmatch(parent+r'_[0-9]+',job):
   subprocess.run(['scontrol','update',f'JobId={job}','TimeLimit=00:30:00'],check=True)
   print(job,'pending walltime set to30min based on completed small-chunk throughput')
PY
