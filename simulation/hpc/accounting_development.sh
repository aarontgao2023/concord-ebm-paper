#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
sacct -X -n -P -j 10254078,10254376,10255861,10258285,10258286,10260447,10260448 --format=JobID,State,ExitCode,ElapsedRaw,AllocCPUS,Start,End > "$root/analysis/accounting_development.psv"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$root/analysis/accounting_development_captured_utc.txt"
