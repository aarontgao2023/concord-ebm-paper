#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
"$HOME/ebmcal-env/bin/python" "$root/toolkit/hpc/v2/reconcile_ledger.py" --ledger "$root/control/submission_ledger.json"
