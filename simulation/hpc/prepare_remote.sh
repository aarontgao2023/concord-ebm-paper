#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
mkdir -p "$root/snapshots" "$root/logs" "$root/runs" "$root/analysis"
test -x "$HOME/ebmcal-env/bin/python"
