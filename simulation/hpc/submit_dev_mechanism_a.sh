#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
snapshot=$root/snapshots/dev_mechanism_a
cd "$snapshot"
sha256sum -c SHA256SUMS
sbatch --parsable --array=0-7%8 --output="$root/logs/dev_mechanism_a_%A_%a.out" --error="$root/logs/dev_mechanism_a_%A_%a.err" hpc/v2/run_array.sbatch "$snapshot" dev_mechanism_a 2
