#!/bin/bash
set -euo pipefail
root=/N/slate/tg11/ebmcal_v2
snapshot=$root/snapshots/dev_paired_bridge_a
cd "$snapshot"
sha256sum -c SHA256SUMS
validation=$(sbatch --parsable --output="$root/logs/paired_validation_%j.out" --error="$root/logs/paired_validation_%j.err" hpc/v2/validate_paired.sbatch "$snapshot")
printf 'validation_job=%s\n' "$validation"
sbatch --parsable --dependency="afterok:$validation" --cpus-per-task=8 --mem=8G --array=0-5%6 --output="$root/logs/dev_paired_bridge_a_%A_%a.out" --error="$root/logs/dev_paired_bridge_a_%A_%a.err" hpc/v2/run_array.sbatch "$snapshot" dev_paired_bridge_a 2
