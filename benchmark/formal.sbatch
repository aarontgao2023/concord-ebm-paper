#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=concord-benchmark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=results/logs/%x_%A_%a.out
# Paired benchmark (Fig. 3b-d, Supplementary Table 2). One array task runs one batch of
# runtime/batch_mapping.json (up to 64 datasets, 16 at a time). Submit from benchmark/, e.g.
#   mkdir -p results/logs
#   CAMPAIGN=power      sbatch --array=0-47%16 formal.sbatch
#   CAMPAIGN=robustness sbatch --array=0-63%16 formal.sbatch
# Edit the account and partition for your cluster, and set CONCORD_PYTHON (Python of the analysis
# environment) and optionally CONCORD_BENCHMARK_OUT (output folder, default results/formal).
set -euo pipefail
umask 077
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 LC_ALL=C
BENCHMARK_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
RUNTIME_ROOT="$BENCHMARK_DIR/runtime"
OUTPUT_DIR="${CONCORD_BENCHMARK_OUT:-$BENCHMARK_DIR/results/formal}"
PYTHON="${CONCORD_PYTHON:-/path/to/analysis-env/bin/python}"
cd "$RUNTIME_ROOT"
: "${CAMPAIGN:?set CAMPAIGN to power or robustness}"
exec "$PYTHON" run_batch.py --campaign "$CAMPAIGN" --batch-index "$SLURM_ARRAY_TASK_ID" --workers 16 --output "$OUTPUT_DIR"
