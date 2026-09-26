#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=common6-formal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-55%16
#SBATCH --output=logs/%x_%A_%a.out
# Six-event ADNI/NACC panel (Figs 1c, 4a-d and 5; Tables 2-3; Supplementary Tables 4-5; Supplementary Fig. 2).
# One array task runs one batch of runtime/batch_mapping.json (56 batches): the controlled composition
# draws (mode n2), the APOE comparisons (mode apoe) and the bootstrap refits (mode bootstrap).
# Required environment:
#   REPO              checkout of this repository; realdata/runtime/inputs/ must hold the prepared inputs
#                     written by realdata/prepare_campaign.py
#   CONCORD_WORK_DIR  output directory (results are written to $CONCORD_WORK_DIR/formal)
#   CONCORD_PYTHON    python of the analysis environment (versions in realdata/runtime/EXECUTION_GATE.json)
# Submit from a directory that contains logs/.
set -euo pipefail
umask 077
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 LC_ALL=C
REPO=${REPO:?set REPO to the repository checkout}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output directory}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
RUNTIME="$REPO/realdata/runtime"
export MPLCONFIGDIR="$WORK/cache/mpl-$SLURM_JOB_ID"
mkdir -p "$MPLCONFIGDIR"
cd "$WORK"
exec "$PYTHON" -B "$RUNTIME/batch_launcher.py" --mapping "$RUNTIME/batch_mapping.json" --index "$SLURM_ARRAY_TASK_ID" \
  --out "$WORK/formal" --gate "$RUNTIME/EXECUTION_GATE.json"
