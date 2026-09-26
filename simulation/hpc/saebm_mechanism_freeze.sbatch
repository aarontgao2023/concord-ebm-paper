#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=saebm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# SA-EBM (pysaebm, conjugate priors, 10,000 iterations, 2,500 burn-in) on the complete-data
# fixed-sequence cells: run method_saebm_a (Fig. 2a,c). Seeds 53400000-53400999, the same seeds as
# method_saebm_reference_a and method_kde_mechanism_a. This run has no JSON configuration; its
# settings are the defaults below.
#
# Usage (from a snapshot made by make_snapshot.py; any configuration name will do):
#   export CONCORD_WORK_DIR=/path/to/work SAEBM_PYTHON=/path/to/saebm-env/bin/python NCHUNKS=50
#   sbatch --array=0-49 SNAPSHOT/hpc/saebm_mechanism_freeze.sbatch SNAPSHOT
# pysaebm runs in its own environment (SAEBM_PYTHON). Output: $CONCORD_WORK_DIR/runs/method_saebm_a,
# one rows_chunkNNN.jsonl per array task; resubmitting a task skips the datasets already fitted.
# Edit the account, partition and resources for your cluster.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
SNAP=$(cd "${1:?usage: sbatch --array=0-(NCHUNKS-1) saebm_mechanism_freeze.sbatch SNAPSHOT}" && pwd)
ROOT=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the directory for run outputs}
PYTHON=${SAEBM_PYTHON:-/path/to/saebm-env/bin/python}
cd "$SNAP" && sha256sum -c --quiet SHA256SUMS
export EBM_SNAPSHOT="$SNAP/scripts"
NCHUNKS=${NCHUNKS:?}; DATASETS=${DATASETS:-1000}; SEED_START=${SEED_START:-53400000}
CELLS=${CELLS:-IID_S,REF_S,BAL_S,REF_S_N4}; N_ITER=${N_ITER:-10000}; BURN_IN=${BURN_IN:-2500}
exec "$PYTHON" "$SNAP/scripts/dev/saebm_mechanism_dev.py" --cells "$CELLS" --n-iter "$N_ITER" --burn-in "$BURN_IN" \
  --seed-start "$SEED_START" --datasets "$DATASETS" --workers 16 --chunk "${SLURM_ARRAY_TASK_ID}" --nchunks "$NCHUNKS" \
  --output "$ROOT/runs/method_saebm_a"
