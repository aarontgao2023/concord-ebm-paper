#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=tk-saebm
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=/N/slate/tg11/ebmcal_v2/logs/saebm_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
export EBM_SNAPSHOT=/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2
NCHUNKS=${NCHUNKS:?}; DATASETS=${DATASETS:?}; CELLS=${CELLS:-IID_S,REF_S,BAL_S,REF_S_N4}; N_ITER=${N_ITER:-10000}; BURN_IN=${BURN_IN:-2500}
exec ~/saebm-env/bin/python /N/slate/tg11/ebmcal_v2/dev/saebm_mechanism_dev.py --cells $CELLS --n-iter $N_ITER --burn-in $BURN_IN \
  --datasets $DATASETS --workers 16 --chunk ${SLURM_ARRAY_TASK_ID} --nchunks $NCHUNKS \
  --output /N/slate/tg11/ebmcal_v2/dev/saebm_mechanism
