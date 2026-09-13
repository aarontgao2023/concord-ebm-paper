#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=tk-stdcons
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=/N/slate/tg11/ebmcal_v2/logs/stdcons_%A_%a.out
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
export EBM_SNAPSHOT=/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2
MODE=${MODE:?}; NCHUNKS=${NCHUNKS:?}; DATASETS=${DATASETS:?}; BPERM=${BPERM:-199}
exec ~/ebmcal-env/bin/python /N/slate/tg11/ebmcal_v2/dev/standardized_consensus_dev.py --mode $MODE \
  --datasets $DATASETS --bperm $BPERM --workers 16 --chunk ${SLURM_ARRAY_TASK_ID} --nchunks $NCHUNKS \
  --output /N/slate/tg11/ebmcal_v2/dev/stdcons_$MODE
