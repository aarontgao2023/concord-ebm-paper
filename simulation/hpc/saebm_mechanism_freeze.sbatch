#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=tk-saebm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/N/slate/tg11/ebmcal_v2/logs/saebm_frz_%A_%a.out
# Second confirmation protocol, SA-EBM arm (method_saebm_a): pysaebm conjugate priors on the
# complete-data fixed-ordering cells, fresh seeds 53.4M, matched seed-for-seed to
# method_saebm_reference_a (run_v2). The script is the frozen copy in the snapshot directory;
# design_v2 comes from the confirm_core_a snapshot as in the development run.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
export EBM_SNAPSHOT=/N/slate/tg11/ebmcal_v2/snapshots/confirm_core_a/scripts/v2
SNAP=/N/slate/tg11/ebmcal_v2/snapshots/method_saebm_a
cd $SNAP && sha256sum -c --quiet SHA256SUMS
NCHUNKS=${NCHUNKS:?}; DATASETS=${DATASETS:-1000}; SEED_START=${SEED_START:-53400000}
CELLS=${CELLS:-IID_S,REF_S,BAL_S,REF_S_N4}; N_ITER=${N_ITER:-10000}; BURN_IN=${BURN_IN:-2500}
exec ~/saebm-env/bin/python $SNAP/saebm_mechanism_dev.py --cells $CELLS --n-iter $N_ITER --burn-in $BURN_IN \
  --seed-start $SEED_START --datasets $DATASETS --workers 16 --chunk ${SLURM_ARRAY_TASK_ID} --nchunks $NCHUNKS \
  --output /N/slate/tg11/ebmcal_v2/runs/method_saebm_a
