#!/bin/bash
#SBATCH -J pkg-dryrun
#SBATCH -A r00579
#SBATCH -p general
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH -o /N/slate/tg11/ebmcal_v2/dev/package/logs/dryrun_%j.out
# Packaged comparison dry run on a simulated REF cohort: invariant_min and standard,
# B=99, 20 stability resamples, one worker per allocated core. Reads its code from the package staging
# directory, not from a confirmation snapshot.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PKG=/N/slate/tg11/ebmcal_v2/dev/package
cd $PKG/scripts_v2
exec ~/ebmcal-env/bin/python dev/ebm_compare_dryrun.py --B ${B:-99} --stability ${STAB:-20} \
    --workers ${SLURM_CPUS_PER_TASK} --estimators ${ESTIMATORS:-invariant_min,standard} --seed ${SEED:-31900000} \
    --output $PKG/dryrun_${SLURM_JOB_ID}
