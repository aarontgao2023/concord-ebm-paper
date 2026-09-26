#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=sharedpanel-apoe
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-1%2
#SBATCH --output=logs/%x_%A_%a.out
# MRI-cognition populations, step 3 of 3 on the cluster: APOE comparisons on the 8-event panel
# (not reported in the paper, but part of the campaign; aggregate_completed_results.py requires
# these two result files). Task 0 ADNI, task 1 NACC.
# Required environment: REPO, CONCORD_WORK_DIR and CONCORD_PYTHON as in preflight.sbatch.
# Submit from a directory that contains logs/.
set -euo pipefail
umask 077
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0 LC_ALL=C
REPO=${REPO:?set REPO to the repository checkout}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
RUNNER="$REPO/realdata/earlier_panels/mri_cognition/run_shared_panel.py"
export MPLCONFIGDIR="$WORK/cache/mpl-$SLURM_JOB_ID"
mkdir -p "$MPLCONFIGDIR"
COHORTS=(adni nacc)
cd "$WORK"
exec "$PYTHON" "$RUNNER" \
  --mode apoe --cohort "${COHORTS[$SLURM_ARRAY_TASK_ID]}" \
  --panel MRIcog8 --estimator invariant_min --workers "$SLURM_CPUS_PER_TASK" \
  --preflight-root "$WORK/mri_cognition/preflight" --out "$WORK/mri_cognition/results"
