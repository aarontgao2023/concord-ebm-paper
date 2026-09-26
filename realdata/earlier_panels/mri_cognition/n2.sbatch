#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=sharedpanel-n2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=0-399%8
#SBATCH --output=logs/%x_%A_%a.out
# MRI-cognition populations, step 2 of 3 on the cluster: two groups drawn from e3/e3 carriers
# (Fig. 4e, ADNI and NACC MRI-cognition rows). 400 tasks = 2 cohorts x 2 panels (4 and 8 events)
# x 2 models (standard = separately fitted DEBM, invariant_min = CONCORD) x 50 blocks of 4 draws.
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
SLOT=$((SLURM_ARRAY_TASK_ID % 8))
FIRST=$((SLURM_ARRAY_TASK_ID / 8 * 4))
LAST=$((FIRST + 4))
COHORTS=(adni nacc)
PANELS=(MRIcog4 MRIcog8)
ESTIMATORS=(standard invariant_min)
cd "$WORK"
exec "$PYTHON" "$RUNNER" \
  --mode n2 --cohort "${COHORTS[$((SLOT / 4))]}" \
  --panel "${PANELS[$(((SLOT % 4) / 2))]}" --estimator "${ESTIMATORS[$((SLOT % 2))]}" \
  --start "$FIRST" --stop "$LAST" --workers "$SLURM_CPUS_PER_TASK" \
  --preflight-root "$WORK/mri_cognition/preflight" --out "$WORK/mri_cognition/results"
