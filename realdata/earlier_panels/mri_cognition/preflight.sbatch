#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=sharedpanel-check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=00:45:00
#SBATCH --array=0-1%2
#SBATCH --output=logs/%x_%A_%a.out
# MRI-cognition populations, step 1 of 3 on the cluster: input check and B = 8 serial/parallel
# comparison (task 0 ADNI, task 1 NACC). Required before n2.sbatch and apoe.sbatch.
# Required environment:
#   REPO              checkout of this repository; realdata/earlier_panels/mri_cognition/ must hold
#                     inputs/ and campaign.json written by freeze_campaign.py
#   CONCORD_WORK_DIR  output folder (writes $CONCORD_WORK_DIR/mri_cognition/preflight)
#   CONCORD_PYTHON    python of the analysis environment (pyebm 2.0.3; concord is vendored)
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
  --mode preflight --cohort "${COHORTS[$SLURM_ARRAY_TASK_ID]}" \
  --workers "$SLURM_CPUS_PER_TASK" --out "$WORK/mri_cognition/preflight"
