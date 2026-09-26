#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=nacc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-199
#SBATCH --output=logs/%x_%A_%a.out
# NACC amyloid PET population: two groups drawn from e3/e3 carriers (Fig. 4e, NACC amyloid
# PET-cognition setting). One draw per array task; results go to $CONCORD_WORK_DIR/earlier/nacc/n2/.
# Required environment:
#   REPO              checkout of this repository
#   TABLE             nacc_p5_deid.csv written by build_cohort_nacc.py
#   CONCORD_WORK_DIR  output folder
#   CONCORD_PYTHON    python of an environment with concord-ebm (commit 40fc046) and pyebm 2.0.3
#   STAGE             n2 (default) or a0 (descriptive counts, no fits; submit without --array)
# Submit from a directory that contains logs/.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
REPO=${REPO:?set REPO to the repository checkout}
TABLE=${TABLE:?set TABLE to nacc_p5_deid.csv}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
STAGE=${STAGE:-n2}
SCRIPTS="$REPO/realdata/earlier_panels/nacc"
sha256sum "$TABLE" "$SCRIPTS/run_nacc.py"
ARGS=""
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then ARGS="--replicates ${SLURM_ARRAY_TASK_ID}:$((SLURM_ARRAY_TASK_ID+1))"; fi
exec "$PYTHON" "$SCRIPTS/run_nacc.py" --table "$TABLE" --out "$WORK/earlier/nacc" \
  --stage $STAGE --workers ${SLURM_CPUS_PER_TASK:-1} $ARGS
