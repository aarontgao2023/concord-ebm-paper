#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=adni-null
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-199
#SBATCH --output=logs/%x_%A_%a.out
# ADNI 12-event population: two groups drawn from e3/e3 carriers, 12 events (Fig. 4e, ADNI
# CSF-MRI-cognition 12-event setting). One draw per array task.
# Required environment:
#   REPO              checkout of this repository
#   TABLE             adni_p12_deid.csv written by build_cohort.py
#   CONCORD_WORK_DIR  output folder; draws go to $CONCORD_WORK_DIR/earlier/adni/null/n2/
#   CONCORD_PYTHON    python of an environment with concord-ebm (commit 40fc046) and pyebm 2.0.3
#   PER               draws per array task (default 1)
# Submit from a directory that contains logs/.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
REPO=${REPO:?set REPO to the repository checkout}
TABLE=${TABLE:?set TABLE to adni_p12_deid.csv}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
PER=${PER:-1}
A=$(( SLURM_ARRAY_TASK_ID * PER )); B=$(( A + PER ))
SCRIPTS="$REPO/realdata/earlier_panels/adni"
sha256sum "$TABLE" "$SCRIPTS/run_null.py" "$SCRIPTS/run_observed.py"
exec "$PYTHON" "$SCRIPTS/run_null.py" --table "$TABLE" --out "$WORK/earlier/adni/null" \
  --stage n2 --replicates $A:$B --workers ${SLURM_CPUS_PER_TASK:-1}
