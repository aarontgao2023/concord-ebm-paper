#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=nacc-stability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x_%A_%a.out
# NACC amyloid PET population: 500 bootstrap refits of CONCORD (Supplementary Fig. 3, NACC panel).
# Array task 0 uses the smallest-proportion common proportions (reported), task 1 pooled
# proportions (not reported). Results go to $CONCORD_WORK_DIR/earlier/nacc_stability/<reference>/.
# Required environment:
#   REPO              checkout of this repository
#   TABLE             nacc_p5_deid.csv written by build_cohort_nacc.py
#   CONCORD_WORK_DIR  output folder
#   CONCORD_PYTHON    python of an environment with concord-ebm (commit 40fc046) and pyebm 2.0.3
# Submit from a directory that contains logs/.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONHASHSEED=0
REPO=${REPO:?set REPO to the repository checkout}
TABLE=${TABLE:?set TABLE to nacc_p5_deid.csv}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
REFERENCES=(min pooled)
REFERENCE=${REFERENCES[$SLURM_ARRAY_TASK_ID]}
SCRIPTS="$REPO/realdata/earlier_panels/nacc"
exec "$PYTHON" "$SCRIPTS/run_nacc_stability.py" --reference "$REFERENCE" --table "$TABLE" \
  --workers "$SLURM_CPUS_PER_TASK" --out "$WORK/earlier/nacc_stability/$REFERENCE"
