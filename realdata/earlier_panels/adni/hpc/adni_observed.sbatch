#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=adni-obs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
# ADNI 12-event population: 500 bootstrap refits within genotype-by-diagnosis cells
# (Supplementary Fig. 3, ADNI panel). Writes $CONCORD_WORK_DIR/earlier/adni/main/a2s_core_P12_w180_DXSUM.json.
# Required environment:
#   REPO, TABLE, CONCORD_WORK_DIR, CONCORD_PYTHON   as in adni_null.sbatch
#   STAGE   a2s (default) or a0 (descriptive counts, no fits)
#   EXTRA   further run_observed.py options (default "--R 500")
# Submit from a directory that contains logs/.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
REPO=${REPO:?set REPO to the repository checkout}
TABLE=${TABLE:?set TABLE to adni_p12_deid.csv}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
STAGE=${STAGE:-a2s}; EXTRA=${EXTRA:---R 500}
SCRIPTS="$REPO/realdata/earlier_panels/adni"
sha256sum "$TABLE" "$SCRIPTS/run_observed.py"
exec "$PYTHON" "$SCRIPTS/run_observed.py" --table "$TABLE" --out "$WORK/earlier/adni/main" \
  --stage $STAGE --workers ${SLURM_CPUS_PER_TASK:-1} $EXTRA
