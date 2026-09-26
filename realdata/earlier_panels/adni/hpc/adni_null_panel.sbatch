#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=adni-nullp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-39
#SBATCH --output=logs/%x_%A_%a.out
# ADNI 12-event population, reduced 8- and 5-event panels (Fig. 4e, ADNI CSF-MRI-cognition 8- and
# 5-event settings): same participants, same adjustment and same draws as the 12-event setting.
# Required environment:
#   REPO, TABLE, CONCORD_WORK_DIR, CONCORD_PYTHON   as in adni_null.sbatch
#   PANEL             P8 or P5; draws go to $CONCORD_WORK_DIR/earlier/adni/null_$PANEL/n2/
#   PER               draws per array task (default 5; 40 tasks x 5 = 200 draws)
# Submit from a directory that contains logs/.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
REPO=${REPO:?set REPO to the repository checkout}
TABLE=${TABLE:?set TABLE to adni_p12_deid.csv}
WORK=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the output folder}
PYTHON=${CONCORD_PYTHON:-YOUR_ENV/bin/python}
PANEL=${PANEL:?set PANEL to P8 or P5}
PER=${PER:-5}
A=$(( SLURM_ARRAY_TASK_ID * PER )); B=$(( A + PER ))
SCRIPTS="$REPO/realdata/earlier_panels/adni"
sha256sum "$TABLE" "$SCRIPTS/run_null.py" "$SCRIPTS/run_observed.py"
exec "$PYTHON" "$SCRIPTS/run_null.py" --table "$TABLE" --out "$WORK/earlier/adni/null_$PANEL" \
  --stage n2 --panel $PANEL --replicates $A:$B --workers ${SLURM_CPUS_PER_TASK:-1}
