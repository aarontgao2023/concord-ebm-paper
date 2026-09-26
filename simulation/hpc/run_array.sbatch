#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=general
#SBATCH --job-name=concord-sim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --signal=B:USR1@1500
# One array task runs one chunk of one cell of a run configuration with run_v2.py.
# Array index = chunk * ncells + cell, so submit --array=0-$((NCELLS * NCHUNKS - 1)).
#
# Usage (from a snapshot made by make_snapshot.py):
#   export CONCORD_WORK_DIR=/path/to/work CONCORD_PYTHON=/path/to/env/bin/python
#   sbatch --array=0-249 SNAPSHOT/hpc/run_array.sbatch SNAPSHOT CONFIG_NAME NCHUNKS
# Output: $CONCORD_WORK_DIR/runs/CONFIG_NAME/cell_<cell>/chunk_<chunk>. The runner is resumable:
# resubmitting a task continues its chunk. Edit the account, partition and resources for your cluster.
set -euo pipefail
snapshot=$1
config_name=$2
nchunks=$3
root=${CONCORD_WORK_DIR:?set CONCORD_WORK_DIR to the directory for run outputs}
python=${CONCORD_PYTHON:-/path/to/concord-env/bin/python}
index=${SLURM_ARRAY_TASK_ID:?}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg PYTHONHASHSEED=0
cd "$snapshot"
sha256sum -c SHA256SUMS
ncells=$("$python" -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["cells"]))' "configs/${config_name}.json")
cell=$(( index % ncells ))
chunk=$(( index / ncells ))
test "$chunk" -lt "$nchunks"
exec "$python" scripts/run_v2.py \
  --config "configs/${config_name}.json" --cell-index "$cell" \
  --chunk "$chunk" --nchunks "$nchunks" --workers "${SLURM_CPUS_PER_TASK}" \
  --time-budget-s 4000 --output "$root/runs/$config_name/cell_${cell}/chunk_${chunk}"
