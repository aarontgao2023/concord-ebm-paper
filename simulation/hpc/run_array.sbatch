#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=ebmv2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --signal=B:USR1@1500
set -euo pipefail
snapshot=$1
config_name=$2
nchunks=$3
root=/N/slate/tg11/ebmcal_v2
index=${SLURM_ARRAY_TASK_ID:?}

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg PYTHONHASHSEED=0
cd "$snapshot"
sha256sum -c SHA256SUMS
ncells=$("$HOME/ebmcal-env/bin/python" -c 'import json,sys; print(len(json.load(open(sys.argv[1]))["cells"]))' "configs/v2/${config_name}.json")
cell=$(( index % ncells ))
chunk=$(( index / ncells ))
test "$chunk" -lt "$nchunks"
exec "$HOME/ebmcal-env/bin/python" scripts/v2/run_v2.py \
  --config "configs/v2/${config_name}.json" --cell-index "$cell" \
  --chunk "$chunk" --nchunks "$nchunks" --workers "${SLURM_CPUS_PER_TASK}" \
  --time-budget-s 4000 --output "$root/runs/$config_name/cell_${cell}/chunk_${chunk}"
