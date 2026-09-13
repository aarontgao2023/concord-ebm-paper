#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=ebmv2-fastval
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
set -euo pipefail
snapshot=$1
root=/N/slate/tg11/ebmcal_v2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg PYTHONHASHSEED=0
cd "$snapshot"
sha256sum -c SHA256SUMS
"$HOME/ebmcal-env/bin/python" -m unittest discover -s scripts/v2 -p test_fast_likelihood_v2.py -q
"$HOME/ebmcal-env/bin/python" scripts/v2/benchmark_fast_likelihood_v2.py --full-fit --events 4 14 --modes original repaired --micro-calls 300 --output "$root/analysis/fast_likelihood_hpc_observed.json"
"$HOME/ebmcal-env/bin/python" scripts/v2/benchmark_fast_likelihood_v2.py --full-fit --events 14 --modes original repaired --micro-calls 100 --permuted --output "$root/analysis/fast_likelihood_hpc_permuted.json"
