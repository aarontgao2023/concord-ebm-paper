#!/bin/bash
#SBATCH --account=r00579
#SBATCH --partition=general
#SBATCH --job-name=ebmv2-validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -euo pipefail
snapshot=$1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
cd "$snapshot"
sha256sum -c SHA256SUMS
"$HOME/ebmcal-env/bin/python" -m unittest discover -s scripts/v2 -p test_paired_engine_v2.py -q
"$HOME/ebmcal-env/bin/python" -m unittest discover -s scripts/v2 -p test_paired_runner_v2.py -q
"$HOME/ebmcal-env/bin/python" scripts/v2/verify_paired_14.py --output /N/slate/tg11/ebmcal_v2/analysis/paired_validation_14.json
