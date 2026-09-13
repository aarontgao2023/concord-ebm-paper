# concord-ebm-paper

Simulations, analysis code and results behind **CONCORD** (COmposition-Normalised Consensus for ORDering
comparison), a composition-invariant way to compare event-based model orderings between groups. The method
itself is the Python package [`concord-ebm`](https://github.com/aarontgao2023/concord-ebm); this repository
holds everything needed to reproduce the paper.

```
simulation/   the ADNI-shaped simulator, the resumable runner, the audited engines as frozen for the
              pre-registered protocols, the analysis scripts, and every configuration (cells, seeds, budgets)
results/      result tables of the frozen protocols (rates with exact decisions, mechanism tables, tallies)
figures/      Figure 1 (mechanism schematic), editable source included
adni/         the ADNI reanalysis (scripts to follow once data access is in place; no data is stored here)
```

## Reproducing the simulations

Every cell is defined by a JSON configuration in `simulation/configs/` (design name and overrides, seeds,
engines, permutation schemes, budget `B=599`, decision rule). Data sets are regenerated exactly from the
seed (`design_v2.simulate`), so only configurations and summary tables are versioned; the raw per-fit records
(several GB) are archived off-repository.

```
cd simulation/scripts
python -m unittest discover -s . -p "test_*.py"            # unit tests of design, engines, runner
python run_v2.py --config ../configs/method_calibration_a.json --cell-index 0 --chunk 0 --nchunks 250 \
       --workers 16 --output /path/to/out                   # one chunk of one cell
```

`simulation/hpc/` contains the Slurm scripts used on IU Quartz (`run_array.sbatch`, `make_snapshot.py`, the
top-up driver). Protocols were frozen as source snapshots before execution; `make_snapshot.py` reproduces
that step.

Two protocols:

| protocol | seeds | content |
|---|---|---|
| first (standard estimator) | 42.xM | IID/REF/STAGE calibration, complete- and partial-null pair tests, power at K=14/27/55 |
| second (method cells) | 53.xM | fixed-ordering mechanism for four estimators, SA-EBM comparison, calibration / power / boundary of the composition-invariant estimator, pair tests, boundary repair |

`results/method_protocol_20260912/README.md` lists the headline numbers and the exact-decision convention.

## Requirements

Python ≥ 3.10, `pyebm==2.0.3` (unmodified; verified by hash before every fit), numpy, pandas, scipy.
The SA-EBM comparison additionally needs `pysaebm` in a separate environment (see
`simulation/scripts/dev/saebm_mechanism_dev.py`).

## Licence

GPL-3.0-or-later (the engines patch and adapt code from the GPL-3 `pyebm`).
