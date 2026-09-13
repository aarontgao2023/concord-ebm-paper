# simulation

| path | content |
|---|---|
| `scripts/design_v2.py` | ADNI-shaped data-generating process; named cells (`REF_H0`, `IID_H0`, `BALANCED_H0`, `STAGE_H0`, `PDF_H0`, `MISS_*_H0`, `PWR_*`) |
| `scripts/engine_v2.py`, `invariant_engine_v2.py`, `fast_likelihood_v2.py`, `dependency_v2.py` | the audited engines exactly as frozen in the protocol snapshots (the package `concord-ebm` carries the same code with package imports) |
| `scripts/run_v2.py` | resumable runner: append-only fit records, exact decisions under a fixed budget, early stopping, permutation operators (`unrestricted`, `diagnosis`, pair-restricted, `diagnosis_observed_count`) |
| `scripts/oracle_*.py`, `paired_engine_v2.py` | mechanism diagnostics used in the first protocol |
| `scripts/ebm_compare.py` | research copy of the packaged procedure |
| `scripts/analyze_*.py`, `scripts/dev/` | analysis and development scripts; `dev/analyze_method_protocol.py`, `dev/rates_exact_all.py`, `dev/core_rates_exact.py`, `dev/tally_*_exact.py` produce the tables in `../results/` |
| `scripts/test_*.py` | unit tests |
| `configs/` | every run configuration; `confirmation_protocol_20260908.json` is the first protocol's freeze record |
| `hpc/` | Slurm array script, snapshot tool, submission/top-up drivers and their audits |

Array index convention: `index = chunk * ncells + cell`; each chunk owns one output directory and is
resumable. A stopped or failed permutation stays *unknown* inside the planned budget and never shrinks the
denominator.
