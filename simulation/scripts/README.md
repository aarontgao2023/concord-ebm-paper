# v2 simulation workflow

This directory is separate from the v1 scripts and results. ADNI is outside the current stage.

- `design_v2.py`: versioned independent simulation streams; exact inversion-count alternatives; explicit nuisance and missingness scope; complete latent truth.
- `engine_v2.py`: original pyebm 2.0.3 and a scoped repair to consensus stopping. Mixture fitting is held fixed. Source fingerprints are checked before fitting.
- `oracle_engine_v2.py`: equal-prior and true-order-informed known-prior diagnostics. Neither is an inferential candidate.
- `oracle_estimated_prior_v2.py`: known component densities, but group mixture priors estimated solely from observed measurements.
- `run_v2.py`: append each fit, retain failures, derive exact decisions relative to the planned fixed budget, and resume only under matching source/config/numerical versions.
- `analyze_v2.py`: dataset-level denominators and Monte Carlo uncertainty, paired comparisons, partial-null endpoints and matched pairwise families.

Run with the pinned HPC environment `~/ebmcal-env/bin/python`. Locally, the review's unpacked wheel may be placed on `PYTHONPATH` without installation. Run necessary tests with `python -m unittest discover -s scripts/v2 -p 'test_*.py'`.

Create a uniquely named JSON config under `configs/v2`, then run `python hpc/v2/make_snapshot.py NAME`. Never edit a frozen snapshot. The snapshot includes checksums, and the HPC submission verifies them before execution. Use `hpc/v2/run_array.sbatch` from that snapshot: array mapping can differ between snapshots and is part of the frozen runner. Current submissions alternate cells across array indices; chunk seeds are deterministic strides of the full planned seed set.

A completed process is not necessarily a completed statistical decision. Do not infer experiment completion from output filenames or a subset of finished seeds. Preserve the full planned denominator and report unknown outcomes when a failure prevents certification. Never use the distribution of p values among completion-selected rows after early stopping. Only a preselected full-budget subset may enter that distribution summary.

A pair-restricted scheme declares `tested_pair_index` and `family`; only its matching pair is an applicable test. Three matching members can be assembled into a Bonferroni family. A global-max test does not localize which pairs differ. The oracle stratification for the present estimator is `(Diagnosis, true stage)`, because diagnosis is part of the fitted input.

Development and confirmation use disjoint seeds. A change prompted by confirmation outcomes requires a newly named candidate and fresh confirmation seeds. Check `docs/STATE_v2.md`, `docs/ITERATIONS_v2.md`, and `runs/v2/job_registry.json` before continuing any work. The shared SSH ControlMaster belongs to another task: reuse it only and never close or reconfigure it.
