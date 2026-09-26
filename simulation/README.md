# Simulations

Code for the simulations of the paper: the data generator, the runner, the models, the run
configurations, the cluster scripts and the summaries that turn per-dataset records into the aggregate
inputs of Fig. 1d, Fig. 2a–d, Fig. 3a and Supplementary Fig. 1. The paired benchmark (Fig. 3b–d and
Supplementary Table 2) has its own runtime in `../benchmark`, which uses the model code in `scripts/`.
Simulated data are not stored; every dataset is regenerated from its seed.

## Contents

| Path | Content |
|---|---|
| `scripts/design_v2.py` | Data generator (`resolve(name, **overrides)`, `simulate(config, seed)`); design version 2.1.0 |
| `scripts/run_v2.py` | Runner: simulates each dataset of one cell, fits the models, runs the permutation tests and writes one record per dataset and model |
| `scripts/engine_v2.py` | Separately fitted DEBM on pyebm 2.0.3 (labels `original` and `repaired`) |
| `scripts/invariant_engine_v2.py` | Pooled-score DEBM (`shared`) and CONCORD (`invariant_min`; `invariant_pooled`) |
| `scripts/kde_engine_v2.py` | Likelihood EBM with the kde_ebm package (`kde_gmm`; `kde_gmm_pooled`) |
| `scripts/paired_engine_v2.py` | Fits `original` and `repaired` in one job from the same mixture fit |
| `scripts/fast_likelihood_v2.py` | Faster evaluation of pyebm's mixture likelihood with identical results |
| `scripts/dependency_v2.py` | Checks the pyebm version and the hashes of three pyebm source files before every run |
| `scripts/oracle_engine_v2.py`, `scripts/oracle_estimated_prior_v2.py` | Diagnostic models with known event distributions; imported by the runner, not used in the paper |
| `scripts/benchmark_fast_likelihood_v2.py` | Checks that `fast_likelihood_v2.py` reproduces pyebm exactly; its helpers are used by the tests |
| `scripts/test_*.py` | Unit tests |
| `scripts/dev/saebm_mechanism_dev.py` | SA-EBM with pysaebm (run `method_saebm_a`) |
| `scripts/dev/toy_target_shift.py` | Three-event example of Fig. 1d |
| `scripts/dev/core_rates_exact.py`, `scripts/dev/rates_exact_all.py` | Rejection rates of the permutation tests (Fig. 3a, Supplementary Fig. 1) |
| `configs/*.json` | The 11 run configurations used in the paper |
| `hpc/make_snapshot.py`, `hpc/run_array.sbatch`, `hpc/saebm_mechanism_freeze.sbatch` | Snapshot tool and Slurm scripts |
| `summaries/*.py` | Summaries of the run records for Fig. 2 |

The Python files directly in `scripts/`, except the tests and `benchmark_fast_likelihood_v2.py`, are
byte-identical to the copies whose SHA-256 hashes the paired benchmark records (the real-data runtime
records those of the model files), so they are left unchanged. Comments inside these engine modules
keep their development wording because the files are kept byte-identical to the code that ran; the
section on code labels gives the terms of the paper.

## Simulated data

`design_v2.simulate` generates one dataset from a design and a seed.

**Groups.** Three groups whose sizes and CN/MCI/AD counts are those of a published comparison of *APOE*
genotypes in ADNI (Venkatraghavan et al. 2021, NeuroImage 227:117646; Fig. 1c):

| Group | Code | Participants | CN | MCI | AD |
|---|---|---|---|---|---|
| CN-heavy | `e2` (g1) | 75 | 57 | 6 | 12 |
| Intermediate | `e33` (g2) | 411 | 244 | 66 | 101 |
| AD-heavy | `e4` (g3) | 485 | 110 | 156 | 219 |

The group codes are labels of simulated groups, not genotypes. The simulated data frame stores the group
in a column named `APOE`, the column that the pyebm calls use for groups.

**Stage.** Each participant's disease stage is Z = round(14 U) with U ~ Beta(0.6, 9) for CN, Beta(2, 2)
for MCI and Beta(9, 0.6) for AD participants, the same in every group. The event at position r of a
group's true sequence has occurred when Z ≥ r.

**Biomarkers 1–14.** Each biomarker defines one event. Measurements are Gaussian with unit variance, mean
0 before the event and δ = √2 Φ⁻¹(AUC) after it, where AUC is the component area under the curve. The
component AUCs were set by modality, and the observation probabilities follow the biomarker availability
reported for the same ADNI comparison (Supplementary Table 1). `design_v2.BIOMARKERS` lists the
biomarkers in this order; its names are code labels only.

| Biomarker | Code label | Modeled on | Component AUC | Observation probability |
|---|---|---|---|---|
| 1 | `ABETA` | CSF | 0.88 | 0.735 |
| 2 | `PTAU` | CSF | 0.85 | 0.735 |
| 3 | `TAU` | CSF | 0.80 | 0.725 |
| 4 | `NG` | CSF | 0.68 | 0.275 |
| 5 | `NFL` | CSF | 0.72 | 0.288 |
| 6 | `ADAS13` | cognitive score | 0.92 | 0.990 |
| 7 | `MMSE` | cognitive score | 0.86 | 1.000 |
| 8 | `Hippocampus` | imaging | 0.82 | 0.990 |
| 9 | `Entorhinal` | imaging | 0.80 | 0.990 |
| 10 | `MidTemp` | imaging | 0.78 | 0.990 |
| 11 | `Fusiform` | imaging | 0.76 | 0.990 |
| 12 | `WholeBrain` | imaging | 0.74 | 0.990 |
| 13 | `Ventricles` | imaging | 0.72 | 0.990 |
| 14 | `Precuneus` | imaging | 0.70 | 0.990 |

**Missing measurements.** Each measurement is present with its biomarker's observation probability,
independently for each participant and biomarker (missing completely at random; the probabilities do not
depend on stage, diagnosis or group). With the override `missing: false` all measurements are present
(cells whose id ends in `_CD`, and the SA-EBM run). Other missingness options of the generator are not used in the paper.

**True sequences.** Without `fixed_base_order`, each dataset draws a new uniformly random true sequence,
shared by the three groups unless the design changes one group's sequence. Fig. 2a,c,d use one fixed
true sequence in every dataset, `fixed_base_order = [5, 9, 1, 12, 11, 10, 3, 8, 13, 0, 6, 7, 2, 4]`
(0-based indices into the table above), that is, biomarkers 6, 10, 2, 13, 12, 11, 4, 9, 14, 1, 7, 8, 3, 5.

**Designs.**
- `REF_H0`, different composition: the fixed counts above; diagnoses are assigned in random order within
  each group.
- `IID_H0`, same composition: the same group sizes; every participant's diagnosis is drawn with the pooled
  probabilities (411/228/332 of 971 = 42.3% CN, 23.5% MCI, 34.2% AD), so the groups share one expected
  composition.
- `sample_scale: 4` or `16` multiplies every group-by-diagnosis count (n = 3,884 or 15,536).
- Changed orderings: `PWR_<group>_K27` reverses 27 of the 91 event pairs in one group (drawn uniformly
  among the orderings with exactly 27 reversed pairs); `K6_SINGLE_<group>` moves one event by 6
  positions; `PWR_<group>_K13` with `geometry: single_displacement` moves one event by 13 positions (the
  first event to last or the last to first); `K6_DISJOINT_<group>` makes six disjoint adjacent swaps.
- `BALANCED_H0` (margin-matched counts) and the stress designs of `resolve` are not used in the paper.

**Seeds.** The datasets of a configuration have seeds `base_seed` to `base_seed + datasets - 1`. Each seed
is split into seven independent random streams (true sequence, changed sequence, group assignment,
diagnosis, stage, measurement noise, missingness), so, for example, a complete-data dataset differs from
the dataset with the same seed only in its missing measurements. The relabelings of a permutation test
use the seed `SeedSequence([dataset seed, first 4 bytes of SHA-256(scheme name) as a little-endian
integer, permutation index, 61000000])`, so all models see the same relabelings of a dataset.

## Models and code labels

| Label in code | Model in the paper |
|---|---|
| `repaired` | Separately fitted DEBM: pyebm 2.0.3 with group-specific mixtures and orderings, with the adjacent-swap search continued until no swap lowers the loss |
| `original` | Separately fitted DEBM with pyebm's own search, which stops after the first accepted swap (in `confirm_core_a`; not shown) |
| `shared` | Pooled-score DEBM: one mixture per biomarker fitted to all participants, unweighted group orderings |
| `invariant_min` | CONCORD: pooled mixtures and weighted group orderings at the common proportions given by the smallest proportion of each diagnosis across groups, rescaled to sum to one |
| `invariant_pooled` | CONCORD with the pooled proportions as common proportions (not in the paper) |
| `kde_gmm` | Likelihood EBM: kde_ebm Gaussian-mixture event models fitted in each group to its CN and AD participants; greedy search for the maximum-likelihood sequence (5 starts of 500 iterations) |
| `kde_gmm_pooled` | Likelihood EBM with mixtures fitted to all participants (not in the paper) |
| `saebm_conjugate` | SA-EBM: pysaebm with conjugate priors, 10,000 iterations, 2,500 burn-in, each group fitted separately with CN participants as non-diseased; the ordering with the highest likelihood |
| `oracle_*` | Diagnostics with known event distributions (not in the paper) |

In these simulations the mixture fitting of separately fitted DEBM uses the pyebm 2.0.3 default stopping
rule, which checks the mixing proportions of the first group only (all of Fig. 2 and Fig. 3a). The
paired benchmark in `../benchmark` requires every group's mixtures to converge. Pooled-score DEBM and
CONCORD use the same continued adjacent-swap search as `repaired`.

| Label in code | Meaning |
|---|---|
| scheme `unrestricted` (`stratify: none`) | Unrestricted permutation of the three group labels |
| scheme `diagnosis` | Within-diagnosis permutation of the three group labels |
| schemes `diagnosis_pair_e2_e33`, `diagnosis_pair_e2_e4`, `diagnosis_pair_e33_e4` | Pairwise within-diagnosis permutation of the two named groups; the third group keeps its labels |
| `e2`, `e33`, `e4`; g1, g2, g3 | CN-heavy, intermediate and AD-heavy simulated groups |
| cell ids `REF`, `IID`, `BAL` | Different composition, same composition, margin-matched counts (not used) |
| cell id suffixes `_H0`, `_S`, `_CD`, `_N4`, `_N16` | No changed sequence; one fixed true sequence; complete data; 4 or 16 times the sample size |
| `bperm`, `rule: le`, `alpha` | B = 599 relabelings; a comparison is rejected when P = (1 + E)/(B + 1) ≤ alpha/3, with E the number of relabelings whose distance is at least the observed distance |
| `full_p_first_n` | The first N datasets use all B relabelings; the others stop once no further relabeling can change the decisions |
| `paired_standard` | `original` and `repaired` are fitted in one job from the same mixture fit |
| `fast_likelihood` | Use `fast_likelihood_v2.py` (identical results) |

## Run configurations and figures

Each configuration in `configs/` names its models (`engines`), cells, permutation schemes, number of
datasets and first seed; its `purpose` field describes it. Models and cells not listed under "Used in the
paper" are part of the same run but are not shown.

| Configuration | First seed | Datasets | Models | Missing data | Used in the paper |
|---|---|---|---|---|---|
| `method_saebm_reference_a` | 53400000 | 1,000 | `repaired`, `shared`, `invariant_min` | no | Fig. 2a,c: all three models, `REF_S_CD` and `IID_S_CD` |
| `method_kde_mechanism_a` | 53400000 | 1,000 | `kde_gmm`, `kde_gmm_pooled` | no | Fig. 2a,c: `kde_gmm`, `REF_S_CD` and `IID_S_CD` |
| SA-EBM run `method_saebm_a` (no JSON file; `hpc/saebm_mechanism_freeze.sbatch`) | 53400000 | 1,000 | `saebm_conjugate` | no | Fig. 2a,c: `REF_S` and `IID_S` |
| `confirm_core_a` | 42100000 | 1,000 | `original`, `repaired` | yes | Fig. 2b and Fig. 3a (separately fitted DEBM): `repaired`, `REF_H0` |
| `method_calibration_a` | 53100000 | 1,000 | `shared`, `invariant_min` | yes | Fig. 2b and Fig. 3a (pooled-score DEBM and CONCORD): `REF_H0` |
| `method_kde_calibration_a` | 53100000 | 300 | `kde_gmm` | no | Fig. 2b (likelihood EBM): `REF_H0_CD` |
| `method_mechanism_a` | 53000000 | 1,000 | `repaired`, `shared`, `invariant_min`, `invariant_pooled` | yes | Fig. 2d, n and 4n: `repaired` and `invariant_min`, `REF_S`, `IID_S`, `REF_S_N4`, `IID_S_N4` |
| `method_mechanism_16n_a` | 53700000 | 200 | `repaired`, `invariant_min` | yes | Fig. 2d, 16n: `REF_S_N16`, `IID_S_N16` |
| `method_pair_complete_a` | 53100000 | 1,000 | `invariant_min` | yes | Supplementary Fig. 1a |
| `method_pair_partial_a` | 53200000 | 500 | `invariant_min` | yes | Supplementary Fig. 1b: comparison between the two unchanged groups in `PWR_E2_K27`, `PWR_E33_K27`, `PWR_E4_K27` |
| `confirm_power_global_a` | 42200000 | 500 | `repaired` | yes | Datasets of the benchmark conditions with 27 of 91 pairs reversed (Fig. 3c,d; Supplementary Table 2) |
| `method_alt_geometry_a` | 53800000 | 500 | `repaired`, `invariant_min` | yes | Datasets and CONCORD results of the benchmark conditions with one event moved 13 positions in the AD-heavy group (`K13_SINGLE_E4`) and six adjacent swaps in the CN-heavy group (`K6_DISJOINT_E2`) |

Notes:
- `method_saebm_reference_a`, `method_kde_mechanism_a` and the SA-EBM run use the same seeds, so all
  models in Fig. 2a,c analyze datasets generated from the same seeds with complete measurements.
- `method_calibration_a` and `method_pair_complete_a` use the same 1,000 datasets, which are also the
  independent Gaussian benchmark condition without a true difference (Fig. 3b). Separately fitted DEBM in
  Fig. 2b and Fig. 3a analyzed another 1,000 datasets of the same design (`confirm_core_a`).
- `method_kde_calibration_a` uses the first 300 seeds of `method_calibration_a` without missing
  measurements, so its datasets are the complete-data versions of those datasets.
- Fig. 1d comes from `scripts/dev/toy_target_shift.py`, which computes by Monte Carlo (400,000 draws per
  stage, NumPy seed 20260914) the expected ranking loss of the six orderings of three events at the
  CN-heavy and AD-heavy compositions (code labels `e2` and `e4`).

## Running the simulations

**Software.** Python 3.12, NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3, scikit-learn 1.5.2, statsmodels
0.14.6 and pyebm 2.0.3 (unmodified; `dependency_v2.py` stops the run if the installed files differ).
The likelihood EBM needs kde_ebm 0.0.3. SA-EBM runs in a separate environment with pysaebm 7.7.7.

**Tests.**

```
python -m unittest discover -s simulation/scripts -p 'test_*.py'
```

Some tests fit small DEBMs and need pyebm. `test_design_v2.py` also compares the generator with the
earlier version 2.0.0, which is in the repository at
`benchmark/runtime/vendor/method_calibration_a/design_v2.py`; the test uses that file by default, and
`CONCORD_FROZEN_DESIGN` can point to another copy.

**One chunk.** A run is split into cells and chunks; chunk c of `nchunks` fits the datasets
`seeds[c::nchunks]` of one cell. Results do not depend on the chunking.

```
export CONCORD_WORK_DIR=/path/to/work
python simulation/scripts/run_v2.py --config simulation/configs/method_calibration_a.json \
    --cell-index 0 --chunk 0 --nchunks 250 --workers 16 \
    --output "$CONCORD_WORK_DIR/runs/method_calibration_a/cell_0/chunk_0"
```

The chunk directory receives `manifest.json` (configuration, source hashes and software versions),
`fit_records.jsonl` (every fit, appended as it finishes), `rows.jsonl` (one record per dataset and model:
fitted orderings, distances, true sequences and permutation counts) and `progress.json`. The runner
resumes an interrupted chunk, refuses to resume if the code or configuration changed, and stops cleanly
when a file named `STOP` appears next to the chunk directories or after `--time-budget-s` seconds.

**On a Slurm cluster.** Copy the code and one configuration into a checksummed snapshot, then submit one
array task per chunk and cell (array index = chunk × number of cells + cell):

```
python simulation/hpc/make_snapshot.py method_calibration_a       # writes $CONCORD_WORK_DIR/snapshots/method_calibration_a
export CONCORD_PYTHON=/path/to/env/bin/python
sbatch --array=0-249 "$CONCORD_WORK_DIR/snapshots/method_calibration_a/hpc/run_array.sbatch" \
    "$CONCORD_WORK_DIR/snapshots/method_calibration_a" method_calibration_a 250
```

The SA-EBM run (`method_saebm_a`) uses `hpc/saebm_mechanism_freeze.sbatch` with `SAEBM_PYTHON` set to the
pysaebm environment and `NCHUNKS=50`; its defaults are the settings of the paper (cells `IID_S`, `REF_S`,
`BAL_S`, `REF_S_N4`; seeds 53400000–53400999; 10,000 iterations, 2,500 burn-in). The Slurm scripts were
written for one university cluster: set the account, partition and resources for yours.

## Summaries

The summaries read the run records (one directory per run under `$CONCORD_WORK_DIR/runs`) and write the
aggregate files that the figure scripts read. Nothing is refitted. On the archived run records, these
commands reproduce the published files exactly (`rates_exact_final.txt` apart from its one-line
header; the counts of `fig2_calibration.csv` are the `REF_H0` rows of `rates_core.txt` and
`rates_calibration.txt`).

```
export CONCORD_RUNS_DIR=$CONCORD_WORK_DIR/runs
python simulation/summaries/build_mechanism_source_data.py --out mechanism_all_cells.csv
python simulation/summaries/mechanism_gap_intervals.py --plan mechanism_all_cells.csv \
    --out v2_mechanism_gap_intervals.csv
python simulation/summaries/runs_eventwise.py --out runs_random_truth_eventwise.csv
python simulation/scripts/dev/core_rates_exact.py confirm_core_a > rates_core.txt
python simulation/scripts/dev/rates_exact_all.py method_calibration_a > rates_calibration.txt
python simulation/scripts/dev/rates_exact_all.py method_pair_complete_a method_pair_partial_a > rates_exact_final.txt
```

| Script | Output | Figure |
|---|---|---|
| `summaries/build_mechanism_source_data.py` | `mechanism_all_cells.csv`: for each fixed-sequence run, cell and model, the mean-position gap between the CN-heavy and AD-heavy groups, systematically misordered event pairs (reversed in more than half of the datasets) per group, error against the true sequence and distances between groups | Fig. 2a,c,d |
| `summaries/mechanism_gap_intervals.py` | `v2_mechanism_gap_intervals.csv`: the mean-position gap with 95% percentile intervals from 2,000 bootstrap resamples of datasets (seed 20260919) | Fig. 2a,d |
| `summaries/runs_eventwise.py` | `runs_random_truth_eventwise.csv`: event-position shift (CN-heavy minus AD-heavy, mean over datasets) and its standard error for each biomarker | Fig. 2b |
| `scripts/dev/core_rates_exact.py`, `scripts/dev/rates_exact_all.py` | Rejection counts and rates with Wilson 95% intervals for each cell, model and permutation scheme. The published `rates_exact_final.txt` keeps only the two blocks of `rates_exact_all.py method_pair_complete_a method_pair_partial_a`; the counts in `fig2_calibration.csv` are the `REF_H0` rows of `core_rates_exact.py confirm_core_a` (model `repaired`) and of `rates_exact_all.py method_calibration_a` (models `shared`, `invariant_min`) | Fig. 3a (`fig2_calibration.csv`); Supplementary Fig. 1 (`rates_exact_final.txt`, `Fig3f_partial_null.csv`) |
| `scripts/dev/toy_target_shift.py` (printed output) | `toy_target_shift.txt` | Fig. 1d |

Notes on the published aggregate files:
- The aggregate input files `mechanism_all_cells.csv`, `v2_mechanism_gap_intervals.csv` and
  `runs_random_truth_eventwise.csv` also keep rows of additional runs that the figures do not read:
  two diagnostic runs with known event distributions (`method_oracle_a`, `method_oracle_estprior_a`)
  in the first two, and `method_calibration_pooled_a` and `method_kde_precision_a` in the third. Their
  configurations are not included here. By default the summaries include these runs; runs without
  records are skipped, and `--runs` selects runs.
- The bootstrap intervals share one random stream that runs through the run/cell/model groups in sorted
  order, so an interval depends on the groups listed before it. To reproduce the published intervals,
  pass the published `mechanism_all_cells.csv` as `--plan`; resamples of groups without records are
  still drawn, so all other intervals are unchanged.
- `runs_eventwise.py --datahash FILE` also checks that `method_calibration_a` generated the same datasets
  as the benchmark's independent Gaussian condition without a true difference (FILE maps
  `GAUSSIAN_H0|<seed>` to the data SHA-256).

## Differences from the executed runs

- **Configurations.** The published files differ from the executed ones only in the removed `stop_file`
  entry (a cluster path; the runner then uses `STOP` next to the chunk directories) and in descriptive
  text: `purpose` was rewritten and development notes were removed. Every key that `run_v2.py` reads is
  unchanged. Because the runner records the configuration in each chunk's manifest, a chunk started with
  an executed configuration cannot be resumed with the published one.
- **Code.** Each run was executed from a snapshot made before it started. The model files are the same.
  Most runs used `design_v2.py` 2.0.0; version 2.1.0 only adds missingness options that are off by
  default, so the datasets are unchanged. Later versions of `run_v2.py` only add model labels and a
  permutation option.
- **SA-EBM.** The run imported `design_v2.py` 2.0.0. In `saebm_mechanism_dev.py` only the docstring and
  the default location of `design_v2.py` changed.
