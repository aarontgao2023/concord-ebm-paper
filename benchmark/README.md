# Paired simulation benchmark (Fig. 3b–d, Supplementary Table 2)

This directory holds the code for the paired benchmark. The benchmark compares separately fitted DEBM, pooled-score DEBM and CONCORD on the same simulated datasets. Every method is tested by within-diagnosis permutation with the same relabelings. The directory also holds the script that builds `rates.csv` and `paired_differences.csv`, which feed Fig. 3b–d and Supplementary Table 2. Fig. 3a is not part of the benchmark (see `../simulation/`).

```
benchmark/
  README.md
  formal.sbatch                   Slurm array script: one batch of datasets per task
  summarize_campaign_audited.py   builds rates.csv, paired_differences.csv and the other summary tables
  runtime/
    campaign.json                 the 11 conditions: seeds, dataset counts, methods fitted, B = 599
    batch_mapping.json            datasets of each Slurm array task (48 + 64 batches)
    run_dataset.py                fits one dataset: observed data plus relabelings, with exact early stopping
    run_batch.py                  runs one batch of batch_mapping.json in parallel
    robust_design.py              correlated and heavy-tailed noise conditions
    allgroup_engine.py            all-groups stopping rule for the pyebm mixture fits
    audit_reuse.py                checks the datasets reused from the simulation runs (no fitting)
    EXECUTION_GATE.json           source, input and environment hashes checked before a formal run
    FROZEN_PACKAGE.json           hashes of every file in runtime/, checked by the summary script
    inputs/historical_reused_rows.jsonl   earlier simulation results that the benchmark reuses (15 MB)
    vendor/confirm_power_global_a/, vendor/method_alt_geometry_a/, vendor/method_calibration_a/
                                  generator and runner snapshots of the three earlier simulation runs
    vendor/w3/                    instrumented pyebm mixture loop used by allgroup_engine.py
```

`vendor/w3/run_stop_isolation.py` is used as a module only: `allgroup_engine.py` imports its class `AllGroupState` (the all-groups stopping rule). Its standalone entry point (`__main__`) needs a script that is not in the repository.

The engine modules (`engine_v2.py`, `invariant_engine_v2.py`, `fast_likelihood_v2.py`, `dependency_v2.py`, `run_v2.py`, `design_v2.py` and the modules `run_v2.py` imports) are loaded from `../simulation/scripts`. The benchmark ran with vendored copies of these modules, whose hashes are recorded in `runtime/EXECUTION_GATE.json` (see "Hash checks" below).

## Design

### Conditions

All conditions use 14 biomarkers and three simulated groups: CN-heavy (75 participants: 57 CN, 6 MCI, 12 AD), intermediate (411: 244/66/101) and AD-heavy (485: 110/156/219), with n = 971. Biomarkers 1–14 follow the input order of `design_v2.BIOMARKERS`: biomarkers 1–5 are modeled on CSF markers, 6–7 on cognitive scores and 8–14 on imaging measures. In every condition, measurements are missing completely at random, with the observation probabilities of Supplementary Table 1.

| Condition (code) | Noise | Changed ordering | Datasets | Dataset seeds | Datasets from |
|---|---|---|---|---|---|
| `GAUSSIAN_H0` | independent Gaussian | none | 1,000 | 53100000–53100999 | simulation run `method_calibration_a`, cell `REF_H0` |
| `EXACT27_G1` | independent Gaussian | 27 of 91 pairs reversed, CN-heavy group | 500 | 42200000–42200499 | `confirm_power_global_a`, cell `PWR_E2_K27` |
| `EXACT27_G3` | independent Gaussian | 27 of 91 pairs reversed, AD-heavy group | 500 | 42200000–42200499 | `confirm_power_global_a`, cell `PWR_E4_K27` |
| `SINGLE13_G3` | independent Gaussian | one event moved 13 positions, AD-heavy group | 500 | 53800000–53800499 | `method_alt_geometry_a`, cell `K13_SINGLE_E4` |
| `DISJOINT6_G1` | independent Gaussian | six adjacent swaps, CN-heavy group | 500 | 53800000–53800499 | `method_alt_geometry_a`, cell `K6_DISJOINT_E2` |
| `CORR_H0` | correlated within modality | none | 1,000 | 92110000–92110999 | new |
| `CORR_SINGLE6_G1` | correlated within modality | one event moved 6 positions, CN-heavy group | 500 | 92120000–92120499 | new |
| `CORR_SINGLE13_G3` | correlated within modality | one event moved 13 positions, AD-heavy group | 500 | 92130000–92130499 | new |
| `T5_H0` | heavy-tailed (t, 5 d.f.) | none | 1,000 | 92210000–92210999 | new |
| `T5_SINGLE6_G1` | heavy-tailed (t, 5 d.f.) | one event moved 6 positions, CN-heavy group | 500 | 92220000–92220499 | new |
| `T5_SINGLE13_G3` | heavy-tailed (t, 5 d.f.) | one event moved 13 positions, AD-heavy group | 500 | 92230000–92230499 | new |

Together these are 7,000 datasets: three conditions without a change (1,000 datasets each) and eight with a changed ordering (500 each). Dataset seeds run from `base_seed` to `base_seed + datasets − 1` in `campaign.json`. A dataset seed fixes the data and the true orderings. The two 27-pair conditions share seeds, and so do `SINGLE13_G3` and `DISJOINT6_G1`; within each of these pairs, the change is placed in a different group or has a different form.

- **True sequences.** In the independent Gaussian conditions, a new random base sequence is drawn for each dataset. The correlated and heavy-tailed conditions use the order of Supplementary Table 1 (biomarkers 1–14 in input order). Without a change, all three groups share the base sequence. With a change, only the named group's sequence differs.
- **Changes.** For "27 of 91 pairs reversed", the changed ordering is drawn uniformly among the orderings with exactly 27 pairs reversed relative to the base (normalized Kendall distance 0.30). "One event moved m positions" moves one event, with the start and end positions drawn uniformly among those m positions apart. With 14 events, a move of 13 places the first event last or the last event first. "Six adjacent swaps" reverses six non-overlapping adjacent pairs. The generator checks the realized number of reversed pairs for every dataset.
- **Noise.** Independent Gaussian noise has unit variance. Correlated noise has correlation 0.55 among biomarkers 1–5, 0.75 between biomarkers 6 and 7 and 0.60 among biomarkers 8–14, and none across modalities (generator option `correlated=True`). Heavy-tailed noise is √(3/5)·t₅, which has unit variance. For each biomarker, the shift after the event is solved by numerical integration and root finding, so that the component AUC equals that of the Gaussian design (`robust_design.py`). In the heavy-tailed conditions, stage, diagnosis, true orderings and the missingness mask are those of the independent Gaussian dataset with the same seed.
- **Models.** All three methods fit Gaussian mixtures in every condition.

### Methods and stopping rule

| Method | Code label |
|---|---|
| Separately fitted DEBM (pyebm 2.0.3 mixtures fitted per group) | `repaired` |
| Pooled-score DEBM (one mixture per biomarker fitted to all three groups; unweighted consensus) | `shared` |
| CONCORD (the same pooled mixture; consensus weighted to the common proportions) | `invariant_min` |

- **Stopping rule.** Separately fitted DEBM used the all-groups stopping rule. The pyebm mixture loop stops when the mean absolute change in the mixing proportions falls below 0.01 in *every* group; pyebm's default rule checks only the first group. The objective, initialization, bounds and inner SLSQP optimizer are unchanged. A fit that reaches 100 outer iterations is counted as failed. The rule is implemented in `allgroup_engine.py`. Observed fits use the instrumented loop in `vendor/w3/`, which also records the trajectory. Relabeled fits use a lighter version that replaces only the stopping test in pyebm's `do_mixturemodel`.
- **Pooled mixture.** Pooled-score DEBM and CONCORD fit one pooled group, so the all-groups rule and the default rule coincide. The pooled mixture does not depend on group labels. It is therefore fitted once per dataset, with the all-groups rule, and reused for every relabeling. Weights and orderings are recomputed at every relabeling.
- **Consensus search.** All three methods use the same adjacent-swap search, continued until no swap lowers the loss (engine mode `repaired`).

### Test

- **Statistic.** For each of the three group comparisons, the statistic is the normalized Kendall distance between the two groups' fitted orderings.
- **Permutation.** Within-diagnosis permutation exchanges all three group labels separately among CN, MCI and AD participants. Relabeling `pid` of the dataset with seed `s` uses `numpy.random.SeedSequence([s, tag, pid, 61000000])`, where `tag` is the first four bytes (little-endian) of SHA-256 of `"diagnosis"` (function `permute` in `simulation/scripts/run_v2.py`). All methods therefore see the same relabelings. Before fitting, `run_dataset.py` records a hash of all 599 relabelings.
- **Decision.** With B = 599 relabelings, P = (1 + E)/600, where E counts relabelings whose distance is at least the observed distance (ties count). A comparison is rejected when P ≤ 0.05/3, that is, when E ≤ 9.
- **Early stopping.** A method stops relabeling a dataset once all three decisions are the same for every possible result of the remaining relabelings. After s successful relabelings with e exceedances, the final P lies in [(1 + e)/600, (1 + e + 599 − s)/600]. Failed relabelings count as unknown slots within the 599. For a dataset stopped early, the summary tables give bounds on the adjusted P value and no point value.

### Datasets and results reused from the simulation runs

The five independent Gaussian conditions reuse datasets from three earlier simulation runs (configurations in `../simulation/configs/`). They also reuse some fitted results, which had been computed with the same relabelings, the same decision rule and B = 599:

- `GAUSSIAN_H0`: the 1,000 datasets of `method_calibration_a` (`REF_H0`). These are the datasets analyzed by pooled-score DEBM and CONCORD in Fig. 2b and Fig. 3a. Their pooled-score DEBM and CONCORD results are reused. Separately fitted DEBM was fitted anew with the all-groups rule.
- `EXACT27_G1`, `EXACT27_G3`: datasets of `confirm_power_global_a`. All three methods were fitted anew.
- `SINGLE13_G3`, `DISJOINT6_G1`: datasets of `method_alt_geometry_a`. The CONCORD results are reused. Separately fitted and pooled-score DEBM were fitted anew.
- The correlated and heavy-tailed conditions are new datasets, and all three methods were fitted.

In total, 18,000 method-dataset results were fitted for the benchmark and 3,000 were reused. In the published run, the new analyses comprised 2,447,908 fits (18,000 observed and 2,429,908 relabeled). No fit failed, and none reached the 100-iteration safeguard.

`inputs/historical_reused_rows.jsonl` holds 5,000 compact records of the earlier runs, one per dataset and method. Each record has the run, cell, seed and engine label, the saved generator configuration and data hash (`truth`), the observed orderings and distances, and the exceedance counts and decisions of the within-diagnosis test. It contains simulation results only. Reused datasets are regenerated with the generator snapshot of their run (`vendor/<run>/design_v2.py`). `run_dataset.py` asserts that each regenerated dataset has the saved data hash before fitting. `audit_reuse.py` performed the full check before the benchmark ran. For all 3,000 reused datasets, it regenerated the data and truth, compared the data, configuration, generator and truth hashes with the per-fit records of the earlier runs, and confirmed that the relabeling function was unchanged. It also wrote `historical_reused_rows.jsonl`.

The file also holds the earlier separately fitted DEBM results for the four Gaussian conditions with a change. These results used the pyebm default stopping rule (2,000 records). They appear only in `historical_first_group_stopping.csv` and are not used in `rates.csv` or `paired_differences.csv`.

## Running

### Environment

The benchmark ran under Python 3.12.14 with NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3, scikit-learn 1.5.2, statsmodels 0.14.6 and pyebm 2.0.3. A formal run refuses to start unless these versions match `EXECUTION_GATE.json`, and `dependency_v2.verify_pyebm` checks the pyebm source files against their recorded hashes. The summary step needs only NumPy. Re-running it on the returned results with NumPy 2.5 reproduced every table byte for byte.

### 1. Check the reused datasets (optional, no fitting)

```
CONCORD_HISTORICAL_ROOT=/path/to/simulation/runs python runtime/audit_reuse.py --workers 8
```

`CONCORD_HISTORICAL_ROOT` is the directory that holds `<configuration>/cell_*/chunk_*/rows.jsonl`, as written by `simulation/scripts/run_v2.py` for `confirm_power_global_a`, `method_alt_geometry_a` and `method_calibration_a`. These per-fit records are regenerated from the seeds and are not versioned. Outputs go to `$CONCORD_WORK_DIR/reuse_audit` (default `benchmark/results/reuse_audit`) or to `$CONCORD_REUSE_OUT`.

### 2. Fit

One dataset:

```
cd runtime
python run_dataset.py --campaign robustness --cell-index 0 --seed-index 0 --output ../results/formal
```

One batch of `batch_mapping.json` (48 batches for `power`, 64 for `robustness`; up to 64 datasets each):

```
python run_batch.py --campaign power --batch-index 0 --workers 16 --output ../results/formal
```

On Slurm, `formal.sbatch` runs one batch per array task. The published run used 16 CPUs, 32 GB and at most 24 h per task. Edit the account and partition for your cluster, set `CONCORD_PYTHON` to the Python of the analysis environment and, optionally, `CONCORD_BENCHMARK_OUT` to the output folder (default `results/formal`).

Each dataset is written to `<output>/<campaign>/<condition>/seed_<seed>/`:

- `manifest.json`: configuration, source and environment hashes and the relabeling hash.
- `truth.json`: the true orderings and generator manifest.
- `fit_records.jsonl`: one line per observed or relabeled fit (append-only; an interrupted run resumes).
- `result.json`: exceedances, decisions and adjusted-P bounds per method.
- `pooled_score_cache.npz`: the cached pooled mixture.

`results/` is ignored by git. In a formal run (no `--benchmark-b`), `run_dataset.py` checks `EXECUTION_GATE.json` first: the hashes of `campaign.json`, of `inputs/historical_reused_rows.jsonl` and of every `.py` file in `runtime/`, and the package versions. `--benchmark-b B` (B ≤ 8) skips these checks, runs B relabelings without early stopping and was used only for timing. Its results are not formal results, and in the correlated and heavy-tailed conditions it uses separate seeds (92310000 + 1000 × condition index + dataset index).

### 3. Summarize

```
python summarize_campaign_audited.py --results results/formal --output results/summary
```

The script checks `runtime/FROZEN_PACKAGE.json` and requires all 7,000 planned result directories. It reads each `result.json` and `truth.json`, and takes the reused results from `inputs/historical_reused_rows.jsonl`. It then writes into a new directory:

| File | Contents |
|---|---|
| `per_dataset_endpoints.csv` | one row per dataset and method: relabeling counts and the endpoints below |
| `pair_decisions_and_bounds.csv` | one row per dataset, method and comparison: exceedances, decision, raw and adjusted P bounds (exact values when all 599 relabelings were run) |
| `truth_error_and_distance.csv` | per dataset and method: normalized Kendall distance of each group's fitted ordering from its true ordering, and between the fitted orderings |
| `rates.csv` | rates with Wilson 95% intervals (see below) |
| `paired_differences.csv` | paired differences between methods with bootstrap 95% intervals (see below) |
| `continuous_summary.csv`, `paired_continuous_differences.csv` | the same summaries for the distances |
| `historical_first_group_stopping.csv` | earlier separately fitted DEBM results with the pyebm default stopping rule, same datasets (four Gaussian conditions with a change) |
| `dataset_provenance.csv` | hashes of each dataset's `result.json`, `manifest.json`, `truth.json` and `fit_records.jsonl` |
| `audited_fit_inventory.csv` | fit counts, warnings and safeguard exits per condition (only with the audit files) |
| `summary_manifest.json` | counts, file hashes, definitions and seeds |

The published tables were produced with four additional options (`--record-audit`, `--transport-audit`, `--historical-audit`, `--native-binding`). These files come from a record-level audit of the results returned from the cluster and are not included in this repository. Without them, the script computes the same tables from the result directories. `rates.csv`, `paired_differences.csv` and all other tables are identical, except the `outer_cap_records` columns, which are left empty and marked "not audited", and `audited_fit_inventory.csv`, which is not written.

## How `rates.csv` and `paired_differences.csv` are built

Endpoints per dataset and method; pairs 1, 2 and 3 are CN-heavy vs intermediate, CN-heavy vs AD-heavy and intermediate vs AD-heavy:

- `any_pair`: at least one of the three comparisons rejected. This is the false-positive endpoint for the three conditions without a change (Fig. 3b; Supplementary Table 2).
- `truly_affected_pair`: at least one comparison involving the changed group rejected. This is the detection endpoint (Fig. 3c,d; Supplementary Table 2). It is empty for conditions without a change.
- `unchanged_pair`: at least one comparison between groups with the same true ordering rejected. For conditions without a change it equals `any_pair`.
- `pair_1`, `pair_2`, `pair_3`: rejection of each comparison.

`rates.csv` has one row per condition, method and endpoint (189 rows). The columns are `campaign`, `cell` (condition), `arm` (method code label), `endpoint`, `rejections`, `planned_R` (number of datasets), `rate`, and `wilson_low`/`wilson_high` (two-sided 95% Wilson interval, z = 1.959963984540054, no adjustment across endpoints). The denominator is always the planned number of datasets.

`paired_differences.csv` has one row per condition, method contrast and endpoint (189 rows). The contrasts are, in order, `shared` − `repaired`, `invariant_min` − `shared` and `invariant_min` − `repaired`. For each dataset, the per-dataset difference is 1, 0 or −1. `difference` is its mean over the datasets of the condition. The interval (`paired_bootstrap_low`, `paired_bootstrap_high`) is the 2.5th and 97.5th percentile, with NumPy's default linear quantile, of 10,000 bootstrap means. Each bootstrap mean resamples the datasets with replacement, so that the three methods stay paired. The resamples are drawn as 100 blocks of 100 from `numpy.random.default_rng(bootstrap_seed)`, with

```
bootstrap_seed = 92600000 + 1000 × condition index + 100 × contrast index + endpoint index
```

The condition index is the position within its campaign in `campaign.json`: `power` holds `EXACT27_G1`, `EXACT27_G3`, `SINGLE13_G3`, `DISJOINT6_G1` and `GAUSSIAN_H0` (indices 0–4), and `robustness` holds `CORR_H0` … `T5_SINGLE13_G3` (indices 0–5). The contrast index is 0–2 in the order above. The endpoint index is 0–5 in the order `any_pair`, `pair_1`, `pair_2`, `pair_3`, `truly_affected_pair`, `unchanged_pair`. The seed is stored in the `bootstrap_seed` column. For example, CONCORD minus separately fitted DEBM for detection with 27 pairs reversed in the AD-heavy group uses seed 92601204, and gives 0.406 [0.354, 0.458]. The distance contrasts in `paired_continuous_differences.csv` use 92700000 (+100000 for `robustness`) + the same terms.

The figures use these rows:

- Fig. 3b: `any_pair` in `GAUSSIAN_H0`, `CORR_H0` and `T5_H0`.
- Fig. 3c: `truly_affected_pair` in `EXACT27_G1`, `EXACT27_G3`, `SINGLE13_G3` and `CORR_SINGLE13_G3`. Its brackets are the `invariant_min` − `shared` rows of `paired_differences.csv`.
- Fig. 3d: the `invariant_min` − `repaired` rows for the same four conditions.
- Supplementary Table 2: `any_pair` for the three conditions without a change and `truly_affected_pair` for the eight with a change.

## Hash checks

- `EXECUTION_GATE.json` lists the SHA-256 of `campaign.json`, `inputs/historical_reused_rows.jsonl` and every `.py` file under `runtime/`, together with the package versions. `run_dataset.py` asserts these in a formal run.
- `FROZEN_PACKAGE.json` lists every file in `runtime/` and the hash of `EXECUTION_GATE.json`.
- `summarize_campaign_audited.py` checks both files and holds the hash of `FROZEN_PACKAGE.json` (`CLOSURE`).
- After editing any file in `runtime/`, recompute `source_sha256` in `EXECUTION_GATE.json`, then `FROZEN_PACKAGE.json`, then `CLOSURE`.

The published copy differs from the copy that ran in three files. `run_dataset.py`, `robust_design.py` and `audit_reuse.py` import the engine modules from `../simulation/scripts` instead of a vendored `vendor/current/` copy, and `audit_reuse.py` also has new default output locations. The hashes were recomputed for this layout. The executed hashes are kept in `EXECUTION_GATE.json` under `published_copy.executed_source_sha256`, and the hash of the executed `FROZEN_PACKAGE.json` is kept as `executed_frozen_package_sha256`. `campaign.json`, `batch_mapping.json`, `inputs/historical_reused_rows.jsonl`, `allgroup_engine.py`, `run_batch.py` and all files under `vendor/` are unchanged. Because the source hashes differ, a re-run writes a different `manifest.json` identity from the published run, and `truth.json` of the correlated and heavy-tailed conditions records the new hash of `robust_design.py`. The data hashes do not change.

To compare the engine modules in `simulation/scripts` with the copies that ran, run this from the repository root:

```
python - <<'EOF'
import hashlib, json, pathlib
gate = json.loads(pathlib.Path("benchmark/runtime/EXECUTION_GATE.json").read_text())
for key, value in gate["published_copy"]["executed_source_sha256"].items():
    if key.startswith("vendor/current/"):
        path = pathlib.Path("simulation/scripts") / key.split("/")[-1]
        print(path, path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == value)
EOF
```

## Code labels

| Label | Meaning |
|---|---|
| `repaired`, `shared`, `invariant_min` | separately fitted DEBM, pooled-score DEBM, CONCORD (column `arm` in the tables; `new_arms` in `campaign.json` lists the methods fitted anew in each condition) |
| `e2`, `e33`, `e4`; `G1`, `G3` in condition names; `g1`–`g3` in column names | CN-heavy, intermediate and AD-heavy simulated groups (code names only; the groups are simulated and are not genotypes). The group label is stored in the column `APOE` of the simulated data |
| `ABETA`, `PTAU`, … in the generator | code names of biomarkers 1–14 (input order of `design_v2.BIOMARKERS`) |
| `diagnosis` (permutation scheme) | within-diagnosis permutation |
| `power`, `robustness` (campaigns) | the five independent Gaussian conditions (datasets from earlier simulation runs) and the six correlated and heavy-tailed conditions (new datasets) |
| `H0`; `EXACT27`; `SINGLE6`, `SINGLE13`; `DISJOINT6` | no change; 27 of 91 pairs reversed; one event moved 6 or 13 positions; six adjacent swaps |
| `CORR`, `T5` | Gaussian noise correlated within modality; heavy-tailed noise (t, 5 d.f.) |
| `allgroup01`, `allgroup_original_tolerance`, `new_allgroup01` | all-groups stopping rule with tolerance 0.01; `new_allgroup01` marks results fitted for the benchmark |
| `certified_historical_reuse` | a reused result from an earlier simulation run |
| `historical_repaired`, `historical_standard_role` | separately fitted DEBM with the pyebm default stopping rule from the earlier runs (used only in `historical_first_group_stopping.csv`) |
| `bperm`, `rule: inclusive_bonferroni` | B = 599 relabelings; reject when P ≤ 0.05/3 |
| `vendor/w3`, `W3`, `W4`, `W5` in comments and version strings | internal work-package names from development |
