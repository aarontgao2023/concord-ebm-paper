# concord-ebm-paper

This repository holds the code for the simulations, the ADNI and NACC analyses, and the figures and
tables of the article

> Tianchuan Gao, Haoyuan Liu and Jingwen Yan. *Comparing biomarker orderings between patient groups
> with different diagnostic composition.*

CONCORD fits one abnormality model to the pooled sample, weights every group to the same CN/MCI/AD
proportions (the common proportions) and tests group differences by permuting group labels within
diagnosis. The method is distributed as the Python package
[concord-ebm](https://github.com/aarontgao2023/concord-ebm). The analyses here call the study
implementation in [`simulation/scripts/`](simulation/scripts/), in which
[`invariant_engine_v2.py`](simulation/scripts/invariant_engine_v2.py) implements pooled-score DEBM
and CONCORD. The scripts for the earlier ADNI and NACC analysis populations use concord-ebm 0.1.1
(see [Software](#software)).

No participant-level data are included (see [Data](#data)). Code is released under the GNU General
Public License v3.0 or later (GPL-3.0-or-later; [`LICENSE`](LICENSE)).

## Contents

| Path | Content | Display items |
|---|---|---|
| [`simulation/`](simulation/) | data generator, runner and model code, the 11 run configurations, Slurm scripts, and the summaries that produce the simulation figure inputs | Figs 1d, 2, 3a; Supplementary Fig. 1; Supplementary Table 1 |
| [`benchmark/`](benchmark/) | paired simulation benchmark: separately fitted DEBM, pooled-score DEBM and CONCORD on the same 7,000 datasets | Fig. 3b–d; Supplementary Table 2 |
| [`realdata/`](realdata/) | six-event ADNI and NACC analyses: cohort selection, controlled composition experiment, *APOE* comparisons, bootstrap refits | Figs 1c, 4a–d, 5; Tables 2, 3; Supplementary Tables 3–5; Supplementary Fig. 2 |
| [`realdata/earlier_panels/`](realdata/earlier_panels/) | earlier ADNI and NACC analysis populations with other biomarker panels | Fig. 4e; Supplementary Fig. 3 |
| [`figures/`](figures/) | figure and table scripts, and the simulation-only aggregate inputs they read | all display items |
| [`requirements/`](requirements/) | Python environments: [`analysis.txt`](requirements/analysis.txt), [`saebm.txt`](requirements/saebm.txt), [`figures.txt`](requirements/figures.txt) | |

Each directory has its own README with the details of its design, commands and outputs:
[`simulation/README.md`](simulation/README.md), [`benchmark/README.md`](benchmark/README.md),
[`realdata/README.md`](realdata/README.md) and [`figures/README.md`](figures/README.md).

## Display items

Figure and table scripts are in [`figures/scripts/`](figures/scripts/) and
[`figures/tables/`](figures/tables/). Inputs with a link are simulation-only aggregates in
[`figures/inputs/simulation/`](figures/inputs/simulation/). Inputs marked † are aggregates of ADNI
or NACC results (counts, orderings, rates and intervals); they are not in the repository and are
read from [`figures/inputs/realdata/`](figures/inputs/realdata/). Run names refer to the
configurations in [`simulation/configs/`](simulation/configs/). Seeds are listed in the next section.

| Item | Script | Inputs | Analysis code and runs |
|---|---|---|---|
| Fig. 1a,b | [`fig1_overview.py`](figures/scripts/fig1_overview.py) | none (schematic) | |
| Fig. 1c | [`fig1_overview.py`](figures/scripts/fig1_overview.py) | CN/MCI/AD counts of Table 2 (written in the script; checked against `reference_counts_weights_ESS.csv`† when present); simulated group sizes | [`realdata/cohort/`](realdata/cohort/), [`prepare_campaign.py`](realdata/prepare_campaign.py), [`summarize_common_panel.py`](realdata/summarize_common_panel.py); [`design_v2.py`](simulation/scripts/design_v2.py) |
| Fig. 1d | [`fig1_overview.py`](figures/scripts/fig1_overview.py) | [`toy_target_shift.txt`](figures/inputs/simulation/toy_target_shift.txt) | [`toy_target_shift.py`](simulation/scripts/dev/toy_target_shift.py) |
| Fig. 2a,c | [`fig2_simulation_composition.py`](figures/scripts/fig2_simulation_composition.py) | [`mechanism_all_cells.csv`](figures/inputs/simulation/mechanism_all_cells.csv), [`v2_mechanism_gap_intervals.csv`](figures/inputs/simulation/v2_mechanism_gap_intervals.csv) | runs `method_saebm_reference_a`, `method_kde_mechanism_a` and the SA-EBM run `method_saebm_a` ([`saebm_mechanism_freeze.sbatch`](simulation/hpc/saebm_mechanism_freeze.sbatch)); [`build_mechanism_source_data.py`](simulation/summaries/build_mechanism_source_data.py), [`mechanism_gap_intervals.py`](simulation/summaries/mechanism_gap_intervals.py) |
| Fig. 2b | [`fig2_simulation_composition.py`](figures/scripts/fig2_simulation_composition.py) | [`runs_random_truth_eventwise.csv`](figures/inputs/simulation/runs_random_truth_eventwise.csv), [`S1_event_parameters.csv`](figures/inputs/simulation/S1_event_parameters.csv) | runs `confirm_core_a`, `method_calibration_a`, `method_kde_calibration_a`; [`runs_eventwise.py`](simulation/summaries/runs_eventwise.py) |
| Fig. 2d | [`fig2_simulation_composition.py`](figures/scripts/fig2_simulation_composition.py) | as Fig. 2a,c | runs `method_mechanism_a` (1× and 4× sample size), `method_mechanism_16n_a` (16×); summaries as Fig. 2a,c |
| Fig. 3a | [`fig3_permutation_testing.py`](figures/scripts/fig3_permutation_testing.py) | [`fig2_calibration.csv`](figures/inputs/simulation/fig2_calibration.csv) | runs `confirm_core_a` (separately fitted DEBM) and `method_calibration_a` (pooled-score DEBM, CONCORD), cell `REF_H0`; [`core_rates_exact.py`](simulation/scripts/dev/core_rates_exact.py), [`rates_exact_all.py`](simulation/scripts/dev/rates_exact_all.py) |
| Fig. 3b–d | [`fig3_permutation_testing.py`](figures/scripts/fig3_permutation_testing.py) | [`rates.csv`](figures/inputs/simulation/rates.csv), [`paired_differences.csv`](figures/inputs/simulation/paired_differences.csv) | [`benchmark/runtime/campaign.json`](benchmark/runtime/campaign.json) (11 conditions), [`run_batch.py`](benchmark/runtime/run_batch.py), [`summarize_campaign_audited.py`](benchmark/summarize_campaign_audited.py) |
| Fig. 4a | [`fig4_real_composition.py`](figures/scripts/fig4_real_composition.py) | none (design constants) | [`prepare_campaign.py`](realdata/prepare_campaign.py) |
| Fig. 4b–d | [`fig4_real_composition.py`](figures/scripts/fig4_real_composition.py) | `real_comp_gap.csv`†, `composition_distance_summary.csv`†, `real_comp_event_shifts.csv`† | [`run_common_panel.py`](realdata/runtime/run_common_panel.py) mode `n2`; [`summarize_common_panel.py`](realdata/summarize_common_panel.py); [`composition_resummary.py`](realdata/composition_resummary.py) |
| Fig. 4e (eight settings) | [`fig4_real_composition.py`](figures/scripts/fig4_real_composition.py) | `fig3_constructed_null_summary.csv`†, `precision_adni_rates.csv`†, `n2_rejection_rates.csv`† | ADNI 12-, 8- and 5-event and NACC amyloid PET settings: [`adni/run_null.py`](realdata/earlier_panels/adni/run_null.py), [`nacc/run_nacc.py`](realdata/earlier_panels/nacc/run_nacc.py), [`aggregate_earlier_panels.py`](realdata/earlier_panels/aggregate_earlier_panels.py); four MRI–cognition settings: [`mri_cognition/run_shared_panel.py`](realdata/earlier_panels/mri_cognition/run_shared_panel.py), [`aggregate_completed_results.py`](realdata/earlier_panels/mri_cognition/aggregate_completed_results.py) |
| Fig. 5a,b | [`fig5_apoe.py`](figures/scripts/fig5_apoe.py) | `APOE_observed_orderings.csv`†, `APOE_pair_tests.csv`† | [`run_common_panel.py`](realdata/runtime/run_common_panel.py) mode `apoe`; [`summarize_common_panel.py`](realdata/summarize_common_panel.py) |
| Fig. 5c–g | [`fig5_apoe.py`](figures/scripts/fig5_apoe.py) | `event_pair_order_summary.csv`†, `bootstrap_orderings.csv`†; optional checks `apoe_cp_e4_prefix_frequencies.csv`†, `apoe_cp_block_by_group.csv`† | [`run_common_panel.py`](realdata/runtime/run_common_panel.py) mode `bootstrap`; [`summarize_common_panel.py`](realdata/summarize_common_panel.py); [`apoe_refit_summaries.py`](realdata/apoe_refit_summaries.py) |
| Results: ε4 and ε3/ε3 orderings from the same refit (70% ADNI, 87% NACC) | no figure | `apoe_cp_block_contrasts.csv`† (statistic `S_csf_before_cog_pairs_0to9`, column `frac_positive`) | [`apoe_refit_summaries.py`](realdata/apoe_refit_summaries.py) |
| Table 1 | [`make_main_tables.py`](figures/tables/make_main_tables.py) | none | |
| Table 2 | [`make_main_tables.py`](figures/tables/make_main_tables.py) | `Table2_demographics_aggregate.csv`†, `reference_counts_weights_ESS.csv`† | [`table2_demographics.py`](realdata/table2_demographics.py); [`summarize_common_panel.py`](realdata/summarize_common_panel.py) |
| Table 3 | [`make_main_tables.py`](figures/tables/make_main_tables.py) | `APOE_pair_tests.csv`† | as Fig. 5a,b |
| Supplementary Table 1 | [`make_si_sim_tables.py`](figures/tables/make_si_sim_tables.py) | [`S1_event_parameters.csv`](figures/inputs/simulation/S1_event_parameters.csv), [`group_stage_parameters.csv`](figures/inputs/simulation/group_stage_parameters.csv) | constants of [`design_v2.py`](simulation/scripts/design_v2.py) |
| Supplementary Table 2 | [`make_si_sim_tables.py`](figures/tables/make_si_sim_tables.py) | [`rates.csv`](figures/inputs/simulation/rates.csv) | as Fig. 3b–d |
| Supplementary Table 3 | [`make_si_cohort_tables.py`](figures/tables/make_si_cohort_tables.py) | `S5_selection_flow.csv`† | [`selection_flow.py`](realdata/cohort/selection_flow.py), from `coverage_audit.json` of [`adni_audit.py`](realdata/cohort/adni_audit.py) and `coverage_summary.json` of [`nacc_audit.py`](realdata/cohort/nacc_audit.py) |
| Supplementary Table 4 | [`make_si_cohort_tables.py`](figures/tables/make_si_cohort_tables.py) | `reference_counts_weights_ESS.csv`† | [`summarize_common_panel.py`](realdata/summarize_common_panel.py) |
| Supplementary Table 5 | [`make_si_cohort_tables.py`](figures/tables/make_si_cohort_tables.py) | `APOE_pair_tests.csv`† (panel `omit_ptau`) | [`run_common_panel.py`](realdata/runtime/run_common_panel.py) mode `apoe`, panel `omit_ptau`; [`summarize_common_panel.py`](realdata/summarize_common_panel.py) |
| Supplementary Fig. 1 | [`suppfig1_pairwise_calibration.py`](figures/scripts/suppfig1_pairwise_calibration.py) | [`rates_exact_final.txt`](figures/inputs/simulation/rates_exact_final.txt), [`Fig3f_partial_null.csv`](figures/inputs/simulation/Fig3f_partial_null.csv) | runs `method_pair_complete_a` (a) and `method_pair_partial_a` (b); [`rates_exact_all.py`](simulation/scripts/dev/rates_exact_all.py) |
| Supplementary Fig. 2 | [`suppfig2_controlled_distance.py`](figures/scripts/suppfig2_controlled_distance.py) | `composition_distance_summary.csv`†, `composition_paired_contrast_summary.csv`† | as Fig. 4b–d |
| Supplementary Fig. 3 | [`suppfig3_earlier_panels.py`](figures/scripts/suppfig3_earlier_panels.py) | `apoe_hist_amyloid_position.csv`† | [`adni/run_observed.py`](realdata/earlier_panels/adni/run_observed.py) (stage `a2s`), [`nacc/run_nacc_stability.py`](realdata/earlier_panels/nacc/run_nacc_stability.py), [`export_stability_source.py`](realdata/earlier_panels/nacc/export_stability_source.py), [`aggregate_earlier_panels.py`](realdata/earlier_panels/aggregate_earlier_panels.py), [`apoe_refit_summaries.py`](realdata/apoe_refit_summaries.py) |

Every figure and table script checks key numbers against the values printed in the manuscript and
stops if one differs.
[`figures/README.md`](figures/README.md) describes each input file and its columns.

## Seeds and permutation budget

| Analysis | Seeds |
|---|---|
| Simulated datasets | dataset *i* (0-based) of a run has seed `base_seed + i`. `base_seed`: 42100000 `confirm_core_a`; 42200000 `confirm_power_global_a`; 53000000 `method_mechanism_a`; 53100000 `method_calibration_a`, `method_pair_complete_a` and `method_kde_calibration_a` (first 300 seeds only); 53200000 `method_pair_partial_a`; 53400000 `method_saebm_reference_a`, `method_kde_mechanism_a` and the SA-EBM run `method_saebm_a`; 53700000 `method_mechanism_16n_a`; 53800000 `method_alt_geometry_a` |
| Three-event example (Fig. 1d) | NumPy seed 20260914, 400,000 draws per stage |
| Mean-position gap intervals (Fig. 2a,d) | 2,000 bootstrap resamples of datasets, seed 20260919 |
| Benchmark datasets | reused datasets keep the seeds of their runs (53100000, 42200000, 53800000 + *i*); new datasets 92110000, 92120000, 92130000 (correlated noise) and 92210000, 92220000, 92230000 (heavy-tailed noise) + *i* |
| Relabelings (simulations, benchmark, *APOE* comparisons) | relabeling *p* of a dataset or comparison with seed *s* uses `numpy.random.SeedSequence([s, tag, p, 61000000])`, where `tag` is the first 4 bytes (little-endian) of SHA-256 of the permutation scheme name; all models therefore see the same relabelings |
| Benchmark paired intervals (Fig. 3d) | 10,000 bootstrap resamples, seed 92600000 + 1000 × condition index + 100 × contrast index + endpoint index ([`benchmark/README.md`](benchmark/README.md)) |
| Controlled composition experiment (Fig. 4a–d) | draws `SeedSequence([2026092300 + cohort, replicate])`, cohort 1 = ADNI, 2 = NACC; fits 2026231000 (ADNI) and 2026232000 (NACC) + replicate |
| Intervals of Fig. 4b–d and Supplementary Fig. 2 | 10,000 resamples of the 200 draws: `random.Random(20260924)` (ADNI) and `random.Random(20260925)` (NACC) for the gap and event shifts; `SeedSequence([20260923, 88000000, cohort])`, cohort 0 = ADNI, 1 = NACC, for the discordant pairs and paired differences |
| *APOE* comparisons (Fig. 5a,b; Table 3; Supplementary Table 5) | 2026239001 (ADNI), 2026239002 (NACC) |
| Bootstrap refits (Fig. 5c–g) | 2026238001 (ADNI), 2026238002 (NACC) + replicate |
| Earlier populations (Fig. 4e; Supplementary Fig. 3) | 20260914 (ADNI 12-event population), 20260917 (NACC amyloid PET population), 20260920 (MRI–cognition populations) |

All permutation tests use B = 599 relabelings and P = (1 + E)/600, where E counts relabelings whose
distance between orderings is at least the observed distance (ties included). Comparisons among three
groups are Bonferroni-corrected for three pairs and rejected when P ≤ 0.05/3 (E ≤ 9). The single
two-group comparisons of Fig. 4e are rejected when P ≤ 0.05 (E ≤ 29). Each directory README gives
the seeds in full.

## Models and code labels

The code and the stored results keep the labels under which the analyses ran. Several of these
files are checked by hash, so their labels are explained here instead of being renamed.

| Label in code or results | Meaning in the article |
|---|---|
| `repaired`, `separate`, `standard` | separately fitted DEBM (pyebm 2.0.3 with group-specific mixtures and orderings); `repaired` marks the continued adjacent-swap search |
| `original`, `standard_original` | separately fitted DEBM with the unmodified pyebm search |
| `shared` | pooled-score DEBM (one mixture per biomarker fitted to all participants; unweighted group orderings) |
| `invariant_min`, `concord`, `concord_min`, `min` | CONCORD. In the simulations and the earlier populations the common proportions are the smallest proportion of each diagnosis across groups, rescaled to sum to one; in the six-event panel they are computed by [`prepare_campaign.py`](realdata/prepare_campaign.py) and fixed in [`realdata/runtime/campaign.json`](realdata/runtime/campaign.json) |
| `invariant_pooled`, `pooled`, `kde_gmm_pooled` | variants with pooled proportions or pooled mixtures (not reported) |
| `kde_gmm`, `likelihood_ebm` | likelihood EBM (kde_ebm Gaussian-mixture event models fitted in each group; greedy search for the maximum-likelihood sequence) |
| `saebm_conjugate` | stage-aware EBM (SA-EBM; pysaebm with conjugate priors) |
| scheme `unrestricted`; suffix `_U`; `U-all` | unrestricted permutation |
| scheme `diagnosis`; suffix `_D`; `D-all` | within-diagnosis permutation |
| schemes `diagnosis_pair_*`; `D-pair` | pairwise within-diagnosis permutation (the third group keeps its labels) |
| `standard_U`, `concord_D` | separately fitted DEBM with unrestricted permutation and CONCORD with within-diagnosis permutation (Fig. 4e) |
| `e2`, `e33`, `e4`; `g1`, `g2`, `g3` in simulation files | CN-heavy, intermediate and AD-heavy simulated groups (the simulated data store the group in a column named `APOE`) |
| `e2`, `e33`, `e4` in ADNI and NACC files | *APOE* ε2 carriers without ε4, ε3/ε3 carriers, ε4 carriers without ε2 |
| `ABETA`, `PTAU`, … `Precuneus` in the generator | biomarkers 1–14 in the input order of `design_v2.BIOMARKERS` (1–5 modeled on CSF markers, 6–7 on cognitive scores, 8–14 on imaging measures) |
| `arm`, `arms` | a fitted model (column `arm` of the benchmark tables) or a combination of model and permutation scheme (real-data files) |
| `Standard` | separately fitted DEBM (in [`benchmark/runtime/campaign.json`](benchmark/runtime/campaign.json) and in older figure input files) |
| `reference`; `target` | the common proportions (real-data files); in the output of the three-event example, `target` is the best-fitting ordering |

**Stopping rule of the mixture fits.** Separately fitted DEBM used the pyebm 2.0.3 default rule,
which stops when the mean change in mixing proportions falls below 0.01 in the first group, in
Fig. 2, Fig. 3a and Fig. 4e. All other analyses required every group's mixtures to meet this
criterion (the all-groups rule of [`benchmark/runtime/allgroup_engine.py`](benchmark/runtime/allgroup_engine.py)
and [`realdata/runtime/allgroup_engine.py`](realdata/runtime/allgroup_engine.py)). Pooled-score DEBM
and CONCORD fit one pooled mixture, for which the two rules coincide.

**Consensus search.** All DEBM variants used the pyebm 2.0.3 adjacent-swap search continued until
no swap lowers the loss (pyebm 2.0.3 stops after the first accepted swap). The exception is
separately fitted DEBM in the ADNI 12-, 8- and 5-event and NACC amyloid PET settings of Fig. 4e,
which used the unmodified pyebm search.

## Software

| File | Used for | Contents |
|---|---|---|
| [`requirements/analysis.txt`](requirements/analysis.txt) | [`simulation/`](simulation/), [`benchmark/`](benchmark/), [`realdata/`](realdata/) | Python 3.12.14, NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3, scikit-learn 1.5.2, statsmodels 0.14.6, pyebm 2.0.3, kde_ebm 0.0.3, concord-ebm 0.1.1 (commit `40fc046`, for the ADNI 12-event and NACC amyloid PET scripts), openpyxl |
| [`requirements/saebm.txt`](requirements/saebm.txt) | SA-EBM ([`saebm_mechanism_dev.py`](simulation/scripts/dev/saebm_mechanism_dev.py)), in a separate environment | pysaebm 7.7.7 |
| [`requirements/figures.txt`](requirements/figures.txt) | [`figures/`](figures/) | Python 3.12, matplotlib 3.11.2, NumPy 2.5.3 |

```bash
python3.12 -m venv .venv-analysis && .venv-analysis/bin/pip install -r requirements/analysis.txt
python3.12 -m venv .venv-saebm    && .venv-saebm/bin/pip install -r requirements/saebm.txt
python3.12 -m venv .venv-figures  && .venv-figures/bin/pip install -r requirements/figures.txt
```

pyebm is used unmodified. [`dependency_v2.py`](simulation/scripts/dependency_v2.py) checks its
version and the hashes of three of its source files before every run, and the formal runs of the
benchmark and of the six-event panel refuse to start unless the package versions equal those
recorded in [`benchmark/runtime/EXECUTION_GATE.json`](benchmark/runtime/EXECUTION_GATE.json) and
[`realdata/runtime/EXECUTION_GATE.json`](realdata/runtime/EXECUTION_GATE.json). The figures use the
Arial font (see [`figures/README.md`](figures/README.md)).

## Reproducing the results

Commands are run from the repository root, except those in
[`realdata/README.md`](realdata/README.md), which are run from `realdata/`. The Slurm scripts were
written for one university cluster; set the account, partition and resources for yours, and the
Python executable through `CONCORD_PYTHON` (`SAEBM_PYTHON` for the SA-EBM environment).

**Figures and tables from the aggregate inputs.** No model is fitted.

```bash
for s in figures/scripts/fig*.py figures/scripts/suppfig*.py; do .venv-figures/bin/python "$s"; done
for s in figures/tables/make_*.py; do .venv-figures/bin/python "$s"; done
```

Figs 1–3, Supplementary Fig. 1 and Supplementary Tables 1 and 2 run from the files in the
repository. The other scripts need the ADNI and NACC aggregates in
[`figures/inputs/realdata/`](figures/inputs/realdata/) and otherwise stop with a list of the missing
files.

**Simulations.** Simulated datasets and per-fit records are not stored; every dataset is regenerated
from its seed. Run the unit tests, fit the runs chunk by chunk with
[`run_v2.py`](simulation/scripts/run_v2.py), then build the figure inputs with the scripts in
[`simulation/summaries/`](simulation/summaries/) and
[`simulation/scripts/dev/`](simulation/scripts/dev/):

```bash
.venv-analysis/bin/python -m unittest discover -s simulation/scripts -p 'test_*.py'
.venv-analysis/bin/python simulation/scripts/run_v2.py --config simulation/configs/method_calibration_a.json \
    --cell-index 0 --chunk 0 --nchunks 250 --workers 16 --output "$CONCORD_WORK_DIR/runs/method_calibration_a/cell_0/chunk_0"
```

[`simulation/README.md`](simulation/README.md) lists every run, the Slurm workflow
([`make_snapshot.py`](simulation/hpc/make_snapshot.py), [`run_array.sbatch`](simulation/hpc/run_array.sbatch))
and the summary commands.

**Paired benchmark.** [`run_batch.py`](benchmark/runtime/run_batch.py) fits one batch of datasets
(one array task of [`formal.sbatch`](benchmark/formal.sbatch)), and
[`summarize_campaign_audited.py`](benchmark/summarize_campaign_audited.py) builds `rates.csv` and
`paired_differences.csv` from the results. The figure scripts read copies of these two files from
[`figures/inputs/simulation/`](figures/inputs/simulation/). The benchmark reuses datasets and fitted
results of three earlier simulation runs; the records it needs are in
[`benchmark/runtime/inputs/historical_reused_rows.jsonl`](benchmark/runtime/inputs/historical_reused_rows.jsonl)
(simulation results only). Details are in [`benchmark/README.md`](benchmark/README.md).

**ADNI and NACC analyses.** These need access to ADNI and NACC data. Six-event panel: the cohort
scripts in [`realdata/cohort/`](realdata/cohort/), then
[`prepare_campaign.py`](realdata/prepare_campaign.py), the model fits with
[`run_common_panel.py`](realdata/runtime/run_common_panel.py) (modes `n2`, `apoe` and `bootstrap`;
[`realdata/formal.sbatch`](realdata/formal.sbatch)), then
[`summarize_common_panel.py`](realdata/summarize_common_panel.py),
[`composition_resummary.py`](realdata/composition_resummary.py),
[`apoe_refit_summaries.py`](realdata/apoe_refit_summaries.py) and
[`table2_demographics.py`](realdata/table2_demographics.py). The earlier populations follow the same
pattern in [`realdata/earlier_panels/`](realdata/earlier_panels/). Scripts that read
participant-level data write to a folder outside the repository, except the prepared inputs that the
two runtimes read next to their code; these are excluded by [`.gitignore`](.gitignore) and
[`realdata/.gitignore`](realdata/.gitignore).
[`realdata/README.md`](realdata/README.md) gives the input files, the order of the steps and the
arguments.

**Hash checks.** The model files in [`simulation/scripts/`](simulation/scripts/) are byte-identical
to the copies that ran in the benchmark and in the six-event panel. Their SHA-256 hashes are recorded
in [`benchmark/runtime/EXECUTION_GATE.json`](benchmark/runtime/EXECUTION_GATE.json) and
[`realdata/runtime/EXECUTION_GATE.json`](realdata/runtime/EXECUTION_GATE.json). The six-event
runtime checks them before a formal run, and [`benchmark/README.md`](benchmark/README.md) shows how
to compare them with the benchmark record.

## Data

ADNI data are available to qualified researchers at <https://adni.loni.usc.edu> after approval of a
data use application. NACC data are available at <https://naccdata.org> after approval of a data
request. No participant-level data are stored in this repository.

- The simulated data are regenerated from the seeds. The simulation aggregates read by the figure
  scripts are in [`figures/inputs/simulation/`](figures/inputs/simulation/).
- Aggregate results derived from ADNI and NACC (counts, orderings, rates and intervals) will be
  provided with the Source Data of the published article. The real-data figure and table scripts
  expect them in [`figures/inputs/realdata/`](figures/inputs/realdata/).
- [`.gitignore`](.gitignore) blocks the cohort data folders, prepared runtime inputs, de-identified
  tables (`*_deid*`), draw-index files (`*.npz`), spreadsheets, and all CSV and JSONL files except
  the simulation-only inputs in [`figures/inputs/simulation/`](figures/inputs/simulation/) and
  [`benchmark/runtime/inputs/`](benchmark/runtime/inputs/).

## Licence and citation

GNU General Public License v3.0 or later (GPL-3.0-or-later); see [`LICENSE`](LICENSE). Two files
adapted from pyebm, [`fast_likelihood_v2.py`](simulation/scripts/fast_likelihood_v2.py) and the
vendored `realdata/earlier_panels/mri_cognition/vendor/concord/likelihood.py`, keep their
GPL-3.0-only notice. If you use this code, please cite the article above.
