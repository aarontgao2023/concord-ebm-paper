# Figures and tables

Code that draws the five figures, the eight tables and the three supplementary figures of the
CONCORD paper from stored aggregate results. No script fits a model, runs a permutation or reads
participant-level data. Every script checks the key numbers against the values printed in the
manuscript (Python `assert`) and stops if one differs.

```
figures/
  scripts/            one script per figure; nc_style.py holds the shared style and paths
  tables/             make_main_tables.py (Tables 1-3), make_si_sim_tables.py (Supplementary
                      Tables 1-2), make_si_cohort_tables.py (Supplementary Tables 3-5)
  inputs/simulation/  simulation-only aggregate inputs (in the repository)
  inputs/realdata/    ADNI/NACC aggregate inputs (not in the repository; see below)
  output/             figures (PDF, PNG, SVG) and tables (TeX, CSV); generated, not versioned
  source_data/        one CSV per figure panel with the plotted numbers; generated, not versioned
```

## Display items

Inputs marked *sim* are in `inputs/simulation/`. Inputs marked *real* are aggregate ADNI or NACC
results that must be placed in `inputs/realdata/`; without them the script stops and lists the
missing files.

| Display item | Script | Inputs |
|---|---|---|
| Fig. 1a, b (schematics) | `scripts/fig1_overview.py` | none; the illustrative weights are computed in the script |
| Fig. 1c (diagnostic composition) | `scripts/fig1_overview.py` | CN/MCI/AD counts of Table 2, written in the script; simulated group sizes (design constants). If present, *real* `reference_counts_weights_ESS.csv` is used to check the counts |
| Fig. 1d (three-event example) | `scripts/fig1_overview.py` | *sim* `toy_target_shift.txt`; `simulation/scripts/dev/toy_target_shift.py` (read for the two compositions) |
| Fig. 2a, c (complete data, fixed true sequence) | `scripts/fig2_simulation_composition.py` | *sim* `mechanism_all_cells.csv`, `v2_mechanism_gap_intervals.csv` |
| Fig. 2b (event shift vs separability) | `scripts/fig2_simulation_composition.py` | *sim* `runs_random_truth_eventwise.csv`, `S1_event_parameters.csv` |
| Fig. 2d (sample size) | `scripts/fig2_simulation_composition.py` | *sim* `mechanism_all_cells.csv`, `v2_mechanism_gap_intervals.csv` |
| Fig. 3a (groups differ only in composition) | `scripts/fig3_permutation_testing.py` | *sim* `fig2_calibration.csv` |
| Fig. 3b-d (paired benchmark) | `scripts/fig3_permutation_testing.py` | *sim* `rates.csv`, `paired_differences.csv` |
| Fig. 4a (design) | `scripts/fig4_real_composition.py` | none (design constants); the script needs the inputs of 4b-e |
| Fig. 4b (mean-position gap) | `scripts/fig4_real_composition.py` | *real* `real_comp_gap.csv` |
| Fig. 4c (discordant event pairs) | `scripts/fig4_real_composition.py` | *real* `composition_distance_summary.csv` |
| Fig. 4d (event shifts) | `scripts/fig4_real_composition.py` | *real* `real_comp_event_shifts.csv` |
| Fig. 4e (groups drawn from ε3/ε3 carriers) | `scripts/fig4_real_composition.py` | *real* `fig3_constructed_null_summary.csv`, `precision_adni_rates.csv`, `n2_rejection_rates.csv` |
| Fig. 5a, b (APOE orderings) | `scripts/fig5_apoe.py` | *real* `APOE_observed_orderings.csv`, `APOE_pair_tests.csv` |
| Fig. 5c, d (pairwise precedence) | `scripts/fig5_apoe.py` | *real* `event_pair_order_summary.csv` |
| Fig. 5e-g (bootstrap refits) | `scripts/fig5_apoe.py` | *real* `bootstrap_orderings.csv`; optional cross-checks *real* `apoe_cp_e4_prefix_frequencies.csv`, `apoe_cp_block_by_group.csv` |
| Table 1 (models compared) | `tables/make_main_tables.py` | none |
| Table 2 (participants) | `tables/make_main_tables.py` | *real* `Table2_demographics_aggregate.csv`, `reference_counts_weights_ESS.csv` |
| Table 3 (APOE comparisons) | `tables/make_main_tables.py` | *real* `APOE_pair_tests.csv` |
| Supplementary Fig. 1a (pairwise test, no group differs) | `scripts/suppfig1_pairwise_calibration.py` | *sim* `rates_exact_final.txt` |
| Supplementary Fig. 1b (pairwise test, one group changed) | `scripts/suppfig1_pairwise_calibration.py` | *sim* `Fig3f_partial_null.csv`, checked against `rates_exact_final.txt` |
| Supplementary Fig. 2 (distance between group orderings) | `scripts/suppfig2_controlled_distance.py` | *real* `composition_distance_summary.csv`, `composition_paired_contrast_summary.csv` |
| Supplementary Fig. 3 (amyloid placed first, other panels) | `scripts/suppfig3_earlier_panels.py` | *real* `apoe_hist_amyloid_position.csv` |
| Supplementary Table 1 (simulated biomarkers) | `tables/make_si_sim_tables.py` | *sim* `S1_event_parameters.csv`, `group_stage_parameters.csv` |
| Supplementary Table 2 (paired benchmark) | `tables/make_si_sim_tables.py` | *sim* `rates.csv` |
| Supplementary Table 3 (selection of participants) | `tables/make_si_cohort_tables.py` | *real* `S5_selection_flow.csv` |
| Supplementary Table 4 (weights, effective sample sizes) | `tables/make_si_cohort_tables.py` | *real* `reference_counts_weights_ESS.csv` |
| Supplementary Table 5 (APOE comparisons without p-tau) | `tables/make_si_cohort_tables.py` | *real* `APOE_pair_tests.csv` |

The table scripts write the table bodies exactly as in the manuscript. Captions and table notes
are edited in the manuscript; the scripts carry a copy of the submitted text.

## Inputs

### Simulation (`inputs/simulation/`, in the repository)

| File | Content | Produced by |
|---|---|---|
| `mechanism_all_cells.csv` | per run, cell and model: mean-position gap, misordered pairs, errors (fixed true sequence) | `simulation/summaries/build_mechanism_source_data.py` |
| `v2_mechanism_gap_intervals.csv` | 95% bootstrap intervals of the mean-position gap | `simulation/summaries/mechanism_gap_intervals.py` |
| `runs_random_truth_eventwise.csv` | per-event position shift, CN-heavy minus AD-heavy group (random true sequence) | `simulation/summaries/runs_eventwise.py` |
| `S1_event_parameters.csv` | component AUC and observation probability of biomarkers 1-14 | constants of `simulation/scripts/design_v2.py` (`BIOMARKERS`) |
| `group_stage_parameters.csv` | CN/MCI/AD counts and stage distributions of the three simulated groups | constants of `design_v2.py` (`REFERENCE_COUNTS`, `BASE_STAGE_BETA`) |
| `fig2_calibration.csv` | rejections with unrestricted and within-diagnosis permutation, one common ordering (1,000 datasets), with Wilson 95% intervals | cell `REF_H0`: `simulation/scripts/dev/core_rates_exact.py confirm_core_a` (separately fitted DEBM) and `simulation/scripts/dev/rates_exact_all.py method_calibration_a` (pooled-score DEBM, CONCORD); column `source` names the run |
| `rates_exact_final.txt` | rejection counts and rates of pairwise within-diagnosis permutation in runs `method_pair_complete_a` and `method_pair_partial_a` | the two blocks printed by `simulation/scripts/dev/rates_exact_all.py method_pair_complete_a method_pair_partial_a`, under a one-line header |
| `Fig3f_partial_null.csv` | CONCORD rejections of the comparison between the two unchanged groups when one group is changed (27 of 91 pairs reversed; 500 datasets per changed group), with Wilson 95% intervals | run `method_pair_partial_a`, cells `PWR_E2_K27`, `PWR_E33_K27` and `PWR_E4_K27`: the counts of `simulation/scripts/dev/rates_exact_all.py` for the two unchanged groups, written to this file by an earlier version of the figure script; `suppfig1_pairwise_calibration.py` checks them against the `method_pair_partial_a` block of `rates_exact_final.txt` |
| `rates.csv`, `paired_differences.csv` | paired benchmark: rates with Wilson intervals, paired differences with paired-bootstrap intervals | `benchmark/summarize_campaign_audited.py` |
| `toy_target_shift.txt` | expected ranking loss of the six orderings in the three-event example | `simulation/scripts/dev/toy_target_shift.py` |

`mechanism_all_cells.csv`, `v2_mechanism_gap_intervals.csv` and `runs_random_truth_eventwise.csv`
also keep rows of additional runs that the figures do not read (see `simulation/README.md`).

### ADNI and NACC (`inputs/realdata/`, not in the repository)

These are aggregate results (counts, orderings, rates and intervals) of the analyses in
`realdata/`; they contain no participant-level data. They are not distributed with the code (see
the Data availability statement of the article).

| File | Used by | Produced by |
|---|---|---|
| `reference_counts_weights_ESS.csv` | Fig. 1c (check), Table 2, Supplementary Table 4 | `realdata/summarize_common_panel.py` |
| `APOE_observed_orderings.csv` | Fig. 5a, b | `realdata/summarize_common_panel.py` |
| `APOE_pair_tests.csv` | Fig. 5a, b, Table 3, Supplementary Table 5 | `realdata/summarize_common_panel.py` |
| `event_pair_order_summary.csv` | Fig. 5c, d | `realdata/summarize_common_panel.py` |
| `bootstrap_orderings.csv` | Fig. 5e-g | `realdata/summarize_common_panel.py` |
| `composition_distance_summary.csv` | Fig. 4c, Supplementary Fig. 2 | `realdata/summarize_common_panel.py` |
| `composition_paired_contrast_summary.csv` | Supplementary Fig. 2 | `realdata/summarize_common_panel.py` |
| `real_comp_gap.csv`, `real_comp_event_shifts.csv` | Fig. 4b, d | `realdata/composition_resummary.py` |
| `apoe_hist_amyloid_position.csv` | Supplementary Fig. 3 | `realdata/apoe_refit_summaries.py` |
| `apoe_cp_e4_prefix_frequencies.csv`, `apoe_cp_block_by_group.csv` | Fig. 5e-g (optional cross-check) | `realdata/apoe_refit_summaries.py` |
| `Table2_demographics_aggregate.csv` | Table 2 | `realdata/table2_demographics.py` |
| `S5_selection_flow.csv` | Supplementary Table 3 | `realdata/cohort/selection_flow.py`, from `coverage_audit.json` (`adni_audit.py`) and `coverage_summary.json` (`nacc_audit.py`) |
| `fig3_constructed_null_summary.csv`, `precision_adni_rates.csv` | Fig. 4e | summaries of the earlier ADNI and NACC panels (`realdata/earlier_panels/`) |
| `n2_rejection_rates.csv` | Fig. 4e | summary of the MRI-cognition panels (`realdata/earlier_panels/mri_cognition/`) |

### Labels in the input files

The input files use the labels of the analysis code:

- Models: `kde_gmm` or `likelihood_ebm` = likelihood EBM; `saebm_conjugate` = stage-aware EBM
  (SA-EBM); `repaired` or `separate` = separately fitted DEBM; `shared` = pooled-score DEBM;
  `invariant_min` or `concord` = CONCORD.
- Permutation schemes: `unrestricted` = unrestricted permutation; `diagnosis` = within-diagnosis
  permutation; `diagnosis_pair_*` = pairwise within-diagnosis permutation.
- `fig2_calibration.csv`: `Standard` = separately fitted DEBM,
  `Pooled only` = pooled-score DEBM; `U-all` and `D-all` = unrestricted and within-diagnosis
  permutation.
- Simulated groups: `g1`, `g2`, `g3` (and, in simulation files, `e2`, `e33`, `e4`) = CN-heavy,
  intermediate and AD-heavy group. In ADNI and NACC files, `e2`, `e33` and `e4` are the APOE
  groups (ε2 carriers without ε4, ε3/ε3, ε4 carriers without ε2).
- Simulated biomarkers: the generator's code names (`ABETA` ... `Precuneus`) are biomarkers 1-14
  in the order of `design_v2.BIOMARKERS` (1-5 modeled on CSF markers, 6-7 on cognitive scores,
  8-14 on imaging measures).

## Running

From the repository root, with the environment below:

```bash
for s in figures/scripts/fig*.py figures/scripts/suppfig*.py; do python "$s"; done
for s in figures/tables/make_*.py; do python "$s"; done
```

Figs 1-3, Supplementary Fig. 1 and Supplementary Tables 1-2 run from the files in the
repository. The other scripts need the ADNI and NACC files in `inputs/realdata/`;
`make_main_tables.py` writes Table 1 and then stops if they are missing.

Each figure script prints `[QA] ... 0 text issue(s)` after checking that all text is 5-7 pt
(panel letters 8 pt), at least 1.5 mm from the edge and not overlapping other text.

Paths can be changed with environment variables:

| Variable | Default | Meaning |
|---|---|---|
| `CONCORD_FIG_INPUTS` | `figures/inputs` | folder containing `simulation/` and `realdata/` |
| `CONCORD_OUTPUT_DIR` | `figures/output` | figures; tables go to its `tables/` subfolder |
| `CONCORD_SOURCE_DATA_DIR` | `figures/source_data` | per-panel source-data CSVs |
| `CONCORD_FONT_DIR` | system font folders | folder with `Arial.ttf`, `Arial Bold.ttf`, `Arial Italic.ttf` |

## Environment

- Python 3.12, matplotlib 3.11 and NumPy 2.5 (the published figures were made with matplotlib
  3.11.2 and NumPy 2.5.3). The table scripts need NumPy only.
- Arial Regular, Bold and Italic. The scripts look for them in the system font folders or in
  `CONCORD_FONT_DIR`. Without Arial, matplotlib uses Helvetica or DejaVu Sans; the figures are
  still drawn, but text widths and the layout checks differ.

With this environment and the Arial fonts supplied with macOS, the PNG files are pixel-identical
to the published figures, and the table scripts reproduce the manuscript's table files exactly.
