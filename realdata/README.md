# ADNI and NACC analyses

Code for the ADNI and NACC analyses of the CONCORD paper: Figs 1c, 4 and 5, Tables 2 and 3,
Supplementary Figs 2 and 3 and Supplementary Tables 3 to 5. The repository holds code only.
ADNI data are available from https://adni.loni.usc.edu and NACC data from https://naccdata.org
under their data use agreements; no participant-level data and no results computed from them are
included, apart from the published genotype-by-diagnosis counts and common proportions in
`runtime/campaign.json`. Scripts that read participant-level data write to a folder you choose outside the
repository; the only exception is the prepared input of the two runtimes, which must sit next to
the runtime code and is excluded by `realdata/.gitignore`. The files that the figure and table
scripts read are aggregates (counts, orderings, rates, intervals).

```
realdata/
  cohort/                      six-event panel: anchor visits, coverage, de-identified inputs, selection flow
  prepare_campaign.py          covariate adjustment, common proportions, draws, runtime/campaign.json
  runtime/                     model fits: composition draws, APOE tests, bootstrap refits
  formal.sbatch                Slurm array for the runtime (56 batches)
  summarize_common_panel.py    aggregate tables from the stored fit records
  composition_resummary.py     mean-position gap and event shifts (Fig. 4b,d)
  apoe_refit_summaries.py      bootstrap-refit summaries (Fig. 5e-g, Results, Supplementary Fig. 3)
  table2_demographics.py       Table 2
  earlier_panels/              the eight settings of Fig. 4e and the two panels of Supplementary Fig. 3
    adni/                      ADNI 12-event population (12-, 8- and 5-event settings; bootstrap refits)
    nacc/                      NACC amyloid PET population (setting and bootstrap refits)
    mri_cognition/             ADNI and NACC MRI-cognition populations (8- and 4-event settings)
    aggregate_earlier_panels.py
```

Paths and commands in this README are relative to `realdata/`: run the commands from that folder
(`cd realdata`), or prefix the script paths with `realdata/` when running from the repository root.

## Software

| Part | Requirements |
|---|---|
| `cohort/`, `prepare_campaign.py`, `table2_demographics.py` | Python 3.12, NumPy, pandas |
| `runtime/` | Python 3.12.14, NumPy 1.26.4, SciPy 1.13.1, pandas 2.2.3, scikit-learn 1.5.2, statsmodels 0.14.6, pyebm 2.0.3, kde_ebm 0.0.3; the engines in `../simulation/scripts` |
| `summarize_common_panel.py` | Python 3.12, NumPy |
| `composition_resummary.py`, `apoe_refit_summaries.py`, `aggregate_earlier_panels.py` | Python standard library only |
| `earlier_panels/adni/`, `earlier_panels/nacc/` | concord-ebm 0.1.1 at commit `40fc046` (`pip install "git+https://github.com/aarontgao2023/concord-ebm@40fc046f1c11e65228e9cf62d99a3d94fd4937bc"`), pyebm 2.0.3, NumPy, SciPy, pandas; `openpyxl` for `build_cohort.py` |
| `earlier_panels/mri_cognition/` | pyebm 2.0.3, NumPy, SciPy, pandas; the concord 0.1.1 code that ran is vendored in `vendor/concord` |

pyebm is verified against the hashes of version 2.0.3 before every fit. The runtime refuses its
production modes unless the Python, NumPy, SciPy, pandas, scikit-learn, statsmodels and pyebm
versions and the kde_ebm source hashes equal those recorded in `runtime/EXECUTION_GATE.json`.

## Data

Paths are passed as arguments or through these environment variables:

| Variable | Meaning |
|---|---|
| `CONCORD_ADNI_DIR` | folder with the ADNI files below |
| `CONCORD_NACC_DIR` | the NACC ADSP-PHC release folder (`ADSP-PHC-122024-investigator_NACC`) |
| `CONCORD_NACC_UDS` | a CSV export of the NACC Uniform Data Set (UDS) investigator file |
| `CONCORD_WORK_DIR` | working folder outside the repository for every output |
| `CONCORD_ENGINE_DIR` | optional; location of the engines if `runtime/` is copied elsewhere (default `../simulation/scripts`) |
| `CONCORD_PYTHON` | Python executable used by the Slurm scripts (`formal.sbatch` and the `.sbatch` files in `earlier_panels/`) |

ADNI files (names as downloaded; the dates in the names identify the releases used):

| Used by | Files |
|---|---|
| six-event panel | `ADSP_PHC_CSF_14Sep2026.csv`, `ADSP_PHC_COGN_14Sep2026.csv`, `ADSP_PHC_DEMODX_14Sep2026.csv`, `APOERES_14Sep2026.csv`, `ADSP_PHC_MERGED_DATADIC_20260521_14Sep2026.csv`; for the coverage audit also `ADSP_PHC_T1_FS_18Oct2025.csv`, `ADSP_PHC_PET_Amyloid_Simple_22Aug2025.csv`, `ADSP_PHC_PLASMA_14Sep2026.csv` |
| 12-event and MRI-cognition populations | `DXSUM_27Aug2025.csv`, `PTDEMOG_27Aug2025.csv`, `MMSE_18Jan2026.csv`, `ADNI1GO23_FS6_alltimepoints_04282021.xlsx`, `UCSFFSX7_27Feb2026.csv`, `ADSP_PHC_T1_FS_18Oct2025.csv`, `ADSP_PHC_CSF_14Sep2026.csv`, `ADSP_PHC_COGN_14Sep2026.csv`, `ADSP_PHC_DEMODX_14Sep2026.csv`, `APOERES_14Sep2026.csv` |

NACC files: from the ADSP-PHC NACC release, `Biomarker/NACC_ADSP_PHC_Biomarker_2024.csv`,
`Cognition/NACC_ADSP_PHC_Cognition_2024.csv`, `Imaging_T1/NACC_ADSP_PHC_T1_Freesurfer_2024.csv` and
`Imaging_PET/NACC_ADSP_PHC_Amyloid_Simple_2024.csv`; from the UDS, the columns `NACCID`,
`NACCVNUM`, `VISITYR`, `VISITMO`, `VISITDAY`, `NACCAPOE`, `NACCUDSD`, `NACCALZD`, `NACCAGE`, `SEX`,
`EDUC`, `NACCMMSE` and `NACCMOCA`.

Participant-level files stay outside the repository, except the prepared inputs that the two
runtimes read from `runtime/inputs/` and `earlier_panels/mri_cognition/inputs/`; `realdata/.gitignore`
excludes both folders, every `local_only/` folder and all `*_deid*.csv` and `*.npz` files. Only the
de-identified inputs (fresh sequential pseudonyms, no identifiers or dates) are needed on a
compute cluster; the linkage to source identifiers stays in `local_only/`.

## Six-event panel (Figs 1c, 4a-d and 5; Tables 2-3; Supplementary Tables 3-5; Supplementary Fig. 2)

Run in this order (the output folders are those used by the default arguments):

1. `cohort/adni_audit.py`, then `cohort/adni_export_core6.py` (both with `--data-dir` and `--out
   $CONCORD_WORK_DIR/audit/adni`): anchor visits, coverage audit (`coverage_audit.json`) and the
   de-identified input `raw_deid.csv` (1,203 participants).
2. `cohort/nacc_audit.py`, then `cohort/nacc_prepare_core6.py` (with `--phc-dir`, `--uds-file`,
   `--out $CONCORD_WORK_DIR/audit/nacc`): the same for NACC (`coverage_summary.json`;
   `prepared/nacc_csfcog6_raw_deid.csv`, 1,584 participants).
   Then `cohort/selection_flow.py` writes `$CONCORD_WORK_DIR/audit/S5_selection_flow.csv`
   (Supplementary Table 3; 8 rows: cohort, stage, participants) from the aggregate counts in
   `coverage_audit.json` and `coverage_summary.json`; its docstring lists the key of each count.
3. `prepare_campaign.py --work-dir $CONCORD_WORK_DIR`: adjusts the six events for age and sex (CSF)
   or age, sex and education (cognition) by OLS in the CN participants of each cohort, computes the
   common proportions (48.9% CN, 28.6% MCI, 22.5% AD), draws the 200 rosters of the controlled
   composition experiment and writes `runtime/inputs/` and `runtime/campaign.json`. With the same
   data and software it writes the published `runtime/campaign.json` unchanged.
4. Optional checks without model fits or with a small budget:
   `python runtime/validate_runtime_logic.py --config runtime/campaign.json --out <file.json>` and
   `python runtime/run_common_panel.py --config runtime/campaign.json --mode preflight --cohort adni --out <folder>`
   (B = 8, separate seeds; repeat for `nacc`).
5. Model fits: `formal.sbatch` runs the 56 batches of `runtime/batch_mapping.json` through
   `runtime/batch_launcher.py` (16 CPUs each). A single unit can also be run directly, for example
   `python runtime/run_common_panel.py --config runtime/campaign.json --gate runtime/EXECUTION_GATE.json --out $CONCORD_WORK_DIR/formal --mode apoe --cohort adni --panel CSFcog6 --workers 16`.
   Modes: `n2` (controlled composition experiment, 200 draws per cohort, four models, two
   compositions), `apoe` (observed CONCORD orderings and the three pairwise within-diagnosis
   permutation tests, 599 relabelings each; panels `CSFcog6` and `omit_ptau`) and `bootstrap`
   (200 refits per cohort). In total 804 result units and 10,792 fits.
6. `summarize_common_panel.py --formal-root $CONCORD_WORK_DIR/formal --new-output $CONCORD_WORK_DIR/summary`.
7. `composition_resummary.py --summary-dir $CONCORD_WORK_DIR/summary --out $CONCORD_WORK_DIR/resummary`.
8. `apoe_refit_summaries.py --summary-dir $CONCORD_WORK_DIR/summary --out $CONCORD_WORK_DIR/resummary`
   (add `--adni-bootstrap` and `--nacc-bootstrap` for Supplementary Fig. 3; see below).
9. `table2_demographics.py --inputs-dir runtime/inputs --out $CONCORD_WORK_DIR/table2`.

The runtime reproduces the fits exactly only with the recorded software versions: formal modes stop
if the SHA-256 of any runtime or engine file, of `campaign.json` or of the kde_ebm sources differs
from `runtime/EXECUTION_GATE.json`. The engines are read from `../simulation/scripts`; they are
byte-identical to the copies that ran and are recorded under their original names
(`vendor/current/...`). The hashes of `engine_adapter.py` and `run_common_panel.py` in the
published gate were recomputed after the two edits that point the runtime at those engines; all
other recorded values are those of the run. `runtime/vendor/w3/run_stop_isolation.py` is used as a
module only: `allgroup_engine.py` imports its class `AllGroupState` (the all-groups stopping rule).
Its standalone entry point (`__main__`) needs a script that is not in the repository.

## Earlier analysis populations (Fig. 4e; Supplementary Fig. 3)

ADNI 12-event population (1,070 participants; 12-, 8- and 5-event settings of Fig. 4e and the ADNI
panel of Supplementary Fig. 3):

1. `earlier_panels/adni/build_cohort.py --data-dir $CONCORD_ADNI_DIR` writes `adni_p12_deid.csv`.
2. `earlier_panels/adni/run_null.py --stage n2` for draws 0-199 (`hpc/adni_null.sbatch`), and with
   `--panel P8` and `--panel P5` (`hpc/adni_null_panel.sbatch`).
3. `earlier_panels/adni/run_observed.py --stage a2s --R 500` (`hpc/adni_observed.sbatch`).

NACC amyloid PET population (1,200 participants):

1. `earlier_panels/nacc/build_cohort_nacc.py` writes `nacc_p5_deid.csv`.
2. `earlier_panels/nacc/run_nacc.py --stage n2` for draws 0-199 (`hpc/nacc.sbatch`).
3. `earlier_panels/nacc/run_nacc_stability.py --reference min` (and `pooled`, not reported;
   `hpc/nacc_stability.sbatch`), then `earlier_panels/nacc/export_stability_source.py`.

ADNI (1,025 participants) and NACC (1,886 participants) MRI-cognition populations:

1. `earlier_panels/mri_cognition/prepare_adni.py` (uses the outputs of `adni/build_cohort.py`) and
   `prepare_nacc.py` (uses `nacc_cohort_full.csv` of `nacc/build_cohort_nacc.py` only to count the
   overlap).
2. `freeze_campaign.py`: adjustment, 200 paired draws per cohort and `campaign.json`.
3. `run_shared_panel.py` in modes `preflight`, `n2` and `apoe` (`preflight.sbatch`, `n2.sbatch`,
   `apoe.sbatch`); the APOE fits are not reported but the aggregation requires them.
4. `aggregate_completed_results.py` writes `n2_rejection_rates.csv`.

Finally `earlier_panels/aggregate_earlier_panels.py` writes `fig3_constructed_null_summary.csv`,
`precision_adni_rates.csv` and `fig4_ADNI_bootstrap_orderings.csv`, and `apoe_refit_summaries.py`
with `--adni-bootstrap <...>/fig4_ADNI_bootstrap_orderings.csv --nacc-bootstrap
<...>/r2_nacc_bootstrap_orderings.csv` writes `apoe_hist_amyloid_position.csv`.

In every Fig. 4e setting two disjoint groups of different composition are drawn without
replacement, within diagnosis, from the e3/e3 carriers of the population, 200 times; each draw is
tested once (B = 599, rejection if E <= 29). Separately fitted DEBM used the pyebm 2.0.3 default
stopping rule, which checks the first group only, with the unmodified pyebm search in the ADNI
12-, 8- and 5-event and NACC amyloid PET settings and the continued search in the four
MRI-cognition settings. In the six-event panel, every DEBM mixture stopped only when all groups
met the 0.01 rule, and every DEBM search was continued until no adjacent swap improved the fit.

## Which script produces which display item

The figure and table scripts are in `figures/`. The aggregate files they read are produced here:

| Display item | Aggregate file(s) | Produced by |
|---|---|---|
| Fig. 1c | `reference_counts_weights_ESS.csv` (rows `actual_APOE`) | `summarize_common_panel.py` |
| Fig. 4a | design constants (`runtime/campaign.json`: `composition_allocation`) | `prepare_campaign.py` |
| Fig. 4b, 4d | `real_comp_gap.csv`, `real_comp_event_shifts.csv` | `composition_resummary.py` |
| Fig. 4c, Supplementary Fig. 2 | `composition_distance_summary.csv`, `composition_paired_contrast_summary.csv` | `summarize_common_panel.py` |
| Fig. 4e | `fig3_constructed_null_summary.csv`, `precision_adni_rates.csv`; `n2_rejection_rates.csv` | `earlier_panels/aggregate_earlier_panels.py`; `earlier_panels/mri_cognition/aggregate_completed_results.py` |
| Fig. 5a,b; Table 3; Supplementary Table 5 | `APOE_observed_orderings.csv`, `APOE_pair_tests.csv` | `summarize_common_panel.py` |
| Fig. 5c,d | `event_pair_order_summary.csv` | `summarize_common_panel.py` |
| Fig. 5e-g | `bootstrap_orderings.csv`; checks `apoe_cp_e4_prefix_frequencies.csv`, `apoe_cp_block_by_group.csv` | `summarize_common_panel.py`; `apoe_refit_summaries.py` |
| Results, paired refits (70% ADNI, 87% NACC) | `apoe_cp_block_contrasts.csv` (statistic `S_csf_before_cog_pairs_0to9`, e33 vs e4, `frac_positive`) | `apoe_refit_summaries.py` |
| Table 2 | `Table2_demographics_aggregate.csv`, `Table2_chi_square.csv` | `table2_demographics.py` |
| Supplementary Table 3 | `S5_selection_flow.csv`, from `coverage_audit.json` (ADNI) and `coverage_summary.json` (NACC) | `cohort/selection_flow.py`, from the outputs of `cohort/adni_audit.py` and `cohort/nacc_audit.py` |
| Supplementary Table 4 | `reference_counts_weights_ESS.csv` | `summarize_common_panel.py` |
| Supplementary Fig. 3 | `apoe_hist_amyloid_position.csv` (rows `ADNI_K12`/`concord_min`, `NACC5_PET`/`min`) | `apoe_refit_summaries.py` |

## Seeds and budgets

| Analysis | Seeds |
|---|---|
| Controlled composition draws | `numpy.random.default_rng(SeedSequence([2026092300 + cohort, replicate]))`, cohort 1 = ADNI, 2 = NACC, replicates 0-199 (`prepare_campaign.py`) |
| Composition fits | `fit_seed_base + replicate`: 2026231000 (ADNI), 2026232000 (NACC) |
| APOE comparisons | `apoe_seed` 2026239001 (ADNI), 2026239002 (NACC); relabeling *p* of a pair uses `SeedSequence([apoe_seed, tag, p, 61000000])`, tag = first 4 bytes (little-endian) of SHA-256 of `diagnosis_pair_<group A>-<group B>` |
| Bootstrap refits | `bootstrap_seed + replicate`, 2026238001 (ADNI), 2026238002 (NACC), through `SeedSequence([seed, 20260923, 62000000])` |
| Paired intervals of the composition experiment | 10,000 resamples of the 200 draws, `SeedSequence([20260923, 88000000, cohort])`, cohort 0 = ADNI, 1 = NACC |
| Gap and event-shift intervals | 10,000 paired resamples, `random.Random(20260924)` (ADNI), `random.Random(20260925)` (NACC) |
| Preflight (optional) | the fit, relabeling and bootstrap seeds above plus 1,000,000,000; B = 8 |
| ADNI 12-event population | seed 20260914: draw *i* `default_rng([20260914, 2, i])`, permutations `20260914 * 1000 + i`; bootstrap refits `resample_within_strata` with seed 20260914 |
| NACC amyloid PET population | seed 20260917 (draws, permutations and the 500 bootstrap refits, as for ADNI) |
| MRI-cognition populations | draws `SeedSequence([20260920, cohort, 9001, replicate])`, fits `2026092000 + 10000 * cohort + replicate`, cohort 1 = ADNI, 2 = NACC; paired intervals seeded from 20260920 and a hash of each endpoint |

Budgets: B = 599 permutations or relabelings; P = (1 + E)/600 with ties counted; APOE comparisons
Bonferroni-corrected for three pairs within each cohort and panel (rejection if E <= 9, smallest
adjusted P 0.005); single two-group comparisons of Fig. 4e rejected if E <= 29. 200 draws per
setting, 200 bootstrap refits per cohort (six-event panel) and 500 per earlier panel.

## Code labels

Several files keep the labels under which the analyses ran. The runtime files, `campaign.json`,
`batch_mapping.json`, `run_shared_panel.py`, `summarize_results.py` and `vendor/concord` are
checked by hash or are the code that ran, so their text is unchanged (apart from the two path edits
in the runtime described above).

| Label in code or output | Meaning in the paper |
|---|---|
| `likelihood_ebm` | likelihood EBM (kde_ebm Gaussian mixtures, greedy search, 5 starts x 500 proposals) |
| `separate`; `standard` (estimator name in concord-ebm 0.1.1) | separately fitted DEBM |
| `shared` | pooled-score DEBM |
| `concord`, `invariant_min`, `concord_min`, `min` | CONCORD; in the six-event panel with the fixed common proportions (`reference` in `campaign.json`) |
| `invariant_pooled`, `pooled` | CONCORD with pooled proportions (not reported) |
| `repaired` / `original` | continued search / unmodified pyebm search |
| `unrestricted`, suffix `_U`, `U-all` | unrestricted permutation |
| `diagnosis`, suffix `_D`, `D-all` | within-diagnosis permutation |
| `diagnosis_pair_<a>-<b>`, `D-pair`, `full599` | pairwise within-diagnosis permutation (599 relabelings) |
| `standard_U`, `Standard_U` | separately fitted DEBM with unrestricted permutation (Fig. 4e) |
| `concord_D`, `CONCORD_D` | CONCORD with within-diagnosis permutation (Fig. 4e) |
| `standardrep_D`, `Standard_D`, `concord_U`, `standard_original` | further combinations; not reported except `standard_original`, which is computed alongside `concord_min` in the bootstrap refits |
| `arm`, `arms` (keys) | one fitted model and permutation combination |
| `Standard`, `shared scoring`, `target` (text in frozen files) | separately fitted DEBM; pooled mixture; the fixed common proportions |
| `e2`, `e33`, `e4` | APOE e2 carriers without e4, e3/e3, e4 carriers without e2 |
| `CN`, `MCI`, `MCIc`, `AD` | cognitively normal, MCI, MCI with progression to dementia within 3 years (ADNI earlier populations), AD dementia |
| `n2` (mode or stage) | the draws of two groups: the controlled composition experiment (six-event panel) or the e3/e3 comparisons of Fig. 4e |
| `matched` / `different` (`imbalanced` in the draw file) | matched / different composition; in `different`, group A = CN-heavy (112/32/16), group B = AD-heavy (48/64/48) |
| `pA_e2comp`, `pB_e4comp`, `pA`, `pB` | groups A and B of a Fig. 4e draw |
| `apoe` | APOE comparisons |
| `bootstrap`, `a2s`, `stability` | bootstrap refits within genotype-by-diagnosis cells |
| `CSFcog6`, `omit_ptau` | six-event panel; five events without p-tau (Supplementary Table 5) |
| `ABETA`, `TAU`, `PTAU`, `MEM`, `EXF`, `LAN` | CSF Abeta42, total tau, p-tau; memory, executive function, language composites |
| `P12`, `P8`, `P5`; `ADNI_K12` | ADNI 12-, 8- and 5-event settings; ADNI 12-event panel |
| `MRIcog8`, `MRIcog4` | MRI-cognition 8- and 4-event settings |
| `NACC5_PET` | NACC amyloid PET panel (5 events) |
| `N2` (in `fig3_constructed_null_summary.csv`) | the e3/e3 comparisons of Fig. 4e |
| `fig3_*`, `fig4_*`, `r2_*`, `precision_*` (file names) | names of earlier source-data files, kept so the figure scripts can read them |

## Licence

GNU General Public License v3.0 or later (GPL-3.0-or-later); see `LICENSE` at the repository root.
