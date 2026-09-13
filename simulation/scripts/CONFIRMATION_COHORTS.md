# Frozen confirmation cohorts

`analyze_confirmation_cohorts.py` uses the standard library. It reads only the
explicit protocol, its frozen snapshot files, and each registered run's
`cell_i/chunk_j/{manifest.json,progress.json,rows.jsonl}`. It never opens large
permutation `fit_records.jsonl` files or generates scientific data.

```sh
python scripts/v2/analyze_confirmation_cohorts.py \
  --protocol runs/v2/confirmation_protocol.json \
  --results-root runs/v2/hpc_results \
  --output-prefix results/v2_confirmation/cohorts

python -m unittest discover -s scripts/v2 -p test_confirmation_cohorts.py -v
```

On a remote host, `--project-root` relocates protocol-relative snapshot paths.
Alternatively `--snapshots-root /path/to/snapshots` selects the registered run ID
under that directory; the frozen config path inside each snapshot remains the
protocol's exact path. Keep `hpc/v2/inventory_run.py` beside this script's repository
layout. `--require-ready` returns status 3 when strict final readiness is false.
Normal incomplete monitoring returns status 0, while an invalid frozen protocol,
snapshot or cohort specification returns status 2 and writes no new outputs.
Outputs are exactly `PREFIX.json` and `PREFIX.md`, outside the raw results root.

## Scope and counting

- K27 e2 and e4 partial-null endpoints use only their declared target-matched
  local scheme: the first500 from `confirm_pair_power_k27_a` and the next500 from
  that direction's topup. Each has R1000. E33 uses its own independent R1000 run.
  These are local false-rejection probabilities with reference bound alpha/3.
- Complete local-family endpoints use only `confirm_pair_complete_ref_a` (REF,
  R1000) and `confirm_pair_power_k27_a` (e2/e4, R500 each). The latter reports any
  true-alternative discovery and both true alternatives discovered, along with
  any-pair and any-true-null rejection. Topups add no family observations.
- The six global power points are the registered e2/e4 × K14/27/55 runs, R500
  each. The endpoint is the single `diagnosis` global-max test at alpha=.05.
  Pairwise Bonferroni decisions are never substituted for that test.
- Source/local observed duplicates are joined by explicit run/cell/engine/seed.
  Resolved design, data hash, biomarkers, and available successful observed
  orderings/taus/truth distances must match exactly. A source failure remains a
  failure. Missing truth on a failed source leaves the check unavailable; it does
  not replace that failure or invalidate an otherwise valid local success.
  A confirmed discrepancy marks both related endpoint records unknown.

Every endpoint keeps its full planned denominator and every unknown identity.
Known rejections k with u unknown among R have bounds `[k/R,(k+u)/R]`.
Rate, Monte Carlo standard error, and Wilson interval are present only when
u=0. Failed or unattempted permutations retain the planned B599 denominator;
independently derived integer completion bounds certify exact early decisions.
Quality/optimizer flags do not automatically exclude otherwise legal outcomes.

Three-valued family logic is used: any known rejection resolves an `any`
endpoint true; all known nonrejections resolve it false. For `both` discovery,
either known missed alternative resolves it false; both known rejections resolve
it true. Other combinations remain unknown.

Common-seed comparisons use the dataset seed as the Monte Carlo unit and report
second minus first. Unknown pair members retain worst/best difference bounds.
With complete pairs, the script provides empirical paired MCSE, a descriptive
normal interval when nondegenerate, and a conservative bounded-variable
Hoeffding interval. A zero empirical difference variance does not establish
equivalence. Intervals are pointwise, not simultaneous across the nine panel
contrasts. Same-seed different cells are related simulations, not duplicate data
and not a pooled independent R3000 cohort.

## Readiness and provenance

`collection_complete` means the inventory accounts for every planned terminal
record and runnable scheme obligation, including exhausted failed budgets.
`all_reported_endpoint_decisions_complete` separately records whether point
estimates can be formed for all reported endpoints. Strict
`ready_for_final_analysis` additionally requires the full-p obligations and
successful provenance/duplicate checks. An endpoint can be resolved before that
run is finished; it must not be presented as evidence of whole-run completion.

`ready_for_bounded_final_analysis` and `status=terminal_with_unknown` allow an
honest terminal report when all work is accounted for but failures leave endpoint
bounds. `status=incomplete` is used while planned work, identities, or integrity
checks remain unresolved. Failure records are never retried, removed, replaced,
or counted as known nonrejections by this analyzer.

The loader verifies protocol config hashes, the snapshot SHA256SUMS digest and
listed file contents, registered seeds/cells/engines/B/full-p settings, runtime
source signatures, pinned pyebm dependency hashes, numerical environments,
generated-design hashes, exact planned/realized inversions and pair masks, and
observed ordering distances. It verifies the fast flag through the frozen config
and runtime sources. It cannot independently recompute a data hash without the
raw generated dataset, or inspect each low-level permutation fit without the
large fit logs. Any inventory integrity error conservatively makes that entire
run unavailable for endpoint estimates until reviewed; no duplicate winner is
selected. All planned missing seeds and failures remain in JSON.

## Small reusable APIs

Sibling analysis wrappers can import this module without triggering analysis:

```python
run = frozen_run(protocol_entry, project_root, snapshots_root=None)
run = load_run(run, results_root)
row, reason = get_record(run, cell, engine, seed)
summary = summarize([(identity, decision_or_none, reason), ...])
contrast = paired_summary(first_seed_outcomes, second_seed_outcomes)
```

`run['records'][(cell,engine,seed)]` remains a list so duplicates are explicit.
The cohort CLI calls `load_run(..., compact=True)` to discard validated latent
arrays and optimizer histories from memory; the reusable loader defaults to
retaining the original row for sibling wrappers.
`get_record` returns a row only when it is unique, successful, and passes the
loader's integrity checks. Wrappers must preserve returned reasons within their
own full planned denominators and compare numerical environments across runs.
No scientific numerical packages are imported.
