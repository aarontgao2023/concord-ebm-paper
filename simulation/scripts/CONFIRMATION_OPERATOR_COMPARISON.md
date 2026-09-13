# Supplementary comparison of pair Bonferroni operators

`analyze_confirmation_operator_comparison.py` reads the four runs named in the explicit post-freeze contract. It compares diagnosis-stratified full-three-group **pair Bonferroni** decisions with the matched-pair diagnosis family. It never substitutes the source global-max rejection, mixes topups, or pools E2/E4 common-random-number cells as independent samples.

Run from the project root with a fresh output prefix:

```sh
python3 scripts/v2/analyze_confirmation_operator_comparison.py \
  --protocol runs/v2/confirmation_protocol.json \
  --contract docs/review_v2/operator_comparison_contract_20260908.json \
  --results-root runs/v2/hpc_results \
  --output-prefix results/v2_confirmation/NEW_CAPTURE/operator_comparison
```

For relocated HPC snapshots, supply `--snapshots-root /path/to/snapshots`; each registered run must be a child of that directory. This is a read-only analysis command, not a submission command. JSON and Markdown outputs must not exist and must be outside actual source/snapshot/raw input paths. Input and analyzer hashes are captured and checked again before releasing the report; a changing capture raises an error without certifying estimates. The script does not read `fit_records.jsonl` or modify runner/analyzer snapshots.

The seven contracted comparisons retain R1,000 for REF and R500 separately for E2/E4 K27. Each arm's estimates require its complete planned endpoint; paired estimates require complete paired endpoints. Otherwise the report contains null estimates/intervals, completion bounds, and every unknown identity. Status and collection completion cover all four runs' work obligations, including source-core IID/original control rows unused by these comparisons. An endpoint may resolve earlier than the overall collection, which is why these flags are separate.

The independent fullgroup extractor checks recorded scheme definitions, successful counts and both saved exact-decision rules against fixed-B completion bounds. Inclusive ties enter `gt+eq`; `nperm` counts successful fits, while failed and unattempted permutations stay unknown within B599. Pair threshold comparison uses integer arithmetic at `.05/3`. A failure can leave an exact decision resolved only if every allowed completion agrees. Local extraction uses the existing cohort analyzer's matched-pair validation. Failed observed fits and duplicates are never replaced by the other run's copy.

JSON includes complete-cohort Wilson intervals/MCSE and paired MCSE, a nondegenerate normal interval when applicable, and conservative pointwise Hoeffding intervals. These seven supplementary intervals are not simultaneous and cannot establish equivalence. This postprocessing is not a complete reproduction of the 2021 procedure or evidence of universal ordering-only strong FWER.

Pure fixture checks:

```sh
PYTHONPATH=scripts/v2 python3 -m unittest test_confirmation_operator_comparison
```

The fixtures read frozen design metadata but use inert count records and fake loader returns. They generate no datasets and fit no model. Cases cover the B599 tie boundary, failed-versus-successful counts, source max versus pair decisions, exact early stopping, malformed declarations, partial-null target masks, all planned denominators, duplicate/mismatched joins, preserved observed failures, legal quality flags, complete paired direction and nondegenerate uncertainty, changing input capture, and output protection.
