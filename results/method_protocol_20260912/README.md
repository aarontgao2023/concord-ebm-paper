# Second confirmation protocol — saved state, 2026-09-13 00:1x UTC

Status at save time: second protocol 42,494 / 42,500 datasets decided (calibration 1,994 / 2,000; every other run
100 %); first protocol 11,223 / 12,000 (93.5 %), tails still running under the repaired driver. Nothing here is
final until the tally shows 100 %; the files below are the frozen inputs for tomorrow's rewrite of report v4,
the deck and the advisor summary.

## Files
| file | content |
|---|---|
| `mechanism_tables.txt` | frozen mechanism tables (R = 1,000 per cell): `method_mechanism_a` × 4 engines, `method_saebm_reference_a` × 3 engines, `method_saebm_a`; produced by `scripts/v2/dev/analyze_method_protocol.py` |
| `rates_exact_20260913T00.txt` | exact frozen-rule rejection rates with Wilson 95 % CIs for every permutation run of both protocols; `scripts/v2/dev/core_rates_exact.py` (core, stage; also shows the old convention side by side) and `rates_exact_all.py` |
| `rates_exact_20260912.txt` | same, one sync earlier |
| `tally_20260913T00.txt` | decision-level completeness per run; `scripts/v2/dev/tally_exact.py`, `tally_method_exact.py` |
| `SHA256SUMS_completed_rows.txt` | fingerprints of every `rows*.jsonl` of the five 100 %-complete runs (1,275 files) — verify unchanged before citing |

Rows themselves: `runs/v2/hpc_results/<run>/cell_*/chunk_*/rows.jsonl` (run_v2 runs) and
`runs/v2/hpc_results/method_saebm_a/rows_chunk*.jsonl` (SA-EBM dev script). Remote originals under
`/N/slate/tg11/ebmcal_v2/runs/<run>/`; snapshots under `/N/slate/tg11/ebmcal_v2/snapshots/method_*`.

## Convention note (applies to every number below and supersedes report v4 §3.2 / §5.4)
A rejection is *exactly determined* when G+E+(B−nperm) ≤ 9 (pair) or ≤ 29 (max); run_v2 stops there
(nperm 590–598). The earlier tally counted rejections only at nperm = 599, dropping those seeds from numerator
and denominator, which biased every rate down. All files here use exact decisions.

## Headline numbers at save time
First protocol (standard estimator), exact: IID 0.030–0.033 (U-all, D-all, both engines); **REF U-all 0.123
[0.104, 0.145]** (both engines; max 0.116); REF D-all 0.033 / 0.027; STAGE D-all 0.020, oracle 0.015;
complete-null D-pair 0.003 / 0.001 / 0.005; partial-null true-null pairs 0.006 / 0.014 / 0.002;
power K27 D-all family ε4 0.057, ε2 0.020; D-pair ε33–ε4 0.117 (these four still topping up).

Second protocol (frozen, seeds 53.xM): gap ε2–ε4 at n / 4n — standard 2.02 / 2.21, shared 0.82 / 1.16,
invariant_pooled 0.14 / 0.06, **invariant_min 0.12 / 0.06**, IID floor 0.07 / 0.04, SA-EBM 2.39 / 0.87.
Calibration REF: invariant_min D-all **0.038 [0.027, 0.051]**, U-all 0.098 [0.081, 0.119]; shared D-all 0.025,
U-all 0.772. Power invariant_min D-all: ε4 K14/27/55 = 0.190 / 0.574 / 0.968; ε2 = 0.072 / 0.214 / 0.618.
Boundary invariant_min D-all: PDF 0.074 [0.054, 0.100]; MISS_CSF 0.052 [0.036, 0.075]; MISS_ALL 0.222
[0.188, 0.260]; U-all 0.194 / 0.176 / 0.422.
