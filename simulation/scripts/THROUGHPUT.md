> Execution update 2026-09-08: the4h/11800s table below describes older snapshots. All ten current confirmation snapshots use2h/4000s with16workers; do not apply the old timing fields to them. See the frozen protocol and confirmation_freeze_audit for current drain margins.

# Read-only throughput audit

Run `analyze_throughput_v2.py` on the machine holding the large fit logs. It uses
only Python's standard library and `hpc/v2/inventory_run.py`. For a separate HPC
audit-tools directory, retain those two relative paths. Do not modify frozen
experiment snapshots to install this analysis utility.

```bash
python scripts/v2/analyze_throughput_v2.py \
  --config runs/v2/snapshots/dev_calibration_a/configs/v2/dev_calibration_a.json \
  --root runs/v2/hpc_results/dev_calibration_a \
  --output-prefix results/v2_dev_calibration_a/throughput
```

On HPC, `--root` points to the actual run directory containing
`cell_i/chunk_j/{manifest,progress}.json`, `rows.jsonl`, and `fit_records.jsonl`.
Keep any `replay_audit.jsonl` files as well. Fit logs are streamed; large truth and
optimizer payloads are discarded after each line. Only compact identity,
status, duration and shared-job metadata stay in memory. Output is a JSON audit
and a Markdown table; stdout is one short status line.

The planned denominator always includes every seed and engine. “Recorded started”
requires a durable fit record, so it excludes queued and currently executing fits.
“Finished” requires the planned primary decisions and any forced full-p obligation;
exhausted failed budgets have a separate terminal count. Mean/p90/p95 durations
include failed recorded attempts, including timeouts.

An unfinished run is explicitly right censored. Savings and overall forecasting
remain null even if some easy datasets have finished. Once the complete run and
log/row counts pass validation, this optional command produces an ideal planning
scenario for the same cells, engines, schemes and B:

```bash
python scripts/v2/analyze_throughput_v2.py \
  --config FROZEN_CONFIG.json --root RUN_DIRECTORY \
  --output-prefix OUTPUT_PREFIX \
  --forecast-datasets 100 --forecast-workers 256
```

The forecast excludes queueing, startup, I/O and idle allocations. It is withheld
if paired replays or conflicting shared durations prevent reconstructing cost.
Shared fit seconds count once per computational job; per-engine duration columns
must never be summed as CPU cost. Reported fit/job avoidance is a reduction in
logical planned work, not a measured CPU-time saving.

Optional allocation measurements may be supplied with `--sacct FILE`. Use a
pipe-delimited file with headers `JobIDRaw` (or `JobID`), `State`, `CPUTimeRAW`, and
`TotalCPU`, restricted to the relevant jobs. The script excludes `.batch`,
`.extern` and other step rows, and rejects duplicate allocation IDs. Utilization
comes only from measured `TotalCPU / CPUTimeRAW`; fit wall seconds are not used to
guess it.

For 16 workers and the runner's maximum 32 submitted tasks, the final wave may
need two task timeouts to drain. The 11,800-second runner budget and Slurm's
signal 1,500 seconds before a four-hour limit are alternative stop triggers:

| Fit timeout | Pairing | Nominal drain | Margin from signal alone | Margin from runner budget, zero startup |
|---:|---|---:|---:|---:|
| 300 s | independent | 600 s | 900 s | 2,000 s |
| 300 s | paired, 600 s/job | 1,200 s | 300 s | 1,400 s |
| 600 s | paired, 1,200 s/job | 2,400 s | −900 s | 200 s |

These margins exclude serialization, fsync, process scheduling and timeout
delivery. The 600-second paired mechanism setting has insufficient signal-only
margin and only 200 seconds before startup/I/O overhead under the runner budget.
It needs a shorter chunk or an earlier runner stop budget. This audit does not
change the runner or Slurm scripts. Timeout results remain unknown outcomes in
the original planned permutation denominator, never a success-only denominator.

Tests: `python -m unittest discover -s scripts/v2 -p test_throughput_v2.py -v`.
