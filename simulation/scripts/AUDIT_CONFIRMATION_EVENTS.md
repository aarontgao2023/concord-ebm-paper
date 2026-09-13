# Streaming raw-event audit for frozen confirmation chunks

`audit_confirmation_events.py` audits one explicitly selected `cell_i/chunk_j`.
It uses only the Python standard library and the existing cohort/inventory
validation helpers. It does not import the scientific runner, generate a
dataset, fit a model, invoke HPC commands, or call the runner's log-repair loader.

Run from the project root, with an output prefix outside raw data and frozen
snapshots:

```sh
python scripts/v2/audit_confirmation_events.py \
  --protocol runs/v2/confirmation_protocol.json \
  --run-id confirm_core_a \
  --chunk runs/v2/hpc_results/confirm_core_a/cell_0/chunk_0 \
  --output-prefix results/v2_confirmation/raw_audit/core_cell0_chunk0 \
  --require-certified
```

For a relocated copy, supply `--project-root` and/or `--snapshots-root`. The
snapshots root contains a subdirectory named after each registered run. The
protocol still selects exactly one registered run and supplies its frozen
config/snapshot hashes; this command never discovers extra runs by glob.
Copy the auditor, `scripts/v2/analyze_confirmation_cohorts.py` and
`hpc/v2/inventory_run.py` with their relative project layout when using a separate
audit installation. No NumPy, SciPy or pyebm import is needed. The snapshot's
runtime files are read and hashed, never executed.

The command writes JSON and Markdown. Exit 0 means the audit ran, not necessarily
that the evidence was certified. `--require-certified` returns 3 for an invalid
or deferred result. Invalid protocol, frozen signature, manifest allocation or
unreadable filesystem errors fail closed with exit 2. Existing raw, protocol,
actual registered snapshot, config and audit source paths cannot be output
targets. Outputs use exclusive temporary files and atomic replacement.

## What is independently recomputed

- The manifest must match the registered config, cell/chunk directory, signed
  runtime sources, exact strided seed allocation and pinned numerical dependency.
- Every raw `(seed, engine, scheme, perm_id)` is checked against that allocation.
  Observed identities use `scheme=null, perm_id=null`; permutation IDs are zero
  based and must lie in `[0, B-1]`. Duplicate identities prevent certification;
  a successful duplicate never replaces an earlier failure.
- Each successful record must contain three legal integer permutations. Pairwise
  Kendall inversion counts and saved normalized distances are checked exactly.
  Four integer histograms accumulate the three pair distances and their maximum.
  This also works when raw events precede their observed record or IDs are not
  a consecutive prefix.
- Raw observed truth, data hash, environment, observed result and diagnostics are
  compared to the derived row. A failed observation can legitimately lack truth;
  its planned identity and failure are retained, with the missing-truth audit
  detail. Successful observed truth errors prevent certification.
- Raw histograms independently reconstruct `gt`, `eq`, `maxgt`, `maxeq`,
  successful/attempted/failed counts, complete flags and exact decisions under
  both `strict` and `le`. They are compared field by field with each derived
  scheme. Bounds always retain the frozen `B+1` denominator and failed unknown
  outcomes. A local scheme uses only its registered tested pair for stopping;
  it does not accidentally require an untested maximum statistic.
- Progress rows, expected rows, statuses, statistical fits and computational
  jobs are checked. Paired engines count one shared computational job while
  retaining two statistical observations. Shared IDs and any recorded optimizer
  mode must agree with their engine/dataset/scheme/ID labels.

## Certification and continuation

`certified_done` means stable input bytes, consistent raw/derived/progress data,
and no further prespecified work for that chunk. `certified_stopped_checkpoint`
means a stable consistent checkpoint, which can still require continuation.
Neither label establishes scientific validity, removes failed outcomes, or
asserts that the whole run is complete.

`deferred` means a coherent settled copy is needed. This includes a live runner
lock, active progress flags, changed inputs, missing files, an uncommitted final
line, differing raw/progress counts, or monotone count differences between raw
and derived checkpoints. Such differences can result from file-transfer timing.
Any accompanying findings are provisional and must not be called corruption.
Same-count conflicting `gt/eq` or observed content in a stable settled snapshot
remains an integrity finding. A writer lock is probed read-only and released
immediately; certification applies to the hashed bytes read, not a promise that
no process can start after the audit.

The JSON contains every planned seed/engine/scheme with inclusive compressed
ID ranges. `missing_permutation_id_ranges` lists every unrecorded planned ID.
`required_resume_id_ranges` is a subset: exact non-full-p early stopping can
make all remaining IDs unnecessary; a prespecified full-p seed still needs its
unattempted IDs. Failed or invalid recorded IDs are listed separately and are
never automatically proposed for retry. Failed observed fits likewise remain
terminal failures. `continuation_manifest_available=false` makes all such lists
provisional; no scheduler should act on a deferred or invalid audit.

The protocol/config/runtime and helper hashes are tied to the start of the
audit and checked again before certification. Each input has a content hash,
size/mtime/inode information, and raw/row logs include committed-line/tail
counts. The auditor cannot prove the implementation that actually ran from
labels alone, nor regenerate permutation labels without the scientific data;
it certifies internal evidence consistency with the frozen contract.

## Memory and tests

Raw and derived JSONL are read one line at a time. Memory keeps one PID byte per
planned permutation, four small inversion histograms per seed/engine/scheme,
and observed hashes/summaries. It never retains the full permutation records.
Finding examples are capped at 200 per class, with uncapped aggregate counts.
JSON outputs contain compact continuation evidence, not fitted trajectories or
scientific effect estimates.

```sh
python -m unittest discover -s scripts/v2 -p 'test_audit_confirmation_events*.py' -v
```

The fixtures use invented three-event records and temporary signed snapshots.
They test hand-calculated pair/max ties, the B599 equality boundary and failed
unknown completion, nonprefix IDs, failure preservation, paired accounting,
exact stopping versus full-p work, snapshot/protocol stability, sync deferral,
malformed records and output protections. They run no simulation or fitting.
