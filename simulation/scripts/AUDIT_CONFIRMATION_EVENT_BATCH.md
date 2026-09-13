# Bounded serial batches of raw-event audits

This standard-library wrapper calls the already reviewed
`audit_confirmation_events.py` serially for selected settled chunks. It performs
no scheduling, SSH, fitting, data generation or raw-file repair. It maintains no
cache and makes no claim about whole-run completion.

```sh
python scripts/v2/audit_confirmation_event_batch.py \
  --protocol runs/v2/confirmation_protocol.json \
  --run-id confirm_core_a \
  --run-id confirm_stage_a \
  --results-root runs/v2/hpc_results \
  --output-dir results/v2_confirmation/raw_batch_unique_label \
  --offset 0 \
  --max-chunks 10 \
  --walltime-budget-s 1800
```

Each `--run-id` must be uniquely registered in the explicit protocol. Repeated
run IDs and unknown IDs are rejected. `--project-root` and `--snapshots-root`
have the same relocation meaning as the underlying auditor. A separate audit
installation needs this wrapper, the underlying auditor, cohort loader and
inventory helper in their existing relative directory layout.

The output directory must not exist, even if empty. Use a new directory for
each batch. Outputs cannot be placed inside the results root or registered
frozen snapshots. Files are published exclusively and never overwritten;
an interrupted batch preserves its already published evidence. A failed
publication or fatal preparation error exits 2 and never converts unprocessed
work into success. A finished batch exits 0 only when every selected chunk has
a certified checkpoint and sources stayed stable; otherwise it exits 3.

## Fixed selection

1. Validate requested registrations and frozen sources. Enumerate only existing
   `run/cell_*/chunk_*/manifest.json` directories inside each explicit run.
2. Record all discovered directories, including errors, missing/unreadable
   metadata and live chunks. A candidate needs a matching signed allocation,
   boolean progress with `done=true` or `stopped=true`, stable selection metadata
   and no held writer lock.
3. Sort settled candidates by run ID, numeric cell and numeric chunk. Apply
   zero-based `--offset`, then `--max-chunks`. Publish the complete immutable
   `selection.json`, including every discovered directory, its selection
   reason, all candidate paths and the selected paths, before any raw scan.
4. Process only that fixed selection, serially. Recheck eligibility before each
   scan so a restarted writer is recorded as live instead of being scanned.
   The underlying auditor independently checks stability again during its scan.

The original manifest hash also binds each selected allocation. It is checked
before the scan, against the auditor's input hash, and at batch finalization.
A manifest replacement or even a valid re-sharding to different planned seeds
is deferred; it cannot inherit certification under the old fixed selection.

Offsets refer to the current batch's fixed settled list. A later invocation can
have a different candidate list because more chunks settled or resumed. Consult
the saved `selection.json` before choosing a subsequent offset; offsets are not
a persistent cursor or a guarantee against repeating an earlier chunk. Chunks
that appear after selection are outside this batch's scope.

`--walltime-budget-s` measures from batch preparation and is checked before
starting each selected chunk. It is a **between-chunk budget**, not a hard
deadline: an ongoing scan can run beyond the remaining time. Set a suitable
margin in the surrounding job's wall limit. The wrapper never interrupts a
scientific runner or kills a raw scan to satisfy this budget.

## Evidence and accounting

Each started chunk produces `RUN__cell_I__chunk_J.json` and `.md`. A per-chunk
exception produces its own error evidence and does not cause automatic retry.
The JSON records the selection hash, batch source hashes and elapsed audit
time. Live or deferred selection/recheck details remain in selection/summary
without opening large raw logs unnecessarily.

`summary.json` and `.md` record all discovered chunks and mutually exclusive
final categories: certified done, certified stopped checkpoint, deferred,
invalid, error, live, outside the offset/limit slice, or selected but unprocessed
because of budget or source changes. Unselected and unprocessed chunks never
count as certified. No existing manifests is an unassessed scope, not success.
The planned row count for each selected run is shown only as context; absent
manifests and uncreated chunks are not audited or treated as completed.

Protocol, config, runtime, wrapper and helper hashes are fixed before selection
and checked before/after scans and at finalization. If any change is observed,
further selected chunks remain unprocessed, the current result is deferred, and
batch-certified counts become zero. Earlier per-chunk files are retained
unchanged as evidence of the bytes audited then. **The final summary governs
batch source stability.** If a process is externally interrupted before the
summary appears, there is no completed batch certificate.

## Toy tests

```sh
python -m unittest discover -s scripts/v2 -p 'test_audit_confirmation_event_batch*.py' -v
```

Tests use temporary signed snapshots and invented logs only. They cover fixed
offset/limit selection, immutable output, live/missing/invalid metadata,
restarted writers, between-chunk budget limits, source changes, independent
per-chunk error artifacts and explicit registered-run scope. They do not invoke
the simulator, fitter, SSH or scheduler.
