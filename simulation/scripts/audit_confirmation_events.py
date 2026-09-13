#!/usr/bin/env python3
"""Read-only streaming audit of one explicitly selected frozen confirmation chunk.

Never imports the scientific runner, repairs a log, generates data or fits a
model. Memory stores PID bytearrays, distance histograms and observed hashes,
not full permutation records. Certification refers to the input bytes audited.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile

import analyze_confirmation_cohorts as C

I = C.I
STATES = {"ok", "error", "invalid"}
COUNT_FIELDS = ("requested_nperm", "nperm", "attempted", "failed", "complete", "gt", "eq", "maxgt", "maxeq", "exact_decision", "definition")


def ranges(values):
    """Sorted inclusive integer ranges; empty is [] and a singleton is [i,i]."""
    result = []
    for value in sorted(values):
        if result and value == result[-1][1] + 1:
            result[-1][1] = value
        else:
            result.append([value, value])
    return result


def file_state(path):
    if not path.exists():
        return None
    stat = path.stat()
    return {"inode": stat.st_ino, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def lock_probe(path):
    """Probe without creating/changing the lock or holding it during the scan."""
    if not path.exists():
        return "absent"
    with path.open("rb") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_SH | fcntl.LOCK_NB)
        except BlockingIOError:
            return "held_by_writer"
        fcntl.flock(handle, fcntl.LOCK_UN)
    return "not_held_at_probe"


class Audit:
    def __init__(self):
        self.findings = []
        self.deferred = []
        self.finding_counts = Counter()
        self.deferred_counts = Counter()
        self.inputs = {}

    def issue(self, code, **context):
        self.finding_counts[code] += 1
        if len(self.findings) < 200:
            self.findings.append({"code": code, **context})

    def defer(self, code, **context):
        self.deferred_counts[code] += 1
        if len(self.deferred) < 200:
            self.deferred.append({"code": code, **context})

    def json_file(self, path):
        before = file_state(path)
        if before is None:
            self.inputs[path.name] = {"path": str(path), "exists": False}
            self.defer("missing_input", path=str(path))
            return None
        raw = path.read_bytes()
        after = file_state(path)
        self.inputs[path.name] = {"path": str(path), "exists": True, "before": before,
                                 "after": after, "sha256": hashlib.sha256(raw).hexdigest()}
        if before != after:
            self.defer("input_changed_during_read", path=str(path))
        try:
            result = json.loads(raw)
            C.require(isinstance(result, dict), "Expected JSON object")
            return result
        except (ValueError, TypeError) as exc:
            self.defer("metadata_unreadable_pending_consistent_copy", path=str(path), message=str(exc))
            return None

    def lines(self, path):
        before = file_state(path)
        if before is None:
            self.inputs[path.name] = {"path": str(path), "exists": False}
            self.defer("missing_input", path=str(path))
            return
        digest, physical, committed, tail_bytes = hashlib.sha256(), 0, 0, 0
        with path.open("rb") as handle:
            for number, raw in enumerate(handle, 1):
                physical += 1
                digest.update(raw)
                if not raw.endswith(b"\n"):
                    tail_bytes = len(raw)
                    self.defer("uncommitted_final_line_no_repair", path=str(path), line=number, bytes=len(raw))
                    continue
                committed += 1
                if not raw.strip():
                    self.issue("blank_event_or_row_line", path=str(path), line=number)
                    continue
                try:
                    value = json.loads(raw)
                    C.require(isinstance(value, dict), "Expected JSON object")
                except (ValueError, TypeError) as exc:
                    self.issue("invalid_committed_json", path=str(path), line=number, message=str(exc))
                    continue
                yield number, value
        after = file_state(path)
        self.inputs[path.name] = {"path": str(path), "exists": True, "before": before, "after": after,
                                 "sha256": digest.hexdigest(), "physical_lines": physical,
                                 "committed_newline_lines": committed, "uncommitted_tail_bytes": tail_bytes,
                                 "final_newline": tail_bytes == 0}
        if before != after:
            self.defer("input_changed_during_read", path=str(path))


def verify_manifest(run, directory, manifest):
    cfg = run["cfg"]
    cell_match = re.fullmatch(r"cell_(\d+)", directory.parent.name)
    chunk_match = re.fullmatch(r"chunk_(\d+)", directory.name)
    C.require(cell_match is not None and chunk_match is not None, "Expected explicit cell_i/chunk_j path")
    cell_index, chunk = int(cell_match[1]), int(chunk_match[1])
    C.require(0 <= cell_index < len(cfg["cells"]), "Cell index outside frozen plan")
    n = I.count(manifest.get("nchunks"), "nchunks")
    I.count(manifest.get("cell_index"), "manifest.cell_index")
    I.count(manifest.get("chunk"), "manifest.chunk")
    C.require(n > 0 and 0 <= chunk < n, "Invalid chunk layout")
    C.require(manifest.get("config") == cfg and manifest.get("cell_index") == cell_index and manifest.get("chunk") == chunk,
              "Manifest differs from frozen config or directory")
    signed = {k: manifest[k] for k in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
    C.require(manifest.get("signature") == I.digest(signed), "Manifest signature mismatch")
    C.require(manifest.get("source_sha256") == run["sources"], "Runtime source hashes differ from frozen snapshot")
    seeds = sorted(run["seeds"])[chunk::n]
    C.require(manifest.get("seeds") == seeds and bool(seeds), "Manifest strided seed allocation mismatch")
    environment = manifest.get("environment", {})
    dep = run["dependency"]
    C.require(isinstance(environment.get("versions"), dict) and set(I.NUMERICAL_PACKAGES) <= set(environment["versions"]), "Missing numerical environment")
    C.require(all(isinstance(v, str) and v for v in environment["versions"].values()), "Invalid numerical version")
    C.require(environment["versions"]["pyebm"] == dep["PYEBM_VERSION"] and environment.get("source_sha256") == dep["SOURCE_SHA256"]
              and environment.get("wheel_sha256") == dep["PYEBM_WHEEL_SHA256"], "Unpinned numerical dependency")
    return cell_index, chunk, seeds, environment


def pair_counts(orderings, events):
    C.require(isinstance(orderings, list) and len(orderings) == 3, "Expected three group orderings")
    for order in orderings:
        C.require(isinstance(order, list) and len(order) == events and all(type(v) is int for v in order)
                  and sorted(order) == list(range(events)), "Ordering is not a legal biomarker permutation")
    return [C.inversions(orderings[a], orderings[b]) for a, b in C.PAIRS]


def validated_taus(record, events):
    counts = pair_counts(record.get("orderings"), events)
    values = record.get("taus")
    expected = [k / math.comb(events, 2) for k in counts]
    C.require(isinstance(values, list) and len(values) == 3 and all(type(v) in (int, float) and math.isfinite(v) for v in values), "Invalid saved taus")
    C.require(values == expected, "Saved taus disagree with legal orderings")
    return counts


def shared_id(cfg, cell, seed, scheme, pid):
    payload = [cfg["run_id"], cell, seed, scheme, pid]
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def observed_digest(record, cfg, full_selected):
    truth = record.get("truth", {})
    manifest = truth.get("manifest") if isinstance(truth, dict) else None
    data_hash = manifest.get("data_sha256") if isinstance(manifest, dict) else None
    return {"status": record["status"], "truth_sha256": C.canonical_sha(record.get("truth", {})),
            "data_sha256": data_hash if I.is_hash(data_hash) else None,
            "observed_sha256": C.canonical_sha({k: record.get(k) for k in ("orderings", "taus", "distance_to_truth", "fit_seconds")}),
            "diagnostics_sha256": C.canonical_sha(record["diagnostics"] | {
                "sequential_stopping": bool(cfg.get("schemes")), "full_p_selected": full_selected})}


def new_state(budget, events):
    support = math.comb(events, 2) + 1
    return {"ids": bytearray(budget), "histograms": [[0] * support for _ in range(4)],
            "attempted": 0, "successful": 0, "failed": 0, "invalid": 0}


def derive(state, observed, definition, cfg):
    if observed is None or observed["status"] != "ok" or "pair_counts" not in observed:
        return None
    reference = observed["pair_counts"]
    hist = state["histograms"]
    gt = [sum(hist[i][reference[i] + 1:]) for i in range(3)]
    eq = [hist[i][reference[i]] for i in range(3)]
    maxgt, maxeq = sum(hist[3][max(reference) + 1:]), hist[3][max(reference)]
    n, budget, alpha = state["successful"], cfg.get("bperm", 0), cfg.get("alpha", .05)
    exact = {"alpha": alpha}
    for rule in ("strict", "le"):
        exact[rule] = {"pair_reject": [I.bounds(g + e, n, budget, alpha, 3, rule) for g, e in zip(gt, eq)],
                       "max_reject": I.bounds(maxgt + maxeq, n, budget, alpha, 1, rule)}
    return {"requested_nperm": budget, "nperm": n, "attempted": state["attempted"], "failed": state["failed"],
            "complete": n == budget, "gt": gt, "eq": eq, "maxgt": maxgt, "maxeq": maxeq,
            "exact_decision": exact, "definition": definition}


def work_item(seed, engine, scheme, state, observed, definition, cfg, selected):
    budget = cfg.get("bperm", 0)
    absent = [pid for pid, value in enumerate(state["ids"]) if value == 0]
    failed = [pid for pid, value in enumerate(state["ids"]) if value == 2]
    invalid = [pid for pid, value in enumerate(state["ids"]) if value == 3]
    value = derive(state, observed, definition, cfg)
    resolved = False
    if value is not None:
        primary = value["exact_decision"][cfg.get("rule", "strict")]
        indexes = [definition["tested_pair_index"]] if "tested_pair_index" in definition else range(3)
        resolved = all(primary["pair_reject"][i] is not None for i in indexes)
        if definition.get("require_max", False):
            resolved = resolved and primary["max_reject"] is not None
    if observed is None:
        action = "await_observed_fit"
    elif observed["status"] != "ok":
        action = "observed_failure_preserved_no_retry"
    elif value is None or invalid:
        action = "inspect_invalid_raw_evidence_no_automatic_retry"
    elif state["successful"] == budget:
        action = "full_budget_complete"
    elif state["attempted"] == budget:
        action = "budget_exhausted_failures_preserved"
    elif resolved and not selected:
        action = "exact_early_stop_no_more_permutations_required"
    else:
        action = "resume_prespecified_full_p" if selected else "resume_unresolved_decision"
    resumable = action.startswith("resume_")
    return {"seed": seed, "engine": engine, "scheme": scheme, "planned_B": budget, "full_p_selected": selected,
            "observed_status": observed["status"] if observed else "missing", "action": action,
            "attempted": state["attempted"], "successful": state["successful"], "failed": state["failed"],
            "invalid_success_records": state["invalid"], "missing_permutation_count": len(absent),
            "missing_permutation_id_ranges": ranges(absent), "failed_permutation_id_ranges": ranges(failed),
            "invalid_record_id_ranges": ranges(invalid), "required_resume_id_ranges": ranges(absent) if resumable else [],
            "primary_decisions_resolved": resolved, "full_p_available": state["successful"] == budget,
            "work_terminal": action in {"observed_failure_preserved_no_retry", "full_budget_complete",
                "budget_exhausted_failures_preserved", "exact_early_stop_no_more_permutations_required"},
            "raw_recomputed": value}


def computational_jobs(observed, states, seeds, cfg):
    count = 0
    standard = {"original", "repaired"}
    paired = cfg.get("paired_standard", False) and standard <= set(cfg["engines"])
    for seed in seeds:
        for engine in cfg["engines"]:
            if paired and engine in standard:
                continue
            count += int((seed, engine) in observed)
            count += sum(states[seed, engine, name]["attempted"] for name in cfg.get("schemes", {}))
        if paired:
            count += int(any((seed, engine) in observed for engine in standard))
            for scheme in cfg.get("schemes", {}):
                a, b = (states[seed, engine, scheme]["ids"] for engine in ("original", "repaired"))
                count += sum(bool(x or y) for x, y in zip(a, b))
    return count


def audit_chunk(run, directory):
    directory = Path(directory).resolve()
    audit = Audit()
    cfg = run["cfg"]
    inputs = [directory / name for name in ("manifest.json", "progress.json", "rows.jsonl", "fit_records.jsonl")]
    initial = {p.name: file_state(p) for p in inputs}
    lock_before = lock_probe(directory / "runner.lock")
    manifest = audit.json_file(directory / "manifest.json")
    progress = audit.json_file(directory / "progress.json")
    if manifest is None:
        return {"schema_version": "confirmation_event_audit_v1", "run_id": cfg["run_id"], "chunk_path": str(directory),
                "status": "deferred", "certified_done_chunk": False, "checkpoint_certified": False,
                "findings": audit.findings, "deferred_reasons": audit.deferred, "inputs": audit.inputs,
                "continuation_manifest_available": False}
    cell_index, chunk, seeds, environment = verify_manifest(run, directory, manifest)
    cell = run["names"][cell_index]
    design = run["designs"][cell]
    events = len(design["biomarker_names"])
    budget, schemes, engines = cfg.get("bperm", 0), cfg.get("schemes", {}), cfg["engines"]
    states = {(seed, engine, scheme): new_state(budget, events) for seed in seeds for engine in engines for scheme in schemes}
    observed, statuses, duplicates = {}, Counter(), 0
    raw_unique, valid_success = 0, 0
    paired = cfg.get("paired_standard", False)
    fast = C.literal_assignments(run["snapshot"] / "scripts/v2/fast_likelihood_v2.py")["FAST_LIKELIHOOD_VERSION"]
    for number, record in audit.lines(directory / "fit_records.jsonl"):
        identity = [record.get(k) for k in ("seed", "engine", "scheme", "perm_id")]
        try:
            seed, engine, scheme, pid = identity
            C.require(type(seed) is int and seed in seeds and engine in engines, "Raw identity outside planned seed/engine allocation")
            C.require((scheme is None and pid is None) or (scheme in schemes and type(pid) is int and 0 <= pid < budget), "Raw scheme/permutation ID outside plan")
            C.require(record.get("kind") == ("observed" if scheme is None else "permutation"), "Raw kind/identity mismatch")
            C.require(record.get("status") in STATES, "Unknown raw fit status")
            C.require(isinstance(record.get("diagnostics"), dict), "Missing raw diagnostics")
            expected_mode = engine.rsplit("_", 1)[-1]
            mode = record["diagnostics"].get("mode")
            if mode is not None and mode != expected_mode:
                audit.issue("raw_diagnostic_mode_engine_mismatch", line=number, identity=identity)
            consensus = record["diagnostics"].get("consensus", [])
            if isinstance(consensus, list) and any(isinstance(v, dict) and v.get("mode", expected_mode) != expected_mode for v in consensus):
                audit.issue("raw_consensus_mode_engine_mismatch", line=number, identity=identity)
            if scheme is None:
                duplicate = (seed, engine) in observed
            else:
                duplicate = states[seed, engine, scheme]["ids"][pid] != 0
            if duplicate:
                duplicates += 1
                audit.issue("duplicate_raw_identity_no_certification", line=number, identity=identity)
                continue
            raw_unique += 1
            statuses[record["status"]] += 1
            expected_implementation = fast if cfg.get("fast_likelihood") else "upstream-frozen-normal"
            if engine not in ("original", "repaired"):
                expected_implementation = "oracle_not_applicable"
            if record.get("likelihood_implementation") != expected_implementation:
                audit.issue("raw_likelihood_implementation_mismatch", line=number, identity=identity)
            shared = record.get("shared_fit_id")
            if paired and engine in ("original", "repaired"):
                if shared != shared_id(cfg, cfg["cells"][cell_index], seed, scheme, pid) or record["diagnostics"].get("shared_fit_id") != shared:
                    audit.issue("invalid_paired_shared_fit_identity", line=number, identity=identity)
            elif shared is not None:
                audit.issue("unexpected_shared_fit_identity", line=number, identity=identity)
            counts = None
            if record["status"] == "ok":
                try:
                    counts = validated_taus(record, events)
                    valid_success += 1
                except (ValueError, TypeError, KeyError) as exc:
                    audit.issue("invalid_successful_ordering_or_tau", line=number, identity=identity, message=str(exc))
            if scheme is None:
                item = observed_digest(record, cfg, seed in run["full_p"])
                item["raw_line"] = number
                if counts is not None:
                    item["pair_counts"] = counts
                observed[seed, engine] = item
                row_adapter = record | {"observed": {k: record.get(k) for k in ("orderings", "taus", "distance_to_truth", "fit_seconds")}}
                try:
                    C.validate_truth(row_adapter, design, run["sources"]["design_v2.py"], environment)
                except (ValueError, KeyError, TypeError) as exc:
                    item["truth_audit_error"] = str(exc)
                    if record["status"] == "ok":
                        audit.issue("invalid_observed_truth_or_provenance", line=number, identity=identity, message=str(exc))
                    # Failed generation can legitimately have no truth: failure
                    # remains in the planned denominator and is never replaced.
            else:
                state = states[seed, engine, scheme]
                state["attempted"] += 1
                if record["status"] != "ok":
                    state["ids"][pid] = 2
                    state["failed"] += 1
                elif counts is None:
                    state["ids"][pid] = 3
                    state["invalid"] += 1
                else:
                    state["ids"][pid] = 1
                    state["successful"] += 1
                    for i, value in enumerate(counts + [max(counts)]):
                        state["histograms"][i][value] += 1
        except (ValueError, TypeError, KeyError) as exc:
            audit.issue("invalid_raw_record", line=number, identity=identity, message=str(exc))
    continuation = [work_item(seed, engine, scheme, states[seed, engine, scheme], observed.get((seed, engine)), definition, cfg, seed in run["full_p"])
                    for seed in seeds for engine in engines for scheme, definition in schemes.items()]
    for seed in seeds:
        hashes = {observed[seed, engine]["data_sha256"] for engine in engines
                  if (seed, engine) in observed and observed[seed, engine]["data_sha256"] is not None}
        if len(hashes) > 1:
            audit.issue("cross_engine_generated_data_hash_mismatch", seed=seed)
    for item in continuation:
        if item["attempted"] and item["observed_status"] != "ok":
            audit.issue("permutation_events_without_successful_observed", seed=item["seed"], engine=item["engine"], scheme=item["scheme"])
    expected_schemes = {(v["seed"], v["engine"], v["scheme"]): v["raw_recomputed"] for v in continuation}
    row_keys, matched_observed = set(), 0
    for number, row in audit.lines(directory / "rows.jsonl"):
        try:
            seed, engine = row.get("seed"), row.get("engine")
            C.require(type(seed) is int and seed in seeds and engine in engines and row.get("run_id") == cfg["run_id"]
                      and row.get("phase") == cfg["phase"] and row.get("cell") == cell, "Derived row identity outside plan")
            key = (seed, engine)
            C.require(key not in row_keys, "Duplicate derived row identity")
            row_keys.add(key)
            original = observed.get(key)
            if original is None:
                audit.defer("derived_observed_not_in_raw_snapshot", seed=seed, engine=engine)
                continue
            observed_fields = {"status": row.get("status"), "truth_sha256": C.canonical_sha(row.get("truth", {})),
                               "observed_sha256": C.canonical_sha(row.get("observed", {})),
                               "diagnostics_sha256": C.canonical_sha(row.get("diagnostics", {}))}
            for field, value in observed_fields.items():
                if original[field] != value:
                    audit.issue("derived_observed_raw_mismatch", seed=seed, engine=engine, field=field)
            if all(original[k] == v for k, v in observed_fields.items()):
                matched_observed += 1
            entries = row.get("schemes", {})
            C.require(isinstance(entries, dict), "Derived schemes is not an object")
            expected_names = set(schemes) if original["status"] == "ok" else set()
            if set(entries) != expected_names:
                audit.issue("derived_scheme_set_mismatch", seed=seed, engine=engine)
            for scheme in set(entries) & set(schemes):
                if not isinstance(entries[scheme], dict):
                    audit.issue("invalid_derived_scheme_types_or_decisions", seed=seed, engine=engine, scheme=scheme,
                                message="Scheme record must be an object")
                    continue
                try:
                    I.scheme_state(entries[scheme], schemes[scheme], cfg, seed in run["full_p"])
                except (ValueError, TypeError, KeyError, AttributeError) as exc:
                    audit.issue("invalid_derived_scheme_types_or_decisions", seed=seed, engine=engine, scheme=scheme, message=str(exc))
                expected = expected_schemes[seed, engine, scheme]
                if expected is None:
                    audit.issue("scheme_without_valid_observed", seed=seed, engine=engine, scheme=scheme)
                    continue
                # A raw log and atomic derived file can be copied on opposite
                # sides of a checkpoint. Only monotone attempt/success/failure
                # changes qualify; same-count gt/eq changes remain findings.
                count_keys = ("attempted", "nperm", "failed")
                old = [entries[scheme].get(k) for k in count_keys]
                new = [expected[k] for k in count_keys]
                if all(type(v) is int and v >= 0 for v in old) and old != new:
                    if all(a <= b for a, b in zip(old, new)) or all(a >= b for a, b in zip(old, new)):
                        audit.defer("derived_raw_snapshot_count_skew", seed=seed, engine=engine, scheme=scheme,
                                    derived=dict(zip(count_keys, old)), raw_recomputed=dict(zip(count_keys, new)))
                for field in COUNT_FIELDS:
                    if entries[scheme].get(field) != expected[field]:
                        audit.issue("derived_scheme_raw_count_mismatch", seed=seed, engine=engine, scheme=scheme, field=field,
                                    derived=entries[scheme].get(field), raw_recomputed=expected[field])
        except (ValueError, TypeError, KeyError) as exc:
            audit.issue("invalid_derived_row", line=number, message=str(exc))
    if set(observed) - row_keys:
        audit.defer("raw_observed_not_in_derived_snapshot", identities=[list(k) for k in sorted(set(observed) - row_keys)])
    jobs = computational_jobs(observed, states, seeds, cfg)
    recomputed_progress = {"fit_records": raw_unique, "rows": len(observed), "expected_rows": len(seeds) * len(engines),
                           "computational_jobs": jobs, "statuses": dict(statuses)}
    if progress is not None:
        if not isinstance(progress.get("statuses"), dict) or any(k not in STATES or type(v) is not int or v < 0 for k, v in progress.get("statuses", {}).items()):
            audit.issue("invalid_progress_status_counts")
        for field, value in recomputed_progress.items():
            if field != "statuses" and type(progress.get(field)) is not int:
                audit.issue("invalid_progress_count_type", field=field)
            if progress.get(field) != value:
                audit.defer("progress_raw_snapshot_count_skew", field=field, progress=progress.get(field), raw_recomputed=value)
        for flag in ("done", "stopped"):
            if type(progress.get(flag)) is not bool:
                audit.issue("invalid_progress_flag", field=flag)
    missing_observed = [{"seed": seed, "engine": engine, "action": "fit_planned_observed_once"}
                        for seed in seeds for engine in engines if (seed, engine) not in observed]
    failed_observed = [{"seed": seed, "engine": engine, "status": value["status"], "action": "preserve_failure_no_retry",
                        "truth_audit_error": value.get("truth_audit_error")}
                       for (seed, engine), value in sorted(observed.items()) if value["status"] != "ok"]
    work_terminal = not missing_observed and all(v["work_terminal"] for v in continuation)
    if progress and progress.get("done") is True and not work_terminal:
        audit.issue("done_flag_with_unfinished_prespecified_work")
    lock_after = lock_probe(directory / "runner.lock")
    if "held_by_writer" in (lock_before, lock_after):
        audit.defer("live_runner_lock")
    if not progress or not (progress.get("done") is True or progress.get("stopped") is True):
        audit.defer("live_or_unsettled_checkpoint")
    for path in inputs:
        final = file_state(path)
        if initial[path.name] != final:
            audit.defer("input_changed_during_chunk_audit", path=str(path))
    deferred = bool(audit.deferred)
    checkpoint = not deferred and not audit.findings
    done = checkpoint and progress.get("done") is True and work_terminal
    result = {"schema_version": "confirmation_event_audit_v1", "inspected_utc": datetime.now(timezone.utc).isoformat(),
              "run_id": cfg["run_id"], "chunk_path": str(directory), "cell": cell, "cell_index": cell_index, "chunk": chunk,
              "nchunks": manifest["nchunks"], "planned_seeds": seeds, "planned_engines": engines, "planned_B": budget,
              "planned_observed_records": len(seeds) * len(engines), "planned_full_budget_permutation_records": len(states) * budget,
              "status": "deferred" if deferred else "invalid" if audit.findings else "certified_done" if done else "certified_stopped_checkpoint",
              "certified_done_chunk": done, "checkpoint_certified": checkpoint, "work_terminal": work_terminal,
              "findings_classification": "requires_consistent_snapshot_before_corruption_claim" if deferred else "stable_input_integrity_findings",
              "findings": audit.findings, "deferred_reasons": audit.deferred, "inputs": audit.inputs,
              "finding_counts": dict(audit.finding_counts), "deferred_reason_counts": dict(audit.deferred_counts),
              "lock_probes": {"before": lock_before, "after": lock_after}, "raw_recomputed_progress": recomputed_progress,
              "unique_successful_ordering_records_validated": valid_success, "raw_duplicate_identities": duplicates,
              "derived_observed_exact_matches": matched_observed, "missing_observed": missing_observed, "failed_observed": failed_observed,
              "continuation_manifest_available": checkpoint, "continuation": continuation,
              "continuation_action_counts": dict(Counter(v["action"] for v in continuation)),
              "provenance": {"config_path": str(run["config_path"]), "config_sha256": C.sha(run["config_path"]),
                  "runtime_source_sha256": run["sources"], "manifest_signature": manifest["signature"],
                  "numerical_environment": {k: environment[k] for k in I.ENVIRONMENT_FIELDS}},
              "memory_model": {"permutation_id_bytes": sum(len(s["ids"]) for s in states.values()),
                  "histogram_integer_bins": sum(sum(len(h) for h in s["histograms"]) for s in states.values()),
                  "observed_summary_count": len(observed), "full_permutation_records_retained": 0},
              "notes": ["Missing/failed IDs stay in the frozen planned budget; failed or invalid recorded fits are never proposed for retry.",
                        "Counts use legal-order Kendall inversion integers. On this finite support, integer greater/equal comparisons match the runner's1e-9 tie rule.",
                        "Exact early-stop decisions can make unattempted IDs unnecessary. Preselected full-p datasets still require every unattempted ID.",
                        "A deferred result is not a corruption finding or a certified restart list. Refresh coherent metadata/log files after the writer settles.",
                        "A certified stopped checkpoint can require further work. A certified done chunk may retain failed-fit uncertainty and is not a scientific-validity claim.",
                        "Paired computational jobs use planned shared identities once; each engine's statistical attempts/counts remain separate.",
                        "No file is repaired, truncated, appended or otherwise changed by the audit; input hashes identify the exact bytes read."]}
    return result


def analyze(protocol_path, run_id, chunk_path, project=C.PROJECT, snapshots_root=None):
    protocol_path = Path(protocol_path).resolve()
    protocol_bytes = protocol_path.read_bytes()
    protocol_hash = hashlib.sha256(protocol_bytes).hexdigest()
    protocol = json.loads(protocol_bytes)
    C.require(protocol.get("configuration_frozen_before_confirmation_data") is True and protocol.get("confirmation_namespace") == "42xxxxxx", "Expected frozen confirmation protocol")
    entries = [entry for entry in protocol["runs"] if entry["run_id"] == run_id]
    C.require(len(entries) == 1, "Run must be explicitly and uniquely registered in protocol")
    run = C.frozen_run(entries[0], Path(project), snapshots_root)
    audit_sources = {p: C.sha(p) for p in (Path(__file__), Path(C.__file__), Path(I.__file__))}
    frozen_inputs = {run["config_path"]: entries[0]["config_sha256"], run["snapshot"] / "SHA256SUMS": entries[0]["snapshot_sha256"]}
    frozen_inputs.update({run["snapshot"] / "scripts/v2" / name: digest for name, digest in run["sources"].items()})
    result = audit_chunk(run, chunk_path)
    watched = frozen_inputs | audit_sources | {protocol_path: protocol_hash}
    changed = [str(p) for p, digest in watched.items() if not p.is_file() or C.sha(p) != digest]
    if changed:
        result.update(status="deferred", checkpoint_certified=False, certified_done_chunk=False, continuation_manifest_available=False,
                      findings_classification="requires_consistent_snapshot_before_corruption_claim")
        result["deferred_reasons"].append({"code": "protocol_frozen_or_audit_source_changed_during_audit", "paths": changed})
        result.setdefault("deferred_reason_counts", {})["protocol_frozen_or_audit_source_changed_during_audit"] = 1
    result["protocol"] = {"path": str(protocol_path), "sha256": protocol_hash, "snapshot_sha256": entries[0]["snapshot_sha256"]}
    result["frozen_input_paths"] = {"snapshot": str(run["snapshot"].resolve()), "config": str(run["config_path"].resolve())}
    result["audit_source_sha256"] = {p.name: digest for p, digest in audit_sources.items()}
    return result


def markdown(result):
    lines = ["# Confirmation raw-event audit", "", f"Status: **{result['status']}**. Certified done: **{result['certified_done_chunk']}**; certified checkpoint: **{result['checkpoint_certified']}**.", "",
             f"Run: {result['run_id']}; chunk: {result['chunk_path']}.", "",
             "Certification applies to the hashed input snapshot, not scientific validity or global run completion. A deferred result must not be treated as corrupted data or an actionable certified restart list."]
    if "raw_recomputed_progress" in result:
        p = result["raw_recomputed_progress"]
        lines += ["", f"Raw unique statistical fits: {p['fit_records']}; computational jobs: {p['computational_jobs']}; observed rows: {p['rows']}/{p['expected_rows']}.", "",
                  "| Continuation action | Dataset / engine / scheme groups |", "|---|---:|"]
        lines += [f"| {action} | {count} |" for action, count in result["continuation_action_counts"].items()]
        lines += ["", "The JSON contains each planned seed/engine/scheme's inclusive missing/failed ID ranges, required-resume ranges, recomputed counts and exact decisions. Failed recorded events remain preserved."]
    lines += ["", f"Integrity observations: {len(result['findings'])}; deferral reasons: {len(result['deferred_reasons'])}.", ""]
    for kind in ("findings", "deferred_reasons"):
        lines += ["- " + json.dumps(item, ensure_ascii=False) for item in result[kind]]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--chunk", required=True, type=Path, help="Exact run/cell_i/chunk_j path; no recursive run discovery")
    parser.add_argument("--project-root", type=Path, default=C.PROJECT)
    parser.add_argument("--snapshots-root", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--require-certified", action="store_true", help="Exit3 unless a consistent done/stopped checkpoint is certified")
    args = parser.parse_args(argv)
    try:
        result = analyze(args.protocol, args.run_id, args.chunk, args.project_root, args.snapshots_root)
        prefix = args.output_prefix.resolve()
        protected = [args.chunk.resolve().parent.parent, (args.project_root / "runs/v2/snapshots").resolve(), (args.project_root / "configs").resolve()]
        if args.snapshots_root:
            protected.append(args.snapshots_root.resolve())
        if "frozen_input_paths" in result:
            protected.append(Path(result["frozen_input_paths"]["snapshot"]).resolve())
        exact_inputs = {args.protocol.resolve()}
        exact_inputs.update(Path(value["path"]).resolve() for value in result.get("inputs", {}).values())
        exact_inputs.update(Path(value).resolve() for value in result.get("frozen_input_paths", {}).values())
        exact_inputs.update((Path(__file__).resolve(), Path(C.__file__).resolve(), Path(I.__file__).resolve()))
        outputs = [Path(str(prefix) + ext) for ext in (".json", ".md")]
        C.require(all(path.resolve() not in exact_inputs and not any(path.resolve().is_relative_to(root) for root in protected) for path in outputs), "Output would modify raw/frozen inputs")
        for path, text in zip(outputs, (json.dumps(result, indent=2, allow_nan=False) + "\n", markdown(result))):
            path.parent.mkdir(parents=True, exist_ok=True)
            # Exclusive random creation cannot follow a pre-existing temp
            # symlink into the read-only evidence. replace replaces a final
            # symlink itself, without writing through its target.
            with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=path.name + ".tmp-", delete=False) as handle:
                temporary = Path(handle.name)
                try:
                    handle.write(text)
                    handle.flush()
                    os.fsync(handle.fileno())
                except BaseException:
                    temporary.unlink(missing_ok=True)
                    raise
            try:
                os.replace(temporary, path)
            finally:
                temporary.unlink(missing_ok=True)
        print(json.dumps({"status": result["status"], "checkpoint_certified": result["checkpoint_certified"], "outputs": [str(p) for p in outputs]}))
        return 3 if args.require_certified and not result["checkpoint_certified"] else 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"raw-event audit error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
