#!/usr/bin/env python3
"""Serial read-only audits of a fixed list of settled confirmation chunks.

No scheduling, SSH, simulation, fitting, repair, cache or implicit run discovery.
The wall-time budget is checked between chunks; an ongoing chunk is not killed.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
import time

import audit_confirmation_events as A

C = A.C
VERSION = "confirmation-event-batch-v1"


def utc():
    return datetime.now(timezone.utc).isoformat()


def exclusive_write(path, value):
    """Publish complete bytes without overwriting even a dangling symlink."""
    text = value if isinstance(value, str) else json.dumps(value, indent=2, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".publish-", delete=False) as handle:
        temporary = Path(handle.name)
        try:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        except BaseException:
            temporary.unlink(missing_ok=True)
            raise
    try:
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def changed_sources(watched):
    changed = []
    for path, digest in watched.items():
        try:
            equal = C.sha(path) == digest
        except OSError:
            equal = False
        if not equal:
            changed.append(str(path))
    return changed


def snapshot_json(path):
    before = A.file_state(path)
    if before is None:
        raise FileNotFoundError(str(path))
    raw = path.read_bytes()
    after = A.file_state(path)
    value = json.loads(raw)
    C.require(isinstance(value, dict), "Expected metadata object")
    return value, {"path": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
                   "before": before, "after": after, "stable": before == after}


def inspect_selection(directory, run):
    """Lightweight selection only; scientific raw bytes are not read here."""
    item = {"run_id": run["run_id"], "chunk_path": str(directory), "selection_status": "error", "inputs": {}}
    try:
        root = directory.parent.parent
        C.require(directory.resolve().is_relative_to(root.resolve()), "Chunk symlink escapes registered run")
        C.require(re.fullmatch(r"cell_\d+", directory.parent.name) is not None and
                  re.fullmatch(r"chunk_\d+", directory.name) is not None, "Noncanonical chunk directory")
        C.require(directory.parent.name == "cell_" + str(int(directory.parent.name[5:])) and
                  directory.name == "chunk_" + str(int(directory.name[6:])), "Noncanonical numeric directory alias")
        manifest, evidence = snapshot_json(directory / "manifest.json")
        item["inputs"]["manifest.json"] = evidence
        cell, chunk, seeds, _ = A.verify_manifest(run, directory, manifest)
        item.update(cell_index=cell, chunk_index=chunk, planned_seeds=seeds)
        progress, evidence = snapshot_json(directory / "progress.json")
        item["inputs"]["progress.json"] = evidence
        C.require(type(progress.get("done")) is bool and type(progress.get("stopped")) is bool,
                  "Progress done/stopped must be booleans")
        item["progress_flags"] = {k: progress[k] for k in ("done", "stopped")}
        item["lock_probe"] = A.lock_probe(directory / "runner.lock")
        if not all(v["stable"] for v in item["inputs"].values()):
            item.update(selection_status="deferred", selection_reason="metadata_changed_during_selection")
        elif item["lock_probe"] == "held_by_writer" or not (progress["done"] or progress["stopped"]):
            item.update(selection_status="live", selection_reason="writer_lock_or_active_progress")
        else:
            item.update(selection_status="candidate", selection_reason="settled_done_or_stopped_and_no_writer_lock")
    except (FileNotFoundError, json.JSONDecodeError) as exc:
        item.update(selection_status="deferred", selection_reason="metadata_missing_or_unreadable", message=str(exc))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        item.update(selection_status="error", selection_reason="invalid_selection_metadata_or_allocation", message=str(exc))
    return item


def prepare(protocol_path, run_ids, results_root, project, snapshots_root):
    C.require(run_ids and len(set(run_ids)) == len(run_ids), "Specify distinct --run-id values")
    protocol_raw = protocol_path.read_bytes()
    protocol = json.loads(protocol_raw)
    C.require(protocol.get("configuration_frozen_before_confirmation_data") is True and
              protocol.get("confirmation_namespace") == "42xxxxxx", "Expected frozen confirmation protocol")
    sources = [Path(__file__).resolve(), Path(A.__file__).resolve(), Path(C.__file__).resolve(), Path(A.I.__file__).resolve()]
    watched = {p: C.sha(p) for p in sources}
    watched[protocol_path] = hashlib.sha256(protocol_raw).hexdigest()
    runs, items, scopes = {}, [], []
    for run_id in sorted(run_ids):
        entries = [entry for entry in protocol["runs"] if entry["run_id"] == run_id]
        C.require(len(entries) == 1, "Requested run must be uniquely registered: " + run_id)
        entry = entries[0]
        run = C.frozen_run(entry, project, snapshots_root)
        runs[run_id] = run
        watched[run["config_path"]] = entry["config_sha256"]
        watched[run["snapshot"] / "SHA256SUMS"] = entry["snapshot_sha256"]
        watched.update({run["snapshot"] / "scripts/v2" / name: digest for name, digest in run["sources"].items()})
        root = C.safe_child(results_root, run_id)
        manifests = sorted(root.glob("cell_*/chunk_*/manifest.json")) if root.is_dir() else []
        scopes.append({"run_id": run_id, "results_path": str(root), "results_path_exists": root.is_dir(),
                       "existing_manifest_count": len(manifests), "planned_rows": len(run["names"])*len(run["seeds"])*len(run["cfg"]["engines"]),
                       "scope": "existing_manifest_chunks_only; absent/uncreated chunks are not certified"})
        items.extend(inspect_selection(path.parent, run) for path in manifests)
    items.sort(key=lambda v: (v["run_id"], v.get("cell_index", 10**9), v.get("chunk_index", 10**9), v["chunk_path"]))
    return runs, watched, items, scopes


def summary_markdown(result):
    lines = ["# Confirmation raw-event batch audit", "", f"Status: **{result['status']}**. Source stability: **{result['source_stable']}**.", "",
             f"Existing manifest directories: {result['discovered_chunks']}; settled candidates: {result['candidate_chunks']}; selected: {result['selected_chunks']}; audits started: {result['audits_started']}.", "",
             "This batch does not assess whole-run completion. Unselected, unprocessed, live, deferred and error chunks are not certified.", "",
             "| Final category | Chunks |", "|---|---:|"]
    lines += [f"| {name} | {count} |" for name, count in result["final_category_counts"].items()]
    lines += ["", "The fixed manifest/candidate list and selection reasons are in selection.json. Each processed chunk has its own immutable JSON/Markdown audit. The summary is authoritative for batch source stability.", "",
              "Wall time is checked before starting each chunk. One ongoing raw scan can exceed the remaining budget; this is not a hard wall-time guarantee.", ""]
    if result["changed_source_paths"]:
        lines += ["Changed sources prevent batch certification:"] + ["- " + p for p in result["changed_source_paths"]]
    return "\n".join(lines) + "\n"


def run_batch(protocol_path, run_ids, results_root, output_dir, *, project=C.PROJECT, snapshots_root=None,
              offset=0, max_chunks=None, walltime_budget_s=None):
    started = time.monotonic()
    C.require(type(offset) is int and offset >= 0, "offset must be a nonnegative integer")
    C.require(max_chunks is None or type(max_chunks) is int and max_chunks > 0, "max_chunks must be positive")
    C.require(walltime_budget_s is None or type(walltime_budget_s) in (int, float) and
              math.isfinite(walltime_budget_s) and walltime_budget_s > 0, "walltime_budget_s must be positive finite")
    protocol_path, results_root, project = (Path(p).resolve() for p in (protocol_path, results_root, project))
    output_dir = Path(output_dir).absolute()
    C.require(not output_dir.exists() and not output_dir.is_symlink(), "Output directory already exists; no overwrite")
    output_dir = output_dir.resolve()
    snapshots_root = Path(snapshots_root).resolve() if snapshots_root is not None else None
    runs, watched, items, scopes = prepare(protocol_path, run_ids, results_root, project, snapshots_root)
    protected = [results_root, (project / "configs").resolve(), (project / "runs/v2/snapshots").resolve()]
    protected += [run["snapshot"].resolve() for run in runs.values()]
    C.require(not any(output_dir.is_relative_to(root) for root in protected) and output_dir not in watched,
              "Output would modify raw/frozen inputs")
    changed = changed_sources(watched)
    C.require(not changed, "Sources changed while preparing candidate list: " + repr(changed))
    candidates = [i for i, item in enumerate(items) if item["selection_status"] == "candidate"]
    selected_indexes = candidates[offset: None if max_chunks is None else offset+max_chunks]
    selected_set = set(selected_indexes)
    for rank, index in enumerate(candidates):
        item = items[index]
        item["candidate_index"] = rank
        item["selected"] = index in selected_set
        item["slice_reason"] = "selected" if index in selected_set else "not_selected_offset" if rank < offset else "not_selected_max_chunks"
    for item in items:
        item.setdefault("selected", False)
    source_hashes = {str(p): digest for p, digest in watched.items()}
    selection = {"schema_version": VERSION, "fixed_at_utc": utc(), "run_ids": sorted(run_ids), "offset": offset,
                 "max_chunks": max_chunks, "walltime_budget_s": walltime_budget_s, "scopes": scopes,
                 "source_sha256": source_hashes, "items": items,
                 "candidate_order": [items[i]["chunk_path"] for i in candidates],
                 "selected_order": [items[i]["chunk_path"] for i in selected_indexes]}
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir()  # Atomic refusal if another process created this output.
    exclusive_write(output_dir / "selection.json", selection)
    selection_hash = C.sha(output_dir / "selection.json")
    outcomes = []
    stop_reason, audits_started, ever_changed = None, 0, set()
    for item in items:
        outcome = {"run_id": item["run_id"], "chunk_path": item["chunk_path"], "selected": item["selected"],
                   "checkpoint_certified": False, "certified_done_chunk": False}
        if not item["selected"]:
            outcome["category"] = item.get("slice_reason", item["selection_status"])
            outcomes.append(outcome)
            continue
        changed = changed_sources(watched)
        ever_changed.update(changed)
        if changed:
            stop_reason = "source_changed"
        if stop_reason is None and walltime_budget_s is not None and time.monotonic()-started >= walltime_budget_s:
            stop_reason = "walltime_budget"
        if stop_reason:
            outcome["category"] = "unprocessed_" + stop_reason
            outcomes.append(outcome)
            continue
        # The candidate list is fixed, but a selected writer may have restarted.
        current = inspect_selection(Path(item["chunk_path"]), runs[item["run_id"]])
        selected_manifest_hash = item["inputs"]["manifest.json"]["sha256"]
        if current["inputs"].get("manifest.json", {}).get("sha256") != selected_manifest_hash:
            outcome.update(category="deferred", recheck=current, manifest_binding_reason="manifest_changed_after_fixed_selection")
            outcomes.append(outcome)
            continue
        if current["selection_status"] != "candidate":
            outcome.update(category=current["selection_status"], recheck=current)
            outcomes.append(outcome)
            continue
        if walltime_budget_s is not None and time.monotonic()-started >= walltime_budget_s:
            stop_reason = "walltime_budget"
            outcome["category"] = "unprocessed_walltime_budget"
            outcomes.append(outcome)
            continue
        audit_started = time.monotonic()
        audits_started += 1
        try:
            result = A.analyze(protocol_path, item["run_id"], item["chunk_path"], project, snapshots_root)
            if result.get("inputs", {}).get("manifest.json", {}).get("sha256") != selected_manifest_hash:
                result.update(status="deferred", checkpoint_certified=False, certified_done_chunk=False,
                              continuation_manifest_available=False, findings_classification="requires_consistent_snapshot_before_corruption_claim")
                result["deferred_reasons"].append({"code": "audit_manifest_differs_from_fixed_selection",
                                                  "selected_sha256": selected_manifest_hash})
            changed = changed_sources(watched)
            ever_changed.update(changed)
            if changed:
                stop_reason = "source_changed"
                result.update(status="deferred", checkpoint_certified=False, certified_done_chunk=False,
                              continuation_manifest_available=False, findings_classification="requires_consistent_snapshot_before_corruption_claim")
                result["deferred_reasons"].append({"code": "batch_source_changed_during_chunk", "paths": changed})
            outcome.update(category=result["status"], checkpoint_certified=result["checkpoint_certified"],
                           certified_done_chunk=result["certified_done_chunk"])
            md = A.markdown(result)
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            result = {"status": "error", "error_type": type(exc).__name__, "message": str(exc),
                      "run_id": item["run_id"], "chunk_path": item["chunk_path"], "checkpoint_certified": False,
                      "certified_done_chunk": False, "continuation_manifest_available": False}
            outcome["category"] = "error"
            md = "# Confirmation raw-event audit error\n\n" + json.dumps(result, indent=2) + "\n"
        outcome["elapsed_seconds"] = time.monotonic()-audit_started
        stem = f"{item['run_id']}__cell_{item['cell_index']}__chunk_{item['chunk_index']}"
        result["batch_evidence"] = {"selection_sha256": selection_hash, "candidate_index": item["candidate_index"],
                                    "batch_source_sha256": source_hashes, "elapsed_seconds": outcome["elapsed_seconds"]}
        for suffix, value in ((".json", result), (".md", md)):
            path = output_dir / (stem + suffix)
            exclusive_write(path, value)
        outcome["audit_json"] = stem + ".json"
        outcome["audit_json_sha256"] = C.sha(output_dir / outcome["audit_json"])
        outcomes.append(outcome)
    # The final summary binds immutable chunk evidence to the original
    # candidate allocation, even if a manifest was replaced after its scan.
    # Earlier files stay intact; they certify only their own audited bytes.
    manifest_changes = []
    for item, outcome in zip(items, outcomes):
        if not item["selected"]:
            continue
        path = Path(item["chunk_path"]) / "manifest.json"
        initial_hash = item["inputs"]["manifest.json"]["sha256"]
        try:
            same = C.sha(path) == initial_hash
        except OSError:
            same = False
        outcome["manifest_matches_fixed_selection_at_finalize"] = same
        if not same:
            manifest_changes.append(str(path))
            outcome["checkpoint_certified"] = outcome["certified_done_chunk"] = False
            outcome["manifest_binding_reason"] = "manifest_changed_after_fixed_selection"
            if outcome["category"] in ("certified_done", "certified_stopped_checkpoint"):
                outcome["category"] = "deferred"
    ever_changed.update(changed_sources(watched))
    changed = sorted(ever_changed)
    stable = not changed
    counts = dict(sorted(Counter(v["category"] for v in outcomes).items()))
    selected_outcomes = [v for v in outcomes if v["selected"]]
    selected_all_certified = bool(selected_outcomes) and stable and all(v["checkpoint_certified"] for v in selected_outcomes)
    summary = {"schema_version": VERSION, "finished_utc": utc(), "status": "source_changed" if not stable else "selected_chunks_certified" if selected_all_certified else "partial_or_uncertified",
               "source_stable": stable, "source_sha256": source_hashes, "changed_source_paths": changed,
               "changed_selected_manifest_paths": manifest_changes,
               "selection_sha256": selection_hash, "scopes": scopes, "discovered_chunks": len(items), "candidate_chunks": len(candidates),
               "selected_chunks": len(selected_indexes), "audits_started": audits_started, "final_category_counts": counts,
               "selected_all_checkpoint_certified": selected_all_certified, "whole_run_completion_assessed": False,
               "batch_checkpoint_certified_count": sum(v["checkpoint_certified"] for v in outcomes) if stable else 0,
               "batch_done_certified_count": sum(v["certified_done_chunk"] for v in outcomes) if stable else 0,
               "stop_reason": stop_reason, "elapsed_seconds": time.monotonic()-started, "outcomes": outcomes,
               "limitations": ["Only the fixed existing-manifest list is in scope; absent and later-created chunks are unassessed.",
                   "Offsets and max-chunks apply to the sorted settled candidate list. Unselected/unprocessed chunks never count as passed.",
                   "Wall-time budget is checked between chunks, not a hard interruption deadline.",
                   "Per-chunk certificates apply to their hashed bytes; batch source stability must also hold.",
                   "No original input is modified, no failures are retried, and no scheduling or scientific computation is invoked."]}
    exclusive_write(output_dir / "summary.json", summary)
    exclusive_write(output_dir / "summary.md", summary_markdown(summary))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--run-id", action="append", required=True, help="Repeat for each explicitly registered run")
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True, help="New directory; existing paths are refused")
    parser.add_argument("--project-root", type=Path, default=C.PROJECT)
    parser.add_argument("--snapshots-root", type=Path)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--max-chunks", type=int)
    parser.add_argument("--walltime-budget-s", type=float)
    args = parser.parse_args(argv)
    try:
        summary = run_batch(args.protocol, args.run_id, args.results_root, args.output_dir, project=args.project_root,
                            snapshots_root=args.snapshots_root, offset=args.offset, max_chunks=args.max_chunks,
                            walltime_budget_s=args.walltime_budget_s)
        print(json.dumps({"status": summary["status"], "source_stable": summary["source_stable"],
                          "final_category_counts": summary["final_category_counts"], "output_dir": str(args.output_dir.resolve())}))
        return 0 if summary["selected_all_checkpoint_certified"] else 3
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"raw-event batch audit error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
