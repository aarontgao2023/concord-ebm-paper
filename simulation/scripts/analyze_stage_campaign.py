"""Fixed-R1000 stage overlap diagnostics from already committed observed evidence.

This wrapper never fits a model. It reconstructs only seeds discovered in saved
observed records, through the existing immutable-source/hash-checked per-chunk
postprocessor. Missing evidence stays pending, errors stay in the denominator,
and incomplete cohorts have no between-dataset mean or quantile estimates.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

import analyze_stage_exchangeability as S

VERSION = "stage-campaign-diagnostics-v1"
RUN_ID = "confirm_stage_a"
PLANNED_R = 1000
QUANTILES = (.05, .25, .5, .75, .95)


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _plan(protocol_path, workspace):
    protocol = S._read_json(protocol_path)
    entries = [r for r in protocol["runs"] if r["run_id"] == RUN_ID]
    _require(len(entries) == 1, "Exactly one registered confirm_stage_a run is required")
    entry = entries[0]
    snapshot, declared, cfg = S._verify_snapshot(workspace, entry)
    _require(cfg.get("run_id") == RUN_ID and cfg.get("phase") == "confirmation"
             and cfg.get("datasets") == PLANNED_R and cfg.get("engines") == ["repaired"]
             and cfg.get("cells") == [{"name": "STAGE_H0", "overrides": {"stage_eta": 1.0}}]
             and cfg.get("bperm") == 599 and cfg.get("alpha") == .05 and cfg.get("rule") == "le",
             "This wrapper only supports the frozen stage eta1 R1000 plan")
    for key, actual in (("base_seed", cfg["base_seed"]), ("datasets_per_cell", cfg["datasets"]),
                        ("cells", cfg["cells"]), ("engines", cfg["engines"]), ("bperm", cfg["bperm"]),
                        ("full_p_first_n", cfg.get("full_p_first_n", 0)),
                        ("full_p_seeds", cfg.get("full_p_seeds"))):
        _require(entry.get(key) == actual, f"Protocol/config mismatch: {key}")
    _require(cfg.get("schemes") == {"diagnosis": {"stratify": "diagnosis", "require_max": True},
                                      "oracle_dx_stage": {"stratify": "oracle_dx_stage", "require_max": True}},
             "Unexpected stage permutation definitions")
    base = S._integer(cfg["base_seed"], "base_seed")
    return entry, snapshot, declared, cfg, list(range(base, base + PLANNED_R))


def _manifest(path, cfg, declared):
    value = S._read_json(path / "manifest.json")
    _require(isinstance(value, dict) and value.get("config") == cfg, "Chunk config differs from frozen stage config")
    cell = S._integer(value.get("cell_index"), "cell_index")
    chunk = S._integer(value.get("chunk"), "chunk")
    nchunks = S._integer(value.get("nchunks"), "nchunks", 1)
    _require(cell == 0 and path.parent.name == "cell_0" and path.name == f"chunk_{chunk}"
             and chunk < nchunks <= PLANNED_R, "Chunk directory/manifest layout mismatch")
    seeds = list(range(cfg["base_seed"], cfg["base_seed"] + PLANNED_R))[chunk::nchunks]
    _require(value.get("seeds") == seeds, "Chunk seed list differs from fixed strided layout")
    sources = {name: declared["scripts/v2/" + name] for name in S.RUNTIME_NAMES
               if "scripts/v2/" + name in declared}
    _require(value.get("source_sha256") == sources, "Chunk runtime source closure differs from snapshot")
    identity = {"config": cfg, "cell_index": cell, "chunk": chunk, "nchunks": nchunks, "source_sha256": sources}
    _require(value.get("signature") == hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
             "Chunk manifest signature mismatch")
    return value


def _discover(path, planned):
    """Discover terminal observed records in a read-only committed file prefix.

    Discovery authorizes no reconstruction itself: analyze_chunk independently
    checks the selected records and all source/config/data hashes again.
    """
    selected, nonterminal, seen = {}, {}, set()
    size = path.stat().st_size
    remaining, digest, tail = size, hashlib.sha256(), None
    with path.open("rb") as stream:
        while remaining:
            line = stream.readline(remaining)
            _require(bool(line), "Fit log truncated during discovery")
            remaining -= len(line)
            digest.update(line)
            if not line.endswith(b"\n"):
                tail = hashlib.sha256(line).hexdigest()
                continue
            if not line.strip():
                continue
            record = json.loads(line)
            _require(isinstance(record, dict), "Committed fit record must be an object")
            if record.get("kind") != "observed":
                continue
            seed = S._integer(record.get("seed"), "observed seed")
            _require(seed in planned and record.get("engine") == "repaired"
                     and record.get("scheme") is None and record.get("perm_id") is None,
                     "Unplanned observed record or invalid observed identity")
            _require(seed not in seen, "Duplicate committed observed identity")
            seen.add(seed)
            status = record.get("status")
            _require(isinstance(status, str) and bool(status), "Observed record lacks a status")
            (selected if status in S.TERMINAL_OBSERVED_STATUSES else nonterminal)[seed] = status
    return selected, nonterminal, {"prefix_bytes": size, "prefix_sha256": digest.hexdigest(),
                                  "ignored_uncommitted_tail_sha256": tail}


def _flatten(metrics):
    """Finite per-dataset scalars; block counts keep their fixed DX/stage keys."""
    result = {}

    def visit(value, prefix):
        if isinstance(value, bool):
            result[prefix] = int(value)
        elif isinstance(value, (int, float)):
            _require(math.isfinite(value), f"Nonfinite block diagnostic: {prefix}")
            result[prefix] = value
        elif value is None:
            raise ValueError(f"Undefined block diagnostic: {prefix}")
        elif isinstance(value, dict):
            for name, child in value.items():
                if name not in {"note", "blocks"}:
                    visit(child, f"{prefix}.{name}")
            for block in value.get("blocks", []):
                key = f"{prefix}.blocks.{block['diagnosis']}.k{block['stage']}"
                for name, child in block.items():
                    if name not in {"diagnosis", "stage"}:
                        visit(child, f"{key}.{name}")
    for scheme in ("oracle_dx_stage", "diagnosis"):
        visit(metrics[scheme], scheme)
    return result


def cohort_summary(rows, errors):
    """No descriptive estimate from an execution-selected subset."""
    available = len(rows) == PLANNED_R and not errors and all(r["diagnostic_status"] == "verified" for r in rows)
    if not available:
        return {"available": False, "planned_datasets": PLANNED_R, "statistics": None,
                "reason": "Every planned diagnostic must be verified and integrity errors absent; partial execution is not a population sample."}
    flattened = [_flatten(row["metrics"]) for row in rows]
    keys = set(flattened[0])
    _require(all(set(value) == keys for value in flattened), "Verified datasets have different block metric layouts")
    output = {}
    for key in sorted(keys):
        values = sorted(value[key] for value in flattened)
        quantiles = {}
        for q in QUANTILES:
            index = (len(values) - 1) * q
            lower, upper = math.floor(index), math.ceil(index)
            quantiles[str(q)] = values[lower] + (values[upper] - values[lower]) * (index - lower)
        output[key] = {"n_datasets": PLANNED_R, "mean": statistics.mean(values), "minimum": values[0],
                       "maximum": values[-1], "quantiles": quantiles}
    return {"available": True, "planned_datasets": PLANNED_R,
            "statistical_unit": "one complete planned dataset, including originally failed fits with verified generation truth",
            "quantile_definition": "linear interpolation at (R-1)*q; descriptive across datasets, not a confidence interval",
            "statistics": output}


def analyze(protocol_path, raw_root, *, workspace=None):
    workspace = Path(workspace or Path(__file__).resolve().parents[2]).resolve()
    protocol_path, raw_root = Path(protocol_path).resolve(), Path(raw_root).resolve()
    protocol_hash = S._sha(protocol_path)
    entry, snapshot, declared, cfg, seeds = _plan(protocol_path, workspace)
    _require(S._sha(protocol_path) == protocol_hash, "Protocol changed during initial verification")
    run_root = (raw_root / RUN_ID).resolve()
    _require(run_root.is_relative_to(raw_root), "Run root escapes raw input root")
    rows = {seed: {"seed": seed, "cell": "STAGE_H0", "observed_status": None,
                   "diagnostic_status": "pending_no_committed_observed", "reconstructed": False} for seed in seeds}
    errors, chunks, valid = [], [], []

    def failure(path, stage, exc, affected=()):
        message = {"path": str(path), "stage": stage, "exception": type(exc).__name__,
                   "message": str(exc), "affected_planned_seeds": sorted(affected)}
        errors.append(message)
        for seed in affected:
            rows[seed].update(diagnostic_status="error_" + stage, error_index=len(errors) - 1)
            rows[seed].pop("metrics", None)
        return len(errors) - 1

    # Preflight every discovered manifest before the first possible reconstruction.
    for path in sorted(run_root.glob("cell_*/chunk_*")):
        report = {"path": str(path), "status": "pending_manifest"}
        chunks.append(report)
        try:
            _require(path.resolve().is_relative_to(run_root) and path.is_dir(), "Chunk path escapes raw root or is not a directory")
            _require(re.fullmatch(r"chunk_(0|[1-9][0-9]*)", path.name) is not None and path.parent.name == "cell_0",
                     "Unexpected stage chunk directory")
            _require((path / "manifest.json").resolve().is_relative_to(run_root), "Chunk manifest escapes raw run root")
            if not (path / "manifest.json").is_file():
                continue
            manifest_hash = S._sha(path / "manifest.json")
            manifest = _manifest(path, cfg, declared)
            _require(S._sha(path / "manifest.json") == manifest_hash, "Manifest changed during preflight")
            report.update(planned_seeds=manifest["seeds"], nchunks=manifest["nchunks"],
                          chunk_manifest_sha256=manifest_hash, status="manifest_verified")
            valid.append((path, manifest, report))
        except Exception as exc:
            report.update(status="error_manifest", error_index=failure(path, "manifest", exc))
    layouts = {manifest["nchunks"] for _, manifest, _ in valid}
    if len(layouts) > 1:
        affected = {seed for _, manifest, _ in valid for seed in manifest["seeds"]}
        index = failure(run_root, "inconsistent_layout", ValueError("Multiple nchunks values in one frozen run"), affected)
        for _, _, report in valid:
            report.update(status="error_inconsistent_layout", error_index=index)
        valid = []
    for path, manifest, report in valid:
        planned = manifest["seeds"]
        for seed in planned:
            rows[seed]["chunk_path"] = str(path)
        log = path / "fit_records.jsonl"
        if not log.resolve().is_relative_to(run_root):
            report.update(status="error_discovery", error_index=failure(path, "discovery", ValueError("Fit log escapes raw run root"), planned))
            continue
        if not log.is_file():
            report["status"] = "pending_fit_log"
            for seed in planned:
                rows[seed]["diagnostic_status"] = "pending_fit_log"
            continue
        try:
            selected, nonterminal, evidence = _discover(log, set(planned))
            report.update(discovery_evidence=evidence, committed_terminal_observed_seeds=sorted(selected),
                          pending_nonterminal_observed_seeds=sorted(nonterminal))
            for seed, status in nonterminal.items():
                rows[seed].update(observed_status=status, diagnostic_status="pending_nonterminal_observed")
            for seed, status in selected.items():
                rows[seed].update(observed_status=status, committed_terminal_observed=True)
        except Exception as exc:
            report.update(status="error_discovery", error_index=failure(path, "discovery", exc, planned))
            continue
        if not selected:
            report["status"] = "pending_no_committed_observed"
            continue
        try:
            _require(S._sha(protocol_path) == protocol_hash, "Protocol changed since campaign verification")
            _require(S._sha(path / "manifest.json") == report["chunk_manifest_sha256"], "Manifest changed since campaign preflight")
            result = S.analyze_chunk(protocol_path, path, sorted(selected), workspace=workspace)
            returned = result["datasets"]
            _require(result.get("run_id") == RUN_ID and result.get("cell") == "STAGE_H0"
                     and result.get("requested_seeds") == sorted(selected)
                     and [item["seed"] for item in returned] == sorted(selected), "Per-chunk result identity mismatch")
            _require(result.get("chunk_manifest_sha256") == report["chunk_manifest_sha256"]
                     and S._sha(path / "manifest.json") == report["chunk_manifest_sha256"]
                     and result.get("planned_chunk_datasets") == len(planned),
                     "Per-chunk manifest hash or fixed chunk denominator changed")
            _require(result.get("protocol_sha256") == protocol_hash
                     and result.get("snapshot_sha256") == entry["snapshot_sha256"]
                     and result.get("generator_sha256") == declared["scripts/v2/design_v2.py"],
                     "Per-chunk result source provenance mismatch")
            for item in returned:
                _require(item.get("observed_status") == selected[item["seed"]]
                         and item.get("diagnostic_status") in {"verified", "unavailable_saved_generation_hash_missing"},
                         "Per-chunk observed status or diagnostic changed")
                if item["diagnostic_status"] == "verified":
                    _require(item.get("reconstructed") is True, "Verified diagnostic lacks reconstruction")
                    _flatten(item["metrics"])
            for item in returned:
                rows[item["seed"]].update(item)
            report.update(status="processed_existing_observed", selected_observed_count=len(selected),
                          per_chunk_result={k: v for k, v in result.items() if k != "datasets"})
        except Exception as exc:
            index = failure(path, "chunk_diagnostic", exc, selected)
            for seed in selected:
                # A per-chunk error may occur after an earlier reconstruction;
                # no partial return is accepted, and no calls are claimed absent.
                rows[seed].update(reconstruction_attempted=True, reconstructed=None)
            report.update(status="error_chunk_diagnostic", error_index=index)
    # Another chunk's work can take time: recheck every accepted preflight
    # manifest before accepting any campaign-level results.
    for path, manifest, report in valid:
        try:
            _require(S._sha(path / "manifest.json") == report["chunk_manifest_sha256"], "Manifest changed during campaign diagnostics")
        except Exception as exc:
            index = failure(path, "manifest_changed", exc, manifest["seeds"])
            report.update(status="error_manifest_changed", error_index=index)
    if S._sha(protocol_path) != protocol_hash:
        failure(protocol_path, "protocol_changed", ValueError("Protocol changed during diagnostics"), seeds)
    ordered = [rows[seed] for seed in seeds]
    try:
        summary = cohort_summary(ordered, errors)
    except Exception as exc:
        failure(run_root, "cohort_summary", exc, seeds)
        summary = cohort_summary(ordered, errors)
    verified = sum(row["diagnostic_status"] == "verified" for row in ordered)
    pending = sum(row["diagnostic_status"].startswith("pending_") for row in ordered)
    unavailable = sum(row["diagnostic_status"].startswith("unavailable_") for row in ordered)
    return {"schema_version": VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
            "run_id": RUN_ID, "cell": "STAGE_H0", "planned_datasets": PLANNED_R,
            "planned_seed_start": seeds[0], "planned_seed_end_inclusive": seeds[-1],
            "committed_terminal_observed_count": sum(row.get("committed_terminal_observed", False) for row in ordered),
            "verified_datasets": verified, "pending_datasets": pending, "unavailable_datasets": unavailable,
            "error_datasets": PLANNED_R - verified - pending - unavailable,
            "status": "complete" if summary["available"] else "error_with_planned_denominator" if errors
                      else "pending_existing_evidence" if pending else "terminal_with_unavailable_diagnostics",
            "cohort_summary": summary, "diagnostic_status_counts": dict(Counter(row["diagnostic_status"] for row in ordered)),
            "protocol_path": str(protocol_path), "protocol_sha256": protocol_hash,
            "snapshot": str(snapshot), "snapshot_sha256": entry["snapshot_sha256"],
            "config_sha256": entry["config_sha256"], "raw_root": str(raw_root),
            "wrapper_sha256": S._sha(__file__), "per_chunk_postprocessor_sha256": S._sha(S.__file__),
            "observed_evidence_required": True, "fitting_performed": False, "new_experiment": False,
            "partial_summary_policy": "Counts and per-seed diagnostics only; no mean/quantiles of early completers.",
            "interpretation": S.INTERPRETATION, "errors": errors, "chunks": chunks, "datasets": ordered}


def render_markdown(result):
    lines = ["# Stage campaign overlap diagnostics", "", S.INTERPRETATION, "",
             f"Status: **{result['status']}**. Verified **{result['verified_datasets']}/{PLANNED_R}**; "
             f"committed terminal observed records {result['committed_terminal_observed_count']}; "
             f"pending {result['pending_datasets']}, unavailable {result['unavailable_datasets']}, errors {result['error_datasets']}.", "",
             "No mean or quantile from a partially executed cohort is a planned-R1000 summary.", ""]
    summary = result["cohort_summary"]
    if summary["available"]:
        lines += ["| Dataset-level diagnostic | Mean over R1000 | 5% | Median | 95% |", "|---|---:|---:|---:|---:|"]
        for key, value in summary["statistics"].items():
            if ".blocks." not in key:
                q = value["quantiles"]
                lines.append(f"| {key} | {value['mean']:.6g} | {q['0.05']:.6g} | {q['0.5']:.6g} | {q['0.95']:.6g} |")
    else:
        lines += ["Cohort mean/quantiles: **unavailable**. Every missing or failed diagnostic remains in the denominator."]
    lines += ["", "| Planned seed | Original fit status | Diagnostic status |", "|---|---|---|"]
    for row in result["datasets"]:
        lines.append(f"| {row['seed']} | {row['observed_status']} | {row['diagnostic_status']} |")
    if result["errors"]:
        lines += ["", "Errors (processing continued to other eligible chunks):", ""]
        for error in result["errors"]:
            lines.append(f"- {error['stage']} at `{error['path']}`: {error['exception']}: {error['message']}")
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True, help="Parent directory containing confirm_stage_a")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()
    S.validate_output_paths(args.protocol, args.raw_root, [args.output_json, args.output_md], workspace=args.workspace)
    result = analyze(args.protocol, args.raw_root, workspace=args.workspace)
    S.write_outputs(args.protocol, args.raw_root,
                    [(args.output_json, json.dumps(result, indent=2, allow_nan=False) + "\n"),
                     (args.output_md, render_markdown(result))], workspace=args.workspace)


if __name__ == "__main__":
    main()
