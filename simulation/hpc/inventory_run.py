#!/usr/bin/env python3
"""Inventory a frozen v2 run using manifest/progress/rows files only.

    python hpc/v2/inventory_run.py --config FROZEN_CONFIG.json \
      --root RUN_OUTPUT_ROOT --output inventory.json

The root contains cell_<index>/chunk_<index> directories. No fit_records file is
opened or scanned. Planned datasets, engines, cells, budgets and full-p cohorts
come exclusively from the supplied config; a done flag never supplies a missing
dataset or a missing permutation decision.

collection_complete means all planned observed identities have terminal records
and all runnable scheme work is terminal (including exhausted failed budgets).
complete additionally requires successful observed fits and every planned primary
decision/full-p obligation. ready_for_final_analysis additionally requires intact,
consistent provenance, unique identities and a stable snapshot. Failures that
leave a primary decision unresolved always prevent readiness. Resolved failed
permutations are retained and counted, not silently dropped.

Readiness concerns this run's declared phase, not scientific validity or adequate
Monte Carlo precision. The program exits zero after producing a valid inventory,
even when the run is incomplete. --require-ready requests exit status 3 instead.
Use --output inventory.json --quiet to save the inventory without stdout JSON.
Only the standard library is required.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys

RUNTIME_FILES = ("run_v2.py", "design_v2.py", "engine_v2.py", "dependency_v2.py", "oracle_engine_v2.py")
OPTIONAL_RUNTIME_FILES = ("oracle_estimated_prior_v2.py", "paired_engine_v2.py", "fast_likelihood_v2.py")
ENVIRONMENT_FIELDS = ("versions", "source_sha256", "wheel_sha256")
NUMERICAL_PACKAGES = ("python", "pyebm", "numpy", "scipy", "pandas", "scikit-learn", "statsmodels")
PAIR_COUNT = 3


def utc():
    return datetime.now(timezone.utc).isoformat()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def is_hash(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def count(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


class Audit:
    def __init__(self):
        self.issues = []
        self.stable = True

    def issue(self, severity, code, message, path=None):
        item = {"severity": severity, "code": code, "message": message}
        if path is not None:
            item["path"] = str(path)
        self.issues.append(item)

    def error(self, code, message, path=None):
        self.issue("error", code, message, path)

    def warning(self, code, message, path=None):
        self.issue("warning", code, message, path)

    def json_file(self, path, required=True):
        if not path.exists():
            self.issue("error" if required else "warning", "missing_file", "File is absent", path)
            return None
        try:
            result = json.loads(path.read_text())
            if not isinstance(result, dict):
                raise ValueError("Expected a JSON object")
            return result
        except (OSError, ValueError) as exc:
            self.error("invalid_json", str(exc), path)
            return None


def validate_config(config):
    for name in ("run_id", "phase"):
        if not isinstance(config.get(name), str) or not config[name]:
            raise ValueError(f"config.{name} must be a nonempty string")
    base, datasets = count(config.get("base_seed"), "base_seed"), count(config.get("datasets"), "datasets")
    if datasets == 0:
        raise ValueError("datasets must be positive")
    engines = config.get("engines")
    if not isinstance(engines, list) or not engines or any(not isinstance(e, str) or not e for e in engines):
        raise ValueError("engines must be a nonempty list of names")
    if len(set(engines)) != len(engines):
        raise ValueError("Duplicate configured engines")
    cells = config.get("cells")
    if not isinstance(cells, list) or not cells:
        raise ValueError("cells must be a nonempty list")
    names = []
    for cell in cells:
        if not isinstance(cell, dict) or not isinstance(cell.get("name"), str):
            raise ValueError("Each cell must declare its name")
        name = cell.get("id", cell["name"])
        if not isinstance(name, str) or not name:
            raise ValueError("Cell IDs must be nonempty strings")
        names.append(name)
    if len(set(names)) != len(names):
        raise ValueError("Cell IDs must be unique within a run")
    schemes = config.get("schemes", {})
    if not isinstance(schemes, dict):
        raise ValueError("schemes must be an object")
    budget = count(config.get("bperm", 0), "bperm")
    if schemes and budget == 0:
        raise ValueError("Nonempty schemes require a positive bperm")
    alpha = Fraction(str(config.get("alpha", .05)))
    if not 0 < alpha < 1 or config.get("rule", "strict") not in ("strict", "le"):
        raise ValueError("Invalid configured alpha/rule")
    for name, definition in schemes.items():
        if not isinstance(definition, dict):
            raise ValueError(f"Scheme {name} definition must be an object")
        target = definition.get("tested_pair_index")
        if target is not None and (type(target) is not int or target not in range(PAIR_COUNT)):
            raise ValueError(f"Scheme {name} has an invalid tested_pair_index")
        if "require_max" in definition and type(definition["require_max"]) is not bool:
            raise ValueError(f"Scheme {name} require_max must be bool")
    seeds = set(range(base, base + datasets))
    full = config.get("full_p_seeds", [])
    if not isinstance(full, list):
        raise ValueError("full_p_seeds must be a list")
    selected = {count(seed, "full_p_seed") for seed in full}
    first = count(config.get("full_p_first_n", 0), "full_p_first_n")
    if first > datasets:
        raise ValueError("full_p_first_n exceeds planned datasets")
    selected.update(range(base, base + first))
    if not selected <= seeds:
        raise ValueError("full_p_seeds includes unplanned seeds")
    return names, seeds, selected


def bounds(exceedances, successes, budget, alpha, comparisons, rule):
    threshold = Fraction(str(alpha)) / comparisons
    def reject(numerator):
        left, right = numerator * threshold.denominator, (budget + 1) * threshold.numerator
        return left < right if rule == "strict" else left <= right
    optimistic, pessimistic = reject(1 + exceedances), reject(1 + exceedances + budget - successes)
    return optimistic if optimistic == pessimistic else None


def scheme_state(entry, definition, config, force_full):
    """Independently certify recorded counts; never trust a done/exact flag alone."""
    if not isinstance(entry, dict):
        raise ValueError("Scheme record must be an object")
    budget = count(entry.get("requested_nperm"), "requested_nperm")
    if budget != config["bperm"]:
        raise ValueError("requested_nperm differs from frozen bperm")
    successes = count(entry.get("nperm"), "nperm")
    attempted = count(entry.get("attempted"), "attempted")
    if not successes <= attempted <= budget:
        raise ValueError("Expected nperm <= attempted <= bperm")
    failed = count(entry.get("failed", attempted - successes), "failed")
    if failed != attempted - successes:
        raise ValueError("failed differs from attempted - nperm")
    if type(entry.get("complete")) is not bool or entry["complete"] != (successes == budget):
        raise ValueError("complete differs from nperm == bperm")
    if "definition" in entry and entry["definition"] != definition:
        raise ValueError("Scheme definition differs from frozen config")
    gt, eq = entry.get("gt"), entry.get("eq")
    if not isinstance(gt, list) or not isinstance(eq, list) or len(gt) != 3 or len(eq) != 3:
        raise ValueError("gt and eq must contain three integer counts")
    counts = [count(g, "gt") + count(e, "eq") for g, e in zip(gt, eq)]
    if any(value > successes for value in counts):
        raise ValueError("gt + eq exceeds nperm")
    max_count = None
    if "maxgt" in entry or "maxeq" in entry:
        max_count = count(entry.get("maxgt"), "maxgt") + count(entry.get("maxeq"), "maxeq")
        if max_count > successes:
            raise ValueError("maxgt + maxeq exceeds nperm")
    if definition.get("require_max", False) and max_count is None:
        raise ValueError("Required max counts are missing")
    alpha, rule = config.get("alpha", .05), config.get("rule", "strict")
    exact = entry.get("exact_decision", {})
    if not isinstance(exact, dict):
        raise ValueError("exact_decision must be an object")
    if exact and Fraction(str(exact.get("alpha"))) != Fraction(str(alpha)):
        raise ValueError("Recorded exact-decision alpha differs from frozen alpha")
    for recorded_rule in ("strict", "le"):
        declaration = exact.get(recorded_rule)
        if declaration is None:
            continue
        values = declaration.get("pair_reject")
        if not isinstance(values, list) or len(values) != 3:
            raise ValueError("Recorded pair decisions must contain three bool/null entries")
        derived = [bounds(value, successes, budget, alpha, 3, recorded_rule) for value in counts]
        if any(value is not None and type(value) is not bool for value in values):
            raise ValueError("Recorded pair decision is not bool/null")
        if values != derived:
            raise ValueError("Recorded exact pair decisions disagree with independent completion bounds")
        if "max_reject" in declaration:
            value = declaration["max_reject"]
            if value is not None and type(value) is not bool:
                raise ValueError("Recorded max decision is not bool/null")
            derived_max = bounds(max_count, successes, budget, alpha, 1, recorded_rule) if max_count is not None else None
            if value != derived_max:
                raise ValueError("Recorded max decision disagrees with independent completion bounds")
    pair_decisions = [bounds(value, successes, budget, alpha, 3, rule) for value in counts]
    targets = [definition["tested_pair_index"]] if "tested_pair_index" in definition else list(range(3))
    resolved = all(pair_decisions[index] is not None for index in targets)
    if definition.get("require_max", False):
        resolved = resolved and bounds(max_count, successes, budget, alpha, 1, rule) is not None
    full = successes == budget
    if full:
        state = "full"
    elif force_full:
        state = "full_p_failed" if attempted == budget else "pending_full_p"
    elif resolved:
        state = "budget_exhausted_resolved" if attempted == budget else "early_stopped_resolved"
    else:
        state = "terminal_unresolved" if attempted == budget else "pending"
    return {"state": state, "nperm": successes, "attempted": attempted, "failed": failed,
            "primary_decisions_resolved": resolved, "full": full,
            "obligation_satisfied": resolved and (full or not force_full),
            "terminal": attempted == budget or (resolved and not force_full)}


def validate_observed(row):
    observed = row.get("observed")
    if not isinstance(observed, dict):
        raise ValueError("Successful observed fit lacks an observed object")
    orders, taus = observed.get("orderings"), observed.get("taus")
    if not isinstance(orders, list) or len(orders) != 3:
        raise ValueError("Successful fit must contain three orderings")
    width = len(orders[0]) if isinstance(orders[0], list) else 0
    if width < 2:
        raise ValueError("Orderings must contain at least two events")
    for order in orders:
        if not isinstance(order, list) or len(order) != width or any(type(v) is not int for v in order) or sorted(order) != list(range(width)):
            raise ValueError("Observed ordering is not an event permutation")
    if not isinstance(taus, list) or len(taus) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) or not 0 <= v <= 1 for v in taus):
        raise ValueError("Observed taus must be three finite values in [0,1]")


def source_expectations(config_path, source_root):
    candidate = Path(source_root) if source_root else config_path.parent.parent.parent / "scripts" / "v2"
    if not candidate.is_dir():
        return {}, None
    # Historical frozen runs have only the required five modules. New modules
    # become mandatory hash matches whenever their frozen source is present.
    existing = {name: hashlib.sha256((candidate / name).read_bytes()).hexdigest()
                for name in RUNTIME_FILES + OPTIONAL_RUNTIME_FILES if (candidate / name).is_file()}
    return existing, str(candidate.resolve())


def inventory(config_path, root, source_root=None, nchunks_override=None):
    started = utc()
    config_path, root = Path(config_path).resolve(), Path(root).resolve()
    audit = Audit()
    config = audit.json_file(config_path)
    if config is None:
        return {"schema_version": "run_inventory_v2.1", "complete": False, "ready_for_final_analysis": False,
                "collection_complete": False, "issues": audit.issues, "inspected_utc": started}
    try:
        cell_names, seeds, full_p = validate_config(config)
        if nchunks_override is not None and count(nchunks_override, "nchunks") == 0:
            raise ValueError("nchunks must be positive")
    except (ValueError, TypeError) as exc:
        audit.error("invalid_config", str(exc), config_path)
        return {"schema_version": "run_inventory_v2.1", "complete": False, "ready_for_final_analysis": False,
                "collection_complete": False, "issues": audit.issues, "inspected_utc": started}
    expected_sources, expected_source_root = source_expectations(config_path, source_root)
    if expected_sources and not set(RUNTIME_FILES) <= set(expected_sources):
        audit.error("incomplete_frozen_sources", "Frozen source directory lacks expected runtime files", expected_source_root)
    if not expected_sources:
        audit.warning("frozen_source_files_unavailable", "Source consistency/signatures will be checked, but no local frozen source copy is available")
    engines, schemes = config["engines"], config.get("schemes", {})
    expected_rows = len(cell_names) * len(seeds) * len(engines)
    records = defaultdict(list)
    chunk_reports, source_profiles, environments = [], defaultdict(list), defaultdict(list)
    source_values, environment_values, environment_paths = {}, {}, set()
    manifests_by_cell = defaultdict(list)
    physical_rows = 0
    for directory in sorted(root.glob("cell_*/chunk_*")):
        if not directory.is_dir():
            continue
        cell_match = re.fullmatch(r"cell_(\d+)", directory.parent.name)
        chunk_match = re.fullmatch(r"chunk_(\d+)", directory.name)
        if cell_match is None or chunk_match is None:
            audit.error("invalid_chunk_path", "Expected cell_i/chunk_j", directory)
            continue
        cell_index, chunk_index = int(cell_match[1]), int(chunk_match[1])
        report = {"path": str(directory.relative_to(root)), "cell_index": cell_index, "chunk": chunk_index}
        manifest = audit.json_file(directory / "manifest.json")
        progress = audit.json_file(directory / "progress.json", required=False)
        report["progress_done"] = progress.get("done") if progress else None
        report["progress_stopped"] = progress.get("stopped") if progress else None
        chunk_seeds = None
        if cell_index >= len(cell_names):
            audit.error("unexpected_cell_directory", "Cell index is not in frozen config", directory)
        if manifest is not None:
            try:
                n_chunks = count(manifest.get("nchunks"), "manifest.nchunks")
                if n_chunks == 0 or not 0 <= chunk_index < n_chunks:
                    raise ValueError("Chunk index is outside manifest.nchunks")
                if nchunks_override is not None and n_chunks != nchunks_override:
                    raise ValueError("manifest.nchunks differs from --nchunks")
                if manifest.get("cell_index") != cell_index or manifest.get("chunk") != chunk_index:
                    raise ValueError("Manifest cell/chunk differs from directory")
                count(manifest.get("cell_index"), "manifest.cell_index")
                count(manifest.get("chunk"), "manifest.chunk")
                if manifest.get("config") != config:
                    raise ValueError("Manifest config differs from supplied frozen config")
                identity = {name: manifest[name] for name in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
                if manifest.get("signature") != digest(identity):
                    raise ValueError("Manifest signature fails recomputation")
                hashes = manifest["source_sha256"]
                if not isinstance(hashes, dict) or not set(RUNTIME_FILES) <= set(hashes) or not all(is_hash(value) for value in hashes.values()):
                    raise ValueError("Incomplete or invalid runtime source hashes")
                if expected_sources and any(hashes.get(name) != value for name, value in expected_sources.items()):
                    audit.error("frozen_source_mismatch", "Runtime source differs from frozen source copy", directory / "manifest.json")
                profile = digest(hashes)
                source_profiles[profile].append(report["path"])
                source_values[profile] = hashes
                environment = manifest.get("environment", {})
                if not isinstance(environment, dict) or any(field not in environment for field in ENVIRONMENT_FIELDS):
                    raise ValueError("Manifest numerical environment is incomplete")
                if not isinstance(environment["versions"], dict) or not set(NUMERICAL_PACKAGES) <= set(environment["versions"]) or any(not isinstance(value, str) or not value for value in environment["versions"].values()):
                    raise ValueError("Missing numerical package versions")
                if not isinstance(environment["source_sha256"], dict) or not environment["source_sha256"] or not all(is_hash(value) for value in environment["source_sha256"].values()):
                    raise ValueError("Invalid numerical dependency source hashes")
                if not is_hash(environment["wheel_sha256"]):
                    raise ValueError("Invalid dependency wheel hash")
                numerical = {field: environment[field] for field in ENVIRONMENT_FIELDS}
                profile = digest(numerical)
                environments[profile].append(report["path"])
                environment_values[profile] = numerical
                if environment.get("pyebm_path"):
                    environment_paths.add(environment["pyebm_path"])
                listed = manifest.get("seeds")
                if not isinstance(listed, list) or any(type(seed) is not int for seed in listed):
                    raise ValueError("Manifest seeds must be an integer list")
                if len(set(listed)) != len(listed):
                    raise ValueError("Duplicate seeds inside manifest")
                intended = sorted(seeds)[chunk_index::n_chunks]
                if listed != intended:
                    raise ValueError("Manifest seeds differ from planned strided chunk allocation")
                chunk_seeds = set(listed)
                manifests_by_cell[cell_index].append((chunk_index, n_chunks, chunk_seeds))
                report.update(manifest_signature_valid=True, n_manifest_seeds=len(listed),
                              expected_rows=len(listed) * len(engines), seeds=list(listed), nchunks=n_chunks)
            except (ValueError, TypeError, KeyError) as exc:
                audit.error("invalid_manifest", str(exc), directory / "manifest.json")
                report["manifest_signature_valid"] = False
        count_in_chunk = 0
        rows_path = directory / "rows.jsonl"
        if not rows_path.exists():
            audit.warning("missing_rows_file", "No dataset rows saved yet", rows_path)
        else:
            before = rows_path.stat()
            try:
                with rows_path.open() as handle:
                    for line_number, text in enumerate(handle, 1):
                        if not text.strip():
                            continue
                        count_in_chunk += 1
                        physical_rows += 1
                        location = f"{rows_path}:{line_number}"
                        try:
                            row = json.loads(text)
                            if not isinstance(row, dict):
                                raise ValueError("Expected a dataset object")
                            run_id, cell, engine, seed = (row.get(name) for name in ("run_id", "cell", "engine", "seed"))
                            if type(seed) is not int or run_id != config["run_id"] or cell not in cell_names or engine not in engines or seed not in seeds:
                                raise ValueError("Unplanned row identity")
                            if "phase" in row and row["phase"] != config["phase"]:
                                raise ValueError("Row phase differs from frozen config")
                            if cell_index >= len(cell_names) or cell != cell_names[cell_index]:
                                raise ValueError("Row belongs to a different cell directory")
                            if chunk_seeds is None or seed not in chunk_seeds:
                                raise ValueError("Row seed is not in this verified chunk manifest")
                            status = row.get("status")
                            if not isinstance(status, str) or not status:
                                raise ValueError("Missing observed status")
                            compact = {"status": status, "source": location, "valid_observed": False,
                                       "schemes": {}, "data_sha256": row.get("truth", {}).get("manifest", {}).get("data_sha256")}
                            if status == "ok":
                                try:
                                    validate_observed(row)
                                    compact["valid_observed"] = True
                                except ValueError as exc:
                                    audit.error("invalid_observed", str(exc), location)
                            diagnostics = row.get("diagnostics", {})
                            if "full_p_selected" in diagnostics and diagnostics["full_p_selected"] is not (seed in full_p):
                                audit.error("full_p_selection_mismatch", "Recorded selection differs from frozen seed plan", location)
                            entries = row.get("schemes", {})
                            if not isinstance(entries, dict):
                                raise ValueError("Row schemes must be an object")
                            for name in set(entries) - set(schemes):
                                audit.error("unexpected_scheme", f"Unplanned scheme {name}", location)
                            for name, entry in entries.items():
                                if name not in schemes:
                                    continue
                                try:
                                    compact["schemes"][name] = scheme_state(entry, schemes[name], config, seed in full_p)
                                except (ValueError, TypeError, KeyError) as exc:
                                    audit.error("invalid_scheme", f"{name}: {exc}", location)
                                    compact["schemes"][name] = {"state": "invalid_data", "obligation_satisfied": False, "terminal": False,
                                                                 "nperm": 0, "attempted": 0, "failed": 0}
                            records[(cell, engine, seed)].append(compact)
                        except (ValueError, TypeError, KeyError, AttributeError) as exc:
                            audit.error("invalid_dataset_row", str(exc), location)
            except (OSError, UnicodeError) as exc:
                audit.error("rows_read_error", str(exc), rows_path)
            after = rows_path.stat()
            if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
                audit.stable = False
                audit.warning("rows_changed_during_read", "Repeat inventory after this checkpoint settles", rows_path)
        report["physical_rows"] = count_in_chunk
        if progress:
            report["progress_rows"] = progress.get("rows")
            if progress.get("rows") != count_in_chunk:
                audit.warning("progress_row_count_mismatch", "Progress and saved rows describe different checkpoints", directory)
            if report.get("expected_rows") is not None and progress.get("expected_rows") != report["expected_rows"]:
                audit.warning("progress_plan_mismatch", "Progress expected_rows differs from independently computed chunk plan", directory)
        chunk_reports.append(report)

    coverage = []
    for index, cell in enumerate(cell_names):
        manifests = manifests_by_cell[index]
        union = set().union(*(entry[2] for entry in manifests)) if manifests else set()
        multiplicity = Counter(seed for _, _, listed in manifests for seed in listed)
        counts = {entry[1] for entry in manifests}
        if len(counts) > 1:
            audit.error("inconsistent_chunk_count", "Manifest nchunks differs within a cell", cell)
        inferred = nchunks_override if nchunks_override is not None else next(iter(counts)) if len(counts) == 1 else None
        missing_chunks = sorted(set(range(inferred)) - {entry[0] for entry in manifests}) if inferred is not None else None
        duplicates = sorted(seed for seed, n in multiplicity.items() if n > 1)
        if duplicates:
            audit.error("overlapping_manifest_seeds", f"{len(duplicates)} planned seeds occur in more than one chunk manifest", cell)
        covered = union == seeds and not duplicates and bool(manifests) and missing_chunks == []
        coverage.append({"cell": cell, "planned_datasets": len(seeds), "manifest_seed_union_count": len(union),
                         "missing_manifest_seeds": sorted(seeds - union), "overlapping_manifest_seeds": duplicates,
                         "inferred_nchunks": inferred, "missing_chunk_indices": missing_chunks, "complete": covered})
    if len(source_profiles) > 1:
        audit.error("source_inconsistency", "Chunks contain different runtime source hashes")
    if len(environments) > 1:
        audit.error("environment_inconsistency", "Chunks contain different numerical environments")
    duplicated = {identity: values for identity, values in records.items() if len(values) > 1}
    for identity, values in duplicated.items():
        audit.error("duplicate_row_identity", f"{identity}: {len(values)} records; none silently selected", [v["source"] for v in values])
    for cell in cell_names:
        for seed in sorted(seeds):
            hashes = {value["data_sha256"] for engine in engines for value in records.get((cell, engine, seed), []) if value.get("data_sha256")}
            if len(hashes) > 1:
                audit.error("paired_data_hash_mismatch", "Engines with one cell/seed contain different generated datasets", f"{cell}:{seed}")

    engine_reports, total_status, overall_schemes = [], Counter(), {}
    missing_total = duplicated_total = observed_ok = observed_failures = observed_pending = 0
    terminal_rows = satisfied_rows = 0
    def empty_scheme(expected):
        return {"planned_rows": expected, "states": Counter(), "successful_permutations": 0,
                "attempted_permutations": 0, "failed_permutations": 0,
                "planned_full_budget_permutations": expected * config.get("bperm", 0)}
    overall_schemes = {name: empty_scheme(expected_rows) for name in schemes}
    for cell in cell_names:
        for engine in engines:
            status_counts, missing, duplicates, failed_seeds, pending_seeds = Counter(), [], [], [], []
            local_schemes = {name: empty_scheme(len(seeds)) for name in schemes}
            for seed in sorted(seeds):
                values = records.get((cell, engine, seed), [])
                row_terminal = row_satisfied = False
                if not values:
                    category = "missing_observed_row"
                    missing.append(seed)
                    missing_total += 1
                    states = {name: {"state": category} for name in schemes}
                elif len(values) > 1:
                    category = "duplicate_identity"
                    duplicates.append(seed)
                    duplicated_total += 1
                    states = {name: {"state": category} for name in schemes}
                else:
                    value = values[0]
                    category = value["status"]
                    if category != "ok" or not value["valid_observed"]:
                        category = category if category != "ok" else "invalid_observed"
                        row_terminal = category not in ("pending", "running")
                        if row_terminal:
                            observed_failures += 1
                            failed_seeds.append(seed)
                        else:
                            observed_pending += 1
                            pending_seeds.append(seed)
                        states = {name: {"state": "observed_failed_or_invalid" if row_terminal else "observed_pending"} for name in schemes}
                    else:
                        observed_ok += 1
                        states = {name: value["schemes"].get(name, {"state": "missing_scheme", "terminal": False,
                                                                   "obligation_satisfied": False}) for name in schemes}
                        row_terminal = all(state.get("terminal", False) for state in states.values())
                        row_satisfied = all(state.get("obligation_satisfied", False) for state in states.values())
                status_counts[category] += 1
                total_status[category] += 1
                terminal_rows += int(row_terminal)
                satisfied_rows += int(row_satisfied)
                for name, state in states.items():
                    for target in (local_schemes[name], overall_schemes[name]):
                        target["states"][state["state"]] += 1
                        target["successful_permutations"] += state.get("nperm", 0)
                        target["attempted_permutations"] += state.get("attempted", 0)
                        target["failed_permutations"] += state.get("failed", 0)
            engine_reports.append({"cell": cell, "engine": engine, "planned_rows": len(seeds),
                                   "present_unique_identities": len(seeds) - len(missing),
                                   "missing_seeds": missing, "duplicate_seeds": duplicates, "failed_or_invalid_observed_seeds": failed_seeds,
                                   "pending_observed_seeds": pending_seeds,
                                   "planned_denominator_status_counts": dict(status_counts), "schemes": local_schemes})
    coverage_complete = all(value["complete"] for value in coverage)
    collection_complete = coverage_complete and terminal_rows == expected_rows and not duplicated_total
    complete = coverage_complete and satisfied_rows == expected_rows and not duplicated_total
    errors = [issue for issue in audit.issues if issue["severity"] == "error"]
    ready = complete and not errors and audit.stable and len(source_profiles) == len(environments) == 1
    for chunk in chunk_reports:
        if chunk["progress_done"] is True and (chunk["physical_rows"] != chunk.get("expected_rows")):
            audit.warning("done_without_row_coverage", "done=true does not establish complete planned rows", chunk["path"])
        if chunk["progress_done"] is True and chunk.get("seeds") is not None and chunk["cell_index"] < len(cell_names):
            unsatisfied = 0
            for engine in engines:
                for seed in chunk["seeds"]:
                    values = records.get((cell_names[chunk["cell_index"]], engine, seed), [])
                    met = len(values) == 1 and values[0]["status"] == "ok" and values[0]["valid_observed"]
                    if met:
                        met = all(values[0]["schemes"].get(name, {}).get("obligation_satisfied", False) for name in schemes)
                    unsatisfied += not met
            if unsatisfied:
                audit.warning("done_without_primary_completion", f"done=true but {unsatisfied} planned rows have unmet obligations", chunk["path"])
    return {"schema_version": "run_inventory_v2.1", "inspection_started_utc": started, "inspection_finished_utc": utc(),
            "run_id": config["run_id"], "phase": config["phase"], "config_path": str(config_path),
            "config_file_sha256": hashlib.sha256(config_path.read_bytes()).hexdigest(), "config_semantic_sha256": digest(config),
            "root": str(root), "complete": complete, "collection_complete": collection_complete,
            "ready_for_final_analysis": ready, "terminal_with_failures": collection_complete and not complete,
            "stable_rows_snapshot": audit.stable,
            "plan": {"cells": cell_names, "engines": engines, "datasets_per_cell": len(seeds),
                     "planned_cell_datasets": len(cell_names) * len(seeds), "planned_rows": expected_rows,
                     "base_seed": config["base_seed"], "last_seed": max(seeds), "primary_rule": config.get("rule", "strict"),
                     "alpha": config.get("alpha", .05), "bperm": config.get("bperm", 0),
                     "full_p_seeds": sorted(full_p), "schemes": schemes},
            "totals": {"physical_rows": physical_rows, "present_expected_identities": len(records),
                       "missing_planned_rows": missing_total, "duplicate_planned_identities": duplicated_total,
                       "observed_ok": observed_ok, "observed_failed_or_invalid": observed_failures,
                       "observed_pending": observed_pending,
                       "terminal_planned_rows": terminal_rows, "satisfied_planned_rows": satisfied_rows,
                       "planned_denominator_status_counts": dict(total_status), "schemes": overall_schemes},
            "manifest_coverage": coverage, "cell_engines": engine_reports, "chunks": chunk_reports,
            "provenance": {"expected_frozen_source_root": expected_source_root, "expected_frozen_sources": expected_sources,
                           "source_profiles": [{"profile": profile, "chunks": chunks, "source_sha256": source_values[profile]} for profile, chunks in source_profiles.items()],
                           "numerical_environment_profiles": [{"profile": profile, "chunks": chunks, "environment": environment_values[profile]} for profile, chunks in environments.items()],
                           "pyebm_path_variants": sorted(environment_paths)},
            "issues": audit.issues,
            "notes": ["All counts use the frozen plan, including missing seeds and failed observed fits.",
                      "Progress flags are advisory; dataset identities, counts, primary decisions and full-p obligations are independently checked.",
                      "Chunk signatures are recomputed individually; different chunk signatures are expected.",
                      "Numerical environments compare versions, dependency source hashes and wheel hash, not installation path.",
                      "No fit_records file was opened; per-permutation event identity and log integrity are outside this inventory.",
                      "Readiness is for analysis of the declared phase, not a claim of validity, power, precision or confirmation."]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, type=Path, help="frozen run configuration")
    parser.add_argument("--root", required=True, type=Path, help="run root containing cell_i/chunk_j")
    parser.add_argument("--source-root", type=Path, help="optional frozen scripts/v2 directory; otherwise inferred next to configs/v2")
    parser.add_argument("--nchunks", type=int, help="optional explicit planned chunk count, including not-yet-started chunks")
    parser.add_argument("--output", type=Path, help="save the same JSON inventory atomically")
    parser.add_argument("--quiet", action="store_true", help="suppress stdout JSON; requires --output")
    parser.add_argument("--require-ready", action="store_true", help="exit 3 when not ready")
    args = parser.parse_args(argv)
    if args.quiet and args.output is None:
        parser.error("--quiet requires --output")
    try:
        result = inventory(args.config, args.root, args.source_root, args.nchunks)
        text = json.dumps(result, indent=2, allow_nan=False) + "\n"
        if args.output:
            if args.output.resolve() == args.config.resolve() or args.output.name in {"manifest.json", "progress.json", "rows.jsonl", "fit_records.jsonl", "runner.lock"}:
                raise ValueError("Inventory output must not overwrite configuration or raw run artifacts")
            args.output.parent.mkdir(parents=True, exist_ok=True)
            temporary = args.output.with_name(args.output.name + f".tmp-{os.getpid()}")
            temporary.write_text(text)
            os.replace(temporary, args.output)
        if not args.quiet:
            sys.stdout.write(text)
        return 3 if args.require_ready and not result["ready_for_final_analysis"] else 0
    except (OSError, ValueError, TypeError) as exc:
        parser.exit(2, f"inventory error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
