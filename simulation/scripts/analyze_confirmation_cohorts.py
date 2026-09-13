#!/usr/bin/env python3
"""Analyze explicitly frozen confirmation cohorts; never pool discovered run globs.

Standard library only. Reads protocol, frozen snapshots, manifests, progress and
dataset rows, never fit_records. Unresolved outcomes retain the planned denominator.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import re
import statistics

PROJECT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("confirmation_inventory", PROJECT / "hpc/v2/inventory_run.py")
I = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(I)
GROUPS = ("e2", "e33", "e4")
PAIRS = ((0, 1), (0, 2), (1, 2))
PAIR_NAMES = ("e2_vs_e33", "e2_vs_e4", "e33_vs_e4")
FAMILY = "diagnosis_within_pair_bonferroni3"
GLOBAL_RUNS = ("confirm_power_global_a", "confirm_power_global_k14_k55_a")
PARTIAL_RUNS = ("confirm_pair_partial_e2_topup_a", "confirm_pair_partial_e4_topup_a", "confirm_pair_partial_e33_a")
LOCAL_RUNS = ("confirm_pair_complete_ref_a", "confirm_pair_power_k27_a")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_child(root, name):
    root = Path(root).resolve()
    path = (root / name).resolve()
    require(path.is_relative_to(root) and path != root, f"Unsafe registered path: {name}")
    return path


def literal_assignments(path):
    result = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                result[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return result


def design_defaults(path):
    """Read dataclass defaults without importing NumPy or running scientific code."""
    tree = ast.parse(Path(path).read_text())
    defaults, biomarkers = {}, None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "BIOMARKERS" for t in node.targets):
            biomarkers = [ast.literal_eval(v.args[0]) for v in node.value.elts]
        if isinstance(node, ast.ClassDef) and node.name == "DesignConfig":
            for item in node.body:
                if isinstance(item, ast.AnnAssign):
                    name = item.target.id
                    if name == "biomarker_names":
                        continue
                    defaults[name] = {} if name == "group_pdf_shift" else ast.literal_eval(item.value)
    require(biomarkers is not None and defaults.get("geometry") == "exact_k", "Unsupported frozen DesignConfig AST")
    defaults.update(biomarker_names=biomarkers, design_version=literal_assignments(path)["DESIGN_VERSION"])
    return defaults


def resolved_design(defaults, cell):
    result = dict(defaults)
    name = cell["name"]
    result["name"] = name
    if name == "IID_H0":
        result.update(diagnosis_mode="iid", composition_lambda=0.0)
    elif name == "STAGE_H0":
        result["stage_eta"] = 1.0
    elif name != "REF_H0":
        match = re.fullmatch(r"PWR_(E2|E33|E4)_K(\d+)", name)
        require(match is not None, f"Unsupported confirmation design {name}")
        result.update(alternative_group=match[1].lower(), target_inversions=int(match[2]))
    require(set(cell.get("overrides", {})) <= set(result), "Unknown design override")
    result.update(cell.get("overrides", {}))
    if result["fixed_base_order"] is not None:
        result["fixed_base_order"] = list(result["fixed_base_order"])
    return result


def planned_null(design):
    changed = design["alternative_group"] if design["target_inversions"] else None
    return [changed not in (GROUPS[a], GROUPS[b]) for a, b in PAIRS]


def inversions(a, b):
    require(len(a) == len(b) and sorted(a) == sorted(b) == list(range(len(a))), "Ordering is not a biomarker permutation")
    rank = {v: i for i, v in enumerate(b)}
    return sum(rank[a[i]] > rank[a[j]] for i in range(len(a)) for j in range(i + 1, len(a)))


def validate_truth(row, design, source_hash, environment):
    truth = row["truth"]
    manifest = truth["manifest"]
    require(truth["config"] == design, "Truth resolved design differs from frozen cell")
    require(truth["seed"] == manifest["seed"] == row["seed"], "Truth seed mismatch")
    require(manifest["config_sha256"] == canonical_sha(design), "Truth config hash mismatch")
    require(manifest["source_sha256"] == source_hash, "Truth generator source hash mismatch")
    require(I.is_hash(manifest.get("data_sha256")), "Missing valid generated-data hash")
    require(truth["biomarker_names"] == design["biomarker_names"] and truth["group_order"] == list(GROUPS), "Biomarker/group order mismatch")
    require(manifest["numpy_version"] == environment["versions"]["numpy"] and manifest["pandas_version"] == environment["versions"]["pandas"], "Generator environment mismatch")
    orders = [truth["group_orderings"][g] for g in GROUPS]
    n = len(design["biomarker_names"])
    require(len(truth["base_order"]) == n, "Wrong truth biomarker count")
    real = {g: inversions(truth["base_order"], order) for g, order in zip(GROUPS, orders)}
    expected = {g: design["target_inversions"] if g == design["alternative_group"] else 0 for g in GROUPS}
    require(real == expected == truth["planned_inversions"] == truth["realized_inversions"], "Planned/realized inversion mismatch")
    null = planned_null(design)
    require(truth["pair_null"] == null, "True-null pair mask mismatch")
    for idx, (a, b) in enumerate(PAIRS):
        k = inversions(orders[a], orders[b])
        value = truth["pair_truth"][PAIR_NAMES[idx]]
        require(value["inversions"] == k and value["null"] is null[idx] and value["normalized_distance"] == k / math.comb(n, 2), "Pair truth distance mismatch")
    require(truth["h1_distance"] == design["target_inversions"] / math.comb(n, 2), "Actual H1 distance mismatch")
    if row["status"] == "ok":
        I.validate_observed(row)
        observed = row["observed"]
        require(all(len(order) == n for order in observed["orderings"]), "Observed biomarker count mismatch")
        distances = [inversions(order, target) / math.comb(n, 2) for order, target in zip(observed["orderings"], orders)]
        taus = [inversions(observed["orderings"][a], observed["orderings"][b]) / math.comb(n, 2) for a, b in PAIRS]
        require(observed["distance_to_truth"] == distances and observed["taus"] == taus, "Observed ordering/distance inconsistency")
        provenance = row.get("diagnostics", {}).get("provenance", {})
        require(all(provenance.get(k) == environment[k] for k in I.ENVIRONMENT_FIELDS), "Row numerical environment differs from manifest")
    return null


def frozen_run(entry, project, snapshots_root=None):
    run_id = entry["run_id"]
    require(re.fullmatch(r"[A-Za-z0-9_-]+", run_id) is not None, "Unsafe run ID")
    snapshot = safe_child(snapshots_root, run_id) if snapshots_root else safe_child(project, entry["snapshot"])
    sums = snapshot / "SHA256SUMS"
    require(sha(sums) == entry["snapshot_sha256"], f"{run_id}: snapshot manifest hash mismatch")
    listed = {}
    for line in sums.read_text().splitlines():
        digest, rel = line.split(maxsplit=1)
        rel = rel.lstrip("*")
        require(I.is_hash(digest) and rel not in listed, "Malformed/duplicate SHA256SUMS entry")
        require(sha(safe_child(snapshot, rel)) == digest, f"{run_id}: frozen file mismatch: {rel}")
        listed[rel] = digest
    config_path = safe_child(snapshot, entry["config"])
    require(entry["config"] in listed and sha(config_path) == entry["config_sha256"], "Frozen config hash mismatch")
    cfg = json.loads(config_path.read_text())
    for key, actual in (("run_id", run_id), ("base_seed", entry["base_seed"]), ("datasets", entry["datasets_per_cell"]),
                        ("cells", entry["cells"]), ("engines", entry["engines"]), ("bperm", entry["bperm"]),
                        ("full_p_first_n", entry["full_p_first_n"])):
        require(cfg.get(key, 0 if key in ("bperm", "full_p_first_n") else None) == actual, f"Registry/config mismatch: {key}")
    require(cfg.get("full_p_seeds") == entry["full_p_seeds"], "Registry full-p seed mismatch")
    for key in ("paired_standard", "fast_likelihood"):
        require(cfg.get(key, False) == entry[key], f"Registry {key} mismatch")
    require(cfg["phase"] == "confirmation" and cfg.get("fast_likelihood") is True, "Expected confirmation with validated fast implementation")
    if cfg.get("bperm", 0):
        require(cfg.get("alpha") == .05 and cfg.get("rule") == "le", "Expected frozen alpha=.05, rule=le")
    names, seeds, full_p = I.validate_config(cfg)
    require(all(str(seed).startswith("42") and len(str(seed)) == 8 for seed in seeds), "Seed outside confirmation namespace42")
    sources = {name: sha(snapshot / "scripts/v2" / name) for name in I.RUNTIME_FILES + I.OPTIONAL_RUNTIME_FILES if (snapshot / "scripts/v2" / name).is_file()}
    require(set(I.RUNTIME_FILES) <= set(sources) and "fast_likelihood_v2.py" in sources, "Missing frozen runtime source")
    require(all(f"scripts/v2/{name}" in listed for name in sources), "Unsigned runtime file")
    defaults = design_defaults(snapshot / "scripts/v2/design_v2.py")
    designs = {name: resolved_design(defaults, spec) for name, spec in zip(names, cfg["cells"])}
    dependency = literal_assignments(snapshot / "scripts/v2/dependency_v2.py")
    return {"run_id": run_id, "cfg": cfg, "config_path": config_path, "snapshot": snapshot,
            "sources": sources, "designs": designs, "names": names, "seeds": seeds, "full_p": full_p,
            "dependency": dependency, "records": {}, "issues": [], "environment": None}


def compact_row(row):
    """Keep audit/decision evidence after validation, omit latent arrays and optimizer histories."""
    result = {key: row[key] for key in ("run_id", "phase", "cell", "engine", "seed", "status", "observed", "schemes") if key in row}
    result["truth"] = {key: value for key, value in row.get("truth", {}).items() if key in {
        "config", "seed", "manifest", "biomarker_names", "group_order", "base_order", "group_orderings",
        "planned_inversions", "realized_inversions", "pair_truth", "pair_null", "h1_distance"}}
    result["diagnostics"] = {key: value for key, value in row.get("diagnostics", {}).items() if key in {
        "provenance", "full_p_selected", "sequential_stopping", "engine_version", "mode", "quality_flags"}}
    return result


def load_run(run, results_root, compact=False):
    root = safe_child(results_root, run["run_id"])
    inv = I.inventory(run["config_path"], root, run["snapshot"] / "scripts/v2")
    run["inventory"] = inv
    issues = run["issues"]
    errors = [v for v in inv["issues"] if v["severity"] == "error"]
    run["trusted"] = not errors and inv["stable_rows_snapshot"]
    if errors:
        issues.append({"code": "inventory_integrity_errors", "details": errors})
    profiles = inv["provenance"]["numerical_environment_profiles"]
    if len(profiles) == 1:
        environment = profiles[0]["environment"]
        dep = run["dependency"]
        if (environment["versions"].get("pyebm") != dep["PYEBM_VERSION"] or
                environment["source_sha256"] != dep["SOURCE_SHA256"] or
                environment["wheel_sha256"] != dep["PYEBM_WHEEL_SHA256"]):
            issues.append({"code": "unpinned_dependency_environment"})
            run["trusted"] = False
        run["environment"] = environment
    records = defaultdict(list)
    for path in sorted(root.glob("cell_*/chunk_*/rows.jsonl")):
        before = path.stat()
        for number, text in enumerate(path.read_text().splitlines(), 1):
            if not text.strip():
                continue
            try:
                row = json.loads(text)
                key = (row["cell"], row["engine"], row["seed"])
                require(row["run_id"] == run["run_id"] and key[0] in run["names"] and key[1] in run["cfg"]["engines"] and key[2] in run["seeds"], "Unplanned row")
                reason = None
                try:
                    validate_truth(row, run["designs"][key[0]], run["sources"]["design_v2.py"], run["environment"])
                except (ValueError, KeyError, TypeError, AttributeError) as exc:
                    reason = "invalid_truth_or_observed" if row.get("status") == "ok" else None
                    issues.append({"code": reason or "failed_row_truth_unavailable_or_invalid", "identity": list(key),
                                   "message": str(exc), "path": f"{path}:{number}",
                                   "failure_preserved": row.get("status") != "ok"})
                records[key].append({"row": compact_row(row) if compact else row, "reason": reason, "path": f"{path}:{number}"})
            except (ValueError, KeyError, TypeError) as exc:
                run["trusted"] = False
                issues.append({"code": "invalid_row", "path": f"{path}:{number}", "message": str(exc)})
        after = path.stat()
        if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
            run["trusted"] = False
            issues.append({"code": "rows_changed_during_cohort_read", "path": str(path)})
    run["records"] = dict(records)
    return run


def get_record(run, cell, engine, seed):
    records = run["records"].get((cell, engine, seed), [])
    if not records:
        return None, "missing_row"
    if len(records) != 1:
        return None, "duplicate_identity"
    item = records[0]
    if not run["trusted"]:
        return None, "run_integrity_invalid"
    if item["reason"]:
        return None, item["reason"]
    row = item["row"]
    if row["status"] != "ok":
        return None, "observed_" + row["status"]
    return row, None


def decision(run, cell, engine, seed, scheme, target=None):
    row, reason = get_record(run, cell, engine, seed)
    if row is None:
        return None, reason
    cfg = run["cfg"]
    definition = cfg["schemes"][scheme]
    entry = row.get("schemes", {}).get(scheme)
    if entry is None:
        return None, "missing_scheme"
    try:
        I.scheme_state(entry, definition, cfg, seed in run["full_p"])
        if target is None:
            require(definition.get("require_max") is True and "tested_pair_index" not in definition, "Not a global-max scheme")
            count = entry["maxgt"] + entry["maxeq"]
            comparisons = 1
        else:
            validate_local_definition(definition, target)
            count = entry["gt"][target] + entry["eq"][target]
            comparisons = 3
        value = I.bounds(count, entry["nperm"], cfg["bperm"], cfg["alpha"], comparisons, cfg["rule"])
        return value, None if value is not None else "permutation_decision_unknown"
    except (ValueError, KeyError, TypeError) as exc:
        return None, "invalid_scheme: " + str(exc)


def any3(values):
    return True if any(v is True for v in values) else False if all(v is False for v in values) else None


def all3(values):
    return False if any(v is False for v in values) else True if all(v is True for v in values) else None


def wilson(k, n):
    z = 1.959963984540054
    p = k / n
    center = (p + z * z / (2 * n)) / (1 + z * z / n)
    radius = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return [max(0.0, center - radius), min(1.0, center + radius)]


def summarize(items):
    """items=(planned identity, bool-or-None, reason); never valid-only estimates."""
    n = len(items)
    require(n > 0, "Empty planned endpoint")
    yes = sum(value is True for _, value, _ in items)
    no = sum(value is False for _, value, _ in items)
    unknown = [{"identity": identity, "reason": reason} for identity, value, reason in items if value is None]
    complete = not unknown
    p = yes / n if complete else None
    return {"planned_R": n, "reject_count": yes, "nonreject_count": no, "unknown_count": len(unknown),
            "decisions_complete": complete, "rate_bounds": [yes / n, (yes + len(unknown)) / n],
            "rate": p, "mcse": math.sqrt(p * (1 - p) / n) if complete else None,
            "wilson_ci95": wilson(yes, n) if complete else None,
            "unknown_reason_counts": dict(Counter(item["reason"] for item in unknown)), "unknown": unknown}


def paired_summary(first, second):
    require(len(first) == len(second) and len(first) > 0, "Paired plans have different sizes")
    require([v[0] for v in first] == [v[0] for v in second], "Paired seed identities differ")
    low = high = 0
    differences, unknown = [], []
    for (seed, a, ar), (_, b, br) in zip(first, second):
        av = [0, 1] if a is None else [int(a)]
        bv = [0, 1] if b is None else [int(b)]
        low += min(bv) - max(av)
        high += max(bv) - min(av)
        if a is None or b is None:
            unknown.append({"seed": seed, "first_reason": ar, "second_reason": br})
        else:
            differences.append(int(b) - int(a))
    n = len(first)
    estimate = statistics.mean(differences) if not unknown else None
    mcse = statistics.stdev(differences) / math.sqrt(n) if not unknown and n > 1 else None
    # The bounded-variable interval remains nondegenerate when all observed
    # differences are zero; the empirical standard error alone cannot prove equivalence.
    radius = math.sqrt(2 * math.log(40) / n)
    return {"planned_pairs": n, "resolved_pairs": len(differences), "unknown_pairs": len(unknown),
            "difference_orientation": "second_minus_first", "difference_bounds": [low / n, high / n],
            "difference": estimate, "paired_mcse": mcse,
            "paired_normal_ci95": [max(-1.0, estimate - 1.959963984540054 * mcse), min(1.0, estimate + 1.959963984540054 * mcse)] if mcse is not None and mcse > 0 else None,
            "paired_hoeffding_ci95": [max(-1.0, estimate - radius), min(1.0, estimate + radius)] if estimate is not None else None,
            "degenerate_empirical_difference": mcse == 0 if mcse is not None else None,
            "interval_scope": "pointwise, not simultaneous across panel; Hoeffding is conservative for independent seed units",
            "unknown": unknown}


def compare_sources(runs):
    reports = []
    for target_id in LOCAL_RUNS:
        if target_id not in runs:
            continue
        target = runs[target_id]
        reuse = target["cfg"]["data_and_fit_reuse"]
        source_id = reuse["source_run_id"]
        require(source_id in runs, f"Unregistered source {source_id}")
        source = runs[source_id]
        engine = reuse["source_engine"] if "source_engine" in reuse else reuse["engine"]
        require(target["seeds"] == set(range(reuse["source_seed_start"], reuse["source_seed_end_inclusive"] + 1))
                and target["seeds"] <= source["seeds"], "Source reuse seed plan mismatch")
        counts, unresolved, mismatches = Counter(), [], []
        for cell in reuse["source_cells"]:
            require(source["designs"][cell] == target["designs"][cell], "Source/local frozen designs differ")
            for seed in sorted(target["seeds"]):
                key = (cell, engine, seed)
                aa, bb = source["records"].get(key, []), target["records"].get(key, [])
                identity = {"cell": cell, "engine": engine, "seed": seed}
                if len(aa) != 1 or len(bb) != 1:
                    counts["missing_or_duplicate"] += 1
                    unresolved.append(identity)
                    continue
                a, b = aa[0]["row"], bb[0]["row"]
                checks, missing_checks = {}, []
                for label, av, bv in (
                    ("resolved_design", a.get("truth", {}).get("config"), b.get("truth", {}).get("config")),
                    ("data_hash", a.get("truth", {}).get("manifest", {}).get("data_sha256"), b.get("truth", {}).get("manifest", {}).get("data_sha256")),
                    ("biomarkers", a.get("truth", {}).get("biomarker_names"), b.get("truth", {}).get("biomarker_names"))):
                    if av is None or bv is None:
                        missing_checks.append(label)
                    else:
                        checks[label] = av == bv
                if missing_checks:
                    counts["unavailable_truth_checks"] += 1
                    unresolved.append(dict(identity, reason="truth_unavailable", checks=missing_checks))
                if a["status"] == b["status"] == "ok":
                    for field in ("orderings", "taus", "distance_to_truth"):
                        checks["observed_" + field] = a.get("observed", {}).get(field) == b.get("observed", {}).get(field)
                    counts["both_observed_ok"] += 1
                else:
                    counts["observed_status_" + a["status"] + "__" + b["status"]] += 1
                if not all(checks.values()):
                    mismatches.append(dict(identity, failed_checks=[k for k, v in checks.items() if not v]))
                    aa[0]["reason"] = bb[0]["reason"] = "source_local_crosscheck_mismatch"
                else:
                    counts["matching_available_checks"] += 1
        planned = len(reuse["source_cells"]) * len(target["seeds"])
        reports.append({"source_run_id": source_id, "local_run_id": target_id, "planned_joins": planned,
                        "counts": dict(counts), "missing_or_duplicate": unresolved, "mismatches": mismatches,
                        "all_planned_observed_exact": counts["both_observed_ok"] == planned and not mismatches and not unresolved,
                        "complete_crosscheck": not unresolved and not mismatches and source["trusted"] and target["trusted"],
                        "source_failures_preserved": True, "independent_replicates_added_by_join": 0})
    return reports


def scientific_definition(definition):
    return {key: definition.get(key) for key in ("stratify", "groups", "tested_pair_index", "require_max")}


def validate_local_definition(definition, index):
    require(type(index) is int and index in range(3), "Invalid local target index")
    require(definition.get("tested_pair_index") == index and definition.get("groups") == list(PAIRS[index])
            and definition.get("stratify") == "diagnosis" and definition.get("require_max") is False,
            "Local definition must match the designated diagnosis-conditioned pair and exclude global max")


def partial_cohorts(runs):
    result = []
    for owner_id in PARTIAL_RUNS:
        if owner_id not in runs:
            continue
        cohort = runs[owner_id]["cfg"]["cohort_definition"]
        cell, engine, scheme, index = (cohort[k] for k in ("cell", "engine", "target_scheme", "target_pair_index"))
        seen, items, profiles = set(), [], []
        for segment in cohort["segments"]:
            run = runs[segment["run_id"]]
            seeds = set(range(segment["seed_start"], segment["seed_end_inclusive"] + 1))
            require(len(seeds) == segment["R"] and seeds <= run["seeds"] and not seen.intersection(seeds), "Invalid/overlapping cohort seed segment")
            seen.update(seeds)
            require(engine in run["cfg"]["engines"] and planned_null(run["designs"][cell])[index], "Cohort target is not a planned true null")
            definition = run["cfg"]["schemes"][scheme]
            validate_local_definition(definition, index)
            profiles.append({"design": run["designs"][cell], "definition": scientific_definition(definition),
                             "B": run["cfg"]["bperm"], "alpha": run["cfg"]["alpha"], "rule": run["cfg"]["rule"], "sources": run["sources"]})
            for seed in sorted(seeds):
                value, reason = decision(run, cell, engine, seed, scheme, index)
                items.append(({"run_id": run["run_id"], "cell": cell, "engine": engine, "seed": seed}, value, reason))
        require(len(items) == cohort["combined_target_R"] and all(p == profiles[0] for p in profiles), "Cohort plan size/scientific method mismatch")
        result.append({"cohort_id": cohort["cohort_id"], "cell": cell, "engine": engine, "target_scheme": scheme,
                       "target_pair_index": index, "target_pair": PAIR_NAMES[index], "null_bound": .05 / 3,
                       "segments": cohort["segments"], "endpoint": "true_null_local_false_rejection", "summary": summarize(items)})
    return result


def global_panel(runs):
    panel, vectors = [], {}
    for run_id in GLOBAL_RUNS:
        if run_id not in runs:
            continue
        run = runs[run_id]
        require(set(run["cfg"]["schemes"]) == {"diagnosis"}, "Global power run includes unexpected schemes")
        require(run["cfg"]["schemes"]["diagnosis"] == {"stratify": "diagnosis", "require_max": True}, "Global power requires the diagnosis-conditioned global-max scheme")
        require(run["cfg"]["engines"] == ["repaired"] and len(run["seeds"]) == 500, "Global power plan is not repaired R500")
        for cell in run["names"]:
            match = re.fullmatch(r"PWR_(E2|E4)_K(14|27|55)", cell)
            require(match is not None and cell not in vectors, "Unexpected/duplicate global power cell")
            items = [(seed, *decision(run, cell, "repaired", seed, "diagnosis")) for seed in sorted(run["seeds"])]
            vectors[cell] = items
            panel.append({"run_id": run_id, "cell": cell, "engine": "repaired", "scheme": "diagnosis", "endpoint": "global_max_rejection",
                          "alternative_group": match[1].lower(), "K": int(match[2]), "actual_distance": int(match[2]) / 91,
                          "summary": summarize(items)})
    comparisons = []
    for a, b in itertools.combinations(sorted(vectors), 2):
        same_direction = a.split("_")[1] == b.split("_")[1]
        same_k = a.split("_")[2] == b.split("_")[2]
        if same_direction or same_k:
            comparisons.append({"first_cell": a, "second_cell": b, "pairing": "common_random_seed_different_cells",
                                "summary": paired_summary(vectors[a], vectors[b])})
    return panel, comparisons


def local_families(runs):
    result = []
    for run_id in LOCAL_RUNS:
        if run_id not in runs:
            continue
        run = runs[run_id]
        family = run["cfg"]["family_definition"]
        require(family["id"] == FAMILY and family["level"] == .05 and family["pair_order"] == list(PAIR_NAMES), "Unexpected local family")
        members = {definition["tested_pair_index"]: name for name, definition in run["cfg"]["schemes"].items() if definition.get("family") == FAMILY}
        require(set(members) == {0, 1, 2} and len(run["cfg"]["schemes"]) == 3, "Incomplete/duplicated local family")
        require(set(family.get("member_schemes", members.values())) == set(members.values()), "Family member declaration mismatch")
        for index, name in members.items():
            validate_local_definition(run["cfg"]["schemes"][name], index)
        for cell in run["names"]:
            null = planned_null(run["designs"][cell])
            endpoints = {"any_pair_rejection": [], "any_true_null_rejection": []}
            if not all(null):
                require(len(run["seeds"]) == 500 and sum(not v for v in null) == 2, "Alternative family must be first500 with two true alternatives")
                endpoints.update(any_true_alternative_discovery=[], both_true_alternatives_discovery=[])
            for seed in sorted(run["seeds"]):
                outcomes = [decision(run, cell, "repaired", seed, members[index], index) for index in range(3)]
                values = [v for v, _ in outcomes]
                true_alt = [v for v, is_null in zip(values, null) if not is_null]
                endpoint_values = {"any_pair_rejection": any3(values), "any_true_null_rejection": any3([v for v, is_null in zip(values, null) if is_null])}
                if true_alt:
                    endpoint_values.update(any_true_alternative_discovery=any3(true_alt), both_true_alternatives_discovery=all3(true_alt))
                reason = "; ".join(sorted({r for v, r in outcomes if v is None and r})) or "unknown_family_decision"
                for endpoint, value in endpoint_values.items():
                    endpoints[endpoint].append((seed, value, None if value is not None else reason))
            result.append({"run_id": run_id, "cell": cell, "engine": "repaired", "family": FAMILY, "pair_null": null,
                           "target_scope": "original_complete_family_only_no_topups", "endpoints": {name: summarize(items) for name, items in endpoints.items()}})
    return result


def analyze(protocol_path, results_root, project=PROJECT, snapshots_root=None):
    protocol_path, project = Path(protocol_path).resolve(), Path(project).resolve()
    protocol = json.loads(protocol_path.read_text())
    require(protocol.get("configuration_frozen_before_confirmation_data") is True and protocol.get("confirmation_namespace") == "42xxxxxx", "Protocol is not frozen namespace42 confirmation")
    entries = protocol["runs"]
    ids = [entry["run_id"] for entry in entries]
    require(len(ids) == len(set(ids)), "Duplicate registered run ID")
    require(set(GLOBAL_RUNS + PARTIAL_RUNS + LOCAL_RUNS) <= set(ids), "Protocol lacks required explicit cohort runs")
    runs = {entry["run_id"]: load_run(frozen_run(entry, project, snapshots_root), results_root, compact=True) for entry in entries}
    environments = {canonical_sha(run["environment"]) for run in runs.values() if run["environment"] is not None}
    sources = {canonical_sha(run["sources"]) for run in runs.values()}
    cross_run_issues = []
    if len(environments) > 1 or len(sources) > 1:
        cross_run_issues.append({"code": "cross_run_environment_or_source_mismatch", "environment_profiles": len(environments), "source_profiles": len(sources)})
        for run in runs.values():
            run["trusted"] = False
    joins = compare_sources(runs)
    partial = partial_cohorts(runs)
    panel, paired = global_panel(runs)
    require({item["cell"] for item in panel} == {f"PWR_{g}_K{k}" for g in ("E2", "E4") for k in (14, 27, 55)}, "Global panel does not contain exactly six frozen points")
    families = local_families(runs)
    reports = []
    for run in runs.values():
        inv = run["inventory"]
        reports.append({"run_id": run["run_id"], "frozen_config": str(run["config_path"]), "config_sha256": sha(run["config_path"]),
                        "runtime_source_sha256": run["sources"], "numerical_environment": run["environment"],
                        "fast_likelihood": run["cfg"]["fast_likelihood"], "paired_standard": run["cfg"].get("paired_standard", False),
                        "inventory_ready": inv["ready_for_final_analysis"], "collection_complete": inv["collection_complete"],
                        "trusted": run["trusted"], "issues": run["issues"], "plan": inv["plan"], "totals": inv["totals"],
                        "planned_missing_and_failed": inv["cell_engines"], "manifest_coverage": inv["manifest_coverage"]})
    ready = all(r["inventory_ready"] and r["trusted"] and not r["issues"] for r in reports) and not cross_run_issues and all(j["complete_crosscheck"] for j in joins)
    collection_complete = all(r["collection_complete"] for r in reports)
    bounded_ready = (collection_complete and not cross_run_issues
                     and all(r["trusted"] and all(issue["code"] == "failed_row_truth_unavailable_or_invalid" for issue in r["issues"]) for r in reports)
                     and all(not j["mismatches"] and all(v.get("reason") == "truth_unavailable" for v in j["missing_or_duplicate"]) for j in joins))
    endpoint_summaries = [item["summary"] for item in partial + panel] + [summary for item in families for summary in item["endpoints"].values()]
    return {"schema_version": "confirmation_cohorts_v2.1", "inspected_utc": datetime.now(timezone.utc).isoformat(),
            "protocol": str(protocol_path), "protocol_sha256": sha(protocol_path), "results_root": str(Path(results_root).resolve()),
            "ready_for_final_analysis": ready, "status": "complete" if ready else "terminal_with_unknown" if bounded_ready else "incomplete",
            "collection_complete": collection_complete, "ready_for_bounded_final_analysis": bounded_ready,
            "all_reported_endpoint_decisions_complete": all(s["decisions_complete"] for s in endpoint_summaries),
            "runs": reports, "cross_run_issues": cross_run_issues, "source_local_crosschecks": joins,
            "partial_null_cohorts": partial, "global_power_panel": panel, "global_power_paired_comparisons": paired,
            "local_family_endpoints": families,
            "notes": ["Only explicitly registered frozen runs and declared seed segments contribute. Other directories are ignored.",
                      "All endpoint denominators are planned R. Unknowns retain bounds; rate, Wilson CI and MCSE are withheld until that endpoint has no unknowns.",
                      "A resolved endpoint can precede run completion; full-p obligations and provenance still govern final readiness.",
                      "Source/local repeats verify the same dataset and observed fit; they add zero independent replicates and never replace failures.",
                      "The three partial-null cohorts concern one local alpha/3 target each. Topups do not contribute complete-family power.",
                      "Paired contrasts use common seed as the Monte Carlo unit, second minus first; differing cells are not identical datasets.",
                      "All confidence intervals are pointwise, not simultaneous across endpoints or contrasts. Zero empirical paired MCSE does not establish equivalence; a conservative Hoeffding interval remains nondegenerate.",
                      "Terminal-with-unknown means every planned record/work obligation is accounted for but failures prevent some exact endpoints. Bounded final analysis can be ready while point-estimate/final-readiness remains false.",
                      "Fast implementation is verified through frozen config and runtime source hashes; rows do not expose every low-level fit call.",
                      "No fit_records were read. This audit cannot independently recompute raw data hashes or inspect every permutation fit.",
                      "Any run-level inventory integrity error conservatively makes that run's endpoint decisions unknown until the error is repaired; no record is silently selected."]}


def markdown(result):
    lines = ["# Frozen confirmation cohort analysis", "", f"Status: **{result['status']}**. Final-analysis ready: **{result['ready_for_final_analysis']}**.", "",
             "All rates use the full planned denominator. Intervals are 95% Wilson intervals when every endpoint decision is known; unfinished endpoints show identification bounds only.", "",
             "| Endpoint | Planned R | Reject | Unknown | Rate or bounds |", "|---|---:|---:|---:|---|"]
    def add(name, summary):
        rate = f"{summary['rate']:.4f} (CI {summary['wilson_ci95'][0]:.4f}–{summary['wilson_ci95'][1]:.4f}; MCSE {summary['mcse']:.4f})" if summary["rate"] is not None else f"[{summary['rate_bounds'][0]:.4f}, {summary['rate_bounds'][1]:.4f}]"
        lines.append(f"| {name} | {summary['planned_R']} | {summary['reject_count']} | {summary['unknown_count']} | {rate} |")
    for item in result["partial_null_cohorts"]:
        add(item["cohort_id"] + " / true-null local error", item["summary"])
    for item in result["global_power_panel"]:
        add(item["cell"] + " / global max", item["summary"])
    for item in result["local_family_endpoints"]:
        for endpoint, summary in item["endpoints"].items():
            add(item["cell"] + " / " + endpoint, summary)
    lines.extend(["", "## Paired global-power contrasts", "", "Differences are second minus first, paired by common seed. Confidence intervals are pointwise; the Hoeffding interval is conservative and remains nondegenerate when every sampled difference is zero.", "",
                  "| First cell | Second cell | Planned pairs | Unknown pairs | Difference or bounds |", "|---|---|---:|---:|---|"])
    for comparison in result["global_power_paired_comparisons"]:
        value = comparison["summary"]
        rendered = (f"{value['difference']:.4f} (Hoeffding CI {value['paired_hoeffding_ci95'][0]:.4f}–{value['paired_hoeffding_ci95'][1]:.4f})"
                    if value["difference"] is not None else f"[{value['difference_bounds'][0]:.4f}, {value['difference_bounds'][1]:.4f}]")
        lines.append(f"| {comparison['first_cell']} | {comparison['second_cell']} | {value['planned_pairs']} | {value['unknown_pairs']} | {rendered} |")
    lines.extend(["", "## Source/local duplicate checks", ""])
    for join in result["source_local_crosschecks"]:
        lines.append(f"- {join['source_run_id']} → {join['local_run_id']}: planned {join['planned_joins']}; both observed successful {join['counts'].get('both_observed_ok', 0)}; missing/duplicate {len(join['missing_or_duplicate'])}; mismatches {len(join['mismatches'])}. Source failures remain source failures.")
    lines.extend(["", "## Run obligations", "", "| Run | Planned rows | Missing rows | Failed observed | Ready |", "|---|---:|---:|---:|---|"])
    for run in result["runs"]:
        total = run["totals"]
        lines.append(f"| {run['run_id']} | {run['plan']['planned_rows']} | {total['missing_planned_rows']} | {total['observed_failed_or_invalid']} | {run['inventory_ready'] and run['trusted'] and not run['issues']} |")
    lines.extend(["", "## Interpretation limits", ""] + ["- " + note for note in result["notes"]])
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=PROJECT)
    parser.add_argument("--snapshots-root", type=Path, help="Relocated frozen snapshots parent; children must equal registered run IDs")
    parser.add_argument("--output-prefix", type=Path, required=True, help="Write PREFIX.json and PREFIX.md")
    parser.add_argument("--require-ready", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = analyze(args.protocol, args.results_root, args.project_root, args.snapshots_root)
        prefix = args.output_prefix.resolve()
        require(not prefix.is_relative_to(args.results_root.resolve()), "Write analysis outside raw results root")
        paths = [Path(str(prefix) + suffix) for suffix in (".json", ".md")]
        require(all(path != args.protocol.resolve() for path in paths), "Output would overwrite protocol")
        protected = [(args.project_root / "configs").resolve()] + [Path(run["frozen_config"]).parents[2] for run in result["runs"]]
        require(not any(path.is_relative_to(directory) for path in paths for directory in protected), "Output would modify frozen snapshot/configuration")
        for path, text in zip(paths, (json.dumps(result, indent=2, allow_nan=False) + "\n", markdown(result))):
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
            temporary.write_text(text)
            os.replace(temporary, path)
        print(json.dumps({"status": result["status"], "ready_for_final_analysis": result["ready_for_final_analysis"], "outputs": [str(p) for p in paths]}))
        return 3 if args.require_ready and not result["ready_for_final_analysis"] else 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(2, f"confirmation cohort error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
