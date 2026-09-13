"""Read-only six-arm confirmation mechanism analysis, JSON/Markdown only.

The public entry point accepts the frozen protocol and an explicitly supplied raw
run directory. No simulation, fitting, HPC access, plots or PDF generation occurs.
All cohort sizes come from the protocol-matched config, never development sizes.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import random
import re
import statistics

import run_record_audit as U

ENGINES = ("original", "repaired", "oracle_equal_original", "oracle_equal_repaired",
           "oracle_estimated_original", "oracle_estimated_repaired")
CONTRASTS = (("original", "repaired"), ("repaired", "oracle_estimated_repaired"),
             ("oracle_equal_original", "oracle_equal_repaired"),
             ("oracle_estimated_original", "oracle_estimated_repaired"))
PAIR_NAMES = ("e2_vs_e33", "e2_vs_e4", "e33_vs_e4")
BOOTSTRAP_SEED = 20260908
LIMITATIONS = [
    "Intervals are pointwise95% dataset-bootstrap intervals for four prespecified within-cell contrasts; there is no simultaneous/multiplicity-adjusted coverage claim.",
    "The primary endpoint is each dataset's mean of three normalized truth distances; groups, biomarkers and optimizer calls are not independent replicates.",
    "REF N and N4 use the same seed as common-random-number coupling. They are not independent or nested samples; no cross-N formal test is performed.",
    "Known-component oracles are mechanism diagnostics, not deployable inference arms or performance bounds. Estimated q uses observed X; equal-prior and estimated-prior initialization differ.",
    "A matching initial true ordering is a descriptive event, not evidence proving or disproving leakage. No-leakage declarations are checked against the frozen implementation design.",
    "An adjacent-neighbor optimum is not a global optimum or a guarantee of improved truth recovery. Recorded optimizer history contains audit fields, not the full numerical objective trajectory.",
    "Quality flags do not silently exclude legal successful orderings. Missing, failed or invalid endpoints remain in the planned denominator, with explicit unknowns and identification bounds.",
    "Bounds for incomplete outcomes are identification bounds, not confidence intervals. A constant observed bootstrap sample does not prove equivalence.",
    "Observed truth error or between-group distance alone establishes neither permutation calibration nor FWER/power; this mechanism run has no permutation tests."
]


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _finite(x):
    return type(x) in (int, float) and math.isfinite(x)


def _legal(order, events):
    return (isinstance(order, list) and len(order) == events and
            all(type(x) is int for x in order) and sorted(order) == list(range(events)))


def _bounded_flag(values):
    yes, no = sum(x is True for x in values), sum(x is False for x in values)
    n = len(values)
    return {"planned_datasets": n, "true": yes, "false": no, "unknown": n - yes - no,
            "all_planned_fraction_bounds": [yes / n, (n - no) / n] if n else None}


def _json_safe(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    return value


def _numbers(values):
    return U.numbers(x for x in values if _finite(x))


def _bounded_mean(values):
    known = [value for value in values if value is not None]
    total = math.fsum(known)
    unknown = len(values) - len(known)
    # Compute the integer missing count first. Adding then subtracting R can
    # otherwise invert a complete-cohort point interval by floating cancellation.
    return [total / len(values), (total + unknown) / len(values)]


def _integer_grid_directions(first, second, events, index):
    """Count second-first directions on the exact Kendall integer grid.

    Each present vector contains three normalized distances and their saved
    dataset mean. Validate each component against k / choose(events,2), without
    an epsilon; the mean's integer numerator is the sum of the three k values.
    Original floating endpoints remain unchanged for means, CIs and bounds.
    """
    if type(events) is not int or events < 2 or index not in (0, 1, 2, 3) or len(first) != len(second):
        raise ValueError("Invalid Kendall-grid direction inputs")
    denominator = math.comb(events, 2)

    def numerator(values):
        if values is None:
            return None
        if len(values) != 4 or not all(_finite(value) and 0 <= value <= 1 for value in values):
            raise ValueError("Direction values must contain three finite distances and their mean")
        integers = [round(value * denominator) for value in values[:3]]
        if any(value != count / denominator for value, count in zip(values[:3], integers)):
            raise ValueError("Direction distance is outside the verified Kendall grid")
        if values[3] != statistics.mean(values[:3]):
            raise ValueError("Direction dataset mean disagrees with its three components")
        return sum(integers) if index == 3 else integers[index]

    pairs = [(numerator(a), numerator(b)) for a, b in zip(first, second)]
    known = [b - a for a, b in pairs if a is not None and b is not None]
    return [sum(value < 0 for value in known), sum(value == 0 for value in known), sum(value > 0 for value in known)]


def _issue(code, cell=None, seed=None, engine=None, **extra):
    return {"code": code, **({"cell": cell} if cell is not None else {}),
            **({"seed": seed} if seed is not None else {}),
            **({"engine": engine} if engine is not None else {}), **extra}



def _literal_assignments(path):
    values = {}
    for node in ast.parse(Path(path).read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                values[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, TypeError):
                pass
    return values


def _design_plans(config_path, cfg):
    """Read the frozen four-cell design without importing/executing the generator."""
    source = Path(config_path).resolve().parents[2] / "scripts/v2/design_v2.py"
    tree = ast.parse(source.read_text())
    defaults, names = {}, None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "BIOMARKERS" for t in node.targets):
            names = [ast.literal_eval(item.args[0]) for item in node.value.elts]
        if isinstance(node, ast.ClassDef) and node.name == "DesignConfig":
            for item in node.body:
                if isinstance(item, ast.AnnAssign):
                    name = item.target.id
                    if name == "biomarker_names":
                        continue
                    defaults[name] = {} if name == "group_pdf_shift" else ast.literal_eval(item.value)
    version = _literal_assignments(source).get("DESIGN_VERSION")
    if names is None or not version or defaults.get("target_inversions") != 0:
        raise ValueError("Unsupported frozen DesignConfig AST")
    defaults.update(biomarker_names=names, design_version=version)
    plans = {}
    for spec in cfg["cells"]:
        values = dict(defaults)
        values["name"] = spec["name"]
        if spec["name"] == "IID_H0":
            values.update(diagnosis_mode="iid", composition_lambda=0.)
        elif spec["name"] == "STAGE_H0":
            values["stage_eta"] = 1.
        elif spec["name"] != "REF_H0":
            raise ValueError("Only the frozen common-order mechanism designs are supported")
        if not set(spec.get("overrides", {})) <= values.keys():
            raise ValueError("Unknown frozen design override")
        values.update(spec.get("overrides", {}))
        if values.get("fixed_base_order") is not None:
            values["fixed_base_order"] = list(values["fixed_base_order"])
        if values.get("target_inversions") != 0 or values.get("alternative_group") is not None:
            raise ValueError("Mechanism plan must have three common true orderings")
        plans[U.cell_id(spec)] = values
    return plans, _sha(source)


def _validate_saved_truth(record, design, generator_hash, environment):
    if record is None:
        return
    truth = record.get("truth", {})
    if record.get("status") != "ok" and not truth.get("manifest", {}).get("data_sha256"):
        return  # an observed failure before generation is retained, not repaired
    if not environment:
        raise ValueError("No verified chunk environment for saved truth")
    manifest = truth["manifest"]
    if truth["config"] != design or truth.get("design_version") != design["design_version"]:
        raise ValueError("Truth DesignConfig/version differs from frozen cell")
    expected_hash = hashlib.sha256(json.dumps(design, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    if manifest.get("config_sha256") != expected_hash or manifest.get("source_sha256") != generator_hash:
        raise ValueError("Truth config/generator hash mismatch")
    if truth.get("seed") != record["seed"] or manifest.get("seed") != record["seed"]:
        raise ValueError("Saved truth seed differs from planned record")
    if not U.I.is_hash(manifest.get("data_sha256")):
        raise ValueError("Missing valid generated-data SHA")
    if (manifest.get("numpy_version") != environment["versions"]["numpy"] or
            manifest.get("pandas_version") != environment["versions"]["pandas"]):
        raise ValueError("Generator NumPy/Pandas differs from chunk environment")
    if record.get("status") == "ok" and any(record.get("diagnostics", {}).get("provenance", {}).get(k) != environment[k]
                                            for k in U.I.ENVIRONMENT_FIELDS):
        raise ValueError("Fit numerical environment differs from verified chunk")
    n = len(design["biomarker_names"])
    if truth["biomarker_names"] != design["biomarker_names"] or truth["group_order"] != list(U.GROUPS):
        raise ValueError("Truth biomarker/group coordinate mismatch")
    base = truth["base_order"]
    if not _legal(base, n) or any(truth["group_orderings"].get(g) != base for g in U.GROUPS):
        raise ValueError("Truth is not the planned common ordering with the frozen event count")
    zeros = dict.fromkeys(U.GROUPS, 0)
    if truth["planned_inversions"] != zeros or truth["realized_inversions"] != zeros:
        raise ValueError("Truth planned/realized inversions differ from mechanism K0")
    if (truth.get("common_order_null") is not True or truth.get("pair_null") != [True] * 3 or
            truth.get("pair_distances") != [0.] * 3 or truth.get("h1_distance") != 0.):
        raise ValueError("Truth null mask/distance differs from three planned null pairs")
    for pair in PAIR_NAMES:
        item = truth["pair_truth"][pair]
        if item.get("null") is not True or item.get("inversions") != 0 or item.get("normalized_distance") != 0.:
            raise ValueError("Pair truth differs from planned K0")


def _validate_record(record, cell, seed, engine):
    """Validate endpoint arithmetic without refitting. Preserve quality warnings."""
    result = {"cell": cell, "seed": seed, "engine": engine,
              "present": record is not None, "status": "missing" if record is None else record.get("status"),
              "endpoint_legal": False, "values": None, "issues": []}
    if record is None:
        return result
    d, truth = record.get("diagnostics", {}), record.get("truth", {})
    result.update(quality_flags=d.get("quality_flags", []), parameters_finite=d.get("parameters_finite"),
                  event_centers_finite=d.get("event_centers_finite"),
                  optimizer_calls=len(d.get("optimizer", [])),
                  optimizer_unsuccessful=sum(x.get("success") is False for x in d.get("optimizer", [])),
                  optimizer_nonfinite=sum(x.get("finite_x") is False or x.get("finite_objective") is False
                                          for x in d.get("optimizer", [])),
                  consensus=d.get("consensus", []))
    if record.get("kind") != "observed" or record.get("scheme") is not None or record.get("perm_id") is not None:
        result["issues"].append(_issue("nonobserved_mechanism_record", cell, seed, engine))
        return result
    if record.get("status") != "ok":
        return result
    try:
        names = truth["biomarker_names"]
        if len(names) < 2 or len(names) != len(set(names)):
            raise ValueError("invalid biomarker names")
        events = len(names)
        orders = record["orderings"]
        if not isinstance(orders, list) or len(orders) != 3 or not all(_legal(x, events) for x in orders):
            raise ValueError("illegal fitted orderings")
        true_orders = [truth["group_orderings"][g] for g in U.GROUPS]
        if not all(_legal(x, events) for x in true_orders):
            raise ValueError("illegal truth orderings")
        distance = [U.kendall(a, b) for a, b in zip(orders, true_orders)]
        taus = [U.kendall(orders[a], orders[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
        for field, actual in (("distance_to_truth", distance), ("taus", taus)):
            saved = record.get(field)
            if not isinstance(saved, list) or len(saved) != 3 or any(
                    not _finite(x) or not 0 <= x <= 1 or abs(x - y) > 1e-12 for x, y in zip(saved, actual)):
                raise ValueError(f"saved {field} disagrees with legal-order recomputation")
        result["values"] = {"distance_to_truth": distance + [statistics.mean(distance)],
                            "taus": taus + [statistics.mean(taus)]}
        result["endpoint_legal"] = True
        if not isinstance(result["quality_flags"], list):
            raise ValueError("quality_flags must be a list")
        # Missing/nonfinite audit fields are reported, not used for success filtering.
        numeric_fields = ("parameters_finite", "event_centers_finite") if engine in ("original", "repaired") else ("event_centers_finite",)
        for field in numeric_fields:
            if type(d.get(field)) is not bool:
                result["issues"].append(_issue("missing_finite_diagnostic", cell, seed, engine, field=field))
        if len(result["consensus"]) != 3:
            result["issues"].append(_issue("incomplete_consensus_audit", cell, seed, engine,
                                           recorded=len(result["consensus"]), expected=3))
        if engine in ("original", "repaired"):
            optimizers = d.get("optimizer")
            if not isinstance(optimizers, list) or not optimizers or d.get("optimizer_calls") != len(optimizers):
                result["issues"].append(_issue("incomplete_optimizer_history", cell, seed, engine))
            elif d.get("optimizer_failures") != sum(x.get("success") is False for x in optimizers):
                result["issues"].append(_issue("optimizer_failure_count_mismatch", cell, seed, engine))
        for item in result["consensus"]:
            if item.get("local_optimal") not in (True, False) or not _finite(item.get("final_score")):
                result["issues"].append(_issue("incomplete_consensus_optimality_audit", cell, seed, engine))
                break
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        result["endpoint_legal"] = False
        result["values"] = None
        result["issues"].append(_issue("invalid_observed_endpoint", cell, seed, engine, message=str(exc)))
    return result


def _oracle_diagnostics(record, cell, seed, engine):
    if record is None or not engine.startswith("oracle_"):
        return None, []
    d, truth = record.get("diagnostics", {}), record.get("truth", {})
    out = {"cell": cell, "seed": seed, "engine": engine, "status": record.get("status"),
           "initial_truth_distances": {}, "initial_matches_truth": {}, "q": {},
           "q_fit_status_counts": {}, "n_boundary": 0, "n_unidentified": 0,
           "n_prior_fits": 0, "expected_prior_fits": None, "provenance_declarations": {}}
    issues = []
    flags = ["inference_eligible", "truth_order_used_in_prior_and_initialization", "individual_latent_stage_used"]
    if engine.startswith("oracle_estimated"):
        flags += ["diagnosis_used_in_prior_estimation", "truth_dx_counts_used"]
    out["provenance_declarations"] = {flag: d.get(flag) for flag in flags}
    if record.get("status") != "ok":
        return out, issues  # partial diagnostics remain unknown; do not demand completed fits
    try:
        if any(d.get(flag) is not False for flag in flags):
            raise ValueError("missing/invalid oracle no-deployability/no-trueS declarations")
        truth_names = truth["biomarker_names"]
        posterior_names = d["posterior_biomarkers"]
        if set(posterior_names) != set(truth_names) or len(posterior_names) != len(truth_names):
            raise ValueError("oracle biomarker coordinate mismatch")
        p = len(truth_names)
        status_counts = Counter()
        for group in U.GROUPS:
            initial = d["initial_orderings"][group]
            if not _legal(initial, p):
                raise ValueError("illegal oracle initialization")
            mapped = [truth_names.index(posterior_names[i]) for i in initial]
            distance = U.kendall(mapped, truth["group_orderings"][group])
            out["initial_truth_distances"][group] = distance
            out["initial_matches_truth"][group] = distance == 0
            q, mixing = d["abnormal_priors"][group], d["normal_mixing_params"][group]
            if len(q) != p or len(mixing) != p or any(not _finite(v) or not 0 <= v <= 1 for v in q):
                raise ValueError("illegal oracle q")
            if any(not _finite(v) or abs(v - (1 - value)) > 1e-12 for v, value in zip(mixing, q)):
                raise ValueError("mixing differs from1-q")
            if any(mixing[initial[i]] > mixing[initial[i + 1]] for i in range(p - 1)):
                raise ValueError("initial ordering does not sort normal mixing; ties may have any order")
            out["q"][group] = dict(zip(posterior_names, q))
            if engine.startswith("oracle_equal"):
                if any(value != .5 for value in q) or d.get("oracle_prior_mode") != "equal_prior":
                    raise ValueError("equal-prior arm is not equal q=.5")
                continue
            if d.get("oracle_prior_mode") != "estimated_group_prior":
                raise ValueError("estimated-prior mode mismatch")
            fits = d["prior_fits"][group]
            if set(fits) != set(posterior_names):
                raise ValueError("incomplete biomarker prior fits")
            for name, value in zip(posterior_names, q):
                fit = fits[name]
                status = fit["status"]
                identified_statuses = {"interior", "lower_boundary", "upper_boundary"}
                unidentified_statuses = {"no_observations", "identical_components", "flat_observed_likelihood"}
                if status not in identified_statuses | unidentified_statuses:
                    raise ValueError("unknown frozen prior-fit status")
                if fit.get("identified") is not (status in identified_statuses):
                    raise ValueError("prior status disagrees with identification flag")
                if status in unidentified_statuses and value != .5:
                    raise ValueError("unidentified prior must use documented q=.5 fallback")
                if not _finite(fit.get("abnormal_prior")) or abs(fit["abnormal_prior"] - value) > 1e-12:
                    raise ValueError("q differs from prior-fit estimate")
                if not 1e-6 <= value <= 1 - 1e-6 or type(fit.get("identified")) is not bool:
                    raise ValueError("prior bounds/identification flag invalid")
                if (status == "lower_boundary" and abs(value - 1e-6) > 1e-12) or (
                        status == "upper_boundary" and abs(value - (1 - 1e-6)) > 1e-12):
                    raise ValueError("boundary status disagrees with q")
                if any(type(fit.get(k)) is not int or fit[k] < 0 for k in ("n_observed", "n_missing")):
                    raise ValueError("invalid prior observed/missing counts")
                status_counts[status] += 1
                out["n_boundary"] += status.endswith("_boundary")
                out["n_unidentified"] += not fit["identified"]
                out["n_prior_fits"] += 1
        if engine.startswith("oracle_estimated"):
            out["expected_prior_fits"] = 3 * p
            if d.get("prior_fit_calls") != out["n_prior_fits"] or out["n_prior_fits"] != 3 * p:
                raise ValueError("prior-fit call count differs from3*events")
            if d.get("prior_boundary_fits") != out["n_boundary"] or d.get("prior_unidentified_fits") != out["n_unidentified"]:
                raise ValueError("prior boundary/unidentified counter mismatch")
        out["q_fit_status_counts"] = dict(status_counts)
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        issues.append(_issue("oracle_diagnostic_integrity", cell, seed, engine, message=str(exc)))
    return out, issues


def _bootstrap_columns(columns, repeats, seed, context):
    """One common dataset-index draw per replicate for ALL endpoint columns."""
    if not columns:
        return {}
    n = len(next(iter(columns.values())))
    if n < 2 or any(len(v) != n for v in columns.values()):
        raise ValueError("Bootstrap columns require equal complete dataset rows, R>=2")
    stream_seed = int.from_bytes(hashlib.sha256(f"{seed}:{context}".encode()).digest()[:8], "big")
    rng = random.Random(stream_seed)
    means = {key: [] for key in columns}
    for _ in range(repeats):
        selected = rng.choices(range(n), k=n)
        for key, values in columns.items():
            means[key].append(math.fsum(values[i] for i in selected) / n)
    return {key: [U.A.percentile(sorted(values), .025), U.A.percentile(sorted(values), .975)]
            for key, values in means.items()}


def _bounded_difference(a, b):
    low = sum(y - x if x is not None and y is not None else y - 1 if y is not None
              else -x if x is not None else -1 for x, y in zip(a, b)) / len(a)
    high = sum(y - x if x is not None and y is not None else y if y is not None
               else 1 - x if x is not None else 1 for x, y in zip(a, b)) / len(a)
    return [low, high]


def _scope(issue, cells):
    if issue.get("cell") in cells:
        return issue["cell"]
    if isinstance(issue.get("key"), list) and issue["key"][0] in cells:
        return issue["key"][0]
    location = str(issue.get("path", ""))
    matches = {int(v) for v in re.findall(r"cell_(\d+)", location)}
    if len(matches) == 1 and next(iter(matches)) < len(cells):
        return cells[next(iter(matches))]
    for cell in cells:
        if location == cell or location.startswith(cell + ":"):
            return cell
    return None



def _audit_persisted_evidence(bundle):
    """Cross-check derived rows against raw observed records; never repair either."""
    issues = []
    cfg = bundle["config"]
    for name, manifest in bundle["manifests"].items():
        directory = Path(name).parent
        cell = U.cell_id(cfg["cells"][manifest["cell_index"]])
        log = directory / "fit_records.jsonl"
        if log.exists() and log.stat().st_size:
            with log.open("rb") as stream:
                stream.seek(-1, 2)
                if stream.read(1) != b"\n":
                    issues.append(_issue("uncommitted_raw_log_tail", cell, path=str(log)))
        rows_path = directory / "rows.jsonl"
        if not rows_path.exists():
            continue  # inventory already records missing derived-row coverage
        before = rows_path.stat()
        try:
            with rows_path.open() as stream:
                for line in stream:
                    if not line.strip():
                        continue
                    row = json.loads(line)
                    record = bundle["records"].get((cell, row["seed"], row["engine"], None, None))
                    if record is None:
                        continue  # missing/duplicate raw identity is already non-estimable
                    fields_match = all(row.get("observed", {}).get(k) == record.get(k)
                                       for k in ("orderings", "taus", "distance_to_truth"))
                    diagnostics_match = all(row.get("diagnostics", {}).get(k) == v
                                            for k, v in record.get("diagnostics", {}).items())
                    if (row.get("status") != record.get("status") or row.get("truth") != record.get("truth")
                            or not fields_match or not diagnostics_match):
                        issues.append(_issue("derived_row_raw_record_mismatch", cell, row["seed"], row["engine"], path=str(rows_path)))
        except (ValueError, KeyError, TypeError, OSError) as exc:
            issues.append(_issue("derived_row_read_error", cell, path=str(rows_path), message=str(exc)))
        after = rows_path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            issues.append(_issue("derived_rows_changed_during_crosscheck", cell, path=str(rows_path)))
    return issues


def _summarize_bundle(bundle, repeats=10000, bootstrap_seed=BOOTSTRAP_SEED):
    """Config-driven core for small fixtures; production entry point verifies R200x4x6."""
    if type(repeats) is not int or repeats < 1 or type(bootstrap_seed) is not int:
        raise ValueError("Positive bootstrap repeats and an integer fixed seed required")
    cfg, records, inventory = bundle["config"], bundle["records"], bundle["inventory"]
    if tuple(cfg["engines"]) != ENGINES or cfg.get("schemes") or cfg.get("bperm", 0) != 0:
        raise ValueError("Expected the six mechanism arms and no permutations")
    seeds = list(range(cfg["base_seed"], cfg["base_seed"] + cfg["datasets"]))
    cells = [U.cell_id(c) for c in cfg["cells"]]
    expected = U.expected_observed(cfg)
    design_plans, generator_hash = _design_plans(bundle["config_path"], cfg)
    environments_by_key = {}
    for manifest in bundle["manifests"].values():
        cell = U.cell_id(cfg["cells"][manifest["cell_index"]])
        for seed in manifest["seeds"]:
            environments_by_key[cell, seed] = manifest.get("environment")
    issues = list(bundle["issues"]) + _audit_persisted_evidence(bundle)
    issues += [{"code": "inventory:" + x["code"], "path": x.get("path"), "message": x.get("message")}
               for x in inventory.get("issues", []) if x.get("severity") == "error"]
    if set(records) - expected:
        issues.append(_issue("unplanned_records", keys=[list(k) for k in sorted(set(records) - expected, key=str)]))
    if not inventory.get("stable_rows_snapshot", False):
        issues.append(_issue("unstable_inventory_snapshot"))
    profiles = U.manifest_profiles(bundle)
    if len(profiles) > 1:
        issues.append(_issue("numerical_environment_mismatch"))
    checks, oracle_rows, pairing, dataset_hashes = {}, [], [], []
    for cell in cells:
        for seed in seeds:
            hashes = {e: U.data_hash(records[cell, seed, e, None, None]) if (cell, seed, e, None, None) in records else None for e in ENGINES}
            existing = [v for v in hashes.values() if isinstance(v, str) and U.I.is_hash(v)]
            same = len(existing) == 6 and len(set(existing)) == 1
            mismatch = len(set(existing)) > 1
            saved_truths = [records[cell, seed, e, None, None].get("truth") for e in ENGINES
                            if (cell, seed, e, None, None) in records and records[cell, seed, e, None, None].get("truth")]
            truth_conflict = bool(saved_truths) and any(t != saved_truths[0] for t in saved_truths[1:])
            dataset_hashes.append({"cell": cell, "seed": seed, "all_six_hashes_present_and_equal": same,
                                   "available_hashes_conflict": mismatch, "hashes": hashes,
                                   "all_six_truths_present_and_equal": len(saved_truths) == 6 and not truth_conflict})
            if mismatch:
                issues.append(_issue("within_dataset_data_hash_mismatch", cell, seed))
            if truth_conflict:
                issues.append(_issue("within_dataset_truth_metadata_mismatch", cell, seed))
            for engine in ENGINES:
                key = cell, seed, engine, None, None
                row = _validate_record(records.get(key), cell, seed, engine)
                try:
                    _validate_saved_truth(records.get(key), design_plans[cell], generator_hash,
                                          environments_by_key.get((cell, seed)))
                except (ValueError, KeyError, TypeError, IndexError) as exc:
                    row["values"], row["endpoint_legal"] = None, False
                    row["issues"].append(_issue("truth_frozen_design_or_provenance_mismatch", cell, seed, engine, message=str(exc)))
                if row["status"] == "ok" and (not U.I.is_hash(hashes[engine]) or mismatch or truth_conflict):
                    row["values"], row["endpoint_legal"] = None, False
                    row["issues"].append(_issue("unverified_endpoint_data_identity", cell, seed, engine))
                checks[key] = row
                issues.extend(row["issues"])
                oracle, extra = _oracle_diagnostics(records.get(key), cell, seed, engine)
                if oracle is not None:
                    oracle_rows.append(oracle)
                issues.extend(extra)
            for first, second in (("original", "repaired"), ("oracle_equal_original", "oracle_equal_repaired"),
                                  ("oracle_estimated_original", "oracle_estimated_repaired")):
                a, b = (records.get((cell, seed, e, None, None)) for e in (first, second))
                success = bool(a and b and a.get("status") == b.get("status") == "ok")
                audit = {"cell": cell, "seed": seed, "first": first, "second": second,
                         "both_successful": success, "checks": {}}
                if success:
                    da, db = a.get("diagnostics", {}), b.get("diagnostics", {})
                    fields = ("optimizer", "optimizer_calls", "optimizer_failures") if first == "original" else (
                        "abnormal_priors", "normal_mixing_params", "initial_orderings", "posterior_biomarkers")
                    if first.startswith("oracle_estimated"):
                        fields += ("prior_fits",)
                    audit["checks"] = {f: f in da and f in db and da[f] == db[f] for f in fields}
                    if first == "original" and cfg.get("paired_standard"):
                        aa, ab = da.get("model_artifacts", {}), db.get("model_artifacts", {})
                        audit["checks"].update({"shared_fit_id": bool(a.get("shared_fit_id")) and a.get("shared_fit_id") == b.get("shared_fit_id"),
                            **{f: f in aa and f in ab and aa[f] == ab[f] for f in
                               ("posterior_sha256", "model_biomarkers", "biomarker_parameters")}})
                    if not all(audit["checks"].values()):
                        issues.append(_issue("consensus_pair_nuisance_audit_mismatch", cell, seed,
                                             first=first, second=second, checks=audit["checks"]))
                pairing.append(audit)
    coverage = {x["cell"]: x.get("complete", False) for x in inventory.get("manifest_coverage", [])}
    inv_engines = {(x["cell"], x["engine"]): x for x in inventory.get("cell_engines", [])}
    cohorts, effects, summaries = [], [], []
    for cell in cells:
        rows = [checks[cell, s, e, None, None] for s in seeds for e in ENGINES]
        missing = [r for r in rows if not r["present"]]
        failures = [r for r in rows if r["present"] and r["status"] != "ok"]
        cell_issues = [i for i in issues if _scope(i, cells) in (None, cell)]
        derived_complete = all((cell, e) in inv_engines and not any(inv_engines[cell, e].get(k)
                                for k in ("missing_seeds", "duplicate_seeds", "pending_observed_seeds")) for e in ENGINES)
        terminal = all(r["status"] in ("ok", "error", "invalid") for r in rows)
        complete = not missing and terminal and coverage.get(cell, False) and derived_complete
        eligible = complete and not cell_issues
        cohorts.append({"cell": cell, "planned_datasets": len(seeds), "planned_records": 6 * len(seeds),
                        "present_records": sum(r["present"] for r in rows), "complete_planned_accounting": complete,
                        "missing_keys": [[r["seed"], r["engine"]] for r in missing],
                        "failed_or_invalid_keys": [[r["seed"], r["engine"], r["status"]] for r in failures],
                        "bootstrap_cohort_gate": eligible, "integrity_issue_count": len(cell_issues)})
        columns, cell_effects = {}, []
        for first, second in CONTRASTS:
            for field, labels in (("distance_to_truth", U.GROUPS + ("dataset_mean",)),
                                  ("taus", PAIR_NAMES + ("dataset_mean",))):
                for index, label in enumerate(labels):
                    values, grid_values = [], []
                    for engine in (first, second):
                        grid_values.append([checks[cell, s, engine, None, None]["values"][field]
                                            if checks[cell, s, engine, None, None]["values"] is not None else None for s in seeds])
                        values.append([checks[cell, s, engine, None, None]["values"][field][index]
                                       if checks[cell, s, engine, None, None]["values"] is not None else None for s in seeds])
                    a, b = values
                    known = [y - x for x, y in zip(a, b) if x is not None and y is not None]
                    estimable = eligible and len(known) == len(seeds)
                    key = first, second, field, label
                    if estimable:
                        columns[key] = known
                    cell_effects.append({"cell": cell, "first": first, "second": second, "field": field,
                        "group_or_pair": label, "contrast_direction": "second_minus_first",
                        "planned_datasets": len(seeds), "known_paired_datasets": len(known),
                        "unknown_paired_datasets": len(seeds) - len(known), "estimable": estimable,
                        "not_estimable_reason": None if estimable else "incomplete_cell_or_integrity_issue" if not eligible else "failed_or_invalid_endpoint_in_planned_pair",
                        "all_planned_difference_identification_bounds": _bounded_difference(a, b),
                        "bounds_trustworthy": not cell_issues,
                        "first_mean": statistics.mean(a) if estimable else None,
                        "second_mean": statistics.mean(b) if estimable else None,
                        "estimate": statistics.mean(known) if estimable else None,
                        "mcse_mean": statistics.stdev(known) / math.sqrt(len(known)) if estimable and len(known) > 1 else None,
                        "bootstrap_ci95": None, "bootstrap_degenerate_sample": estimable and len(set(known)) == 1,
                        "n_lower_equal_higher": _integer_grid_directions(*grid_values, len(design_plans[cell]["biomarker_names"]), index),
                        "direction_count_method": "exact_Kendall_integer_grid; recover and validate three k/choose(m,2) components; dataset mean compares sum(k) with denominator3*choose(m,2), without epsilon",
                        "direction_count_grid_denominator": math.comb(len(design_plans[cell]["biomarker_names"]), 2) * (3 if index == 3 else 1),
                        "count_denominator_note": "Known paired datasets only for direction counts; planned R stays unchanged."})
        intervals = _bootstrap_columns(columns, repeats, bootstrap_seed, cfg["run_id"] + ":" + cell)
        for effect in cell_effects:
            effect["bootstrap_ci95"] = intervals.get((effect["first"], effect["second"], effect["field"], effect["group_or_pair"]))
        effects.extend(cell_effects)
        for engine in ENGINES:
            erows = [checks[cell, s, engine, None, None] for s in seeds]
            present = [r for r in erows if r["present"]]
            consensuses = [c for r in present for c in r.get("consensus", [])]
            finite_flags = lambda key: [r.get(key) if r["present"] else None for r in erows]
            local_flags = [any(c.get("local_optimal") is False for c in r.get("consensus", []))
                           if len(r.get("consensus", [])) == 3 and all(c.get("local_optimal") in (True, False)
                                                                    for c in r["consensus"]) else None for r in erows]
            means = [r["values"]["distance_to_truth"][3] if r["values"] else None for r in erows]
            known = [x for x in means if x is not None]
            summaries.append({"cell": cell, "engine": engine, "planned_datasets": len(seeds),
                "status_counts_planned": dict(Counter(r["status"] for r in erows)),
                "endpoint_legal": _bounded_flag([r["endpoint_legal"] if r["present"] else None for r in erows]),
                "truth_error_mean": statistics.mean(known) if eligible and len(known) == len(seeds) else None,
                "truth_error_all_planned_identification_bounds": _bounded_mean(means),
                "has_quality_flags": _bounded_flag([bool(r.get("quality_flags")) if r["present"] and r["status"] == "ok" else None for r in erows]),
                "quality_flag_counts": dict(Counter(f for r in present for f in r.get("quality_flags", []) if isinstance(f, str))),
                "parameters_finite": _bounded_flag(finite_flags("parameters_finite")) if engine in ("original", "repaired") else {"applicable": False},
                "event_centers_finite": _bounded_flag(finite_flags("event_centers_finite")),
                "any_remaining_improving_neighbor": _bounded_flag(local_flags),
                "recorded_consensus_count": len(consensuses), "planned_group_consensuses": 3 * len(seeds),
                "nonoptimal_recorded_consensuses": sum(c.get("local_optimal") is False for c in consensuses),
                "residual_improvement": _numbers(c.get("residual_improvement") for c in consensuses),
                "optimizer_calls_recorded": sum(r.get("optimizer_calls", 0) for r in present),
                "optimizer_unsuccessful_recorded": sum(r.get("optimizer_unsuccessful", 0) for r in present),
                "optimizer_nonfinite_recorded": sum(r.get("optimizer_nonfinite", 0) for r in present),
                "subunit_counts_are_not_independent_R": True})
    oracle_summaries = []
    for cell in cells:
        for engine in ENGINES[2:]:
            selected = [r for r in oracle_rows if r["cell"] == cell and r["engine"] == engine]
            matches = [v for r in selected for v in r["initial_matches_truth"].values()]
            counts = Counter()
            for row in selected:
                counts.update(row["q_fit_status_counts"])
            oracle_summaries.append({"cell": cell, "engine": engine, "planned_datasets": len(seeds),
                "present_datasets": len(selected), "planned_group_initializations": 3 * len(seeds),
                "recorded_group_initializations": len(matches), "initial_exact_truth_count": sum(matches),
                "initialization_unknown_count": 3 * len(seeds) - len(matches),
                "initial_distance_to_truth_recorded": _numbers(v for r in selected for v in r["initial_truth_distances"].values()),
                "prior_fits_recorded": sum(r["n_prior_fits"] for r in selected),
                "q_boundary_fits_recorded": sum(r["n_boundary"] for r in selected),
                "q_unidentified_fits_recorded": sum(r["n_unidentified"] for r in selected),
                "q_fit_status_counts": dict(counts),
                "subunit_counts_are_not_independent_R": True})
    result = {"schema_version": "confirmation-mechanism-six-arm-v1", "created_utc": U.I.utc(),
            "run_id": cfg["run_id"], "phase": cfg["phase"], "config": cfg,
            "planned_cell_datasets": len(cells) * len(seeds), "planned_fit_records": len(expected),
            "collection_complete": all(c["complete_planned_accounting"] for c in cohorts),
            "ready": all(c["complete_planned_accounting"] for c in cohorts) and not issues,
            "all_prespecified_effects_estimable": all(e["estimable"] for e in effects),
            "bootstrap": {"repeats": repeats, "fixed_seed": bootstrap_seed,
                          "resampling_unit": "entire_dataset_shared_index_draw_across_all_arms_groups_and_endpoints_in_each_cell",
                          "coverage": "pointwise95_not_simultaneous", "cross_N_formal_test": False},
            "inputs": bundle["inputs"], "inventory": inventory, "issues": issues,
            "resolved_frozen_designs": design_plans, "frozen_generator_sha256": generator_hash,
            "cohorts": cohorts, "effects": effects, "engine_summaries": summaries,
            "dataset_hash_checks": dataset_hashes, "nuisance_pair_checks": pairing,
            "dataset_diagnostics": [{k: v for k, v in r.items() if k not in ("values", "issues")} for r in checks.values()],
            "oracle_diagnostics": oracle_rows, "oracle_summaries": oracle_summaries,
            "fitting_or_generation_performed": False,
            "limitations": LIMITATIONS, "nonfinite_json_encoding": "null; associated diagnostic integrity issue retained"}
    return _json_safe(result)


def _protocol_config(protocol_path, workspace):
    protocol = json.loads(Path(protocol_path).read_text())
    runs = [r for r in protocol["runs"] if r["run_id"] == "confirm_mechanism_b"]
    if len(runs) != 1:
        raise ValueError("Exactly one frozen confirm_mechanism_b entry required")
    run = runs[0]
    snapshot = (workspace / run["snapshot"]).resolve()
    if not snapshot.is_relative_to(workspace) or _sha(snapshot / "SHA256SUMS") != run["snapshot_sha256"]:
        raise ValueError("Frozen snapshot manifest does not match protocol")
    names = set()
    for line in (snapshot / "SHA256SUMS").read_text().splitlines():
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            raise ValueError("Malformed snapshot checksum line")
        digest, name = match.groups()
        path = (snapshot / name).resolve()
        if name in names or not path.is_relative_to(snapshot) or _sha(path) != digest:
            raise ValueError("Frozen source path/hash mismatch")
        names.add(name)
    config_path = snapshot / run["config"]
    if run["config"] not in names or _sha(config_path) != run["config_sha256"]:
        raise ValueError("Frozen config hash mismatch")
    cfg = json.loads(config_path.read_text())
    expected_cells = ["IID_H0", "REF_H0", "REF_H0_N4", "STAGE_H0"]
    if (cfg["run_id"] != "confirm_mechanism_b" or cfg["phase"] != "confirmation" or cfg["datasets"] != 200
            or tuple(cfg["engines"]) != ENGINES or [U.cell_id(c) for c in cfg["cells"]] != expected_cells
            or cfg.get("schemes") or cfg.get("bperm", 0) != 0):
        raise ValueError("Unexpected confirmation design: expected frozen R200x4x6 mechanism plan")
    return config_path, run


def analyze(protocol_path, raw_root, *, repeats=10000, bootstrap_seed=BOOTSTRAP_SEED, workspace=None):
    workspace = Path(workspace or Path(__file__).resolve().parents[2]).resolve()
    config_path, run = _protocol_config(protocol_path, workspace)
    bundle = U.load_run(config_path, Path(raw_root))
    result = _summarize_bundle(bundle, repeats, bootstrap_seed)
    result["protocol"] = {"path": str(Path(protocol_path).resolve()), "sha256": _sha(protocol_path),
                          "config_path": str(config_path), "config_sha256": run["config_sha256"],
                          "snapshot_sha256": run["snapshot_sha256"]}
    result["analysis_source_sha256"] = {name: _sha(Path(__file__).with_name(name)) for name in
        ("analyze_confirmation_mechanism.py", "run_record_audit.py", "analyze_v2.py")}
    result["analysis_source_sha256"]["inventory_run.py"] = _sha(U.I.__file__)
    return result


def markdown_lines(result):
    lines = ["# Confirmation mechanism B：六臂机制分析", "",
             f"计划 {result['planned_cell_datasets']} 个cell-dataset、{result['planned_fit_records']} 条拟合记录。"
             f"完整计划记账：{result['collection_complete']}；全部预定对比可估计：{result['all_prespecified_effects_estimable']}。", "",
             "未完成格不提供bootstrap。失败／非法端点保留计划分母；识别区间不是置信区间。差值为second−first，真序误差负差表示恢复改善。", "",
             "| Cell | 已有/计划记录 | 完整记账 | 失败/非法记录 | CI cohort gate |",
             "|---|---:|---|---:|---|"]
    lines += [f"| {c['cell']} | {c['present_records']}/{c['planned_records']} | {c['complete_planned_accounting']} | "
              f"{len(c['failed_or_invalid_keys'])} | {c['bootstrap_cohort_gate']} |" for c in result["cohorts"]]
    lines += ["", "| Cell | 预定对比 | 计划/可用pairs | 平均真序误差差 | Pointwise95% CI | 全计划识别区间 |",
              "|---|---|---:|---:|---|---|"]
    for e in result["effects"]:
        if e["field"] == "distance_to_truth" and e["group_or_pair"] == "dataset_mean":
            lines.append(f"| {e['cell']} | {e['first']} → {e['second']} | {e['planned_datasets']}/{e['known_paired_datasets']} | "
                         f"{U.fmt(e['estimate'])} | {e['bootstrap_ci95']} | {e['all_planned_difference_identification_bounds']} |")
    lines += ["", "逐组真序距离、逐pair距离、完整质量标记、remaining-neighbor审计、q边界及初始化真S匹配见JSON。"
              "标准臂optimizer记录是可见审计历史，不是重新验证完整数值轨迹。", ""]
    lines += ["| Cell / oracle | 初始序列等于真S / 计划组序列 | 记录q拟合数 | q边界 / 未识别 |",
              "|---|---:|---:|---:|"]
    lines += [f"| {x['cell']} / {x['engine']} | {x['initial_exact_truth_count']}/{x['planned_group_initializations']} | "
              f"{x['prior_fits_recorded']} | {x['q_boundary_fits_recorded']}/{x['q_unidentified_fits_recorded']} |"
              for x in result["oracle_summaries"]]
    lines += ["", "未记录初始化列为unknown，不能把它算成不匹配；以上组／biomarker计数不是独立R。", ""]
    lines += ["- " + item for item in result["limitations"]]
    if result["issues"]:
        lines += ["", f"输入／一致性问题共{len(result['issues'])}项，详见JSON；这些问题不能通过丢弃记录消除。"]
    return lines



def _output_directory(out, protocol_path, raw_root, workspace):
    """Resolve final artifact paths before reading results; preserve all raw/frozen inputs."""
    out, protocol_path, raw_root, workspace = (Path(x).resolve() for x in (out, protocol_path, raw_root, workspace))
    protocol = json.loads(protocol_path.read_text())
    protected = [raw_root, workspace / "runs/v2/snapshots"]
    protected += [(workspace / run["snapshot"]).resolve() for run in protocol["runs"]]
    for suffix in (".json", ".md"):
        path = (out / ("confirmation_mechanism_analysis" + suffix)).resolve()
        if path == protocol_path or any(path.is_relative_to(root) for root in protected):
            raise ValueError("Analysis output would overwrite raw results, protocol or a frozen snapshot")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--raw", type=Path, required=True, help="Explicit synchronized raw run root; no automatic HPC retrieval")
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=10000)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    workspace = Path(args.workspace or Path(__file__).resolve().parents[2]).resolve()
    out = _output_directory(args.out, args.protocol, args.raw, workspace)
    result = analyze(args.protocol, args.raw, repeats=args.bootstrap, bootstrap_seed=args.bootstrap_seed, workspace=workspace)
    U.save(out, "confirmation_mechanism_analysis", result, markdown_lines(result))


if __name__ == "__main__":
    main()
