#!/usr/bin/env python3
"""Post-freeze, fixed-cohort comparison of pair Bonferroni permutation operators.

Standard library, read only scientific inputs. No fits, raw-event import, or
automatic discovery/pooling of runs. See the explicitly required JSON contract.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile

import analyze_confirmation_cohorts as C

VERSION = "confirmation_operator_comparison_v1"
RUN_IDS = ("confirm_core_a", "confirm_pair_complete_ref_a", "confirm_power_global_a", "confirm_pair_power_k27_a")
COHORTS = (
    ("REF_H0", RUN_IDS[0], RUN_IDS[1], 42100000, 1000, (True, True, True)),
    ("PWR_E2_K27", RUN_IDS[2], RUN_IDS[3], 42200000, 500, (False, False, True)),
    ("PWR_E4_K27", RUN_IDS[2], RUN_IDS[3], 42200000, 500, (True, False, False)),
)
ENDPOINTS = ("any_true_null_rejection", "any_false_null_rejection", "all_false_null_rejection")


def validate_fullgroup_definition(definition):
    C.require(isinstance(definition, dict)
              and definition.get("stratify") == "diagnosis"
              and definition.get("groups", [0, 1, 2]) == [0, 1, 2]
              and definition.get("require_max") is True
              and "tested_pair_index" not in definition and "family" not in definition,
              "Expected full-three-group diagnosis reference, not a local or max-only replacement")


def fullgroup_pair_decision(run, cell, seed, index):
    """Independent pair extractor: fixed B+1, inclusive ties, three comparisons.

    The reference also stores max counts, but max rejection is never this endpoint.
    I.scheme_state checks all saved declarations against completion bounds first.
    """
    row, reason = C.get_record(run, cell, "repaired", seed)
    if row is None:
        return None, reason
    try:
        C.require(type(index) is int and 0 <= index < 3, "Invalid target pair")
        cfg = run["cfg"]
        C.require((cfg["bperm"], cfg["alpha"], cfg["rule"]) == (599, .05, "le"), "Expected B599 alpha .05 le")
        definition = cfg["schemes"]["diagnosis"]
        validate_fullgroup_definition(definition)
        entry = row.get("schemes", {}).get("diagnosis")
        C.require(isinstance(entry, dict) and entry.get("definition") == definition, "Missing or mismatched fullgroup scheme definition")
        C.I.scheme_state(entry, definition, cfg, seed in run["full_p"])
        result = C.I.bounds(entry["gt"][index] + entry["eq"][index], entry["nperm"], 599, .05, 3, "le")
        return result, None if result is not None else "permutation_decision_unknown"
    except (ValueError, TypeError, KeyError, AttributeError, IndexError) as exc:
        return None, "invalid_fullgroup_scheme: " + str(exc)


def local_pair_decision(run, cell, seed, index, scheme):
    row, reason = C.get_record(run, cell, "repaired", seed)
    if row is None:
        return None, reason
    try:
        entry = row.get("schemes", {}).get(scheme)
        C.require(isinstance(entry, dict) and entry.get("definition") == run["cfg"]["schemes"][scheme], "Missing or mismatched local scheme definition")
        return C.decision(run, cell, "repaired", seed, scheme, target=index)
    except (ValueError, TypeError, KeyError, AttributeError, IndexError) as exc:
        return None, "invalid_local_scheme: " + str(exc)


def local_schemes(run):
    result = {}
    for name, definition in run["cfg"]["schemes"].items():
        index = definition.get("tested_pair_index")
        C.validate_local_definition(definition, index)
        C.require(definition.get("family") == C.FAMILY and index not in result, "Duplicate/wrong local family member")
        result[index] = name
    C.require(set(result) == {0, 1, 2}, "Incomplete local family")
    family = run["cfg"].get("family_definition", {})
    C.require(family.get("id") == C.FAMILY and family.get("level") == .05
              and family.get("pair_order") == list(C.PAIR_NAMES)
              and family.get("member_schemes") == [result[i] for i in range(3)], "Frozen family declaration differs")
    return result


def validate_plan(runs):
    C.require(set(runs) == set(RUN_IDS), "Require exactly the four contracted runs; no topups or other K values")
    for cell, source_id, local_id, base, n, null in COHORTS:
        source, local = runs[source_id], runs[local_id]
        expected = set(range(base, base+n))
        for run in (source, local):
            cfg = run["cfg"]
            C.require(run["seeds"] == expected and cfg["datasets"] == n and cfg["base_seed"] == base, "Frozen R/absolute seed plan differs")
            C.require((cfg["bperm"], cfg["alpha"], cfg["rule"], cfg.get("fast_likelihood")) == (599, .05, "le", True), "Frozen B/alpha/rule/fast implementation differs")
            C.require("repaired" in cfg["engines"] and cell in run["names"], "Missing contracted cell/engine")
            C.require(C.planned_null(run["designs"][cell]) == list(null), "Wrong planned null mask")
        C.require(source["designs"][cell] == local["designs"][cell], "Source/local resolved design differs")
        C.require(source["cfg"].get("paired_standard", False) is (source_id == "confirm_core_a")
                  and local["cfg"].get("paired_standard", False) is False, "Unexpected paired execution policy")
        validate_fullgroup_definition(source["cfg"]["schemes"]["diagnosis"])
        local_schemes(local)
        reuse = local["cfg"]["data_and_fit_reuse"]
        expected_cells = ["REF_H0"] if local_id == RUN_IDS[1] else ["PWR_E2_K27", "PWR_E4_K27"]
        C.require(reuse["source_run_id"] == source_id and reuse.get("source_engine", reuse.get("engine")) == "repaired"
                  and reuse["source_cells"] == expected_cells and reuse["source_seed_start"] == base
                  and reuse["source_seed_end_inclusive"] == base+n-1
                  and reuse.get("source_reference_scheme") == "diagnosis"
                  and reuse.get("same_data_not_new_replicates") is True, "Source reuse declaration differs")
    C.require(runs[RUN_IDS[0]]["names"] == ["IID_H0", "REF_H0"]
              and runs[RUN_IDS[0]]["cfg"]["engines"] == ["original", "repaired"], "Unexpected source core layout")
    C.require(runs[RUN_IDS[1]]["names"] == ["REF_H0"], "Unexpected local REF cells")
    for run_id in RUN_IDS[2:]:
        C.require(runs[run_id]["names"] == ["PWR_E2_K27", "PWR_E4_K27"], "Unexpected K27 cells")
    for run_id in RUN_IDS[1:]:
        C.require(runs[run_id]["cfg"]["engines"] == ["repaired"], "Unexpected local/power engine")


def endpoint_value(outcomes, null, endpoint):
    indexes = [i for i in range(3) if null[i] == (endpoint == "any_true_null_rejection")]
    C.require(indexes and endpoint in ENDPOINTS, "Endpoint has no planned target pairs")
    values = [outcomes[i][0] for i in indexes]
    value = C.all3(values) if endpoint == "all_false_null_rejection" else C.any3(values)
    reasons = sorted({outcomes[i][1] or "unknown_pair_decision" for i in indexes if outcomes[i][0] is None})
    return value, None if value is not None else "; ".join(reasons)


def join_copy(runs):
    """Normalize malformed failed truth containers only for legacy join traversal.

    Validation reasons/status remain intact. Nothing repairs or substitutes a fit.
    Returning copies also avoids modifying caller records when compare_sources
    marks source/local mismatches on the in-memory validation wrappers.
    """
    result = {}
    for run_id, run in runs.items():
        clone = dict(run, records={})
        for key, items in run["records"].items():
            kept = []
            for item in items:
                row = dict(item["row"])
                truth = row.get("truth")
                truth = dict(truth) if isinstance(truth, dict) else {}
                if not isinstance(truth.get("manifest"), dict):
                    truth["manifest"] = {}
                row["truth"] = truth
                kept.append(dict(item, row=row))
            clone["records"][key] = kept
        result[run_id] = clone
    return result


def summarize_comparisons(runs):
    """Pure postprocessing of loader-validated records, preserving planned R."""
    validate_plan(runs)
    runs = join_copy(runs)
    cross_issues = []
    environments = {C.canonical_sha(r["environment"]) for r in runs.values() if r["environment"] is not None}
    sources = {C.canonical_sha(r["sources"]) for r in runs.values()}
    if len(environments) > 1 or len(sources) != 1:
        cross_issues.append({"code": "cross_run_source_or_environment_mismatch"})
        for run in runs.values():
            run["trusted"] = False
    joins = C.compare_sources(runs)
    endpoints = []
    for cell, source_id, local_id, base, n, null in COHORTS:
        source, local = runs[source_id], runs[local_id]
        schemes = local_schemes(local)
        names = ENDPOINTS[:1] if all(null) else ENDPOINTS
        values = {name: {"fullgroup": [], "local": []} for name in names}
        for seed in range(base, base+n):
            full = [fullgroup_pair_decision(source, cell, seed, i) for i in range(3)]
            matched = [local_pair_decision(local, cell, seed, i, schemes[i]) for i in range(3)]
            for name in names:
                values[name]["fullgroup"].append((seed, *endpoint_value(full, null, name)))
                values[name]["local"].append((seed, *endpoint_value(matched, null, name)))
        for name in names:
            a, b = values[name]["fullgroup"], values[name]["local"]
            paired = C.paired_summary(a, b)
            paired["difference_orientation"] = "local_minus_fullgroup"
            paired["paired_table"] = dict(Counter(f"{int(x)}{int(y)}" for (_, x, _), (_, y, _) in zip(a, b) if x is not None and y is not None))
            paired["paired_table_order"] = ["fullgroup", "local"]
            endpoints.append({"endpoint_id": cell + ":" + name, "cell": cell, "endpoint": name,
                              "engine": "repaired", "fullgroup_run": source_id, "local_run": local_id,
                              "planned_R": n, "pair_null": list(null), "fullgroup": C.summarize(a),
                              "local": C.summarize(b), "paired_difference": paired})
    C.require(len(endpoints) == 7, "Expected seven contracted comparisons")
    reports = []
    for run_id in RUN_IDS:
        run = runs[run_id]
        inv = run["inventory"]
        reports.append({"run_id": run_id, "trusted": run["trusted"], "issues": run["issues"],
                        "collection_complete": inv["collection_complete"], "inventory_ready": inv["ready_for_final_analysis"],
                        "plan": inv["plan"], "totals": inv["totals"], "environment": run["environment"],
                        "sources": run["sources"], "frozen_config": str(run["config_path"]),
                        "snapshot": str(run["snapshot"]), "paired_standard": run["cfg"].get("paired_standard", False),
                        "fit_timeout_s": run["cfg"].get("fit_timeout_s")})
    complete = all(e["paired_difference"]["unknown_pairs"] == 0 for e in endpoints)
    collected = all(r["collection_complete"] for r in reports)
    provenance_ok = not cross_issues and all(r["trusted"] for r in reports) and not any(j["mismatches"] for j in joins)
    ready = collected and complete and provenance_ok and all(j["all_planned_observed_exact"] for j in joins) and all(r["inventory_ready"] for r in reports)
    return {"schema_version": VERSION, "analysis_role": "post_freeze_supplementary", "difference_orientation": "local_minus_fullgroup",
            "planned_dataset_cell_units": 2000, "planned_method_observations": 4000,
            "counting_scope": "Combinatorial counts only. Dataset seed is the MC unit within each cohort; E2/E4 share CRN seeds, so cells are not pooled as independent replications.",
            "endpoint_count": 7, "all_endpoint_decisions_complete": complete, "collection_complete": collected,
            "collection_scope": "all obligations of the four source/local runs, including source core IID/original controls not used in endpoints",
            "ready_for_final_analysis": ready, "status": "complete" if ready else "terminal_with_unknown" if collected and provenance_ok else "incomplete",
            "cross_run_issues": cross_issues, "source_local_crosschecks": joins, "run_inventory": reports,
            "endpoints": endpoints, "point_estimate_scope": "each arm requires its own complete planned endpoint; contrasts require complete paired endpoints",
            "interval_scope": "seven pointwise supplementary comparisons, not simultaneous and not equivalence tests",
            "raw_event_audit_performed": False, "data_generation_or_fitting_performed": False}


def validate_contract(contract, protocol_bytes):
    C.require(contract.get("schema_version") == "post_freeze_operator_comparison_contract_v1", "Unsupported contract")
    C.require(contract["protocol_sha256"] == hashlib.sha256(protocol_bytes).hexdigest(), "Protocol differs from frozen supplementary contract")
    C.require(set(contract["registered_runs"]) == set(RUN_IDS) and contract["engine"] == "repaired"
              and (contract["B"], contract["alpha"], contract["rule"]) == (599, .05, "le")
              and contract["difference_orientation"] == "local_minus_fullgroup" and contract["number_of_endpoints"] == 7, "Contract operator/endpoint plan differs")
    expected = [{"cell": c, "fullgroup_run": s, "local_run": l, "base_seed": b, "R": n,
                 "pair_null": list(null), "endpoints": list(ENDPOINTS[:1] if all(null) else ENDPOINTS)} for c, s, l, b, n, null in COHORTS]
    C.require(contract["cohorts"] == expected, "Contract fixed cohort definitions differ")


def input_files(protocol_path, contract_path, project, snapshots, results_root):
    files = {Path(protocol_path).resolve(), Path(contract_path).resolve(), Path(__file__).resolve(), Path(C.__file__).resolve(), Path(C.I.__file__).resolve()}
    for snapshot in snapshots:
        files.update(p.resolve() for p in snapshot.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    for run_id in RUN_IDS:
        root = Path(results_root) / run_id
        for filename in ("manifest.json", "progress.json", "rows.jsonl"):
            files.update(p.resolve() for p in root.glob("cell_*/chunk_*/" + filename))
    return sorted(files)


def load_compact_run(run, results_root):
    # Legacy validation reads the complete row first. Retain only join/decision
    # evidence afterwards, including failed identities, without changing C.
    run = C.load_run(run, results_root)
    for items in run["records"].values():
        for item in items:
            row = item["row"]
            for field, keys in (("truth", ("config", "manifest", "biomarker_names")),
                                ("observed", ("orderings", "taus", "distance_to_truth"))):
                if isinstance(row.get(field), dict):
                    row[field] = {key: row[field][key] for key in keys if key in row[field]}
            row.pop("diagnostics", None)  # Already validated; run inventory records environment.
    return run


def analyze(protocol_path, results_root, contract_path, project=C.PROJECT, snapshots_root=None):
    protocol_path, contract_path = Path(protocol_path).resolve(), Path(contract_path).resolve()
    protocol_bytes = protocol_path.read_bytes()
    contract_bytes = contract_path.read_bytes()
    protocol, contract = json.loads(protocol_bytes), json.loads(contract_bytes)
    validate_contract(contract, protocol_bytes)
    C.require(protocol.get("configuration_frozen_before_confirmation_data") is True and protocol.get("confirmation_namespace") == "42xxxxxx", "Expected frozen confirmation protocol")
    entries = {e["run_id"]: e for e in protocol["runs"]}
    C.require(len(entries) == len(protocol["runs"]) and set(RUN_IDS) <= set(entries), "Missing/duplicate registered runs")
    for run_id in RUN_IDS:
        C.require(all(entries[run_id][field] == value for field, value in contract["registered_runs"][run_id].items()), "Contract snapshot/config pin differs")
    snapshots = [C.safe_child(snapshots_root, run_id) if snapshots_root else C.safe_child(project, entries[run_id]["snapshot"]) for run_id in RUN_IDS]
    files = input_files(protocol_path, contract_path, project, snapshots, results_root)
    hashes = {str(p): C.sha(p) for p in files}
    C.require(hashes[str(protocol_path)] == hashlib.sha256(protocol_bytes).hexdigest()
              and hashes[str(contract_path)] == hashlib.sha256(contract_bytes).hexdigest(), "Protocol/contract changed before loading")
    runs = {run_id: load_compact_run(C.frozen_run(entries[run_id], project, snapshots_root), results_root) for run_id in RUN_IDS}
    result = summarize_comparisons(runs)
    current = input_files(protocol_path, contract_path, project, snapshots, results_root)
    changed = sorted(set(hashes).symmetric_difference(str(p) for p in current) | {str(p) for p in current if str(p) in hashes and C.sha(p) != hashes[str(p)]})
    # Refuse releasing estimates derived from a mixture of input captures.
    C.require(not changed, "Inputs or source changed during read; no report certified: " + ", ".join(changed))
    result.update(inspected_utc=datetime.now(timezone.utc).isoformat(), protocol=str(protocol_path),
                  protocol_sha256=hashes[str(protocol_path)], contract=str(contract_path),
                  contract_sha256=hashes[str(contract_path)], results_root=str(Path(results_root).resolve()),
                  source_sha256={str(p): hashes[str(p)] for p in (Path(__file__).resolve(), Path(C.__file__).resolve(), Path(C.I.__file__).resolve())},
                  input_sha256=hashes, inputs_stable=True, protected_snapshot_roots=[str(p.resolve()) for p in snapshots],
                  scope_limits=contract["scope_limits"], uncertainty_contract=contract["uncertainty"])
    return result


def markdown(result):
    lines = ["# Post-freeze supplementary operator comparison", "", f"Status: **{result['status']}**; collection complete: **{result['collection_complete']}**; endpoint decisions complete: **{result['all_endpoint_decisions_complete']}**.", "",
             "Differences are local minus fullgroup. Fullgroup uses three pair Bonferroni decisions, not global-max rejection. Missing/failed/untrusted outcomes retain planned denominators. Bracketed completion bounds are not confidence intervals.", "",
             "| Endpoint | Planned R | Fullgroup rate or bounds (unknown) | Local rate or bounds (unknown) | Difference or bounds (unknown pairs) |", "|---|---:|---|---|---|"]
    def rate(s):
        point = f"{s['rate']:.5f}" if s["rate"] is not None else f"[{s['rate_bounds'][0]:.5f}, {s['rate_bounds'][1]:.5f}]"
        return point + f" ({s['unknown_count']})"
    for e in result["endpoints"]:
        p = e["paired_difference"]
        value = f"{p['difference']:.5f}" if p["difference"] is not None else f"[{p['difference_bounds'][0]:.5f}, {p['difference_bounds'][1]:.5f}]"
        lines.append(f"| {e['endpoint_id']} | {e['planned_R']} | {rate(e['fullgroup'])} | {rate(e['local'])} | {value} ({p['unknown_pairs']}) |")
    lines.extend(["", "All intervals are pointwise across seven supplementary endpoints; they do not establish equivalence. Complete-cohort Wilson/paired intervals and explicit unknown identities are in JSON. No complete-case population estimate is reported.", "", "## Source/local joins", ""])
    for j in result["source_local_crosschecks"]:
        lines.append(f"- {j['source_run_id']} / {j['local_run_id']}: planned {j['planned_joins']}; observed both successful {j['counts'].get('both_observed_ok', 0)}; unavailable/duplicate {len(j['missing_or_duplicate'])}; mismatches {len(j['mismatches'])}. Repeated observations add no independent replicates.")
    lines.extend(["", "## Interpretation limits", ""] + ["- " + note for note in result.get("scope_limits", [])])
    return "\n".join(lines) + "\n"


def write_outputs(prefix, result, protocol_path, contract_path, results_root, project):
    paths = [Path(str(prefix) + ext).absolute() for ext in (".json", ".md")]
    protected_roots = [Path(results_root).resolve(), (Path(project) / "configs").resolve()] + [Path(p).resolve() for p in result["protected_snapshot_roots"]]
    protected_files = {Path(p).resolve() for p in result["input_sha256"]} | {Path(protocol_path).resolve(), Path(contract_path).resolve()}
    for path in paths:
        C.require(not path.exists() and not path.is_symlink(), "Refuse existing output: " + str(path))
        C.require(path.resolve() not in protected_files and not any(path.resolve().is_relative_to(root) for root in protected_roots), "Output would modify scientific inputs/source")
    texts = [json.dumps(result, indent=2, allow_nan=False) + "\n", markdown(result)]
    for path, text in zip(paths, texts):
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, prefix=".operator-report-", delete=False) as stream:
            temporary = Path(stream.name)
            try:
                stream.write(text); stream.flush(); os.fsync(stream.fileno())
            except BaseException:
                temporary.unlink(missing_ok=True)
                raise
        try:
            # Hard-link publication atomically refuses an output created meanwhile.
            os.link(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
    return [str(p) for p in paths]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--project-root", type=Path, default=C.PROJECT)
    parser.add_argument("--snapshots-root", type=Path)
    parser.add_argument("--output-prefix", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        result = analyze(args.protocol, args.results_root, args.contract, args.project_root, args.snapshots_root)
        outputs = write_outputs(args.output_prefix, result, args.protocol, args.contract, args.results_root, args.project_root)
        print(json.dumps({"status": result["status"], "ready_for_final_analysis": result["ready_for_final_analysis"], "outputs": outputs}))
        return 0
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        parser.exit(2, f"operator comparison error: {exc}\n")


if __name__ == "__main__":
    raise SystemExit(main())
