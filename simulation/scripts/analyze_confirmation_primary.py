"""Explicit-direction, fixed-denominator core/stage confirmation postprocessing.

No generation or fitting. summarize_run is a pure function over validated rows
and the frozen config. The CLI reuses the separate confirmation-cohort provenance
loader; it never treats the subset of discovered rows as the planned cohort.
The frozen analyze_v2.py and its backward-compatible summaries are unchanged.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import statistics

import analyze_v2 as A

VERSION = "confirmation-primary-v1"
SUPPORTED = {"confirm_core_a": ("unrestricted", "diagnosis"),
             "confirm_stage_a": ("diagnosis", "oracle_dx_stage")}
LIMITATIONS = [
    "Core primary Bonferroni FWER contrast is U minus D; its p-shift is D minus U. These signs are intentionally different.",
    "Stage contrast is oracle(DX,k) minus diagnosis: a negative rejection difference means fewer oracle rejections. Oracle stage is unavailable in practice.",
    "Top-level estimates/Monte Carlo intervals require the complete prespecified endpoint cohort. Incomplete-cohort estimates are explicitly nested as complete_case only.",
    "All-planned bounds retain failed, missing and unfinished outcomes; these completion bounds are not confidence intervals.",
    "Collection completion and endpoint completion are separate: terminal_with_unknown means every planned computation is terminal but at least one reported endpoint remains unevaluable. It does not convert failures to nonrejections or permit complete-case primary estimates.",
    "Each dataset seed is one Monte Carlo unit. Three pair p-differences are averaged before between-dataset inference; they are not three independent replicates.",
    "Global-max rejection and p-values are separate from the three-pair Bonferroni family. Positive p-shift does not establish calibration.",
    "Stage full-p10 is a prespecified execution diagnostic, with no population p-shift/distribution claim. Core full-p100 supports only its prespecified mean endpoint at the achieved precision.",
    "Pointwise bootstrap intervals and finite sampled outcomes do not prove equivalence or universal ordering-only FWER control.",
    "A constant observed difference sample has no reported primary percentile interval. Complete cohorts also receive a conservative pointwise Hoeffding interval using the prespecified [-1,1] range and independent dataset units.",
]


def _cell(spec):
    return spec.get("id", spec["name"])


def _selected(config):
    first = A.integer(config.get("full_p_first_n", 0), "full_p_first_n")
    if first > config["datasets"]:
        raise ValueError("full-p first_n exceeds planned datasets")
    return set(range(config["base_seed"], config["base_seed"]+first)) | set(config.get("full_p_seeds", []))


def _rate(values):
    summary = A.binomial_summary(values)
    complete = summary["n_unevaluable"] == 0
    return {"planned_datasets": len(values), "known_datasets": summary["n_evaluable"],
            "unknown_datasets": summary["n_unevaluable"], "reject_count": summary["reject_count"],
            "complete": complete, "endpoint_complete": complete, "estimate": summary["rate"] if complete else None,
            "mcse": summary["mcse"] if complete else None,
            "wilson_ci95": summary["wilson_ci95"] if complete else None,
            "all_planned_rate_bounds": summary["all_attempted_rate_bounds"],
            "complete_case": {"estimate": summary["rate"], "mcse": summary["mcse"],
                              "wilson_ci95": summary["wilson_ci95"],
                              "scope": "evaluable subset; not the planned-cohort estimate when any outcome is unknown"}}


def _mean(values, bounds, *, direction, context, repeats, bootstrap_seed):
    summary = A.continuous_summary(values, repeats, bootstrap_seed, context)
    n, known = len(values), summary["n_datasets"]
    complete = n > 0 and n == known
    radius = math.sqrt(2*math.log(40)/n) if complete else None
    return {"direction": direction, "planned_datasets": n, "known_datasets": known,
            "unknown_datasets": n-known, "complete": complete, "endpoint_complete": complete,
            "estimate": summary["estimate"] if complete else None,
            "mcse_mean": summary["mcse_mean"] if complete else None,
            "bootstrap_ci95": summary["bootstrap_ci95"] if complete and not summary["bootstrap_degenerate_sample"] else None,
            "bootstrap_degenerate_sample": summary["bootstrap_degenerate_sample"],
            "interval_caution": summary["interval_caution"],
            "hoeffding_ci95": [max(-1., summary["estimate"]-radius), min(1., summary["estimate"]+radius)] if complete else None,
            "hoeffding_scope": "conservative pointwise finite-sample interval; range [-1,1], independent dataset units, complete cohort only",
            "bootstrap_unit": "dataset seed",
            "all_planned_mean_bounds": [math.fsum(b[i] for b in bounds)/n for i in (0, 1)] if n else None,
            "bounds_kind": "bounds over allowable completions, not a confidence interval",
            "complete_case": summary | {"scope": "known matched-dataset subset only; not the planned-cohort estimate when incomplete"}}


def _difference(positive, negative, *, positive_name, negative_name, context, repeats, bootstrap_seed):
    values, bounds = [], []
    for pos, neg in zip(positive, negative):
        values.append(int(pos)-int(neg) if pos is not None and neg is not None else None)
        lo_pos, hi_pos = (0, 1) if pos is None else (int(pos), int(pos))
        lo_neg, hi_neg = (0, 1) if neg is None else (int(neg), int(neg))
        bounds.append((lo_pos-hi_neg, hi_pos-lo_neg))
    result = _mean(values, bounds, direction=f"{positive_name} minus {negative_name}",
                   context=context, repeats=repeats, bootstrap_seed=bootstrap_seed)
    result.update(paired_table=dict(Counter(f"{int(p)}{int(n)}" for p, n in zip(positive, negative)
                                           if p is not None and n is not None)),
                  paired_table_order=[positive_name, negative_name])
    return result


def _evaluate(row, name, config):
    # A present invalid/pending/failed row must never supply stale positive or
    # negative decisions, even if it accidentally retains old scheme fields.
    if row is None or row.get("status") != "ok":
        return A.evaluate_scheme(None, config["alpha"], config["rule"]), None
    entry = row.get("schemes", {}).get(name)
    if entry is None:
        return A.evaluate_scheme(None, config["alpha"], config["rule"]), None
    if entry.get("definition") != config["schemes"][name]:
        raise ValueError(f"{name} definition differs from frozen config")
    if entry.get("requested_nperm") != config["bperm"]:
        raise ValueError(f"{name} budget differs from frozen B")
    # Certification comes from counts and the fixed budget. A saved decision
    # alone cannot certify an otherwise unresolved completion interval.
    counts = {k: v for k, v in entry.items() if k != "exact_decision"}
    return A.evaluate_scheme(counts, config["alpha"], config["rule"]), counts


def _p_interval(evaluation, entry, index):
    values = evaluation["pmax"] if index is None else evaluation["pge"]
    value = values if index is None else (values[index] if values is not None else None)
    if value is not None:
        return (value, value)
    if entry is None:
        return (0.0, 1.0)
    if index is None:
        if "maxgt" not in entry or "maxeq" not in entry:
            return (0.0, 1.0)
        count = A.integer(entry["maxgt"], "maxgt")+A.integer(entry["maxeq"], "maxeq")
    else:
        if entry.get("gt") is None or entry.get("eq") is None:
            return (0.0, 1.0)
        count = A.integer(entry["gt"][index], "gt")+A.integer(entry["eq"][index], "eq")
    budget, n = evaluation["budget"], evaluation["nperm"]
    return ((1+count)/(1+budget), (1+count+budget-n)/(1+budget))


def summarize_run(config, rows, *, unavailable_reasons=None, bootstrap_repeats=2000, bootstrap_seed=20260908):
    """Summarize every configured cell/engine/seed; rows may be incomplete.

    Input provenance is a caller obligation for this pure function. CLI inputs
    pass through the confirmation-cohort loader before arriving here. No input
    row or frozen configuration is mutated.
    """
    run_id = config["run_id"]
    if run_id not in SUPPORTED or config.get("rule") != "le" or config.get("alpha") != .05 or config.get("bperm") != 599:
        raise ValueError("Supported frozen primary runs require B599, alpha=.05 and le")
    scheme_a, scheme_b = SUPPORTED[run_id]
    if set(config["schemes"]) != {scheme_a, scheme_b}:
        raise ValueError("Unexpected primary scheme set")
    for name, definition in config["schemes"].items():
        expected = "none" if name == "unrestricted" else name
        if (definition.get("stratify") != expected or definition.get("groups", [0, 1, 2]) != [0, 1, 2]
                or not definition.get("require_max") or definition.get("tested_pair_index") is not None
                or definition.get("family") is not None):
            raise ValueError("Primary contrasts require their joint-label frozen scheme definitions")
    total = A.integer(config["datasets"], "datasets")
    base = A.integer(config["base_seed"], "base_seed")
    if total == 0:
        raise ValueError("The planned cohort must be nonempty")
    seeds = list(range(base, base+total))
    full = _selected(config)
    if not full <= set(seeds):
        raise ValueError("Unplanned full-p seed")
    cells, engines = [_cell(c) for c in config["cells"]], config["engines"]
    if len(set(cells)) != len(cells) or len(set(engines)) != len(engines):
        raise ValueError("Repeated configured cell or engine")
    indexed = {}
    for row in rows:
        key = (row["cell"], row["engine"], row["seed"])
        if (row["run_id"] != run_id or key[0] not in cells or key[1] not in engines
                or type(key[2]) is not int or key[2] not in seeds):
            raise ValueError("Unplanned row identity")
        if key in indexed:
            raise ValueError("Duplicate dataset identity; no winner is selected")
        flag = row.get("diagnostics", {}).get("full_p_selected")
        if type(flag) is not bool or flag != (key[2] in full):
            raise ValueError("Recorded full-p selection differs from the frozen absolute seed set")
        null = A.null_pairs(row)
        if null is not None and not all(null):
            raise ValueError("Core/stage H0 truth mask is not a complete common-order null")
        indexed[key] = row
    unavailable_reasons = unavailable_reasons or {}
    groups = []
    for cell in cells:
        for engine in engines:
            context = f"{VERSION}:{run_id}:{cell}:{engine}"
            selected_values = {"mean_pair_p_shift": [], **{f"pair_p_shift_{p}": [] for p in A.PAIRS}, "max_p_shift": []}
            selected_bounds = {name: [] for name in selected_values}
            outcomes = {name: {"bonferroni": [], "global_max": []} for name in (scheme_a, scheme_b)}
            statuses = Counter()
            for seed in seeds:
                key = (cell, engine, seed)
                row = indexed.get(key)
                statuses[row.get("status", "unknown") if row else unavailable_reasons.get(key, "missing_row")] += 1
                evaluated = {name: _evaluate(row, name, config) for name in (scheme_a, scheme_b)}
                for name, (evaluation, _) in evaluated.items():
                    outcomes[name]["bonferroni"].append(A.any_known(evaluation["pair_reject"]))
                    outcomes[name]["global_max"].append(evaluation["max_reject"])
                if seed not in full:
                    continue
                av, ae = evaluated[scheme_a]
                bv, be = evaluated[scheme_b]
                pair_values, pair_bounds = [], []
                for index, pair in enumerate(A.PAIRS):
                    value = bv["pge"][index]-av["pge"][index] if av["pge"] is not None and bv["pge"] is not None else None
                    al, ah = _p_interval(av, ae, index)
                    bl, bh = _p_interval(bv, be, index)
                    pair_values.append(value)
                    pair_bounds.append((bl-ah, bh-al))
                    selected_values[f"pair_p_shift_{pair}"].append(value)
                    selected_bounds[f"pair_p_shift_{pair}"].append(pair_bounds[-1])
                selected_values["mean_pair_p_shift"].append(statistics.mean(pair_values) if all(v is not None for v in pair_values) else None)
                selected_bounds["mean_pair_p_shift"].append(tuple(statistics.mean(b[i] for b in pair_bounds) for i in (0, 1)))
                selected_values["max_p_shift"].append(bv["pmax"]-av["pmax"] if av["pmax"] is not None and bv["pmax"] is not None else None)
                al, ah = _p_interval(av, ae, None)
                bl, bh = _p_interval(bv, be, None)
                selected_bounds["max_p_shift"].append((bl-ah, bh-al))
            positive, negative = (scheme_a, scheme_b) if run_id == "confirm_core_a" else (scheme_b, scheme_a)
            differences = {kind: _difference(outcomes[positive][kind], outcomes[negative][kind],
                positive_name=positive, negative_name=negative, context=context+":"+kind,
                repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed) for kind in ("bonferroni", "global_max")}
            p_shifts = {name: _mean(selected_values[name], selected_bounds[name],
                direction=f"{scheme_b} minus {scheme_a}", context=context+":"+name,
                repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed) for name in selected_values}
            p_shifts["mean_pair_p_shift"]["estimand"] = "mean across datasets of the within-dataset mean across three pair p-differences"
            groups.append({"run_id": run_id, "cell": cell, "engine": engine,
                "role": "primary_REF_repaired" if run_id == "confirm_core_a" and cell == "REF_H0" and engine == "repaired"
                        else "stage_boundary_primary" if run_id == "confirm_stage_a" else "control_or_engine_sensitivity",
                "planned_datasets": total, "status_or_unavailable_reason_counts": dict(statuses),
                "primary_bonferroni_rejection_difference": differences["bonferroni"],
                "scheme_bonferroni_rejection_rates": {name: _rate(value["bonferroni"]) for name, value in outcomes.items()},
                "global_max_separate": {"difference": differences["global_max"],
                                        "rates": {name: _rate(value["global_max"]) for name, value in outcomes.items()}},
                "preselected_p": {"absolute_seeds": sorted(full), "planned_datasets": len(full),
                    "role": "prespecified_core_mean_endpoint" if run_id == "confirm_core_a" else "execution_diagnostic_only_no_population_p_shift_claim",
                    "unselected_complete_datasets_excluded": True, "metrics": p_shifts}})
    for group in groups:
        group["endpoint_completion"] = {
            "primary_bonferroni_rejection_difference": group["primary_bonferroni_rejection_difference"]["endpoint_complete"],
            "global_max_difference": group["global_max_separate"]["difference"]["endpoint_complete"],
            "preselected_mean_pair_p_shift": group["preselected_p"]["metrics"]["mean_pair_p_shift"]["endpoint_complete"],
            "preselected_max_p_shift": group["preselected_p"]["metrics"]["max_p_shift"]["endpoint_complete"],
        }
        group["endpoint_complete"] = all(group["endpoint_completion"].values())
    return {"schema_version": VERSION, "run_id": run_id,
            "configuration_sha256_canonical": hashlib.sha256(json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "input_row_count": len(indexed), "planned_row_count": total*len(cells)*len(engines),
            "endpoint_complete": bool(groups) and all(group["endpoint_complete"] for group in groups),
            "endpoint_completion_scope": "all reported primary, secondary and diagnostic endpoints; inspect each endpoint separately",
            "statistical_unit": "one simulated dataset seed within cell and engine", "groups": groups,
            "legacy_analyzer_unchanged": True, "limitations": LIMITATIONS,
            "source_analyzer_direction_note": "Frozen comparison is D-U for core and oracle-D for stage. This report computes the core primary U-D directly from the same matched decisions."}


def build_report(protocol_path, results_root, *, project=None, snapshots_root=None, run_ids=None, bootstrap_repeats=2000):
    # Lazy import keeps the pure statistics interface independent of I/O helpers.
    import analyze_confirmation_cohorts as C
    project = Path(project or Path(__file__).resolve().parents[2]).resolve()
    protocol = json.loads(Path(protocol_path).read_text())
    requested = list(run_ids or SUPPORTED)
    if len(set(requested)) != len(requested) or not set(requested) <= SUPPORTED.keys():
        raise ValueError("Select core and/or stage once each")
    entries = {e["run_id"]: e for e in protocol["runs"]}
    if len(entries) != len(protocol["runs"]):
        raise ValueError("Duplicate protocol run IDs")
    runs, reports = [], []
    for run_id in requested:
        run = C.load_run(C.frozen_run(entries[run_id], project, snapshots_root), results_root, compact=True)
        rows, reasons = [], {}
        for cell in run["names"]:
            for engine in run["cfg"]["engines"]:
                for seed in sorted(run["seeds"]):
                    row, reason = C.get_record(run, cell, engine, seed)
                    if row is not None:
                        rows.append(row)
                    else:
                        reasons[(cell, engine, seed)] = reason
        summary = summarize_run(run["cfg"], rows, unavailable_reasons=reasons, bootstrap_repeats=bootstrap_repeats)
        summary["provenance"] = {"snapshot": str(run["snapshot"]), "snapshot_sha256": entries[run_id]["snapshot_sha256"],
            "config_sha256": entries[run_id]["config_sha256"], "trusted_run": run["trusted"],
            "environment": run["environment"], "issues": run["issues"],
            "inventory_ready_for_final_analysis": run["inventory"]["ready_for_final_analysis"],
            "inventory_collection_complete": run["inventory"].get("collection_complete", False)}
        summary["collection_complete"] = run["inventory"].get("collection_complete", False)
        integrity_valid = run["trusted"] and all(issue.get("code") == "failed_row_truth_unavailable_or_invalid"
                                                  for issue in run["issues"])
        summary["integrity_valid"] = integrity_valid
        summary["status"] = ("integrity_invalid" if not integrity_valid else "collecting" if not summary["collection_complete"]
                             else "complete" if summary["endpoint_complete"] else "terminal_with_unknown")
        reports.append(summary)
        runs.append(run)
    environment_status = ("not_applicable_single_run" if len(runs) < 2 else
        "unavailable_profile" if any(r["environment"] is None for r in runs) else
        "match" if len({json.dumps(r["environment"], sort_keys=True) for r in runs}) == 1 else "mismatch")
    return {"schema_version": VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
            "protocol_sha256": C.sha(protocol_path), "postprocessor_sha256": C.sha(__file__),
            "cohort_loader_sha256": C.sha(C.__file__), "analyzer_sha256": C.sha(A.__file__),
            "runs": reports, "fitting_or_generation_performed": False,
            "cross_run_environment_status": environment_status,
            "limitations": LIMITATIONS}


def render_markdown(report):
    descriptions = {"collecting": "Results are still missing or unfinished.",
                    "terminal_with_unknown": "Collection is finished; unknown endpoints retain bounds and have no complete-cohort estimate.",
                    "complete": "Collection is finished and every reported endpoint is evaluable.",
                    "integrity_invalid": "Source validation failed; this is not a successful terminal analysis.",
                    "collection_not_assessed": "Pure statistical summary; collection state was not assessed."}
    lines = ["# Confirmation primary endpoints", "", *[f"- {note}" for note in LIMITATIONS], "",
             "| Run | Collection complete | All reported endpoints complete | Status | Interpretation |",
             "|---|---|---|---|---|"]
    for run in report["runs"]:
        status = run.get("status", "collection_not_assessed")
        lines.append(f"| {run['run_id']} | {run.get('collection_complete', 'not assessed')} | "
                     f"{run['endpoint_complete']} | {status} | {descriptions[status]} |")
    lines += ["", "All-reported endpoint completion includes secondary and diagnostic p summaries; the per-endpoint columns below show which quantities remain unknown.", "",
              "| Run / cell / engine | Rejection direction | Complete pairs / planned | Rejection endpoint complete | Complete-cohort difference | Planned completion bounds | p-shift complete / selected | p endpoint complete |",
              "|---|---|---:|---|---:|---|---:|---|"]
    for run in report["runs"]:
        for group in run["groups"]:
            effect = group["primary_bonferroni_rejection_difference"]
            p = group["preselected_p"]["metrics"]["mean_pair_p_shift"]
            estimate = "unavailable" if effect["estimate"] is None else f"{effect['estimate']:.6g}"
            lines.append(f"| {run['run_id']} / {group['cell']} / {group['engine']} | {effect['direction']} | "
                f"{effect['known_datasets']}/{effect['planned_datasets']} | {effect['endpoint_complete']} | {estimate} | {effect['all_planned_mean_bounds']} | "
                f"{p['known_datasets']}/{p['planned_datasets']} | {p['endpoint_complete']} |")
    return "\n".join(lines)+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True, help="Parent of explicitly registered run directories")
    parser.add_argument("--project", type=Path)
    parser.add_argument("--snapshots-root", type=Path)
    parser.add_argument("--run-id", choices=list(SUPPORTED), action="append")
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()
    import analyze_stage_exchangeability as S
    S.validate_output_paths(args.protocol, args.results_root, [args.output_json, args.output_md], workspace=args.project)
    if args.snapshots_root and any(p is not None and p.resolve().is_relative_to(args.snapshots_root.resolve())
                                  for p in (args.output_json, args.output_md)):
        raise ValueError("Output cannot overwrite snapshots")
    report = build_report(args.protocol, args.results_root, project=args.project, snapshots_root=args.snapshots_root,
                          run_ids=args.run_id, bootstrap_repeats=args.bootstrap_repeats)
    for path, content in [(args.output_json, json.dumps(report, indent=2, allow_nan=False)+"\n"),
                          (args.output_md, render_markdown(report))]:
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)


if __name__ == "__main__":
    main()
