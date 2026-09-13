#!/usr/bin/env python3
"""Core selected ties and supplementary paired engine contrasts, JSON/MD only.

Output contract dated 2026-09-08. The frozen secondary endpoint is mean p-shift;
the fixed full-p sample and >= tie convention also permit descriptive ties.
The added repair contrasts are explicitly supplementary, not retrospectively
declared primary/preregistered contrasts. No simulation, fitting or HPC access.
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
import analyze_confirmation_primary as P

VERSION = "confirmation-ties-repair-v1"
CONTRACT_DATE = "2026-09-08"
RUN_ID = "confirm_core_a"
ENGINES = ("original", "repaired")
SCHEMES = ("unrestricted", "diagnosis")
TIE_METRICS = tuple(A.PAIRS) + ("global_max", "within_dataset_mean_pair")
LIMITATIONS = [
    "Output contract dated 2026-09-08; no confirmation effects were read to define this contract.",
    "The frozen named secondary endpoint is mean p-shift. Ties use its frozen selected100 and >= convention; their descriptive summaries are not relabeled as a frozen primary inferential endpoint.",
    "Original minus repaired rejection contrasts and the difference of U-minus-D gaps are supplementary output, not retrospectively declared exact preregistered contrasts.",
    "Ties are E/599, without the p-value +1. Only all599 successful permutations yield a usable tie fraction; early exact decisions alone do not.",
    "Core tie differences are diagnosis minus unrestricted. Repair contrasts are original minus repaired; gap contrast is (U-D)_original minus (U-D)_repaired.",
    "All planned100 or1000 dataset units are retained. Primary means, distribution summaries and Monte Carlo intervals require that endpoint's complete planned cohort.",
    "Unknown completion bounds are not confidence intervals; failed fits are never nonrejections. Collection complete and endpoint complete are different.",
    "Resampling uses whole matched datasets. Intervals are pointwise and do not establish equivalence, universal ordering-only validity or a causal tie explanation.",
    "Same-input source/truth checks establish pairing; compact rows audit does not inspect every raw permutation fit or independently reconstruct data.",
]


def _summary(values, bounds, seeds, reasons, *, support, context, repeats, bootstrap_seed,
             distribution=False):
    """Fixed-denominator mean and bounds, never a mean of early completers."""
    n = len(seeds)
    if not n or not (len(values) == len(bounds) == len(reasons) == n):
        raise ValueError("Endpoint arrays must match a nonempty planned cohort")
    lo, hi = support
    for value, (lower, upper) in zip(values, bounds):
        if not all(math.isfinite(v) for v in (lower, upper)) or not lo <= lower <= upper <= hi:
            raise ValueError("Invalid endpoint completion bounds")
        if value is not None and (not math.isfinite(value) or not lo <= value <= hi or lower != value or upper != value):
            raise ValueError("Known value must have equal finite bounds")
    known = [float(value) for value in values if value is not None]
    complete = len(known) == n
    estimate = math.fsum(known)/n if complete else None
    degenerate = complete and min(known) == max(known)
    radius = (hi-lo)*math.sqrt(math.log(40)/(2*n)) if complete else None
    ordered = sorted(known) if complete else None
    return {"planned_datasets": n, "known_datasets": len(known), "unknown_datasets": n-len(known),
            "endpoint_complete": complete, "estimate": estimate,
            "all_planned_mean_bounds": [math.fsum(b[i] for b in bounds)/n for i in (0, 1)],
            "bounds_kind": "allowable-completion bounds, not a confidence interval", "support": list(support),
            "mcse_mean": statistics.stdev(known)/math.sqrt(n) if complete and n > 1 else None,
            "bootstrap_ci95": A.bootstrap(known, lambda x: math.fsum(x)/len(x), repeats, bootstrap_seed, context)
                if complete and not degenerate else None,
            "bootstrap_degenerate_sample": degenerate,
            "hoeffding_ci95": [max(lo, estimate-radius), min(hi, estimate+radius)] if complete else None,
            "confidence_scope": "pointwise dataset Monte Carlo mean uncertainty; no equivalence or simultaneous coverage",
            "distribution": {"minimum": ordered[0], "maximum": ordered[-1],
                "quantiles": {str(q): A.percentile(ordered, q) for q in (.05, .25, .5, .75, .95)},
                "quantile_definition": "linear interpolation at (R-1)*q; descriptive, not confidence limits"}
                if complete and distribution else None,
            "unknown_reason_counts": dict(Counter(reason or "unknown" for value, reason in zip(values, reasons) if value is None)),
            "unknown_seeds": [seed for seed, value in zip(seeds, values) if value is None]}


def _linear(parts, coefficients):
    """A per-dataset linear contrast, with sharp bounds for missing components."""
    values, bounds = [], []
    for items in zip(*parts):
        values.append(math.fsum(c*value for c, value in zip(coefficients, items)) if all(v is not None for v in items) else None)
        lower = math.fsum(c*value if value is not None else min(0, c) for c, value in zip(coefficients, items))
        upper = math.fsum(c*value if value is not None else max(0, c) for c, value in zip(coefficients, items))
        bounds.append((lower, upper))
    return values, bounds


def _validate_config(cfg, *, frozen=False):
    if (cfg.get("run_id") != RUN_ID or cfg.get("bperm") != 599 or cfg.get("alpha") != .05
            or cfg.get("rule") != "le" or cfg.get("engines") != list(ENGINES)
            or set(cfg.get("schemes", {})) != set(SCHEMES)):
        raise ValueError("Only the paired original/repaired core B599 le plan is supported")
    for scheme in SCHEMES:
        d = cfg["schemes"][scheme]
        if (d.get("stratify") != ("none" if scheme == "unrestricted" else scheme)
                or d.get("groups", [0, 1, 2]) != [0, 1, 2] or d.get("require_max") is not True
                or d.get("tested_pair_index") is not None or d.get("family") is not None):
            raise ValueError("Only joint three-group U/D scheme definitions are supported")
    n, base = A.integer(cfg["datasets"], "datasets"), A.integer(cfg["base_seed"], "base_seed")
    seeds = list(range(base, base+n))
    selected = sorted(P._selected(cfg))
    cells = [P._cell(c) for c in cfg["cells"]]
    if not n or not selected or not set(selected) <= set(seeds) or len(cells) != len(set(cells)):
        raise ValueError("Invalid planned cells or full-p cohort")
    if frozen and (n != 1000 or base != 42100000 or selected != seeds[:100]
                   or cells != ["IID_H0", "REF_H0"] or cfg.get("phase") != "confirmation"):
        raise ValueError("Actual frozen core must retain its two cells, R1000 and selected100")
    return cells, seeds, selected


def _entry(row, scheme, cfg):
    if row is None or row.get("status") != "ok":
        return A.evaluate_scheme(None, .05, "le"), None
    evaluation, entry = P._evaluate(row, scheme, cfg)
    if entry is not None:
        # In addition to count bounds, require exact accounting in pure fixtures;
        # production input has already passed the inventory's stronger checks.
        n = A.integer(entry["nperm"], "nperm")
        failed = A.integer(entry.get("failed", 0), "failed")
        attempted = A.integer(entry.get("attempted", n+failed), "attempted")
        if attempted != n+failed or attempted > 599:
            raise ValueError("Invalid attempted/successful/failed permutation accounting")
    return evaluation, entry


def summarize_run(cfg, rows, *, unavailable_reasons=None, bootstrap_repeats=2000, bootstrap_seed=20260908):
    """Pure summary of already validated rows; small namespace31 toys allowed.

    Production callers must use build_report for frozen plan/source/truth audit.
    No input objects are modified. Successful paired rows must share truth/data.
    """
    cells, seeds, selected = _validate_config(cfg)
    if type(bootstrap_repeats) is not int or bootstrap_repeats <= 0:
        raise ValueError("bootstrap_repeats must be a positive integer")
    indexed, issues = {}, []
    unavailable_reasons = dict(unavailable_reasons or {})
    for row in rows:
        key = (row["cell"], row["engine"], row["seed"])
        if row.get("run_id") != RUN_ID or key[0] not in cells or key[1] not in ENGINES or type(key[2]) is not int or key[2] not in seeds:
            raise ValueError("Unplanned row identity")
        if key in indexed:
            raise ValueError("Duplicate row identity; no winner is selected")
        if row.get("diagnostics", {}).get("full_p_selected") is not (key[2] in selected):
            raise ValueError("Full-p selection differs from the frozen absolute set")
        if row.get("status") == "ok":
            manifest = row.get("truth", {}).get("manifest", {})
            digest = manifest.get("data_sha256")
            if A.null_pairs(row) != [True]*3 or not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
                raise ValueError("Successful core row requires common-order truth and valid data hash")
        indexed[key] = row
    # Same planned seed alone is insufficient for a paired analysis.
    for cell in cells:
        for seed in seeds:
            keys = [(cell, engine, seed) for engine in ENGINES]
            pair = [indexed.get(key) for key in keys]
            if all(row is not None and row.get("status") == "ok" for row in pair) and pair[0]["truth"] != pair[1]["truth"]:
                issues.append({"code": "paired_truth_or_data_mismatch", "cell": cell, "seed": seed})
                for key in keys:
                    indexed[key] = None
                    unavailable_reasons[key] = "paired_truth_or_data_mismatch"
    tie_groups, tie_differences, repair, dataset_ties, counts = [], [], [], [], {}
    for cell in cells:
        outcomes = {(engine, scheme, kind): [] for engine in ENGINES for scheme in SCHEMES for kind in ("bonferroni", "global_max")}
        ties = {(engine, scheme, metric): [] for engine in ENGINES for scheme in SCHEMES for metric in TIE_METRICS}
        reasons = {(engine, scheme): [] for engine in ENGINES for scheme in SCHEMES}
        for engine in ENGINES:
            statuses = Counter()
            for seed in seeds:
                key = (cell, engine, seed)
                row = indexed.get(key)
                reason = (row.get("status") if row and row.get("status") != "ok" else None) if row else unavailable_reasons.get(key, "missing_row")
                statuses[reason or "ok"] += 1
                for scheme in SCHEMES:
                    evaluation, entry = _entry(row, scheme, cfg)
                    outcomes[engine, scheme, "bonferroni"].append(A.any_known(evaluation["pair_reject"]))
                    outcomes[engine, scheme, "global_max"].append(evaluation["max_reject"])
                    if seed not in selected:
                        continue
                    full = entry is not None and evaluation["complete"] and evaluation["pge"] is not None and evaluation["pmax"] is not None
                    values = {metric: None for metric in TIE_METRICS}
                    if full:
                        fractions = [A.integer(e, "eq")/599 for e in entry["eq"]]
                        values.update(zip(A.PAIRS, fractions))
                        values["global_max"] = A.integer(entry["maxeq"], "maxeq")/599
                        values["within_dataset_mean_pair"] = math.fsum(fractions)/3
                    why = None if full else reason or ("missing_scheme" if entry is None else "incomplete_B599_successful_permutations")
                    reasons[engine, scheme].append(why)
                    for metric, value in values.items():
                        ties[engine, scheme, metric].append(value)
                    dataset_ties.append({"cell": cell, "engine": engine, "scheme": scheme, "seed": seed,
                                         "full_B599_available": full, "unknown_reason": why, "tie_fractions": values})
            counts[f"{cell}:{engine}"] = dict(statuses)
        for engine in ENGINES:
            for scheme in SCHEMES:
                metrics = {}
                for metric in TIE_METRICS:
                    values = ties[engine, scheme, metric]
                    metrics[metric] = _summary(values, [(v, v) if v is not None else (0., 1.) for v in values], selected,
                        reasons[engine, scheme], support=(0., 1.), context=f"{VERSION}:{cell}:{engine}:{scheme}:{metric}",
                        repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed, distribution=True)
                tie_groups.append({"cell": cell, "engine": engine, "scheme": scheme, "metrics": metrics})
            metrics = {}
            for metric in TIE_METRICS:
                values, bounds = _linear([ties[engine, "diagnosis", metric], ties[engine, "unrestricted", metric]], [1, -1])
                missing = [";".join(filter(None, r)) or None for r in zip(reasons[engine, "diagnosis"], reasons[engine, "unrestricted"])]
                metrics[metric] = _summary(values, bounds, selected, missing, support=(-1., 1.),
                    context=f"{VERSION}:{cell}:{engine}:D-U:{metric}", repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed,
                    distribution=True)
            tie_differences.append({"cell": cell, "engine": engine, "direction": "diagnosis minus unrestricted", "metrics": metrics})
        for kind in ("bonferroni", "global_max"):
            comparisons = [(scheme, [("original", scheme, kind), ("repaired", scheme, kind)], [1, -1], (-1., 1.)) for scheme in SCHEMES]
            comparisons.append(("U_minus_D_gap", [("original", "unrestricted", kind), ("original", "diagnosis", kind),
                                                   ("repaired", "unrestricted", kind), ("repaired", "diagnosis", kind)], [1, -1, -1, 1], (-2., 2.)))
            for name, keys, coefficients, support in comparisons:
                parts = [outcomes[key] for key in keys]
                values, bounds = _linear(parts, coefficients)
                result = _summary(values, bounds, seeds, ["one_or_more_paired_decisions_unknown" if v is None else None for v in values],
                    support=support, context=f"{VERSION}:{cell}:{kind}:{name}", repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed)
                result["observed_joint_decision_counts"] = dict(Counter("".join("?" if v is None else str(int(v)) for v in items) for items in zip(*parts)))
                repair.append({"cell": cell, "endpoint": kind, "comparison": name,
                    "direction": "(U-D)_original minus (U-D)_repaired" if name == "U_minus_D_gap" else "original minus repaired",
                    "component_order": [list(key) for key in keys], "coefficients": coefficients,
                    "role": "supplementary_contract_2026_09_08_not_exact_preregistered_contrast", "summary": result})
    endpoints = [metric for group in tie_groups+tie_differences for metric in group["metrics"].values()] + [item["summary"] for item in repair]
    return {"schema_version": VERSION, "output_contract_date": CONTRACT_DATE, "run_id": RUN_ID,
            "planned_datasets_per_cell": len(seeds), "planned_tie_datasets_per_cell": len(selected),
            "planned_rows": len(seeds)*len(cells)*len(ENGINES), "input_rows": len(rows),
            "selected_absolute_seeds": selected, "unselected_naturally_complete_rows_excluded_from_ties": True,
            "bootstrap_seed": bootstrap_seed, "bootstrap_repeats": bootstrap_repeats, "bootstrap_unit": "whole matched dataset seed",
            "status_or_unavailable_reason_counts": counts, "pairing_integrity_issues": issues,
            "endpoint_complete": all(e["endpoint_complete"] for e in endpoints),
            "endpoint_completion_scope": "every reported selected-tie and supplementary rejection contrast",
            "ties_by_scheme": tie_groups, "tie_differences": tie_differences, "repair_rejection_contrasts": repair,
            "selected_dataset_ties": dataset_ties, "limitations": LIMITATIONS}


def build_report(protocol_path, results_root, *, project=None, snapshots_root=None, bootstrap_repeats=2000, bootstrap_seed=20260908):
    import analyze_confirmation_cohorts as C
    project = Path(project or C.PROJECT).resolve()
    before = C.sha(protocol_path)
    protocol = json.loads(Path(protocol_path).read_text())
    entries = {e["run_id"]: e for e in protocol["runs"]}
    if len(entries) != len(protocol["runs"]) or RUN_ID not in entries:
        raise ValueError("Protocol must register exactly one core run")
    entry = entries[RUN_ID]
    run = C.load_run(C.frozen_run(entry, project, snapshots_root), results_root, compact=True)
    _validate_config(run["cfg"], frozen=True)
    rows, reasons = [], {}
    for cell in run["names"]:
        for engine in ENGINES:
            for seed in sorted(run["seeds"]):
                row, reason = C.get_record(run, cell, engine, seed)
                if row is not None:
                    rows.append(row)
                else:
                    reasons[cell, engine, seed] = reason
    result = summarize_run(run["cfg"], rows, unavailable_reasons=reasons,
                           bootstrap_repeats=bootstrap_repeats, bootstrap_seed=bootstrap_seed)
    integrity = run["trusted"] and not result["pairing_integrity_issues"] and all(i.get("code") == "failed_row_truth_unavailable_or_invalid" for i in run["issues"])
    if C.sha(protocol_path) != before:
        raise ValueError("Protocol changed while building report")
    result.update(created_utc=datetime.now(timezone.utc).isoformat(), protocol_sha256=before,
        postprocessor_sha256=C.sha(__file__), primary_helper_sha256=C.sha(P.__file__),
        base_analyzer_sha256=C.sha(A.__file__), cohort_loader_sha256=C.sha(C.__file__), inventory_sha256=C.sha(C.I.__file__),
        fitting_or_generation_performed=False, collection_complete=run["inventory"].get("collection_complete", False),
        integrity_valid=integrity, provenance={"snapshot": str(run["snapshot"]), "snapshot_sha256": entry["snapshot_sha256"],
            "config_sha256": entry["config_sha256"], "environment": run["environment"], "issues": run["issues"],
            "trusted_run": run["trusted"], "compact_rows": True, "raw_fit_events_read": False})
    result["status"] = "integrity_invalid" if not integrity else "collecting" if not result["collection_complete"] else "complete" if result["endpoint_complete"] else "terminal_with_unknown"
    return result


def render_markdown(result):
    lines = ["# Core selected ties and supplementary repair contrasts", "", f"Output contract: {CONTRACT_DATE}; status: **{result.get('status', 'collection_not_assessed')}**.",
             f"Collection complete: {result.get('collection_complete', 'not assessed')}; all reported endpoints complete: {result['endpoint_complete']}.",
             "Terminal-with-unknown means completed accounting with unavailable endpoints; no complete-case primary estimates are substituted.", "",
             *[f"- {note}" for note in LIMITATIONS], "", "## Selected tie fractions", "",
             "| Cell / engine / scheme | Metric | Known / planned | Mean or completion bounds | Descriptive quartiles |", "|---|---|---:|---|---|"]
    for group in result["ties_by_scheme"]+result["tie_differences"]:
        for metric, value in group["metrics"].items():
            estimate = value["estimate"] if value["estimate"] is not None else value["all_planned_mean_bounds"]
            quantiles = value["distribution"]["quantiles"] if value["distribution"] is not None else None
            lines.append(f"| {group['cell']} / {group['engine']} / {group.get('scheme', 'D minus U')} | {metric} | {value['known_datasets']}/{value['planned_datasets']} | {estimate} | {quantiles} |")
    lines += ["", "## Supplementary paired rejection contrasts", "", "| Cell / endpoint / comparison | Direction | Known / planned | Difference or completion bounds | Pointwise Hoeffding CI |", "|---|---|---:|---|---|"]
    for item in result["repair_rejection_contrasts"]:
        value = item["summary"]
        estimate = value["estimate"] if value["estimate"] is not None else value["all_planned_mean_bounds"]
        lines.append(f"| {item['cell']} / {item['endpoint']} / {item['comparison']} | {item['direction']} | {value['known_datasets']}/{value['planned_datasets']} | {estimate} | {value['hoeffding_ci95']} |")
    return "\n".join(lines)+"\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--results-root", type=Path, required=True)
    parser.add_argument("--project", type=Path)
    parser.add_argument("--snapshots-root", type=Path)
    parser.add_argument("--bootstrap-repeats", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260908)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args(argv)
    import analyze_stage_exchangeability as S
    def protect():
        S.validate_output_paths(args.protocol, args.results_root, [args.output_json, args.output_md], workspace=args.project)
        project = Path(args.project or Path(__file__).resolve().parents[2]).resolve()
        snapshot_roots = [(project/"runs/v2/snapshots").resolve(), (project/"snapshots").resolve()]
        if args.snapshots_root:
            snapshot_roots.append(args.snapshots_root.resolve())
        if any(p is not None and any(p.resolve().is_relative_to(root) for root in snapshot_roots)
               for p in (args.output_json, args.output_md)):
            raise ValueError("Output cannot overwrite project or relocated snapshots")
    protect()
    result = build_report(args.protocol, args.results_root, project=args.project, snapshots_root=args.snapshots_root,
                          bootstrap_repeats=args.bootstrap_repeats, bootstrap_seed=args.bootstrap_seed)
    protect()
    S.write_outputs(args.protocol, args.results_root,
        [(args.output_json, json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False)+"\n"),
         (args.output_md, render_markdown(result))], workspace=args.project)


if __name__ == "__main__":
    main()
