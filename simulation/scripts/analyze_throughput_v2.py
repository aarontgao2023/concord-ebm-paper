#!/usr/bin/env python3
"""Stream v2 fit logs into a compact, censoring-aware throughput audit.

Requires only the standard library and hpc/v2/inventory_run.py. This is read-only
with respect to raw runs. A completed-fit timing distribution is descriptive:
unfinished runs never yield average savings or a completion-based overall ETA.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "hpc/v2"))
import inventory_run as I


def percentile(values, probability):
    if not values:
        return None
    ordered = sorted(values)
    index = (len(ordered) - 1) * probability
    lo, hi = math.floor(index), math.ceil(index)
    return ordered[lo] + (ordered[hi] - ordered[lo]) * (index - lo)


def duration_summary(values):
    values = list(values)
    return {"n": len(values), "sum": sum(values), "mean": statistics.mean(values) if values else None,
            "p90": percentile(values, .9), "p95": percentile(values, .95), "max": max(values, default=None)}


def finite_seconds(value):
    return type(value) in (int, float) and math.isfinite(value) and value >= 0


def load_logs(root, cfg):
    """Discard each large fit dictionary after keeping identity/status/duration.

    Memory grows with compact identity metadata, not optimizer histories, truth,
    orderings or data. Duplicate identities are excluded as ambiguous, reported,
    and block completion claims; no arbitrary duplicate is selected.
    """
    cells, seeds, _ = I.validate_config(cfg)
    records, duplicates, inputs, issues = {}, set(), [], []
    physical, replay_count = 0, 0
    for path in sorted(Path(root).glob("cell_*/chunk_*/fit_records.jsonl")):
        try:
            cell = cells[int(path.parent.parent.name.removeprefix("cell_"))]
            manifest = json.loads((path.parent / "manifest.json").read_text())
            assigned = set(manifest["seeds"])
            before, digest = path.stat(), hashlib.sha256()
            with path.open("rb") as stream:
                for number, line in enumerate(stream, 1):
                    digest.update(line)
                    if not line.strip():
                        continue
                    physical += 1
                    r = json.loads(line)
                    seed, engine, scheme, pid = (r.get(n) for n in ("seed", "engine", "scheme", "perm_id"))
                    if type(seed) is not int or seed not in seeds or seed not in assigned or engine not in cfg["engines"]:
                        raise ValueError(f"Unplanned fit identity at line {number}")
                    if (scheme is None and pid is not None) or (scheme is not None and
                        (scheme not in cfg.get("schemes", {}) or type(pid) is not int or not 0 <= pid < cfg["bperm"])):
                        raise ValueError(f"Unplanned scheme/pid at line {number}")
                    key = cell, seed, engine, scheme, pid
                    if key in records or key in duplicates:
                        duplicates.add(key)
                        records.pop(key, None)
                        issues.append({"code": "duplicate_fit_identity", "key": list(key), "path": str(path), "line": number})
                        continue
                    seconds = r.get("fit_seconds")
                    if not finite_seconds(seconds):
                        issues.append({"code": "missing_or_invalid_fit_seconds", "key": list(key)})
                        seconds = None
                    shared = r.get("shared_fit_id")
                    if shared is not None and (not isinstance(shared, str) or not shared):
                        raise ValueError(f"Invalid shared_fit_id at line {number}")
                    if cfg.get("paired_standard") and engine in ("original", "repaired") and shared is None:
                        issues.append({"code": "paired_record_missing_shared_id", "key": list(key)})
                    status = r.get("status")
                    if not isinstance(status, str) or not status:
                        raise ValueError(f"Invalid fit status at line {number}")
                    records[key] = {"status": status, "seconds": seconds, "shared": shared,
                                    "timeout": "timeout" in r.get("diagnostics", {}).get("error", "").lower(),
                                    "source": str(path)}
            after = path.stat()
            if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
                issues.append({"code": "fit_log_changed_during_read", "path": str(path)})
            inputs.append({"path": str(path.resolve()), "sha256": digest.hexdigest(), "bytes": after.st_size})
        except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
            issues.append({"code": "fit_log_read_error", "path": str(path), "message": str(exc)})
    for path in sorted(Path(root).glob("cell_*/chunk_*/replay_audit.jsonl")):
        try:
            with path.open() as stream:
                for line in stream:
                    if line.strip():
                        json.loads(line)
                        replay_count += 1
        except (OSError, ValueError) as exc:
            issues.append({"code": "replay_log_read_error", "path": str(path), "message": str(exc)})
    return records, issues, inputs, physical, len(duplicates), replay_count


def load_rows(root):
    rows, duplicates = {}, set()
    for path in sorted(Path(root).glob("cell_*/chunk_*/rows.jsonl")):
        with path.open() as stream:
            for line in stream:
                if not line.strip():
                    continue
                row = json.loads(line)
                key = row["cell"], row["seed"], row["engine"]
                if key in rows or key in duplicates:
                    duplicates.add(key)
                    rows.pop(key, None)
                else:
                    rows[key] = {"status": row["status"], "schemes": row.get("schemes", {})}
    return rows


def group_reports(cfg, rows, records, issues):
    cells, seeds, full = I.validate_config(cfg)
    reports = []
    indexed = defaultdict(dict)
    for key, value in records.items():
        indexed[key[0], key[2], key[3]][key] = value
    for cell in cells:
        for engine in cfg["engines"]:
            for scheme in [None] + list(cfg.get("schemes", {})):
                group_records = indexed[cell, engine, scheme]
                by_seed = defaultdict(list)
                for key, value in group_records.items():
                    by_seed[key[1]].append(value)
                tally, states = Counter(), Counter()
                for seed in sorted(seeds):
                    row = rows.get((cell, seed, engine))
                    saved = by_seed[seed]
                    if scheme is None:
                        tally["started"] += bool(saved)
                        if row is not None:
                            good = row["status"] == "ok"
                            tally["finished"] += good
                            tally["terminal"] += row["status"] not in ("pending", "running")
                            tally["observed_failed"] += not good
                            states[row["status"]] += 1
                            if len(saved) != 1 or (saved and saved[0]["status"] != row["status"]):
                                issues.append({"code": "observed_row_log_mismatch", "key": [cell, seed, engine]})
                        else:
                            states["missing_row"] += 1
                        continue
                    tally["started"] += bool(saved)
                    if row is None or row["status"] != "ok":
                        states["missing_observed" if row is None else "observed_failed"] += 1
                        tally["blocked_observed_failure"] += row is not None
                        continue
                    entry = row.get("schemes", {}).get(scheme)
                    if entry is None:
                        states["missing_scheme"] += 1
                        continue
                    try:
                        state = I.scheme_state(entry, cfg["schemes"][scheme], cfg, seed in full)
                        tally["decision_resolved"] += state["primary_decisions_resolved"]
                        tally["finished"] += state["obligation_satisfied"]
                        tally["terminal"] += state["terminal"]
                        tally["full_budget"] += state["full"]
                        tally["selected_fullp_finished"] += state["full"] and seed in full
                        states[state["state"]] += 1
                        successes = sum(v["status"] == "ok" for v in saved)
                        if len(saved) != state["attempted"] or successes != state["nperm"]:
                            issues.append({"code": "scheme_row_log_count_mismatch", "key": [cell, seed, engine, scheme],
                                           "row_attempted_successes": [state["attempted"], state["nperm"]],
                                           "log_attempted_successes": [len(saved), successes]})
                    except (ValueError, KeyError, TypeError) as exc:
                        issues.append({"code": "invalid_scheme_row", "key": [cell, seed, engine, scheme], "message": str(exc)})
                durations = [v["seconds"] for v in group_records.values() if v["seconds"] is not None]
                reports.append({"cell": cell, "engine": engine, "scheme": scheme or "observed",
                    "planned_datasets": len(seeds), "recorded_started_datasets": tally["started"],
                    "decision_resolved_datasets": tally["decision_resolved"] if scheme is not None else None,
                    "finished_datasets": tally["finished"], "terminal_datasets": tally["terminal"],
                    "blocked_observed_failure_datasets": tally["blocked_observed_failure"],
                    "planned_selected_fullp_datasets": len(full) if scheme is not None else None,
                    "full_budget_datasets": tally["full_budget"] if scheme is not None else None,
                    "selected_fullp_finished_datasets": tally["selected_fullp_finished"] if scheme is not None else None,
                    "planned_full_budget_statistical_fits": len(seeds) * (1 if scheme is None else cfg["bperm"]),
                    "statistical_fit_records": len(group_records),
                    "successful_fit_records": sum(v["status"] == "ok" for v in group_records.values()),
                    "failed_fit_records": sum(v["status"] != "ok" for v in group_records.values()),
                    "timeout_fit_records": sum(v["timeout"] for v in group_records.values()),
                    "status_counts": dict(Counter(v["status"] for v in group_records.values())),
                    "states": dict(states), "fit_seconds_all_recorded_attempts": duration_summary(durations),
                    "unique_jobs_touching_engine_nonadditive": len({v["shared"] or str(k) for k, v in group_records.items()})})
    return reports


def computation_reports(records, issues):
    jobs = defaultdict(list)
    for key, record in records.items():
        # Prefixes prevent a user-supplied shared id colliding with a fallback key.
        job_id = ("shared", record["shared"]) if record["shared"] else ("independent", key)
        jobs[job_id].append((key, record))
    by_group = defaultdict(list)
    conflicts, half_pairs = 0, 0
    for job_id, members in jobs.items():
        units = {(k[0], k[1], k[3], k[4]) for k, _ in members}
        if len(units) != 1:
            issues.append({"code": "shared_id_crosses_dataset_or_scheme", "id": str(job_id)})
            continue
        unit = next(iter(units))
        durations = {v["seconds"] for _, v in members if v["seconds"] is not None}
        if len(durations) > 1:
            conflicts += 1
        if job_id[0] == "shared" and len(members) != 2:
            half_pairs += 1
        if job_id[0] == "shared" and not {k[2] for k, _ in members} <= {"original", "repaired"}:
            issues.append({"code": "shared_job_has_nonstandard_engine", "id": str(job_id)})
        by_group[unit[0], unit[2] or "observed"].append({
            "seconds": max(durations) if durations else None,
            "failed_engine_results": sum(v["status"] != "ok" for _, v in members),
            "has_shared_id": job_id[0] == "shared"})
    groups = []
    for (cell, scheme), data in sorted(by_group.items()):
        groups.append({"cell": cell, "scheme": scheme, "unique_computational_jobs": len(data),
            "shared_jobs": sum(v["has_shared_id"] for v in data),
            "jobs_with_failed_engine_result": sum(v["failed_engine_results"] > 0 for v in data),
            "seconds_counted_once_per_job": duration_summary(v["seconds"] for v in data if v["seconds"] is not None)})
    return {"unique_computational_jobs": len(jobs), "shared_duration_conflicts": conflicts,
            "partially_logged_shared_jobs": half_pairs, "cell_schemes": groups,
            "recorded_unique_job_seconds": sum(g["seconds_counted_once_per_job"]["sum"] for g in groups)}


def time_seconds(value):
    """Parse Slurm duration fields such as 02:03:04.5 or 1-02:03:04."""
    value = value.strip()
    if not value or value in ("Unknown", "None", "N/A"):
        return None
    days, clock = (value.split("-", 1) if "-" in value else ("0", value))
    parts = [float(p) for p in clock.split(":")]
    if not 1 <= len(parts) <= 3:
        raise ValueError(f"Invalid sacct duration: {value}")
    return float(days) * 86400 + sum(p * 60 ** i for i, p in enumerate(reversed(parts)))


def sacct_summary(path):
    if path is None:
        return None
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream, delimiter="|"))
    allocations, steps, seen = [], [], set()
    for row in rows:
        job = row.get("JobIDRaw") or row.get("JobID")
        if not job:
            raise ValueError("sacct input requires a JobIDRaw or JobID header")
        if "." in job:
            steps.append(job)
            continue
        if job in seen:
            raise ValueError(f"Duplicate sacct allocation row {job}")
        seen.add(job)
        allocated = float(row["CPUTimeRAW"]) if row.get("CPUTimeRAW") else None
        measured = time_seconds(row.get("TotalCPU", ""))
        allocations.append({"job_id": job, "state": row.get("State"), "allocated_cpu_seconds": allocated,
                            "measured_total_cpu_seconds": measured})
    allocated = sum(r["allocated_cpu_seconds"] for r in allocations) if allocations and all(r["allocated_cpu_seconds"] is not None for r in allocations) else None
    measured = sum(r["measured_total_cpu_seconds"] for r in allocations) if allocations and all(r["measured_total_cpu_seconds"] is not None for r in allocations) else None
    return {"path": str(Path(path).resolve()), "allocation_rows": allocations, "excluded_step_ids": steps,
            "allocated_cpu_seconds": allocated, "measured_total_cpu_seconds": measured,
            "measured_cpu_utilization": measured / allocated if measured is not None and allocated else None,
            "scope": "User-supplied allocation rows only; never sum parent allocations with batch/extern steps."}


def drain_budget(timeout, paired=False, workers=16, wall=14400, runtime=11800, signal_lead=1500, startup=0):
    jobs = max(2 * workers, 8)
    job_timeout = timeout * (2 if paired else 1)
    drain = math.ceil(jobs / workers) * job_timeout
    trigger = min(startup + runtime, wall - signal_lead)
    return {"fit_timeout_s": timeout, "paired": paired, "workers": workers,
            "max_wave_jobs": jobs, "max_job_timeout_s": job_timeout, "nominal_drain_s": drain,
            "runtime_trigger_with_startup_s": startup + runtime, "signal_trigger_s": wall - signal_lead,
            "earliest_stop_trigger_s": trigger, "remaining_after_nominal_drain_s": wall - trigger - drain,
            "signal_only_remaining_after_drain_s": signal_lead - drain,
            "runtime_only_remaining_after_drain_s": wall - startup - runtime - drain,
            "overhead_excluded": "Serialization, filesystem/fsync, process scheduling and timeout delivery are not hard-bounded here."}


def analyze(config_path, root, sacct=None, forecast_datasets=None, forecast_workers=None):
    cfg = json.loads(Path(config_path).read_text())
    inventory = I.inventory(config_path, root)
    records, issues, inputs, physical, duplicates, replays = load_logs(root, cfg)
    try:
        rows = load_rows(root)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        rows = {}
        issues.append({"code": "row_read_error", "message": str(exc)})
    groups = group_reports(cfg, rows, records, issues)
    computation = computation_reports(records, issues)
    # Both snapshots must describe the same persisted attempt counts. Checkpoints
    # can lag an active append; this blocks a complete-cohort extrapolation.
    planned_statistical = sum(g["planned_full_budget_statistical_fits"] for g in groups)
    paired_reduction = int(bool(cfg.get("paired_standard")) and {"original", "repaired"} <= set(cfg["engines"]))
    planned_jobs = len(cfg["cells"]) * cfg["datasets"] * (len(cfg["engines"]) - paired_reduction) * (1 + len(cfg.get("schemes", {})) * cfg.get("bperm", 0))
    complete = inventory["ready_for_final_analysis"] and not issues and not duplicates and all(
        g["finished_datasets"] == g["planned_datasets"] for g in groups)
    savings = None
    forecast = None
    if complete:
        savings = {"statistical_fit_fraction_avoided": 1 - len(records) / planned_statistical,
                   "unique_logical_computational_job_fraction_avoided": 1 - computation["unique_computational_jobs"] / planned_jobs,
                   "not_a_cpu_time_savings_estimate": True,
                   "failed_attempts_retained_in_numerator": True}
        for group in groups:
            group["complete_cohort_mean_statistical_fit_fraction_avoided"] = 1 - group["statistical_fit_records"] / group["planned_full_budget_statistical_fits"]
        if forecast_datasets is not None and forecast_workers is not None and not replays and not computation["shared_duration_conflicts"]:
            total_seconds = computation["recorded_unique_job_seconds"] * forecast_datasets / cfg["datasets"]
            forecast = {"datasets_per_cell": forecast_datasets, "workers": forecast_workers,
                "same_cells_engines_schemes_bperm_assumed": True,
                "ideal_worker_seconds": total_seconds, "ideal_wall_seconds": total_seconds / forecast_workers,
                "basis": "All planned datasets in the completed run, including failed permutation attempts.",
                "limitations": "Planning scenario, not a guaranteed ETA; excludes queue time, idle allocation, startup/I/O and future distribution/runtime shifts."}
    else:
        for group in groups:
            group["complete_cohort_mean_statistical_fit_fraction_avoided"] = None
    return {"schema_version": "throughput_audit_v1", "created_utc": I.utc(), "run_id": cfg["run_id"],
        "phase": cfg["phase"], "complete_for_throughput_estimation": bool(complete),
        "execution_collection_terminal": inventory["collection_complete"],
        "right_censored_or_incomplete": not complete, "inventory": inventory, "issues": issues,
        "inputs": inputs, "cell_engine_schemes": groups,
        "physical_fit_log_lines": physical, "unambiguous_statistical_fit_records": len(records),
        "duplicate_fit_identities": duplicates, "failed_statistical_fit_records": sum(v["status"] != "ok" for v in records.values()),
        "planned_full_budget_statistical_fits": planned_statistical,
        "planned_full_budget_computational_jobs_with_declared_pairing": planned_jobs,
        "computation": computation, "paired_replay_audit_entries": replays,
        "complete_cohort_savings": savings, "future_run_planning_scenario": forecast,
        "sacct": sacct_summary(sacct),
        "drain_checks_context": "Illustrative historical development scenarios: 16workers, 4h walltime, runtime11800s and signal lead1500s. These are not inferred settings or safety certification for the analyzed run; use its frozen batch script and execution amendments.",
        "drain_checks": [drain_budget(300), drain_budget(300, True), drain_budget(600, True)],
        "notes": ["Started means at least one durable fit record; queued, active and torn/unlogged attempts are not observable from this log.",
            "Finished means successful observed fit or primary decisions plus planned full-p obligation satisfied; terminal failed budgets are reported separately.",
            "All duration quantiles include successful and failed recorded attempts. They omit still-running attempts and are descriptive under right censoring.",
            "Shared paired durations appear in both engine descriptions; engine columns are not additive costs. Computational summaries count each shared id once.",
            "Persisted unique-job seconds omit unlogged work and retries of an already saved half-pair. Replay audits have no duration; conflicting shared durations use their maximum as a lower bound, not an exact total.",
            "No savings mean or whole-run ETA is inferred from completion-selected subsets of an unfinished run.",
            "Fit wall seconds do not establish CPU utilization. Only optional sacct TotalCPU/CPUTimeRAW supplies measured allocation utilization.",
            "Timeouts remain failed/unknown planned outcomes in inference. Never recompute denominators using successful permutations only.",
            "In the illustrative historical4h/runtime11800 scenario, paired timeout600 implies2400s nominal drain: signal1500 alone is insufficient and runtime11800 leaves200s. This is not the frozen confirmation run's execution budget."]}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path, help="frozen config")
    ap.add_argument("--root", required=True, type=Path, help="run directory with complete local or HPC fit logs")
    ap.add_argument("--output-prefix", required=True, type=Path, help="write PREFIX.json and PREFIX.md")
    ap.add_argument("--sacct", type=Path, help="pipe-separated sacct -P output with headers, CPUTimeRAW and TotalCPU")
    ap.add_argument("--forecast-datasets", type=int, help="future datasets per same cell; used only after complete-cohort validation")
    ap.add_argument("--forecast-workers", type=int, help="workers for ideal same-configuration planning scenario")
    args = ap.parse_args()
    if (args.forecast_datasets is None) != (args.forecast_workers is None):
        ap.error("Both --forecast-datasets and --forecast-workers are required together")
    if any(v is not None and v < 1 for v in (args.forecast_datasets, args.forecast_workers)):
        ap.error("Forecast datasets and workers must be positive")
    output_json = Path(str(args.output_prefix) + ".json")
    output_md = Path(str(args.output_prefix) + ".md")
    if output_json.resolve() == args.config.resolve() or output_json.name in {"manifest.json", "progress.json", "rows.jsonl", "fit_records.jsonl"}:
        ap.error("Output must not overwrite configuration or raw run artifacts")
    result = analyze(args.config, args.root, args.sacct, args.forecast_datasets, args.forecast_workers)
    status = ("Complete cohort." if result["complete_for_throughput_estimation"] else
              "TERMINAL BUT INCOMPLETE / INVALID: no complete-cohort savings or forecast." if result["execution_collection_terminal"] else
              "RIGHT CENSORED / INCOMPLETE: no complete-cohort savings mean or whole-run ETA is reported.")
    lines = [f"# {result['run_id']} throughput audit", "", status,
        "", "| Cell / engine / scheme | Planned | Recorded started | Decision resolved | Finished | Full B | Selected full-p finished/planned | Fits / failed | Seconds mean / p90 / p95 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    def f(v):
        return "NA" if v is None else f"{v:.4g}"
    for g in result["cell_engine_schemes"]:
        d = g["fit_seconds_all_recorded_attempts"]
        lines.append(f"| {g['cell']} / {g['engine']} / {g['scheme']} | {g['planned_datasets']} | {g['recorded_started_datasets']} | {g['decision_resolved_datasets']} | {g['finished_datasets']} | {g['full_budget_datasets']} | {g['selected_fullp_finished_datasets']}/{g['planned_selected_fullp_datasets']} | {g['statistical_fit_records']} / {g['failed_fit_records']} | {f(d['mean'])} / {f(d['p90'])} / {f(d['p95'])} |")
    lines += ["", f"Statistical fit records: {result['unambiguous_statistical_fit_records']}; unique computational jobs: {result['computation']['unique_computational_jobs']}; recorded unique-job seconds: {f(result['computation']['recorded_unique_job_seconds'])}.",
        f"Duplicate fit identities: {result['duplicate_fit_identities']}; replay audit entries: {result['paired_replay_audit_entries']}; shared-duration conflicts: {result['computation']['shared_duration_conflicts']}.",
        "", "Complete-cohort savings and optional future planning:", "```json",
        json.dumps({"savings": result["complete_cohort_savings"], "future": result["future_run_planning_scenario"]}, indent=2), "```", "",
        "Illustrative historical development scenarios only: 16workers,4h wall,runtime11800s,signal1500s. Not the analyzed run's execution certification.", "",
        "| Fit timeout | Paired | Nominal drain | Signal-only remaining | Runtime-only remaining (startup=0) |",
        "|---:|---|---:|---:|---:|"]
    lines += [f"| {d['fit_timeout_s']} | {d['paired']} | {d['nominal_drain_s']} | {d['signal_only_remaining_after_drain_s']} | {d['runtime_only_remaining_after_drain_s']} |" for d in result["drain_checks"]]
    lines += ["", "All time-budget margins are seconds and exclude startup/fsync/serialization/scheduling overhead.", ""]
    lines += ["- " + note for note in result["notes"]]
    if result["issues"]:
        lines += ["", "Input/checkpoint issues:", "```json", json.dumps(result["issues"], indent=2), "```"]
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    output_md.write_text("\n".join(lines) + "\n")
    print(json.dumps({"prefix": str(args.output_prefix), "complete": result["complete_for_throughput_estimation"],
                      "fit_records": result["unambiguous_statistical_fit_records"], "issues": len(result["issues"])}))


if __name__ == "__main__":
    main()
