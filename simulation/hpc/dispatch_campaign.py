#!/usr/bin/env python3
"""Idempotent, narrowly scoped sbatch dispatch of a fresh campaign plan.

Default: validate only. --apply writes a durable intent before each subprocess.
This command has no cancellation, requeue, throttle or accounting-mutation API.
A fresh normalized queue file is reread under the project lock; collection itself
is external. See DISPATCH_CAMPAIGN.md for the scheduler-race boundary.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import subprocess

import inventory_run as I
import manage_campaign as M

PROJECT = "ebmcal_v2"
DEFAULT_ROOT = Path("/N/slate/tg11/ebmcal_v2")
IDENTITY_FIELDS = ("run_id", "snapshot_sha256", "array_task_ids", "nchunks", "mapping",
                   "cpus_per_task", "max_concurrent", "walltime_s", "dispatch_generation")
EXECUTION_IDENTITY_FIELDS = ("execution_walltime_policy", "execution_evidence")


def submission_identity(submission):
    """Preserve ordinary receipt hashes; bind both fields of any new exception."""
    identity = {name: submission[name] for name in IDENTITY_FIELDS}
    present = [name in submission for name in EXECUTION_IDENTITY_FIELDS]
    if any(present):
        if not all(present) or submission["execution_walltime_policy"] is None:
            raise ValueError("Execution override identity requires both policy and evidence")
        identity.update({name: submission[name] for name in EXECUTION_IDENTITY_FIELDS})
    return identity


def utc():
    return datetime.now(timezone.utc)


def durable_json(path, value):
    """Persist bytes and the rename before any submission is allowed."""
    path = Path(path)
    temporary = path.with_name(path.name + f".tmp-{os.getpid()}")
    with temporary.open("w") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    directory_fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def contained(path, root, name):
    path = Path(path)
    if not path.is_absolute():
        raise ValueError(f"{name} must be an absolute HPC path")
    resolved = path.resolve()
    if resolved == root or not resolved.is_relative_to(root):
        raise ValueError(f"{name} is outside the allowed project root")
    return resolved


def fresh(timestamp, now, max_age_s, name):
    age = (now-M.parse_utc(timestamp)).total_seconds()
    if not 0 <= age <= max_age_s:
        raise ValueError(f"{name} is stale or future-dated ({age:.1f}s; limit {max_age_s}s)")


def captured_text(value):
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value or "")[:10000]


def check_batch_topology(batch):
    values = {"nodes": [], "ntasks": []}
    for line in batch.splitlines():
        directive = re.match(r"^\s*#SBATCH\s+(.+)$", line)
        if not directive:
            continue
        tokens = shlex.split(directive[1])
        for index, token in enumerate(tokens):
            for name, short in (("nodes", "N"), ("ntasks", "n")):
                if token in {f"--{name}", f"-{short}"}:
                    values[name].append(tokens[index+1] if index+1 < len(tokens) else None)
                elif token.startswith(f"--{name}="):
                    values[name].append(token.split("=", 1)[1])
                elif re.fullmatch(rf"-{short}[0-9]+", token):
                    values[name].append(token[2:])
    for name, declarations in values.items():
        if declarations != ["1"]:
            raise ValueError(f"Frozen batch script must declare exactly one --{name}=1")


def validate_plan_header(plan, now, plan_max_age_s, weekly_used):
    if plan.get("schema_version") != "campaign_plan_v1" or plan.get("project") != PROJECT:
        raise ValueError("Only an ebmcal_v2 campaign_plan_v1 is accepted")
    fresh(plan["created_utc"], now, plan_max_age_s, "Plan")
    fresh(plan["queue_captured_utc"], now, plan_max_age_s, "Plan's queue snapshot")
    if plan.get("pause_new_submissions"):
        raise ValueError("Plan pauses new submissions")
    used = weekly_used if weekly_used is not None else plan.get("weekly_used_percent")
    if type(used) not in (int, float) or not math.isfinite(used) or not 0 <= used < 90:
        raise ValueError("An explicit current weekly-used percentage below 90 is required")
    if not isinstance(plan.get("submission_plans"), list):
        raise ValueError("Plan lacks a submission list")
    ids = [s.get("plan_id") for s in plan["submission_plans"]]
    if len(ids) != len(set(ids)) or any(not I.is_hash(value) for value in ids):
        raise ValueError("Submission plan IDs must be unique SHA256 identities")
    M.positive(plan.get("cpu_budget", {}).get("cap"), "project CPU cap")
    slots = plan.get("user_job_budget", {})
    if not isinstance(slots, dict):
        raise ValueError("User submitted-job budget must be an object")
    if slots.get("cap") is not None:
        M.positive(slots["cap"], "max_user_submitted_jobs")


def validate_submission(submission, run_report, root):
    identity = submission_identity(submission)
    if I.digest(identity) != submission.get("plan_id"):
        raise ValueError("Submission identity hash does not match plan_id")
    run_id = identity["run_id"]
    if not isinstance(run_id, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", run_id):
        raise ValueError("Invalid run ID")
    if run_report.get("run_id") != run_id or run_report.get("snapshot_sha256") != identity["snapshot_sha256"]:
        raise ValueError("Run report and submission snapshot identity disagree")
    argv = submission.get("argv")
    core_signal_override = identity.get("execution_walltime_policy") == M.CORE_PAIRED_1H_POLICY
    shift = 1 if core_signal_override else 0
    if not isinstance(argv, list) or len(argv) != 12 + shift or any(not isinstance(a, str) for a in argv):
        raise ValueError(f"Expected the exact {12 + shift}-argument generated sbatch command")
    if argv[0] != "sbatch":
        raise ValueError("Only sbatch is allowed")
    snapshot_text, name = argv[9 + shift], argv[10 + shift]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", name):
        raise ValueError("Invalid frozen configuration name")
    if name != run_id:
        raise ValueError("config_name must equal run_id; aliases cannot share an output directory")
    snapshot = contained(snapshot_text, root, "Snapshot")
    if not snapshot.is_relative_to(root/"snapshots"):
        raise ValueError("Snapshot must belong to the project's snapshots directory")
    script = contained(argv[8 + shift], root, "Batch script")
    if script != snapshot/"hpc/v2/run_array.sbatch":
        raise ValueError("Unexpected batch script")
    frozen = M.snapshot_info(snapshot, name, identity["snapshot_sha256"])
    check_batch_topology(script.read_text())
    if Path(frozen["remote_root"]).resolve() != root:
        raise ValueError("Frozen batch script targets a different project root")
    config = M.read_json(frozen["config_path"])
    cells, seeds, full_p = I.validate_config(config)
    if config["run_id"] != run_id or config["phase"] not in {"development", "confirmation", "diagnostic"}:
        raise ValueError("Frozen run identity or phase differs from the planned run")
    nchunks = M.positive(identity["nchunks"], "nchunks")
    if nchunks > config["datasets"] or nchunks*len(cells) > 1000:
        raise ValueError("Chunk allocation is empty or exceeds checked MaxArraySize=1000")
    if identity["mapping"] != frozen["mapping"] or run_report.get("mapping") != frozen["mapping"]:
        raise ValueError("Frozen array mapping changed")
    indexes = identity["array_task_ids"]
    if (not isinstance(indexes, list) or not indexes
            or any(type(i) is not int or not 0 <= i < len(cells)*nchunks for i in indexes)
            or indexes != sorted(set(indexes))):
        raise ValueError("Array task IDs must be sorted, unique members of the frozen layout")
    cpus = M.positive(identity["cpus_per_task"], "cpus_per_task")
    concurrency = M.positive(identity["max_concurrent"], "max_concurrent")
    if concurrency > len(indexes):
        raise ValueError("Array throttle exceeds the number of submitted members")
    wall = M.positive(identity["walltime_s"], "walltime_s")
    M.positive(identity["dispatch_generation"], "dispatch_generation")
    evidence = M.validate_execution_walltime(frozen, config, cpus, wall,
                                            identity.get("execution_walltime_policy"), nchunks=nchunks)
    if evidence != identity.get("execution_evidence"):
        raise ValueError("Execution evidence differs from reverified frozen source and timing")
    output_root = root/"runs"/name
    if Path(run_report.get("output_root", "")).resolve() != output_root:
        raise ValueError("Plan was not built against this HPC output root")
    contained(str(output_root), root, "Output root")
    hh, rem = divmod(wall, 3600)
    mm, ss = divmod(rem, 60)
    array = ",".join(map(str, indexes))+f"%{concurrency}"
    logs = f"{frozen['remote_root']}/logs/{name}"
    expected = ["sbatch", "--parsable", f"--array={array}", f"--cpus-per-task={cpus}",
                f"--time={hh:02d}:{mm:02d}:{ss:02d}", f"--comment={PROJECT}:{submission['plan_id']}",
                f"--output={logs}_%A_%a.out", f"--error={logs}_%A_%a.err",
                f"{snapshot_text}/hpc/v2/run_array.sbatch", snapshot_text, name, str(nchunks)]
    if core_signal_override:
        expected.insert(8, M.CORE_PAIRED_1H_SIGNAL_OPTION)
    if argv != expected:
        raise ValueError("Command differs from the strict sbatch argument allowlist")
    intent = submission.get("ledger_intent", {})
    if any(intent.get(k) != v for k, v in identity.items()) or intent.get("plan_id") != submission["plan_id"] or intent.get("state") != "intent":
        raise ValueError("Ledger intent differs from validated submission identity")
    if any((name in intent) != (name in identity) for name in EXECUTION_IDENTITY_FIELDS):
        raise ValueError("Ledger intent contains an unbound execution override")
    return {"identity": identity, "snapshot": snapshot, "config": config, "cells": cells,
            "seeds": seeds, "full_p": full_p, "output_root": output_root,
            "argv": expected, "frozen": frozen}


def check_stop(validated):
    config, root = validated["config"], validated["output_root"]
    identity = validated["identity"]
    project_root = root.parent.parent
    stop_files = [project_root/"STOP"]
    if config.get("stop_file"):
        stop = Path(config["stop_file"])
        if not stop.is_absolute():
            stop = validated["snapshot"]/stop
        stop_files.append(stop)
    if any(p.exists() for p in stop_files):
        raise ValueError("Project/configuration STOP file is present")
    if not config.get("stop_file"):
        for index in identity["array_task_ids"]:
            cell = index//identity["nchunks"] if identity["mapping"] == "cell_major" else index%len(validated["cells"])
            if (root/f"cell_{cell}"/"STOP").exists():
                raise ValueError("Cell STOP file is present")


def recheck_work(validated):
    config, root = validated["config"], validated["output_root"]
    identity = validated["identity"]
    check_stop(validated)
    inventory = I.inventory(validated["frozen"]["config_path"], root,
                            source_root=validated["snapshot"]/"scripts/v2",
                            nchunks_override=identity["nchunks"])
    errors = [e for e in inventory["issues"] if e["severity"] == "error"]
    if errors or not inventory["stable_rows_snapshot"]:
        raise ValueError(f"Fresh inventory is not stable/intact: {errors}")
    for index in identity["array_task_ids"]:
        if identity["mapping"] == "cell_major":
            cell, chunk = divmod(index, identity["nchunks"])
        else:
            chunk, cell = divmod(index, len(validated["cells"]))
        directory = root/f"cell_{cell}"/f"chunk_{chunk}"
        if not config.get("stop_file") and (directory.parent/"STOP").exists():
            raise ValueError("Cell STOP file is present")
        seeds = sorted(validated["seeds"])[chunk::identity["nchunks"]]
        work = M.chunk_work(directory, seeds, config, validated["full_p"])
        if work["state"] != "needs_resume":
            raise ValueError(f"Array index {index} is now {work['state']}; regenerate the plan")


def recheck_capacity(plan, submissions, queue, ledger, now, max_age_s):
    q = M.queue_view(queue, PROJECT, now, max_age_s)
    if M.parse_utc(queue["captured_utc"]) < M.parse_utc(plan["queue_captured_utc"]):
        raise ValueError("Queue file predates the plan's queue evidence")
    if ledger.get("project") != PROJECT:
        raise ValueError("Wrong project submission ledger")
    pending, unknown_reserved = M.ledger_view(ledger, q["jobs_by_id"], now)
    existing_ids = {entry["plan_id"] for entry in ledger["submissions"]}
    new_reserved = 0
    targets = set()
    for submission in submissions:
        if submission["plan_id"] in existing_ids:
            raise ValueError("Plan ID already has a ledger record; never replay the same dispatch")
        for index in submission["array_task_ids"]:
            key = (submission["run_id"], index)
            if key in targets or q["active"].get(key) or pending.get(key):
                raise ValueError(f"Array target {key} is duplicated, active or awaiting submission reconciliation")
            targets.add(key)
        new_reserved += submission["cpus_per_task"]*submission["max_concurrent"]
    used = q["project_reserved_cpus"]+unknown_reserved
    cap = plan["cpu_budget"]["cap"]
    if used+new_reserved > cap:
        raise ValueError(f"Fresh project reservations {used}+{new_reserved} exceed cap {cap}")
    slots = M.user_job_budget(ledger, q, plan.get("user_job_budget", {}).get("cap"))
    new_members = sum(len(s["array_task_ids"]) for s in submissions)
    if slots["cap"] is not None and slots["used_before_plan"]+new_members > slots["cap"]:
        raise ValueError(f"Fresh user submitted-job members {slots['used_before_plan']}+{new_members} "
                         f"exceed cap {slots['cap']}")
    return {"project_reserved": q["project_reserved_cpus"], "unknown_reserved": unknown_reserved,
            "new_reserved": new_reserved, "cap": cap,
            "user_job_budget": slots | {"new_reserved": new_members,
                "remaining_after_plan": None if slots["cap"] is None else slots["cap"]-slots["used_before_plan"]-new_members}}


def dispatch(plan_path, ledger_path, queue_path, project_root=DEFAULT_ROOT, *, apply=False,
             timeout_s=30, max_age_s=60, plan_max_age_s=None, weekly_used=None, now_func=utc, submitter=None):
    """The injectable submitter is for fixtures; CLI always uses subprocess.run."""
    submitter = submitter or subprocess.run
    root = Path(project_root).resolve()
    if not root.is_dir():
        raise ValueError("Project root must already exist on the dispatch host")
    ledger_path = contained(str(Path(ledger_path).absolute()), root, "Submission ledger")
    queue_path = Path(queue_path).resolve()
    plan_path = Path(plan_path).resolve()
    if not 0 < timeout_s <= 60 or not 0 < max_age_s <= 60:
        raise ValueError("Submission timeout and freshness limit must be in (0,60] seconds")
    plan_age = max_age_s if plan_max_age_s is None else plan_max_age_s
    if type(plan_age) not in (int, float) or not math.isfinite(plan_age) or not 0 < plan_age <= 300:
        raise ValueError("Plan freshness limit must be in (0,300] seconds")
    limits = {"plan_max_age_s": plan_age, "live_queue_max_age_s": max_age_s}
    result = {"schema_version": "campaign_dispatch_v1", "created_utc": now_func().isoformat(),
              "project": PROJECT, "mode": "apply" if apply else "dry_run",
              "dispatcher_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "accounting_validator_sha256": hashlib.sha256(Path(M.A.__file__).read_bytes()).hexdigest(),
              "freshness_limits_s": limits,
              "outcomes": [], "status": "blocked", "notes": [
                  "No cancellation, throttle change or automatic submission retry is supported.",
                  "Queue collection is external; all project submission entrances must share this lock and ledger."]}
    lock_path = root/"campaign.dispatch.lock"
    binding_path = root/"campaign.dispatch.binding.json"
    with lock_path.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plan = M.read_json(plan_path)
        validate_plan_header(plan, now_func(), plan_age, weekly_used)
        binding = {"project": PROJECT, "ledger_path": str(ledger_path)}
        if binding_path.exists():
            existing = M.read_json(binding_path)
            if any(existing.get(k) != v for k, v in binding.items()):
                raise ValueError("Project dispatch root is already bound to a different ledger")
        ledger = M.read_json(ledger_path)
        queue = M.read_json(queue_path)
        submissions = plan["submission_plans"]
        run_reports = {run["run_id"]: run for run in plan.get("runs", [])}
        validated = [validate_submission(s, run_reports[s["run_id"]], root) for s in submissions]
        result["capacity"] = recheck_capacity(plan, submissions, queue, ledger, now_func(), max_age_s)
        for value in validated:
            recheck_work(value)
        if not apply:
            result.update(status="validated_dry_run", validated_plan_ids=[s["plan_id"] for s in submissions])
            return result
        if not binding_path.exists():
            durable_json(binding_path, binding | {"created_utc": now_func().isoformat()})
        for offset, (submission, value) in enumerate(zip(submissions, validated)):
            # Time spent validating/dispatching earlier arrays cannot make a later
            # array exempt from freshness, STOP, ledger or capacity checks.
            try:
                validate_plan_header(plan, now_func(), plan_age, weekly_used)
                ledger = M.read_json(ledger_path)
                queue = M.read_json(queue_path)
                recheck_capacity(plan, submissions[offset:], queue, ledger, now_func(), max_age_s)
                value = validate_submission(submission, run_reports[submission["run_id"]], root)
                recheck_work(value)
                # Filesystem inventory/hash checks can take time. Refresh again
                # after those checks, immediately before committing the intent.
                ledger = M.read_json(ledger_path)
                queue = M.read_json(queue_path)
                validate_plan_header(plan, now_func(), plan_age, weekly_used)
                recheck_capacity(plan, submissions[offset:], queue, ledger, now_func(), max_age_s)
            except (ValueError, KeyError, TypeError, OSError) as exc:
                result["outcomes"].append({"plan_id": submission["plan_id"], "state": "not_dispatched", "reason": str(exc)})
                result["status"] = "stopped_before_next_submission"
                break
            intent = dict(submission["ledger_intent"])
            intent.update(created_utc=now_func().isoformat(), argv=value["argv"],
                          dispatcher_sha256=result["dispatcher_sha256"], freshness_limits_s=limits,
                          plan_created_utc=plan["created_utc"],
                          planning_queue_captured_utc=plan["queue_captured_utc"],
                          dispatch_queue_captured_utc=queue["captured_utc"])
            intent["history"] = [{"state": "intent", "utc": intent["created_utc"]}]
            ledger["submissions"].append(intent)
            durable_json(ledger_path, ledger)
            try:
                check_stop(value)
                validate_plan_header(plan, now_func(), plan_age, weekly_used)
                fresh(queue["captured_utc"], now_func(), max_age_s, "Queue before subprocess")
            except (ValueError, OSError) as exc:
                # The program has not entered submitter yet; this is affirmative
                # no-submit evidence, unlike absence from squeue after a timeout.
                outcome = {"state": "not_submitted_verified", "reason": str(exc),
                           "reconciliation_evidence": {"kind": "dispatcher_pre_invocation_guard",
                                                       "subprocess_invoked": False,
                                                       "utc": now_func().isoformat()}}
                intent.update(outcome, updated_utc=now_func().isoformat())
                intent["history"].append(outcome | {"utc": intent["updated_utc"]})
                durable_json(ledger_path, ledger)
                result["outcomes"].append({"plan_id": submission["plan_id"], **outcome})
                result["status"] = "stopped_before_subprocess"
                break
            # The only external operation: exact validated argv, shell=False.
            environment = {k:v for k,v in os.environ.items()
                           if not k.startswith("SBATCH_") and k != "SLURM_CLUSTERS"}
            try:
                response = submitter(value["argv"], capture_output=True, text=True,
                                     timeout=timeout_s, shell=False, check=False,
                                     stdin=subprocess.DEVNULL, cwd=str(root), env=environment)
                stdout, stderr = response.stdout or "", response.stderr or ""
                job_id = stdout.strip()
                success = response.returncode == 0 and re.fullmatch(r"[1-9][0-9]*", job_id) is not None
                if success and any(e is not intent and e.get("job_id") == job_id and e.get("state") in M.LEDGER_OPEN
                                   for e in ledger["submissions"]):
                    success = False
                outcome = {"state": "submitted" if success else "submission_unknown",
                           "returncode": response.returncode, "stdout": stdout[:10000], "stderr": stderr[:10000]}
                if success:
                    outcome["job_id"] = job_id
                else:
                    outcome["reason"] = "Nonzero exit, non-pure job ID or repeated open receipt; reconcile before any retry"
            except Exception as exc:
                outcome = {"state": "submission_unknown", "reason": str(exc),
                           "exception_type": type(exc).__name__,
                           "stdout": captured_text(getattr(exc, "stdout", "")),
                           "stderr": captured_text(getattr(exc, "stderr", ""))}
            intent.update(outcome, updated_utc=now_func().isoformat())
            intent["history"].append(outcome | {"utc": intent["updated_utc"]})
            # If this save fails, the prior durable intent still blocks a retry.
            durable_json(ledger_path, ledger)
            result["outcomes"].append({"plan_id": submission["plan_id"], **outcome})
            if outcome["state"] != "submitted":
                result["status"] = "submission_unknown_reconciliation_required"
                break
        else:
            result["status"] = "submitted" if submissions else "nothing_to_submit"
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--weekly-used-percent", type=float)
    parser.add_argument("--timeout-s", type=float, default=30)
    parser.add_argument("--max-age-s", type=float, default=60)
    parser.add_argument("--plan-max-age-s", type=float,
                        help="Optional proposal age limit up to300s; live queue remains limited by --max-age-s (at most60s)")
    args = parser.parse_args(argv)
    protected = {p.resolve() for p in (args.plan, args.ledger, args.queue,
                 args.project_root/"campaign.dispatch.lock", args.project_root/"campaign.dispatch.binding.json")}
    target = args.output.resolve()
    if (target in protected or target.is_relative_to(args.project_root.resolve()/"snapshots")
            or target.is_relative_to(args.project_root.resolve()/"runs")):
        parser.error("--output cannot overwrite inputs, dispatch state, frozen sources or simulation results")
    try:
        result = dispatch(args.plan, args.ledger, args.queue, args.project_root, apply=args.apply,
                          timeout_s=args.timeout_s, max_age_s=args.max_age_s,
                          plan_max_age_s=args.plan_max_age_s,
                          weekly_used=args.weekly_used_percent)
    except (ValueError, KeyError, TypeError, OSError) as exc:
        result = {"schema_version": "campaign_dispatch_v1", "created_utc": utc().isoformat(),
                  "mode": "apply" if args.apply else "dry_run", "status": "blocked",
                  "outcomes": [], "error": str(exc)}
    M.write_json(args.output, result)
    return 0 if result["status"] in {"validated_dry_run", "submitted", "nothing_to_submit"} else 3


if __name__ == "__main__":
    raise SystemExit(main())
