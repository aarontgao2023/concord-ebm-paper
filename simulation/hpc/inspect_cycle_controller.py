#!/usr/bin/env python3
"""Observe accepted cycle jobs without submitting or changing any HPC state.

Example on the HPC host (also works as Python source supplied through stdin):
    python - --cycle /N/slate/tg11/ebmcal_v2/control/cycles/CYCLE < this_file.py

Only direct, real cycle directories below the fixed control/cycles root are
accepted. All inputs are checked before any scheduler query. Output is one JSON
document on stdout; no files are written. Command success is not job completion
or an independent verification of a receipt's identity. Queries are sequential;
each has a 30-second timeout, including any accounting fallback.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


ROOT = Path("/N/slate/tg11/ebmcal_v2")
MAX_JOBS = 64
TIMEOUT_SECONDS = 30
SACCT_FORMAT = (
    "--format=JobID%80,JobName%120,State%40,ExitCode,ElapsedRaw,"
    "AllocCPUS,ReqCPUS,Submit,Start,End,TimelimitRaw,Comment%160"
)


def utc():
    return datetime.now(timezone.utc).isoformat()


def as_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value or ""


def query(argv, *, utc_environment=False):
    started = time.monotonic()
    record = {"argv": argv, "started_utc": utc(), "timeout_s": TIMEOUT_SECONDS,
              "environment_override": {"TZ": "UTC"} if utc_environment else {},
              "stdout": "", "stderr": "", "returncode": None,
              "timed_out": False, "exception": None}
    env = os.environ.copy()
    if utc_environment:
        env["TZ"] = "UTC"
    try:
        result = subprocess.run(argv, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", env=env,
                                timeout=TIMEOUT_SECONDS, check=False, shell=False)
        record.update(stdout=result.stdout, stderr=result.stderr,
                      returncode=result.returncode)
    except subprocess.TimeoutExpired as exc:
        record.update(stdout=as_text(exc.stdout), stderr=as_text(exc.stderr),
                      timed_out=True, exception="TimeoutExpired")
    except OSError as exc:
        record["exception"] = f"{type(exc).__name__}: {exc}"
    record.update(finished_utc=utc(), elapsed_s=time.monotonic() - started)
    return record


def collect_inputs(cycles):
    fixed_cycles = ROOT / "control/cycles"
    # Reject a redirected control root, including a symlink into another tree.
    if fixed_cycles.resolve(strict=True) != fixed_cycles:
        raise ValueError("The fixed control/cycles root must be a real directory")
    sources, jobs, seen_cycles = [], {}, set()
    for supplied in cycles:
        path = Path(supplied)
        if not path.is_absolute():
            raise ValueError("--cycle requires an absolute cycle directory path")
        real = path.resolve(strict=True)
        if real.parent != fixed_cycles or not real.is_dir() or path != real:
            raise ValueError("Cycle must be a direct real directory under the fixed control/cycles root")
        if real in seen_cycles:
            continue
        seen_cycles.add(real)
        dispatch_path = real / "dispatch.json"
        if dispatch_path.resolve(strict=True) != dispatch_path or not dispatch_path.is_file():
            raise ValueError("dispatch.json must be a real file in the selected cycle")
        before = dispatch_path.stat()
        raw = dispatch_path.read_bytes()
        after = dispatch_path.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
                after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
            raise ValueError("dispatch.json changed during capture; retry after it settles")
        data = json.loads(raw)
        if (not isinstance(data, dict) or data.get("schema_version") != "campaign_dispatch_v1"
                or data.get("project") != "ebmcal_v2" or not isinstance(data.get("outcomes"), list)):
            raise ValueError("Expected an ebmcal_v2 campaign_dispatch_v1 with outcomes")
        sources.append({"cycle": str(real), "dispatch_path": str(dispatch_path),
                        "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw),
                        "captured_utc": utc()})
        for index, outcome in enumerate(data["outcomes"]):
            if not isinstance(outcome, dict):
                raise ValueError("Malformed dispatch outcome")
            if outcome.get("state") != "submitted":
                continue
            job_id, plan_id = outcome.get("job_id"), outcome.get("plan_id")
            if not isinstance(job_id, str) or re.fullmatch(r"[1-9][0-9]*", job_id) is None:
                raise ValueError("Submitted outcome requires a positive decimal parent job ID")
            if not isinstance(plan_id, str) or re.fullmatch(r"[0-9a-f]{64}", plan_id) is None:
                raise ValueError("Submitted outcome requires a SHA256-shaped plan ID")
            refs = jobs.setdefault(job_id, [])
            refs.append({"cycle": str(real), "outcome_index": index, "plan_id": plan_id})
            if len(jobs) > MAX_JOBS:
                raise ValueError(f"At most {MAX_JOBS} distinct submitted jobs may be observed")
    return sources, jobs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycle", action="append", required=True, metavar="ABSOLUTE_PATH",
                        help="Repeat for each real directory under the fixed control/cycles root")
    args = parser.parse_args(argv)
    report = {"schema_version": "cycle_controller_observation_v1", "started_utc": utc(),
              "root": str(ROOT), "read_only": True, "sources": [], "jobs": [],
              "identity_verification_performed": False,
              "job_completion_assessed": False, "errors": [],
              "limitations": [
                  "Observations are sequential, not an atomic queue snapshot.",
                  "Missing controller output or accounting rows prove neither failure nor completion.",
                  "Successful queries and SHA256-shaped plan IDs do not verify receipt/controller identity.",
                  "Accounting comment availability depends on cluster accounting configuration.",
                  "Timeouts preserve captured partial output; no failed query is automatically retried."]}
    try:
        sources, jobs = collect_inputs(args.cycle)
        report["sources"] = sources
    except (OSError, ValueError, TypeError) as exc:
        report.update(status="input_rejected_no_scheduler_queries", finished_utc=utc())
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))
        return 2
    for job_id, references in jobs.items():
        control = query(["scontrol", "-o", "show", "job", job_id])
        item = {"job_id": job_id, "submitted_references": references,
                "conflicting_plan_ids": len({r["plan_id"] for r in references}) > 1,
                "controller": control, "accounting": None, "fallback_reason": None}
        if control["returncode"] != 0 or not control["stdout"].strip():
            item["fallback_reason"] = "controller_query_unsuccessful_or_empty_not_a_job_state"
            item["accounting"] = query(
                ["sacct", "-X", "-n", "-P", "-j", job_id, "--array", SACCT_FORMAT],
                utc_environment=True)
        report["jobs"].append(item)
    report.update(status="observation_finished", distinct_jobs=len(jobs), finished_utc=utc())
    print(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
