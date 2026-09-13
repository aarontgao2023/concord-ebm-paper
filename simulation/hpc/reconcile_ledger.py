"""Record member accounting under the project lock; close only complete parents.

Each sacct snapshot replaces earlier member evidence. Missing records never
inherit an earlier terminal state. Failed queries are recorded without releasable
allocations. Original receipt identities remain intact. This tool never submits,
cancels or retries a job and never generates scientific data or performs fitting.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import re
import subprocess

from accounted_members import (FIELD, SCHEMA, TERMINAL_STATES, receipt_identity,
                               releasable_members, sacct_command, utc_timestamp)
from manage_campaign import write_json


def parse_sacct(text, entry):
    """Retain only exact planned allocations; a parent summary releases nothing."""
    identity = receipt_identity(entry)
    parent = identity["job_id"]
    expected = {f"{parent}_{task}" for task in identity["array_task_ids"]}
    seen, records = set(), []
    for line in text.splitlines():
        if not line.strip():
            continue
        parts = line.strip().split("|")
        if len(parts) != 7:
            raise ValueError("Malformed sacct accounting row")
        job, state, code, elapsed, cpus, start, end = parts
        if job in seen:
            raise ValueError("Duplicate sacct allocation identity")
        seen.add(job)
        if job == parent:
            continue
        if job not in expected:
            raise ValueError(f"Unexpected sacct allocation identity {job!r}")
        match = re.fullmatch(r"([A-Z_]+)(?: by [0-9]+)?", state)
        if match is None:
            raise ValueError("Malformed/unexpanded sacct state")
        if not re.fullmatch(r"[0-9]+:[0-9]+", code):
            raise ValueError("Malformed sacct exit code")
        if not re.fullmatch(r"[0-9]+", elapsed) or not re.fullmatch(r"[0-9]+", cpus):
            raise ValueError("Malformed sacct allocation counts")
        records.append({"job_id": job, "state": match[1], "raw_state": state,
                        "exit_code": code, "elapsed_s": int(elapsed), "allocated_cpus": int(cpus),
                        "start": start, "end": end})
    return records


def reconcile(ledger, query, now):
    """Transform a copy, using an injected one-query-per-job accounting reader."""
    if ledger.get("schema_version") != "campaign_ledger_v1" or not isinstance(ledger.get("submissions"), list):
        raise ValueError("Explicit campaign_ledger_v1 required")
    result = deepcopy(ledger)
    closed, unresolved, recorded, errors = [], [], [], []
    for entry in result["submissions"]:
        if entry.get("state") not in ("intent", "submission_unknown", "submitted"):
            continue
        job = entry.get("job_id")
        if not isinstance(job, str) or not re.fullmatch(r"[1-9][0-9]*", job):
            unresolved.append({"plan_id": entry.get("plan_id"), "reason": "requires explicit submission reconciliation"})
            continue
        identity = receipt_identity(entry)
        command = sacct_command(job)
        try:
            allocations = parse_sacct(query(command), entry)
            status, error = "ok", None
        except (OSError, subprocess.SubprocessError, ValueError) as exc:
            # Failed/malformed queries invalidate old release evidence. No retry.
            allocations, status, error = [], "failed", f"{type(exc).__name__}: {exc}"
            errors.append({"job_id": job, "error": error})
        captured = now()
        utc_timestamp(captured)
        evidence = {"schema_version": SCHEMA, "source": "sacct", "command": command,
                    "captured_utc": captured, "receipt": identity, "query_status": status,
                    "allocations": allocations}
        if error:
            evidence["error"] = error
        entry[FIELD] = evidence
        released = releasable_members(entry)
        recorded.append({"job_id": job, "planned_members": len(entry["array_task_ids"]),
                         "accounted_members": len(allocations), "terminal_ended_members": len(released),
                         "query_status": status})
        expected = set(entry["array_task_ids"])
        if released == expected:
            entry["state"] = "terminal_verified"
            entry["updated_utc"] = captured
            entry["reconciliation_evidence"] = {
                "source": "sacct -X allocations; every originally submitted array member explicitly terminal with End",
                "captured_utc": captured,
                "allocations": {r["job_id"]: {k: v for k, v in r.items() if k != "job_id"} for r in allocations}}
            closed.append(job)
        else:
            present = {r["job_id"] for r in allocations}
            unresolved.append({"job_id": job, "missing_accounting_members": len(expected) - len(present),
                               "known_nonterminal_members": sum(r["state"] not in TERMINAL_STATES for r in allocations),
                               "unverified_end_or_nonterminal_members": len(allocations) - len(released),
                               "query_status": status})
        if receipt_identity(entry) != identity:
            raise AssertionError("Reconciliation mutated the original receipt identity")
    return result, {"closed": closed, "unresolved": unresolved, "member_evidence_recorded": recorded, "errors": errors}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, default=Path("/N/slate/tg11/ebmcal_v2"))
    args = parser.parse_args(argv)
    root, path = args.project_root.resolve(), args.ledger.resolve()
    if not path.is_relative_to(root):
        raise ValueError("Ledger must remain inside the project root")
    with (root / "campaign.dispatch.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        binding = root / "campaign.dispatch.binding.json"
        if binding.exists():
            bound = json.loads(binding.read_text())
            if bound.get("ledger_path") and Path(bound["ledger_path"]).resolve() != path:
                raise ValueError("Project is bound to another ledger")
        ledger = json.loads(path.read_text())
        def query(command):
            response = subprocess.run(command, check=True, text=True, capture_output=True, timeout=30)
            return response.stdout
        updated, report = reconcile(ledger, query, lambda: datetime.now(timezone.utc).isoformat())
        write_json(path, updated)
    print(json.dumps(report))
    return 2 if report["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
