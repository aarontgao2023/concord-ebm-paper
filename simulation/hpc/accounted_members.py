"""Pure validation of sacct evidence for individual array submission slots.

This only certifies a member set; it does not mutate receipt identities. The
caller controls capacity and continuation policy. Active queue members always
win over earlier terminal accounting evidence.
"""
from datetime import datetime, timedelta, timezone
import re

SCHEMA = "sacct_member_evidence_v1"
FIELD = "member_accounting_evidence"
TERMINAL_STATES = frozenset({"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY",
                             "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"})
RECEIPT_FIELDS = ("job_id", "plan_id", "run_id", "snapshot_sha256", "array_task_ids")
SACCT_FORMAT = "--format=JobID%80,State%40,ExitCode,ElapsedRaw,AllocCPUS,Start,End"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sacct_command(job_id):
    require(isinstance(job_id, str) and re.fullmatch(r"[1-9][0-9]*", job_id) is not None, "Invalid accepted parent job ID")
    return ["sacct", "-X", "-n", "-P", "-j", job_id, "--array", SACCT_FORMAT]


def receipt_identity(entry):
    """Validate and copy the original identity; never rewrite a receipt field."""
    sacct_command(entry.get("job_id"))
    for key in ("plan_id", "snapshot_sha256"):
        require(isinstance(entry.get(key), str) and re.fullmatch(r"[0-9a-f]{64}", entry[key]) is not None, f"Invalid receipt {key}")
    require(isinstance(entry.get("run_id"), str) and bool(entry["run_id"]), "Missing receipt run ID")
    tasks = entry.get("array_task_ids")
    require(isinstance(tasks, list) and bool(tasks) and all(type(v) is int and v >= 0 for v in tasks)
            and len(tasks) == len(set(tasks)), "Receipt needs unique nonnegative integer member IDs")
    return {key: list(entry[key]) if key == "array_task_ids" else entry[key] for key in RECEIPT_FIELDS}


def utc_timestamp(value):
    require(isinstance(value, str) and bool(value), "Missing sacct capture UTC")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("Malformed sacct capture UTC") from exc
    require(parsed.tzinfo is not None and parsed.utcoffset() == timedelta(0), "Accounting capture timestamp must explicitly be UTC")
    return parsed


def ended(value):
    """Slurm End may be local ISO time; missing/Unknown/None is not evidence."""
    if not isinstance(value, str) or not value or value.strip() != value:
        return False
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d+)?(?:Z|[+-]\d{2}:\d{2})?", value):
        return False
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
        return True
    except ValueError:
        return False


def releasable_members(entry, active_members=(), now=None):
    """Certified off-queue task indices. Malformed present evidence raises.

    No evidence, nonterminal state, or missing End releases no such member.
    Each evidence object is one query snapshot: missing members never inherit a
    previous terminal state. Optional ``now`` checks future capture timestamps.
    This function reads no clock, filesystem, queue, or accounting service.
    """
    evidence = entry.get(FIELD)
    if evidence is None:
        return set()
    identity = receipt_identity(entry)
    require(isinstance(evidence, dict) and evidence.get("schema_version") == SCHEMA, "Unsupported member accounting evidence")
    require(isinstance(evidence.get("receipt"), dict), "Missing accounting receipt binding")
    receipt_identity(evidence["receipt"])
    require(evidence.get("receipt") == identity, "Accounting evidence belongs to a different original receipt")
    require(evidence.get("source") == "sacct" and evidence.get("command") == sacct_command(entry["job_id"]), "Accounting evidence must come from the exact sacct allocation query")
    captured = utc_timestamp(evidence.get("captured_utc"))
    if now is not None:
        current = utc_timestamp(now) if isinstance(now, str) else now
        require(isinstance(current, datetime) and current.tzinfo is not None and captured <= current.astimezone(timezone.utc), "Accounting evidence is future-dated")
    active = set(active_members)
    require(all(type(v) is int and v >= 0 for v in active), "Active members must be nonnegative task integers")
    allocations = evidence.get("allocations")
    require(isinstance(allocations, list), "Accounting allocations must be a list")
    status = evidence.get("query_status")
    require(status in {"ok", "failed"}, "Missing accounting query status")
    if status == "failed":
        require(not allocations and isinstance(evidence.get("error"), str) and bool(evidence["error"]), "Failed query cannot provide releasable allocations")
        return set()
    seen, released = set(), set()
    expected = {f"{entry['job_id']}_{task}": task for task in entry["array_task_ids"]}
    for allocation in allocations:
        require(isinstance(allocation, dict), "Malformed allocation record")
        job = allocation.get("job_id")
        require(isinstance(job, str) and job in expected and re.fullmatch(r"[1-9][0-9]*_(?:0|[1-9][0-9]*)", job) is not None,
                "Allocation is not an exact planned parent_task identity")
        require(job not in seen, "Duplicate accounting member identity")
        seen.add(job)
        state = allocation.get("state")
        require(isinstance(state, str) and re.fullmatch(r"[A-Z_]+", state) is not None, "Malformed normalized accounting state")
        if state in TERMINAL_STATES and ended(allocation.get("end")):
            released.add(expected[job])
    return released - active
