"""Invented accounting evidence only; no scheduler or remote calls."""
from copy import deepcopy
import unittest

import accounted_members as A


def receipt():
    return {"state": "submitted", "job_id": "123", "plan_id": "a" * 64, "run_id": "toy",
            "snapshot_sha256": "b" * 64, "array_task_ids": [0, 2, 7]}


def evidence(entry, allocations):
    return {"schema_version": A.SCHEMA, "source": "sacct", "command": A.sacct_command(entry["job_id"]),
            "captured_utc": "2026-09-08T17:00:00+00:00", "receipt": A.receipt_identity(entry),
            "query_status": "ok", "allocations": allocations}


def allocation(task=0, state="COMPLETED", end="2026-09-08T12:59:00"):
    return {"job_id": f"123_{task}", "state": state, "end": end}


class AccountedMemberTests(unittest.TestCase):
    def test_sacct_expands_large_arrays_without_changing_parent_argument_position(self):
        command = A.sacct_command("10262276")
        self.assertEqual(command[4:7], ["-j", "10262276", "--array"])
        self.assertEqual(command.count("--array"), 1)
        self.assertEqual(command[-1], A.SACCT_FORMAT)

    def entry(self, allocations):
        entry = receipt()
        entry[A.FIELD] = evidence(entry, allocations)
        return entry

    def test_absence_or_missing_member_never_releases(self):
        self.assertEqual(A.releasable_members(receipt()), set())
        self.assertEqual(A.releasable_members(self.entry([])), set())
        self.assertEqual(A.releasable_members(self.entry([allocation(2)])), {2})

    def test_all_explicit_terminal_states_including_failure_release(self):
        for state in A.TERMINAL_STATES:
            with self.subTest(state=state):
                self.assertEqual(A.releasable_members(self.entry([allocation(state=state)])), {0})

    def test_nonterminal_unknown_and_missing_end_stay_reserved(self):
        for state in ("RUNNING", "PENDING", "COMPLETING", "REQUEUED", "UNKNOWN"):
            self.assertEqual(A.releasable_members(self.entry([allocation(state=state)])), set())
        for end in (None, "", "Unknown", "None", "-", "2026-09-08", "2026-99-99T12:00:00"):
            self.assertEqual(A.releasable_members(self.entry([allocation(end=end)])), set())

    def test_queue_active_member_overrides_previous_terminal(self):
        entry = self.entry([allocation(0), allocation(2)])
        self.assertEqual(A.releasable_members(entry, active_members=[0]), {2})
        self.assertEqual(A.releasable_members(entry, active_members=[0, 2]), set())

    def test_exact_parent_and_task_identity_required(self):
        for job in ("124_0", "123_1", "123", "123_0.batch", "123_[0-2]", "123_00", "0123_0", 123):
            entry = self.entry([allocation()])
            entry[A.FIELD]["allocations"][0]["job_id"] = job
            with self.subTest(job=job), self.assertRaises(ValueError):
                A.releasable_members(entry)

    def test_duplicate_member_and_receipt_identity_changes_are_rejected(self):
        with self.assertRaises(ValueError):
            A.releasable_members(self.entry([allocation(), allocation()]))
        for key, value in (("job_id", "124"), ("plan_id", "c" * 64), ("run_id", "other"),
                           ("snapshot_sha256", "c" * 64), ("array_task_ids", [0, 2])):
            entry = self.entry([allocation()])
            entry[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                A.releasable_members(entry)

    def test_source_command_schema_and_utc_are_mandatory(self):
        for key, value in (("source", "squeue"), ("command", ["sacct", "-j", "123"]), ("schema_version", "other"),
                           ("captured_utc", "2026-09-08T17:00:00"), ("captured_utc", "2026-09-08T13:00:00-04:00"),
                           ("captured_utc", "garbage"), ("query_status", None), ("allocations", {})):
            entry = self.entry([allocation()])
            entry[A.FIELD][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                A.releasable_members(entry)
        entry = self.entry([allocation()])
        with self.assertRaises(ValueError):
            A.releasable_members(entry, now="2026-09-08T16:59:59Z")

    def test_failed_query_drops_all_old_evidence(self):
        entry = self.entry([])
        entry[A.FIELD].update(query_status="failed", error="TimeoutExpired")
        self.assertEqual(A.releasable_members(entry), set())
        entry[A.FIELD]["allocations"] = [allocation()]
        with self.assertRaises(ValueError):
            A.releasable_members(entry)

    def test_helper_does_not_mutate_inputs(self):
        entry = self.entry([allocation()])
        before = deepcopy(entry)
        self.assertEqual(A.releasable_members(entry), {0})
        self.assertEqual(entry, before)


if __name__ == "__main__":
    unittest.main()
