"""Fake sacct and local temporary ledgers only; no scheduler operations."""
from copy import deepcopy
import fcntl
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import accounted_members as A
import reconcile_ledger as R
from test_accounted_members import receipt


def line(task=0, state="COMPLETED", end="2026-09-08T12:59:00", job=None):
    return f"{job or '123_' + str(task)}|{state}|0:0|60|16|2026-09-08T12:58:00|{end}\n"


class ReconcileLedgerTests(unittest.TestCase):
    def ledger(self):
        return {"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": [receipt()]}

    def apply(self, output, ledger=None):
        commands = []
        def query(command):
            commands.append(command)
            return output
        updated, report = R.reconcile(ledger or self.ledger(), query, lambda: "2026-09-08T17:00:00Z")
        self.assertEqual(commands, [A.sacct_command("123")])
        return updated, report

    def test_partial_terminal_evidence_saved_without_parent_close(self):
        original = self.ledger()
        before = deepcopy(original)
        updated, report = self.apply(line(0) + line(2, "RUNNING", "Unknown"), original)
        entry = updated["submissions"][0]
        self.assertEqual(entry["state"], "submitted")
        self.assertEqual(A.releasable_members(entry), {0})
        self.assertEqual(report["closed"], [])
        self.assertEqual(report["unresolved"][0]["missing_accounting_members"], 1)
        self.assertEqual(original, before)
        self.assertEqual(A.receipt_identity(entry), A.receipt_identity(before["submissions"][0]))

    def test_complete_terminal_parent_retains_existing_close_semantics(self):
        updated, report = self.apply(line(0) + line(2, "FAILED") + line(7, "TIMEOUT"))
        entry = updated["submissions"][0]
        self.assertEqual(entry["state"], "terminal_verified")
        self.assertEqual(report["closed"], ["123"])
        self.assertEqual(set(entry["reconciliation_evidence"]["allocations"]), {"123_0", "123_2", "123_7"})
        self.assertEqual(A.releasable_members(entry), {0, 2, 7})

    def test_missing_end_does_not_close_parent(self):
        updated, report = self.apply(line(0) + line(2) + line(7, end="Unknown"))
        self.assertEqual(updated["submissions"][0]["state"], "submitted")
        self.assertEqual(A.releasable_members(updated["submissions"][0]), {0, 2})
        self.assertEqual(report["closed"], [])

    def test_new_empty_or_nonterminal_snapshot_does_not_inherit_old_terminal(self):
        first, _ = self.apply(line(0))
        for output in ("", line(0, "PENDING", "Unknown")):
            updated, _ = self.apply(output, first)
            self.assertEqual(A.releasable_members(updated["submissions"][0]), set())

    def test_bad_query_format_duplicate_and_wrong_job_clear_release_evidence(self):
        first, _ = self.apply(line(0))
        for output in ("malformed\n", line(0) + line(0), line(job="999_0"), line(job="123_0.batch"),
                       line(job="123_[0-2]"), line(state="COMPLETED+"), line().replace("|60|", "|-1|")):
            updated, report = self.apply(output, first)
            self.assertEqual(A.releasable_members(updated["submissions"][0]), set())
            self.assertEqual(updated["submissions"][0][A.FIELD]["query_status"], "failed")
            self.assertEqual(len(report["errors"]), 1)

    def test_cancelled_uid_normalization_and_parent_summary(self):
        updated, report = self.apply(line(job="123") + line(0, "CANCELLED by 1234"))
        self.assertEqual(A.releasable_members(updated["submissions"][0]), {0})
        self.assertEqual(report["errors"], [])

    def test_timeout_is_not_retried_and_clears_old_release_evidence(self):
        first, _ = self.apply(line(0))
        calls = []
        def timeout(command):
            calls.append(command)
            raise subprocess.TimeoutExpired(command, 30)
        updated, report = R.reconcile(first, timeout, lambda: "2026-09-08T17:00:00Z")
        self.assertEqual(len(calls), 1)
        self.assertEqual(A.releasable_members(updated["submissions"][0]), set())
        self.assertEqual(len(report["errors"]), 1)

    def test_unknown_submission_without_parent_id_is_not_queried(self):
        ledger = self.ledger()
        ledger["submissions"][0].update(state="submission_unknown", job_id=None)
        def query(command):
            self.fail("Unknown submission must not guess an accounting identity")
        updated, report = R.reconcile(ledger, query, lambda: "2026-09-08T17:00:00Z")
        self.assertNotIn(A.FIELD, updated["submissions"][0])
        self.assertEqual(len(report["unresolved"]), 1)

    def test_main_queries_and_writes_inside_original_project_lock(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "ledger.json"
            ledger.write_text(json.dumps(self.ledger()))
            (root / "campaign.dispatch.binding.json").write_text(json.dumps({"ledger_path": str(ledger)}))
            def query(command, **kwargs):
                self.assertEqual(command, A.sacct_command("123"))
                with (root / "campaign.dispatch.lock").open("a+") as other:
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(other, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return subprocess.CompletedProcess(command, 0, line(0), "")
            original_writer = R.write_json
            def writer(path, value):
                with (root / "campaign.dispatch.lock").open("a+") as other:
                    with self.assertRaises(BlockingIOError):
                        fcntl.flock(other, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return original_writer(path, value)
            with patch.object(R.subprocess, "run", side_effect=query), patch.object(R, "write_json", side_effect=writer):
                self.assertEqual(R.main(["--ledger", str(ledger), "--project-root", str(root)]), 0)
            saved = json.loads(ledger.read_text())
            self.assertEqual(A.releasable_members(saved["submissions"][0]), {0})
            self.assertEqual(A.receipt_identity(saved["submissions"][0]), A.receipt_identity(self.ledger()["submissions"][0]))


if __name__ == "__main__":
    unittest.main()
