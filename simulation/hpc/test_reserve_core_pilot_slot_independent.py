"""Independent bounded helper regressions; subprocess is entirely mocked."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M

PROJECT = Path(__file__).resolve().parents[2]
PARENT = "10265062"
PLAN = "0661eee6714148115c1ae9789e84ec63590e425416276a0a0302017ab83005a6"
MEMBERS = [25, *range(430, 445)]


class IndependentReserveTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="core_throttle_independent_")
        self.addCleanup(temporary.cleanup)
        self.root = (Path(temporary.name) / "ebmcal_v2").resolve()
        self.ledger = self.root / "control/submission_ledger.json"
        self.ledger.parent.mkdir(parents=True)
        self.evidence = self.root / "control/campaign_versions/confirmation_20260908_core_paired_1h_pilot/reserve_core_pilot_slot.json"
        self.snapshot = self.root / "snapshots/confirm_core_a"
        self.snapshot.mkdir(parents=True)
        self.cfg = self.snapshot / "config.json"
        source_snapshot = PROJECT / "runs/v2/snapshots/confirm_core_a"
        self.cfg.write_bytes((source_snapshot / "configs/v2/confirm_core_a.json").read_bytes())
        (self.snapshot / "SHA256SUMS").write_bytes((source_snapshot / "SHA256SUMS").read_bytes())
        receipt = {"run_id": "confirm_core_a", "snapshot_sha256": hashlib.sha256((self.snapshot / "SHA256SUMS").read_bytes()).hexdigest(),
                   "array_task_ids": MEMBERS, "nchunks": 500, "mapping": "chunk_major", "cpus_per_task": 16,
                   "max_concurrent": 3, "walltime_s": 7200, "dispatch_generation": 5,
                   "state": "submitted", "plan_id": PLAN, "job_id": PARENT}
        self.assertEqual(M.I.digest(D.submission_identity(receipt)), PLAN)
        self.ledger.write_text(json.dumps({"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": [receipt]}))
        (self.root / "campaign.dispatch.binding.json").write_text(json.dumps({"project": "ebmcal_v2", "ledger_path": str(self.ledger)}))
        source = (PROJECT / "hpc/v2/reserve_core_pilot_slot_20260908.sh").read_text()
        self.helper = self.root / "helper.sh"
        self.helper.write_text(source)
        body = source.split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0]
        body = body.replace("Path('/N/slate/tg11/ebmcal_v2').resolve()", f"Path({str(self.root)!r}).resolve()", 1)
        self.code = compile(body, "independent_inert_throttle_helper", "exec", optimize=2)
        self.calls = []
        self.control_reads = 0
        self.queue = "".join(f"{PARENT}_{member}|PENDING|2:00:00|ebmcal_v2:{PLAN}\n" for member in MEMBERS)
        self.control_mutation = lambda body: body
        self.post_timeout = False
        self.after_intent = None
        self.did_mutate = False

    def process(self, argv, **kwargs):
        self.calls.append(argv[:])
        if argv == ["squeue", "-r", "-j", PARENT, "-h", "-o", "%i|%T|%l|%k"]:
            return subprocess.CompletedProcess(argv, 0, self.queue, "")
        if argv == ["scontrol", "show", "job", PARENT]:
            self.control_reads += 1
            if self.control_reads > 1 and self.post_timeout:
                raise subprocess.TimeoutExpired(argv, 30)
            throttle = 3 if self.control_reads == 1 else 2
            body = (f"JobId={PARENT} ArrayJobId={PARENT} ArrayTaskId=25,430-444%{throttle} "
                    f"ArrayTaskThrottle={throttle} JobState=PENDING Comment=ebmcal_v2:{PLAN} "
                    "TimeLimit=2:00:00 NumCPUs=16")
            return subprocess.CompletedProcess(argv, 0, self.control_mutation(body), "")
        self.assertEqual(argv, ["scontrol", "update", "JobId=10265062", "ArrayTaskThrottle=2"])
        self.assertEqual(json.loads(self.evidence.read_text())["status"], "update_unknown")
        self.assertIsNone(json.loads(self.evidence.read_text())["subprocess_invoked"])
        return subprocess.CompletedProcess(argv, 0, "", "")

    def execute(self):
        real_fsync = os.fsync
        def fsync(fd):
            real_fsync(fd)
            if (self.after_intent and not self.did_mutate and self.evidence.exists()
                    and json.loads(self.evidence.read_text())["status"] == "update_unknown"):
                self.did_mutate = True
                self.after_intent()
        with mock.patch.object(M, "snapshot_info", return_value={"config_path": self.cfg}), \
                mock.patch.object(M, "validate_execution_walltime", return_value=None), \
                mock.patch.object(sys, "argv", ["-", str(self.helper)]), \
                mock.patch.object(sys, "path", sys.path[:]), \
                mock.patch("subprocess.run", side_effect=self.process), \
                mock.patch("os.fsync", side_effect=fsync), mock.patch("builtins.print"):
            exec(self.code, {"__name__": "__main__"})

    def updates(self):
        return [argv for argv in self.calls if argv[:2] == ["scontrol", "update"]]

    def test_wrong_controller_parent_and_unplanned_queue_member_never_mutate(self):
        self.control_mutation = lambda body: body.replace(f"JobId={PARENT} ", "JobId=99999999 ", 1)
        original = self.ledger.read_bytes()
        with self.assertRaises(SystemExit):
            self.execute()
        self.assertEqual(self.updates(), [])
        self.assertEqual(self.ledger.read_bytes(), original)
        self.assertEqual(json.loads(self.evidence.read_text())["status"], "skipped_no_change")
        # A separate fresh attempt fixture, without any prior scheduler mutation.
        self.evidence.unlink(); self.calls.clear(); self.control_reads = 0
        self.control_mutation = lambda body: body
        self.queue += f"{PARENT}_445|PENDING|2:00:00|ebmcal_v2:{PLAN}\n"
        with self.assertRaises(SystemExit):
            self.execute()
        self.assertEqual(self.updates(), [])
        self.assertEqual(self.ledger.read_bytes(), original)

    def test_helper_source_changes_after_intent_block_before_scheduler_update(self):
        original = self.ledger.read_bytes()
        self.after_intent = lambda: self.helper.write_text(self.helper.read_text() + "\n# inert race\n")
        with self.assertRaises(SystemExit):
            self.execute()
        value = json.loads(self.evidence.read_text())
        self.assertEqual(value["status"], "skipped_stale_or_changed_no_change")
        self.assertIs(value["sources_and_ledger_unchanged"], False)
        self.assertIs(value["subprocess_invoked"], False)
        self.assertEqual(self.updates(), [])
        self.assertEqual(self.ledger.read_bytes(), original)

    def test_ledger_changes_after_intent_are_preserved_and_block_update(self):
        changed_bytes = self.ledger.read_bytes() + b"\n"
        self.after_intent = lambda: self.ledger.write_bytes(changed_bytes)
        with self.assertRaises(SystemExit):
            self.execute()
        self.assertEqual(json.loads(self.evidence.read_text())["status"], "skipped_stale_or_changed_no_change")
        self.assertEqual(self.ledger.read_bytes(), changed_bytes)
        self.assertEqual(self.updates(), [])

    def test_successful_update_with_postread_timeout_is_uncertain_and_not_retried(self):
        identity = D.submission_identity(json.loads(self.ledger.read_text())["submissions"][0])
        self.post_timeout = True
        with self.assertRaisesRegex(ValueError, "no automatic retry"):
            self.execute()
        value = json.loads(self.evidence.read_text())
        self.assertEqual((value["returncode"], value["status"]), (0, "needs_reconciliation"))
        self.assertEqual(value["exception_type"], "TimeoutExpired")
        saved = json.loads(self.ledger.read_text())["submissions"][0]
        self.assertEqual(D.submission_identity(saved), identity)
        self.assertEqual(saved["execution_adjustments"][0]["state"], "needs_reconciliation")
        count = len(self.calls)
        with self.assertRaisesRegex(ValueError, "One-shot evidence exists"):
            self.execute()
        self.assertEqual(len(self.calls), count)
        self.assertEqual(len(self.updates()), 1)


if __name__ == "__main__":
    unittest.main()
