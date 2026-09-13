"""Pure-mock orchestration fixtures; no scheduler, simulator or real HPC calls."""
from contextlib import redirect_stdout
from datetime import datetime, timedelta, timezone
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import campaign_cycle as C


TOOLS = ("campaign_cycle.py", "reconcile_ledger.py", "accounted_members.py", "collect_queue.py",
         "manage_campaign.py", "dispatch_campaign.py", "inventory_run.py")


class CampaignCycleTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = (Path(self.temp.name) / "toy_hpc").resolve()
        self.proposal = self.root / "control/campaign_versions/proposal.json"
        self.proposal.parent.mkdir(parents=True)
        self.proposal.write_bytes(b'{ "project": "toy", "runs": [] }\n')
        self.ledger = self.root / "control/submission_ledger.json"
        self.ledger.write_bytes(b'{ "existing_authoritative_ledger": true }\n')
        self.code = self.root / "toolkit/hpc/v2"
        self.code.mkdir(parents=True)
        for name in TOOLS:
            (self.code / name).write_text(f"# Never executed mock tool: {name}\n")
        self.now = datetime(2026, 9, 8, 17, 0, tzinfo=timezone.utc)
        self.planning_time = (self.now - timedelta(seconds=120)).isoformat()
        self.plan_time = (self.now - timedelta(seconds=90)).isoformat()
        self.dispatch_time = (self.now - timedelta(seconds=5)).isoformat()
        # Deliberate formatting/sentinel bytes expose timestamp rewriting or
        # a freshly serialized replacement of the original planning evidence.
        self.queue_bytes = ('{\n "captured_utc": "' + self.planning_time + '", "sentinel": "planning queue"\n}\n').encode()
        self.new_queue_bytes = ('{\n "captured_utc": "' + self.dispatch_time + '", "sentinel": "dispatch queue"\n}\n').encode()
        self.plan_bytes = ('{\n "created_utc": "' + self.plan_time + '", "queue_captured_utc": "' + self.planning_time
                           + '", "submission_plans": [], "runs": [], "issues": []\n}\n').encode()
        self.calls = []
        self.fail_step = None
        self.timeout_step = None

    def fake_subprocess(self, argv, **kwargs):
        self.assertIsInstance(argv, list)
        self.assertEqual(kwargs, {"capture_output": True, "text": True, "timeout": 180, "check": False})
        self.assertEqual(Path(argv[1]).parent, self.code)
        tool = Path(argv[1]).name
        output = Path(argv[argv.index("--output") + 1]) if "--output" in argv else None
        if tool == "collect_queue.py":
            step = "queue" if output.name == "queue.json" else "queue_dispatch"
            self.assertEqual(output.name, "queue.json" if step == "queue" else "queue_dispatch.json")
        else:
            step = {"reconcile_ledger.py": "reconcile", "manage_campaign.py": "plan", "dispatch_campaign.py": "dispatch"}[tool]
        self.calls.append((step, argv[:]))
        if self.timeout_step == step:
            raise subprocess.TimeoutExpired(argv, 180)
        if self.fail_step == step:
            return subprocess.CompletedProcess(argv, 7, "mock failure output\n", "mock failure stderr\n")
        if step == "queue":
            output.write_bytes(self.queue_bytes)
        elif step == "queue_dispatch":
            output.write_bytes(self.new_queue_bytes)
        elif step == "plan":
            queue = Path(argv[argv.index("--queue") + 1])
            self.assertEqual(queue.name, "queue.json")
            self.assertEqual(queue.read_bytes(), self.queue_bytes)
            output.write_bytes(self.plan_bytes)
        elif step == "dispatch":
            plan = Path(argv[argv.index("--plan") + 1])
            queue = Path(argv[argv.index("--queue") + 1])
            self.assertEqual(queue.name, "queue_dispatch.json")
            self.assertEqual(queue.read_bytes(), self.new_queue_bytes)
            self.assertEqual(plan.read_bytes(), self.plan_bytes)
            self.assertEqual((plan.parent / "queue.json").read_bytes(), self.queue_bytes)
            self.assertEqual(argv[argv.index("--plan-max-age-s") + 1], "300")
            self.assertNotIn("--max-age-s", argv)  # live queue keeps dispatcher default freshness
            output.write_text(json.dumps({"status": "mock_submitted" if "--apply" in argv else "validated_dry_run", "outcomes": []}))
        return subprocess.CompletedProcess(argv, 0, f"mock {step} stdout\n", "")

    def run_cycle(self, apply=False):
        argv = ["--campaign", str(self.proposal), "--weekly-used-percent", "17",
                "--usage-observed-utc", self.now.isoformat()] + (["--apply"] if apply else [])
        captured = io.StringIO()
        with patch.object(C, "ROOT", self.root), patch.object(C, "datetime", wraps=datetime) as clock, \
                patch.object(C.subprocess, "run", side_effect=self.fake_subprocess), redirect_stdout(captured):
            clock.now.return_value = self.now
            exit_code = C.main(argv)
        result = json.loads(captured.getvalue())
        return exit_code, result, Path(result["cycle"])

    def test_fresh_dispatch_queue_follows_plan_without_rewriting_old_evidence(self):
        ledger_before, proposal_before = self.ledger.read_bytes(), self.proposal.read_bytes()
        exit_code, result, cycle = self.run_cycle()
        self.assertEqual(exit_code, 0)
        self.assertEqual([name for name, _ in self.calls], ["reconcile", "queue", "plan", "queue_dispatch", "dispatch"])
        self.assertEqual(result["status"], "cycle_finished")
        self.assertEqual(result["dispatch_status"], "validated_dry_run")
        self.assertEqual((cycle / "queue.json").read_bytes(), self.queue_bytes)
        self.assertEqual((cycle / "queue_dispatch.json").read_bytes(), self.new_queue_bytes)
        self.assertEqual((cycle / "plan.json").read_bytes(), self.plan_bytes)
        self.assertEqual(self.ledger.read_bytes(), ledger_before)
        self.assertEqual(self.proposal.read_bytes(), proposal_before)
        saved_input = json.loads((cycle / "input.json").read_text())
        self.assertEqual(saved_input["tools_sha256"], {name: hashlib.sha256((self.code / name).read_bytes()).hexdigest() for name in TOOLS})
        self.assertFalse(saved_input["apply_requested"])
        self.assertEqual(json.loads((cycle / "cycle.json").read_text()), result)
        for step, _ in self.calls:
            self.assertEqual((cycle / (step + ".stdout")).read_text(), f"mock {step} stdout\n")

    def test_apply_is_forwarded_only_when_explicit_and_only_to_dispatch(self):
        self.run_cycle()
        self.assertFalse(any("--apply" in argv for _, argv in self.calls))
        self.calls.clear()
        _, _, cycle = self.run_cycle(apply=True)
        self.assertEqual([name for name, argv in self.calls if "--apply" in argv], ["dispatch"])
        self.assertTrue(json.loads((cycle / "input.json").read_text())["apply_requested"])

    def test_dispatch_queue_collection_failure_never_dispatches(self):
        self.fail_step = "queue_dispatch"
        exit_code, result, cycle = self.run_cycle(apply=True)
        self.assertEqual(exit_code, 3)
        self.assertEqual([name for name, _ in self.calls], ["reconcile", "queue", "plan", "queue_dispatch"])
        self.assertEqual(result["status"], "incomplete")
        self.assertIn("queue_dispatch returned7", result["error"])
        self.assertEqual((cycle / "plan.json").read_bytes(), self.plan_bytes)
        self.assertEqual((cycle / "queue.json").read_bytes(), self.queue_bytes)
        self.assertFalse((cycle / "dispatch.json").exists())
        self.assertEqual((cycle / "queue_dispatch.stderr").read_text(), "mock failure stderr\n")

    def test_dispatch_queue_collection_timeout_never_dispatches_or_retries(self):
        self.timeout_step = "queue_dispatch"
        exit_code, result, cycle = self.run_cycle(apply=True)
        self.assertEqual(exit_code, 3)
        self.assertEqual([name for name, _ in self.calls], ["reconcile", "queue", "plan", "queue_dispatch"])
        self.assertEqual(result["status"], "incomplete")
        self.assertFalse((cycle / "dispatch.json").exists())
        self.assertEqual((cycle / "plan.json").read_bytes(), self.plan_bytes)

    def test_initial_queue_failure_stops_before_plan_and_dispatch(self):
        self.fail_step = "queue"
        exit_code, result, cycle = self.run_cycle(apply=True)
        self.assertEqual(exit_code, 3)
        self.assertEqual([name for name, _ in self.calls], ["reconcile", "queue"])
        self.assertEqual(result["status"], "incomplete")
        self.assertFalse((cycle / "plan.json").exists())
        self.assertFalse((cycle / "dispatch.json").exists())

    def test_plan_failure_does_not_refresh_queue_or_dispatch(self):
        self.fail_step = "plan"
        exit_code, result, cycle = self.run_cycle(apply=True)
        self.assertEqual(exit_code, 3)
        self.assertEqual([name for name, _ in self.calls], ["reconcile", "queue", "plan"])
        self.assertEqual(result["status"], "incomplete")
        self.assertFalse((cycle / "queue_dispatch.json").exists())
        self.assertFalse((cycle / "dispatch.json").exists())


if __name__ == "__main__":
    unittest.main()
