"""Fake-submit fixtures: no real sbatch or scheduler operation is performed."""
import copy
from datetime import datetime, timedelta, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M
import test_inventory_run as fixtures


class DispatchTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.InventoryFixtures()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        self.root = (self.fixture.base/"ebmcal_v2").resolve()
        (self.root/"snapshots").mkdir(parents=True)
        target = self.root/"snapshots/fixture"
        self.fixture.frozen.rename(target)
        self.fixture.frozen = target
        self.fixture.config_path = target/"configs/v2/fixture.json"
        self.fixture.source_root = target/"scripts/v2"
        self.fixture.root = self.root/"runs/fixture"
        script = target/"hpc/v2/run_array.sbatch"
        script.parent.mkdir(parents=True)
        script.write_text("\n".join(["#!/bin/bash", "#SBATCH --nodes=1", "#SBATCH --ntasks=1", "#SBATCH --signal=B:USR1@1500",
            f"root={self.root}", "cell=$(( index / nchunks ))", "chunk=$(( index % nchunks ))",
            "python scripts/v2/run_v2.py --time-budget-s 11800"]))
        files = sorted(p for p in target.rglob("*") if p.is_file())
        sha = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(target)}\n" for p in files)
        (target/"SHA256SUMS").write_text(sha)
        self.now = datetime(2026, 9, 8, tzinfo=timezone.utc)
        self.queue = {"schema_version": "campaign_queue_v1", "captured_utc": self.now.isoformat(),
                      "scope": "all_user_jobs", "complete": True, "jobs": []}
        self.ledger = {"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": []}
        campaign = {"project": "ebmcal_v2", "max_project_cpus": 256,
                    "runs": [{"run_id": "fixture", "snapshot": str(target),
                              "snapshot_sha256": hashlib.sha256(sha.encode()).hexdigest(),
                              "nchunks": 2, "output_root": str(self.fixture.root)}]}
        self.plan = M.make_plan(campaign, self.queue, self.ledger, self.now, weekly_used=12)
        self.assertEqual(len(self.plan["submission_plans"]), 1, self.plan)
        self.plan_path = self.root/"plan.json"
        self.queue_path = self.root/"queue.json"
        self.ledger_path = self.root/"ledger.json"
        self.save_inputs()
        self.calls = []

    def save_inputs(self):
        for path, value in ((self.plan_path, self.plan), (self.queue_path, self.queue), (self.ledger_path, self.ledger)):
            path.write_text(json.dumps(value))

    def success(self, argv, **kwargs):
        saved = json.loads(self.ledger_path.read_text())
        self.assertEqual(saved["submissions"][-1]["state"], "intent")
        self.calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0, "12345\n", "")

    def dispatch(self, apply=False, submitter=None, **kwargs):
        return D.dispatch(self.plan_path, self.ledger_path, self.queue_path, self.root,
                          apply=apply, now_func=lambda: self.now,
                          submitter=submitter or self.success, **kwargs)

    def test_default_dry_run_validates_without_sbatch_or_ledger_mutation(self):
        before = self.ledger_path.read_bytes()
        result = self.dispatch()
        self.assertEqual(result["status"], "validated_dry_run")
        self.assertEqual(self.calls, [])
        self.assertEqual(self.ledger_path.read_bytes(), before)
        self.assertFalse((self.root/"campaign.dispatch.binding.json").exists())

    def test_success_persists_intent_first_and_repeated_plan_is_blocked(self):
        result = self.dispatch(apply=True)
        self.assertEqual(result["status"], "submitted")
        receipt = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(receipt["job_id"], "12345")
        self.assertEqual([v["state"] for v in receipt["history"]], ["intent", "submitted"])
        with self.assertRaisesRegex(ValueError, "already has a ledger"):
            self.dispatch(apply=True)
        self.assertEqual(len(self.calls), 1)

    def test_timeout_records_unknown_and_retains_partial_stdout(self):
        def timeout(argv, **kwargs):
            raise subprocess.TimeoutExpired(argv, 30, output=b"12345\n")
        result = self.dispatch(apply=True, submitter=timeout)
        self.assertEqual(result["status"], "submission_unknown_reconciliation_required")
        saved = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(saved["state"], "submission_unknown")
        self.assertEqual(saved["stdout"], "12345\n")
        with self.assertRaisesRegex(ValueError, "already has a ledger"):
            self.dispatch(apply=True)

    def test_non_pure_jobid_and_nonzero_exit_are_uncertain(self):
        for returncode, stdout in ((0, "12345;cluster\n"), (0, "Submitted batch job 12345"), (1, "")):
            with self.subTest(returncode=returncode, stdout=stdout):
                self.ledger_path.write_text(json.dumps(self.ledger))
                fake = lambda argv, **kw: subprocess.CompletedProcess(argv, returncode, stdout, "diagnostic")
                result = self.dispatch(apply=True, submitter=fake)
                self.assertEqual(result["outcomes"][0]["state"], "submission_unknown")

    def test_environment_and_shell_do_not_add_sbatch_options(self):
        with mock.patch.dict("os.environ", {"SBATCH_PARTITION": "elsewhere", "SLURM_CLUSTERS": "other"}):
            self.dispatch(apply=True)
        _, kwargs = self.calls[0]
        self.assertIs(kwargs["shell"], False)
        self.assertNotIn("SBATCH_PARTITION", kwargs["env"])
        self.assertNotIn("SLURM_CLUSTERS", kwargs["env"])

    def test_changed_command_is_rejected_before_intent(self):
        self.plan["submission_plans"][0]["argv"][0] = "sh"
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "Only sbatch"):
            self.dispatch(apply=True)
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])

    def test_extra_argument_cannot_escape_the_allowlist(self):
        self.plan["submission_plans"][0]["argv"].insert(1, "--wrap=echo unsafe")
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "12-argument"):
            self.dispatch(apply=True)

    def test_snapshot_outside_project_is_rejected(self):
        self.plan["submission_plans"][0]["argv"][9] = "/other/snapshot"
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "outside"):
            self.dispatch(apply=True)

    def test_snapshot_changed_after_plan_is_rejected(self):
        (self.fixture.source_root/"engine_v2.py").write_text("# changed\n")
        with self.assertRaisesRegex(ValueError, "byte mismatch"):
            self.dispatch(apply=True)

    def test_stop_file_and_usage_threshold_block_before_intent(self):
        (self.root/"STOP").touch()
        with self.assertRaisesRegex(ValueError, "STOP"):
            self.dispatch(apply=True)
        (self.root/"STOP").unlink()
        with self.assertRaisesRegex(ValueError, "below 90"):
            self.dispatch(apply=True, weekly_used=90)

    def test_stale_plan_or_queue_is_rejected(self):
        self.now += timedelta(seconds=61)
        with self.assertRaisesRegex(ValueError, "stale"):
            self.dispatch(apply=True)

    def test_ledger_is_reread_and_unknown_submission_consumes_capacity(self):
        intent = copy.deepcopy(self.plan["submission_plans"][0]["ledger_intent"])
        intent.update(run_id="another_run", array_task_ids=list(range(16)), max_concurrent=16,
                      plan_id="f"*64, state="submission_unknown")
        self.ledger["submissions"].append(intent)
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "exceed cap"):
            self.dispatch(apply=True)

    def test_fresh_other_project_job_slots_block_before_intent(self):
        self.plan["user_job_budget"]["cap"] = 2
        self.queue["jobs"] = [{"job_id": "900", "state": "PENDING", "project": "phase-test",
                                "array_task_ids": None, "cpus_per_task": 1}]
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "user submitted-job members"):
            self.dispatch(apply=True)
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])

    def test_submit_slot_recheck_uses_array_cardinality_not_cpu_throttle(self):
        self.plan["user_job_budget"]["cap"] = 601
        self.queue["jobs"] = [{"job_id": "900", "state": "PENDING", "project": "phase-test",
            "array_task_ids": list(range(600)), "cpus_per_task": 1,
            "max_concurrent": 1, "running_tasks": 0}]
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "600\\+2 exceed cap 601"):
            self.dispatch(apply=True)
        self.assertEqual(self.calls, [])

    def test_fresh_unknown_receipt_consumes_user_slots(self):
        self.plan["user_job_budget"]["cap"] = 2
        intent = copy.deepcopy(self.plan["submission_plans"][0]["ledger_intent"])
        intent.update(run_id="another_run", array_task_ids=[0], max_concurrent=1,
                      plan_id="f"*64, state="submission_unknown")
        self.ledger["submissions"] = [intent]
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "user submitted-job members"):
            self.dispatch(apply=True)
        self.assertEqual(len(json.loads(self.ledger_path.read_text())["submissions"]), 1)

    def test_user_slot_limit_is_rechecked_after_expensive_inventory(self):
        self.plan["user_job_budget"]["cap"] = 2
        self.save_inputs()
        original = D.recheck_work
        count = 0
        def inspect(value):
            nonlocal count
            original(value)
            count += 1
            if count == 2:
                self.queue["jobs"] = [{"job_id": "900", "state": "PENDING", "project": "phase-test",
                    "array_task_ids": None, "cpus_per_task": 1}]
                self.queue_path.write_text(json.dumps(self.queue))
        with mock.patch.object(D, "recheck_work", side_effect=inspect):
            result = self.dispatch(apply=True)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertIn("user submitted-job members", result["outcomes"][0]["reason"])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])
        self.assertEqual(self.calls, [])

    def test_user_slot_cap_must_be_positive_integer(self):
        self.plan["user_job_budget"]["cap"] = True
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "max_user_submitted_jobs"):
            self.dispatch(apply=True)

    def test_active_target_in_fresh_queue_blocks_duplicate(self):
        self.queue["jobs"] = [{"job_id": "123", "project": "ebmcal_v2", "run_id": "fixture",
            "state": "PENDING", "array_task_ids": [0], "cpus_per_task": 16,
            "max_concurrent": 1, "running_tasks": 0,
            "snapshot_sha256": self.plan["submission_plans"][0]["snapshot_sha256"]}]
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "active"):
            self.dispatch(apply=True)

    def test_newly_complete_chunk_requires_a_fresh_plan(self):
        self.fixture.write_chunk(0, 2)
        with self.assertRaisesRegex(ValueError, "now complete"):
            self.dispatch(apply=True)

    def test_project_lock_excludes_concurrent_dispatcher(self):
        with (self.root/"campaign.dispatch.lock").open("a") as stream:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaises(BlockingIOError):
                self.dispatch(apply=True)
        self.assertEqual(self.calls, [])

    def test_project_binding_prevents_switching_to_empty_ledger(self):
        self.dispatch(apply=True)
        other = self.root/"other_ledger.json"
        other.write_text(json.dumps(self.ledger))
        with self.assertRaisesRegex(ValueError, "different ledger"):
            D.dispatch(self.plan_path, other, self.queue_path, self.root,
                       apply=True, now_func=lambda: self.now, submitter=self.success)

    def test_intent_survives_exception_before_returned_receipt(self):
        def unexpected(argv, **kwargs):
            raise RuntimeError("unexpected transport wrapper failure")
        result = self.dispatch(apply=True, submitter=unexpected)
        self.assertEqual(result["outcomes"][0]["state"], "submission_unknown")

    def two_plans(self):
        base = self.plan["submission_plans"][0]
        submissions = []
        for index in (0, 1):
            submission = copy.deepcopy(base)
            submission.update(array_task_ids=[index], max_concurrent=1)
            identity = {name: submission[name] for name in D.IDENTITY_FIELDS}
            submission["plan_id"] = M.I.digest(identity)
            submission["argv"][2] = f"--array={index}%1"
            submission["argv"][5] = f"--comment=ebmcal_v2:{submission['plan_id']}"
            submission["ledger_intent"] = identity | {"plan_id": submission["plan_id"], "state": "intent",
                "created_utc": self.now.isoformat()}
            submissions.append(submission)
        self.plan["submission_plans"] = submissions
        self.save_inputs()

    def test_new_other_project_job_between_arrays_stops_second_intent(self):
        self.plan["user_job_budget"]["cap"] = 2
        self.two_plans()
        def first_success(argv, **kwargs):
            response = self.success(argv, **kwargs)
            self.queue["jobs"] = [{"job_id": "900", "state": "PENDING", "project": "phase-test",
                "array_task_ids": None, "cpus_per_task": 1}]
            self.queue_path.write_text(json.dumps(self.queue))
            return response
        result = self.dispatch(apply=True, submitter=first_success)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertEqual(len(self.calls), 1)
        receipts = json.loads(self.ledger_path.read_text())["submissions"]
        self.assertEqual(len(receipts), 1)
        self.assertEqual(receipts[0]["state"], "submitted")
        self.assertIn("user submitted-job members", result["outcomes"][1]["reason"])

    def test_uncertain_first_submission_stops_without_intenting_second(self):
        self.two_plans()
        def uncertain(argv, **kwargs):
            self.calls.append(argv)
            return subprocess.CompletedProcess(argv, 1, "", "unknown")
        result = self.dispatch(apply=True, submitter=uncertain)
        self.assertEqual(result["status"], "submission_unknown_reconciliation_required")
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(len(json.loads(self.ledger_path.read_text())["submissions"]), 1)

    def test_first_submit_can_age_plan_out_before_second(self):
        self.two_plans()
        def first(argv, **kwargs):
            self.now += timedelta(seconds=61)
            return subprocess.CompletedProcess(argv, 0, "12345", "")
        result = self.dispatch(apply=True, submitter=first)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertEqual([v["state"] for v in result["outcomes"]], ["submitted", "not_dispatched"])
        self.assertEqual(len(json.loads(self.ledger_path.read_text())["submissions"]), 1)

    def test_duplicate_open_job_id_is_uncertain(self):
        self.two_plans()
        result = self.dispatch(apply=True)
        self.assertEqual(result["status"], "submission_unknown_reconciliation_required")
        self.assertEqual([v["state"] for v in result["outcomes"]], ["submitted", "submission_unknown"])

    def test_long_inventory_cannot_submit_an_expired_plan(self):
        original = D.recheck_work
        calls = 0
        def slow(value):
            nonlocal calls
            original(value)
            calls += 1
            if calls == 2:
                self.now += timedelta(seconds=61)
        with mock.patch.object(D, "recheck_work", side_effect=slow):
            result = self.dispatch(apply=True)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])

    def test_slow_durable_intent_is_closed_with_no_submit_evidence(self):
        original = D.durable_json
        def slow(path, value):
            original(path, value)
            if path == self.ledger_path and value["submissions"][-1]["state"] == "intent":
                self.now += timedelta(seconds=61)
        with mock.patch.object(D, "durable_json", side_effect=slow):
            result = self.dispatch(apply=True)
        self.assertEqual(result["status"], "stopped_before_subprocess")
        self.assertEqual(self.calls, [])
        receipt = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(receipt["state"], "not_submitted_verified")
        self.assertIs(receipt["reconciliation_evidence"]["subprocess_invoked"], False)

    def test_config_alias_cannot_change_the_logical_target_key(self):
        value = copy.deepcopy(self.plan["submission_plans"][0])
        value["run_id"] = "alias"
        identity = {name: value[name] for name in D.IDENTITY_FIELDS}
        value["plan_id"] = M.I.digest(identity)
        report = copy.deepcopy(self.plan["runs"][0])
        report["run_id"] = "alias"
        with self.assertRaisesRegex(ValueError, "config_name must equal run_id"):
            D.validate_submission(value, report, self.root)

    def test_short_topology_override_cannot_increase_hidden_cpu_reservation(self):
        batch = "#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --mem=12G -n8\n"
        with self.assertRaisesRegex(ValueError, "ntasks=1"):
            D.check_batch_topology(batch)

    def test_cli_output_cannot_overwrite_ledger(self):
        args = ["--plan", str(self.plan_path), "--queue", str(self.queue_path),
                "--ledger", str(self.ledger_path), "--project-root", str(self.root),
                "--output", str(self.ledger_path)]
        with self.assertRaises(SystemExit), mock.patch("sys.stderr"):
            D.main(args)

    def age_proposal(self, seconds=90):
        old = (self.now-timedelta(seconds=seconds)).isoformat()
        self.plan.update(created_utc=old, queue_captured_utc=old)
        self.save_inputs()

    def test_explicit_older_proposal_uses_fresh_queue_without_rewriting_evidence(self):
        self.age_proposal()
        plan_before = self.plan_path.read_bytes()
        result = self.dispatch(apply=True, plan_max_age_s=300)
        self.assertEqual(result["status"], "submitted")
        self.assertEqual(self.plan_path.read_bytes(), plan_before)
        receipt = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(receipt["plan_created_utc"], self.plan["created_utc"])
        self.assertEqual(receipt["planning_queue_captured_utc"], self.plan["queue_captured_utc"])
        self.assertEqual(receipt["dispatch_queue_captured_utc"], self.queue["captured_utc"])
        self.assertEqual(receipt["freshness_limits_s"], {"plan_max_age_s": 300, "live_queue_max_age_s": 60})

    def test_default_proposal_age_remains_sixty_seconds(self):
        self.age_proposal(seconds=61)
        with self.assertRaisesRegex(ValueError, "Plan is stale"):
            self.dispatch(apply=True)
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])

    def test_old_proposal_exception_does_not_accept_stale_live_queue(self):
        self.age_proposal()
        self.queue["captured_utc"] = (self.now-timedelta(seconds=61)).isoformat()
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "Queue snapshot is stale"):
            self.dispatch(apply=True, plan_max_age_s=300)
        self.assertEqual(self.calls, [])

    def test_explicit_limit_rejects_future_proposal_or_original_queue_time(self):
        for key in ("created_utc", "queue_captured_utc"):
            with self.subTest(key=key):
                self.age_proposal()
                self.plan[key] = (self.now+timedelta(seconds=1)).isoformat()
                self.save_inputs()
                with self.assertRaisesRegex(ValueError, "future-dated"):
                    self.dispatch(plan_max_age_s=300)

    def test_proposal_and_original_planning_queue_each_have_the_300s_cap(self):
        for key in ("created_utc", "queue_captured_utc"):
            with self.subTest(key=key):
                self.age_proposal()
                self.plan[key] = (self.now-timedelta(seconds=301)).isoformat()
                self.save_inputs()
                with self.assertRaisesRegex(ValueError, "stale"):
                    self.dispatch(plan_max_age_s=300)

    def test_new_queue_must_not_predate_original_planning_evidence(self):
        self.age_proposal(seconds=10)
        self.queue["captured_utc"] = (self.now-timedelta(seconds=11)).isoformat()
        self.save_inputs()
        with self.assertRaisesRegex(ValueError, "predates"):
            self.dispatch(plan_max_age_s=300)

    def test_plan_limit_cannot_exceed_300_or_expand_live_queue_limit(self):
        for limit in (0, -1, 301, float("nan"), float("inf"), True):
            with self.subTest(limit=limit), self.assertRaisesRegex(ValueError, "Plan freshness"):
                self.dispatch(plan_max_age_s=limit)
        with self.assertRaisesRegex(ValueError, "freshness limit"):
            self.dispatch(plan_max_age_s=300, max_age_s=61)

    def test_old_proposal_live_queue_rechecked_after_inventory_before_intent(self):
        self.age_proposal()
        original, count = D.recheck_work, 0
        def slow(value):
            nonlocal count
            original(value)
            count += 1
            if count == 2:
                self.now += timedelta(seconds=61)
        with mock.patch.object(D, "recheck_work", side_effect=slow):
            result = self.dispatch(apply=True, plan_max_age_s=300)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertIn("Queue snapshot is stale", result["outcomes"][0]["reason"])
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])

    def test_old_proposal_stale_queue_after_intent_fsync_never_invokes_sbatch(self):
        self.age_proposal()
        original = D.durable_json
        def slow(path, value):
            original(path, value)
            if path == self.ledger_path and value["submissions"][-1]["state"] == "intent":
                self.now += timedelta(seconds=61)
        with mock.patch.object(D, "durable_json", side_effect=slow):
            result = self.dispatch(apply=True, plan_max_age_s=300)
        self.assertEqual(result["status"], "stopped_before_subprocess")
        self.assertEqual(self.calls, [])
        receipt = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(receipt["state"], "not_submitted_verified")
        self.assertIs(receipt["reconciliation_evidence"]["subprocess_invoked"], False)

    def test_proposal_age_rechecked_after_inventory_even_when_queue_is_refreshed(self):
        self.age_proposal(seconds=250)
        original, count = D.recheck_work, 0
        def slow(value):
            nonlocal count
            original(value)
            count += 1
            if count == 2:
                self.now += timedelta(seconds=51)
                self.queue["captured_utc"] = self.now.isoformat()
                self.queue_path.write_text(json.dumps(self.queue))
        with mock.patch.object(D, "recheck_work", side_effect=slow):
            result = self.dispatch(apply=True, plan_max_age_s=300)
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertIn("Plan is stale", result["outcomes"][0]["reason"])
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])


if __name__ == "__main__":
    unittest.main()
