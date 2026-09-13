"""Small continuation, duplicate-dispatch, layout and capacity fixtures."""
import copy
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import inventory_run as I
import manage_campaign as M
import test_inventory_run as fixtures


class CampaignTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.InventoryFixtures()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        self.now = datetime(2026, 9, 8, tzinfo=timezone.utc)
        self.queue = {"schema_version": "campaign_queue_v1", "captured_utc": self.now.isoformat(),
                      "scope": "all_user_jobs", "complete": True, "jobs": []}
        self.ledger = {"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": []}
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 256,
                         "runs": [{"run_id": "fixture", "snapshot": str(self.fixture.frozen),
                                   "remote_snapshot": "/remote/snapshots/fixture", "nchunks": 2,
                                   "output_root": str(self.fixture.root),
                                   "remote_output_root": "/remote/ebmcal_v2/runs/fixture"}]}
        self.script = self.fixture.frozen/"hpc/v2/run_array.sbatch"
        self.script.parent.mkdir(parents=True)
        self.script.write_text("\n".join(["#!/bin/bash", "#SBATCH --signal=B:USR1@1500",
            "root=/remote/ebmcal_v2", "cell=$(( index / nchunks ))", "chunk=$(( index % nchunks ))",
            "python scripts/v2/run_v2.py --time-budget-s 11800"]))
        self.seal()

    def seal(self):
        paths = sorted(p for p in self.fixture.frozen.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        manifest = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(self.fixture.frozen)}\n" for p in paths)
        (self.fixture.frozen/"SHA256SUMS").write_text(manifest)
        self.campaign["runs"][0]["snapshot_sha256"] = hashlib.sha256(manifest.encode()).hexdigest()

    def plan(self, **kwargs):
        return M.make_plan(self.campaign, self.queue, self.ledger, self.now, **kwargs)

    def active_array(self, task_ids=(0,), cpus=16, throttle=1):
        return {"job_id": "1234", "state": "MIXED", "project": "ebmcal_v2", "run_id": "fixture",
                "array_task_ids": list(task_ids), "cpus_per_task": cpus,
                "max_concurrent": throttle, "running_tasks": 0,
                "snapshot_sha256": self.campaign["runs"][0]["snapshot_sha256"]}

    def test_unstarted_chunks_are_planned_and_original_layout_is_detected(self):
        p = self.plan()
        self.assertEqual(p["runs"][0]["mapping"], "cell_major")
        self.assertEqual(p["submission_plans"][0]["array_task_ids"], [0, 1])
        self.assertEqual(p["runs"][0]["chunks"][0]["seeds"], [100, 102])
        self.assertEqual(p["cpu_budget"]["new_reserved"], 32)
        self.assertFalse(self.fixture.root.exists())

    def test_new_mapping_and_frozen_registration_mismatch(self):
        self.script.write_text(self.script.read_text().replace("index / nchunks", "index % ncells").replace("index % nchunks", "index / ncells"))
        self.seal()
        self.assertEqual(self.plan()["runs"][0]["mapping"], "chunk_major")
        self.assertEqual(M.array_index(1, 2, 3, 4, "chunk_major"), 7)
        self.assertEqual(M.array_index(1, 2, 3, 4, "cell_major"), 6)
        self.campaign["runs"][0]["mapping"] = "cell_major"
        self.assertEqual(self.plan()["submission_plans"], [])

    def test_snapshot_byte_change_blocks_resume(self):
        self.script.write_text(self.script.read_text()+"\n# changed\n")
        p = self.plan()
        self.assertEqual(p["runs"][0]["state"], "blocked_validation")
        self.assertEqual(p["submission_plans"], [])

    def test_complete_chunks_are_not_replayed(self):
        self.fixture.write_chunk(0, 2)
        self.fixture.write_chunk(1, 2)
        p = self.plan()
        self.assertTrue(p["all_registered_runs_complete"])
        self.assertEqual(p["submission_plans"], [])

    def test_done_flag_does_not_hide_missing_engine_rows(self):
        rows = [self.fixture.make_row(seed, "original") for seed in (100, 102)]
        self.fixture.write_chunk(0, 2, rows=rows)
        self.fixture.write_chunk(1, 2)
        p = self.plan()
        self.assertEqual(p["submission_plans"][0]["array_task_ids"], [0])

    def test_terminal_failures_are_retained_without_perpetual_retries(self):
        rows = [self.fixture.make_row(seed, engine) for seed in (100, 102) for engine in self.fixture.config["engines"]]
        rows[0].update(status="error", observed={})
        self.fixture.write_chunk(0, 2, rows=rows)
        self.fixture.write_chunk(1, 2)
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertEqual(p["runs"][0]["chunks"][0]["state"], "terminal_with_failures")
        self.assertEqual(p["runs"][0]["chunks"][0]["failed_observed"], 1)
        self.assertFalse(p["all_registered_runs_complete"])

    def test_pending_job_protects_target_and_only_project_capacity_is_limited(self):
        self.queue["jobs"] = [self.active_array(), {"job_id": "2222", "state": "RUNNING",
            "project": "phase-test", "cpus_per_task": 1024, "array_task_ids": None}]
        p = self.plan()
        self.assertEqual(p["submission_plans"][0]["array_task_ids"], [1])
        self.assertEqual(p["cpu_budget"]["project_reserved"], 16)
        self.assertEqual(p["cpu_budget"]["all_user_reserved"], 1040)

    def test_array_reservation_uses_throttle_not_total_array_size(self):
        self.fixture.config.update(datasets=1000, cells=[{"name": "REF"}, {"name": "IID"}])
        self.fixture.save_config()
        self.campaign.update(max_project_cpus=512)
        self.campaign["runs"][0].update(nchunks=500, max_parallel_tasks=32)
        self.seal()
        self.queue["jobs"] = [self.active_array(range(1000), cpus=16, throttle=32)]
        p = self.plan()
        self.assertEqual(p["cpu_budget"]["project_reserved"], 512)
        self.assertEqual(p["submission_plans"], [])

    def test_user_submit_slots_count_all_array_members_and_other_project_jobs(self):
        self.fixture.config.update(datasets=1000, cells=[{"name": "REF"}, {"name": "IID"}])
        self.fixture.save_config()
        self.campaign.update(max_project_cpus=512, max_user_submitted_jobs=1000,
                             max_tasks_per_submission=1000)
        self.campaign["runs"][0].update(nchunks=500, max_parallel_tasks=32)
        self.seal()
        other = self.active_array(range(900), cpus=1, throttle=2)
        other.update(project="phase-test", run_id="other")
        self.queue["jobs"] = [other] + [{"job_id": str(9000+i), "project": "phase-test",
            "state": "RUNNING", "cpus_per_task": 1, "array_task_ids": None} for i in range(3)]
        p = self.plan()
        s = p["submission_plans"][0]
        self.assertEqual(s["array_task_ids"], list(range(97)))
        self.assertEqual((s["nchunks"], s["max_concurrent"]), (500, 32))
        self.assertEqual(p["cpu_budget"]["project_reserved"], 0)
        self.assertEqual(p["cpu_budget"]["new_reserved"], 512)
        self.assertEqual(p["user_job_budget"]["queue_active_members"], 903)
        self.assertEqual(p["user_job_budget"]["remaining_after_plan"], 0)

    def test_optional_user_cap_absent_preserves_plans_and_still_reports_counts(self):
        self.queue["jobs"] = [{"job_id": "99", "state": "PENDING", "project": "phase-test",
                                "array_task_ids": None, "cpus_per_task": 1000}]
        p = self.plan()
        self.assertEqual(p["submission_plans"][0]["array_task_ids"], [0, 1])
        self.assertIsNone(p["user_job_budget"]["cap"])
        self.assertEqual(p["user_job_budget"]["queue_active_members"], 1)

    def test_exhausted_user_cap_does_not_change_or_cancel_other_jobs(self):
        self.campaign["max_user_submitted_jobs"] = 1
        self.queue["jobs"] = [{"job_id": "99", "state": "PENDING", "project": "phase-test",
                                "array_task_ids": None, "cpus_per_task": 1}]
        before = copy.deepcopy(self.queue)
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertEqual(p["cpu_budget"]["project_reserved"], 0)
        self.assertEqual(self.queue, before)
        self.assertEqual(p["runs"][0]["state_counts"], {"eligible_resume": 2})

    def test_slot_limited_batches_keep_layout_and_skip_open_first_batch(self):
        self.campaign["max_user_submitted_jobs"] = 1
        first = self.plan()["submission_plans"][0]
        self.assertEqual(first["array_task_ids"], [0])
        self.ledger["submissions"] = [first["ledger_intent"]]
        self.assertEqual(self.plan()["submission_plans"], [])
        self.campaign["max_user_submitted_jobs"] = 2
        p = self.plan()
        second = p["submission_plans"][0]
        self.assertEqual(second["array_task_ids"], [1])
        self.assertEqual(second["nchunks"], first["nchunks"])
        self.assertEqual(p["runs"][0]["chunks"][1]["seeds"], [101, 103])
        self.assertNotEqual(second["plan_id"], first["plan_id"])
        self.assertEqual(p["user_job_budget"]["unresolved_submission_reserved"], 1)

    def test_partial_visible_receipt_reserves_only_unseen_members_in_addition_to_queue(self):
        first = self.plan()["submission_plans"][0]
        self.ledger["submissions"] = [first["ledger_intent"] | {"state": "submitted", "job_id": "1234"}]
        self.queue["jobs"] = [self.active_array([1], throttle=2)]
        self.campaign["max_user_submitted_jobs"] = 3
        p = self.plan()
        budget = p["user_job_budget"]
        self.assertEqual(budget["queue_active_members"], 1)
        self.assertEqual(budget["unresolved_submission_reserved"], 1)
        self.assertEqual(budget["remaining_before_plan"], 1)
        self.queue["jobs"][0]["array_task_ids"] = [0, 1]
        self.assertEqual(self.plan()["user_job_budget"]["unresolved_submission_reserved"], 0)

    def test_unknown_receipt_absent_queue_reserves_slots_even_for_another_run(self):
        entry = self.plan()["submission_plans"][0]["ledger_intent"]
        entry.update(state="submission_unknown", run_id="unregistered_run")
        self.ledger["submissions"] = [entry]
        self.campaign["max_user_submitted_jobs"] = 3
        p = self.plan()
        self.assertEqual(p["user_job_budget"]["unresolved_submission_reserved"], 2)
        self.assertEqual(p["submission_plans"][0]["array_task_ids"], [0])
        entry.update(state="not_submitted_verified", reconciliation_evidence={"subprocess_invoked": False})
        self.assertEqual(self.plan()["submission_plans"][0]["array_task_ids"], [0, 1])

    def test_invalid_user_submit_cap_fails_closed(self):
        for value in (0, -1, True, 1.5, "1000"):
            with self.subTest(value=value):
                self.campaign["max_user_submitted_jobs"] = value
                with self.assertRaisesRegex(ValueError, "max_user_submitted_jobs"):
                    self.plan()

    def test_terminal_queue_rows_do_not_consume_user_slots(self):
        self.campaign["max_user_submitted_jobs"] = 2
        self.queue["jobs"] = [{"job_id": "99", "state": "COMPLETED"}]
        self.assertEqual(self.plan()["user_job_budget"]["queue_active_members"], 0)

    def test_running_tasks_can_exceed_a_recently_reduced_throttle(self):
        job = self.active_array((0, 1), throttle=1)
        job["running_tasks"] = 2
        self.queue["jobs"] = [job]
        self.assertEqual(self.plan()["cpu_budget"]["project_reserved"], 32)

    def test_unknown_submission_reserves_cpus_and_blocks_repeat(self):
        p = self.plan()
        self.ledger["submissions"] = [p["submission_plans"][0]["ledger_intent"]]
        repeated = self.plan()
        self.assertEqual(repeated["submission_plans"], [])
        self.assertEqual(repeated["cpu_budget"]["unresolved_submission_reserved"], 32)
        self.assertEqual(repeated["runs"][0]["chunks"][0]["state"], "submission_reconciliation_required")

    def test_visible_receipt_does_not_double_reserve_capacity(self):
        p = self.plan()
        receipt = p["submission_plans"][0]["ledger_intent"] | {"state": "submitted", "job_id": "1234"}
        self.ledger["submissions"] = [receipt]
        self.queue["jobs"] = [self.active_array((0, 1), throttle=2)]
        repeated = self.plan()
        self.assertEqual(repeated["cpu_budget"]["unresolved_submission_reserved"], 0)
        self.assertEqual(repeated["cpu_budget"]["project_reserved"], 32)

    def test_wrong_visible_job_identity_cannot_release_an_uncertain_reservation(self):
        p = self.plan()
        self.ledger["submissions"] = [p["submission_plans"][0]["ledger_intent"] |
            {"state": "submitted", "job_id": "1234"}]
        job = self.active_array((0, 1), throttle=2)
        job["project"] = "phase-test"
        self.queue["jobs"] = [job]
        with self.assertRaisesRegex(ValueError, "identity"):
            self.plan()

    def test_verified_terminal_receipt_allows_new_generation(self):
        first = self.plan()["submission_plans"][0]
        self.ledger["submissions"] = [first["ledger_intent"] | {"state": "terminal_verified",
            "reconciliation_evidence": {"sacct_state": "TIMEOUT", "job_id": "1234"}}]
        second = self.plan()["submission_plans"][0]
        self.assertNotEqual(first["plan_id"], second["plan_id"])
        self.assertEqual(second["dispatch_generation"], 2)

    def test_stale_or_incomplete_queue_is_not_permission_to_submit(self):
        self.queue["captured_utc"] = (self.now-timedelta(minutes=4)).isoformat()
        with self.assertRaisesRegex(ValueError, "stale"):
            self.plan()
        self.queue["captured_utc"] = self.now.isoformat()
        self.queue["complete"] = False
        with self.assertRaisesRegex(ValueError, "complete"):
            self.plan()

    def test_short_walltime_and_too_large_array_block_before_dispatch(self):
        self.campaign["runs"][0]["walltime_s"] = 3600
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertIn("Walltime", p["runs"][0]["issues"][0])
        self.campaign["runs"][0]["walltime_s"] = 14400
        self.campaign["max_array_size"] = 1
        self.assertEqual(self.plan()["submission_plans"], [])

    def test_usage_threshold_and_stop_file_block_only_new_work(self):
        self.assertEqual(self.plan(weekly_used=90)["submission_plans"], [])
        stop = self.fixture.base/"STOP"
        stop.touch()
        self.fixture.config["stop_file"] = str(stop)
        self.fixture.save_config()
        self.seal()
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertEqual(p["runs"][0]["chunks"][0]["state"], "stop_file_present")

    def test_json_cli_has_no_shell_side_effect_and_saves_errors(self):
        base = self.fixture.base
        for name, obj in (("campaign", self.campaign), ("queue", self.queue), ("ledger", self.ledger)):
            (base/f"{name}.json").write_text(json.dumps(obj))
        args = ["--campaign", str(base/"campaign.json"), "--queue", str(base/"queue.json"),
                "--ledger", str(base/"ledger.json"), "--output", str(base/"plan.json")]
        # The fixture timestamp is old relative to actual time, so fail closed.
        self.assertEqual(M.main(args), 2)
        self.assertEqual(json.loads((base/"plan.json").read_text())["submission_plans"], [])


if __name__ == "__main__":
    unittest.main()
