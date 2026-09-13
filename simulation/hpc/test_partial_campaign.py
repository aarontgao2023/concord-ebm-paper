"""Toy per-member continuation/capacity tests; no accounting or submission I/O."""
import copy
from datetime import timedelta
import json
import unittest

import accounted_members as A
import dispatch_campaign as D
import manage_campaign as M
import test_manage_campaign as manager_fixtures
import test_dispatch_campaign as dispatch_fixtures


def evidence(entry, now, tasks=(0,), state="COMPLETED", end="2026-09-07T23:59:00"):
    return {"schema_version": A.SCHEMA, "receipt": A.receipt_identity(entry),
            "source": "sacct", "command": A.sacct_command(entry["job_id"]),
            "captured_utc": now.isoformat(), "query_status": "ok",
            "allocations": [{"job_id": f"{entry['job_id']}_{i}", "state": state, "end": end}
                            for i in tasks]}


class PartialCampaignTests(unittest.TestCase):
    def setUp(self):
        self.h = manager_fixtures.CampaignTests()
        self.h.setUp()
        self.addCleanup(self.h.doCleanups)
        self.old = self.h.plan()["submission_plans"][0]["ledger_intent"] | {"state": "submitted", "job_id": "1234"}
        self.h.ledger["submissions"] = [self.old]

    def release(self, **kwargs):
        self.old[A.FIELD] = evidence(self.old, self.h.now, **kwargs)

    def views(self):
        q = M.queue_view(self.h.queue, "ebmcal_v2", self.h.now, 60)
        return M.ledger_view(self.h.ledger, q["jobs_by_id"], self.h.now), M.user_job_budget(self.h.ledger, q, 2)

    def test_partial_parent_releases_only_certified_member_for_same_snapshot_resume(self):
        self.h.queue["jobs"] = [self.h.active_array(task_ids=(1,))]
        self.release()
        before = copy.deepcopy(self.old)
        self.h.campaign["max_user_submitted_jobs"] = 2
        p = self.h.plan()
        s = p["submission_plans"][0]
        self.assertEqual(s["array_task_ids"], [0])
        self.assertEqual(s["dispatch_generation"], 2)
        self.assertEqual(s["snapshot_sha256"], self.old["snapshot_sha256"])
        self.assertEqual(s["nchunks"], self.old["nchunks"])
        self.assertEqual(p["user_job_budget"]["queue_active_members"], 1)
        self.assertEqual(p["user_job_budget"]["unresolved_submission_reserved"], 0)
        self.assertEqual(p["user_job_budget"]["certified_terminal_members_released"], 1)
        self.assertEqual(self.old, before)
        (pending, unseen_cpu), _ = self.views()
        self.assertNotIn(("fixture", 0), pending)
        self.assertIn(("fixture", 1), pending)
        self.assertEqual(unseen_cpu, 0)

    def test_absent_parent_reserves_only_uncertified_remaining_members_cpu_and_slots(self):
        self.release()
        (pending, cpu), slots = self.views()
        self.assertEqual(set(pending), {("fixture", 1)})
        self.assertEqual(cpu, 16)
        self.assertEqual(slots["unresolved_submission_reserved"], 1)
        self.assertEqual(self.h.plan()["submission_plans"][0]["array_task_ids"], [0])

    def test_all_members_certified_can_release_zero_cpu_without_rewriting_parent(self):
        self.release(tasks=(0, 1))
        (pending, cpu), slots = self.views()
        self.assertEqual(dict(pending), {})
        self.assertEqual(cpu, 0)
        self.assertEqual(slots["used_before_plan"], 0)
        self.assertEqual(self.h.plan()["submission_plans"][0]["array_task_ids"], [0, 1])
        self.assertEqual(self.old["state"], "submitted")

    def test_queue_active_member_wins_over_earlier_terminal_sacct(self):
        self.release()
        self.h.queue["jobs"] = [self.h.active_array(task_ids=(0,))]
        (pending, cpu), slots = self.views()
        self.assertEqual(set(pending), {("fixture", 0), ("fixture", 1)})
        self.assertEqual(slots["certified_terminal_members_released"], 0)
        self.assertEqual(slots["used_before_plan"], 2)
        self.assertEqual(self.h.plan()["submission_plans"], [])

    def test_unknown_and_intent_keep_all_targets_cpus_and_slots_even_with_evidence(self):
        self.release(tasks=(0, 1))
        for state in ("intent", "submission_unknown"):
            with self.subTest(state=state):
                self.old["state"] = state
                (pending, cpu), slots = self.views()
                self.assertEqual(len(pending), 2)
                self.assertEqual(cpu, 32)
                self.assertEqual(slots["used_before_plan"], 2)
                self.assertEqual(slots["certified_terminal_members_released"], 0)
                self.assertEqual(self.h.plan()["submission_plans"], [])

    def test_wrong_receipt_duplicate_and_future_evidence_fail_closed(self):
        self.release()
        valid = copy.deepcopy(self.old[A.FIELD])
        malformed = []
        bad = copy.deepcopy(valid); bad["receipt"]["snapshot_sha256"] = "f"*64; malformed.append(bad)
        bad = copy.deepcopy(valid); bad["allocations"].append(bad["allocations"][0]); malformed.append(bad)
        bad = copy.deepcopy(valid); bad["captured_utc"] = (self.h.now+timedelta(seconds=1)).isoformat(); malformed.append(bad)
        for value in malformed:
            with self.subTest(value=value):
                self.old[A.FIELD] = value
                with self.assertRaises(ValueError):
                    self.h.plan()

    def test_failed_query_missing_end_and_nonterminal_cannot_release_member(self):
        self.release(end="Unknown")
        self.assertEqual(len(self.views()[0][0]), 2)
        self.release(state="RUNNING")
        self.assertEqual(len(self.views()[0][0]), 2)
        self.old[A.FIELD].update(query_status="failed", allocations=[], error="temporary accounting failure")
        self.assertEqual(len(self.views()[0][0]), 2)

    def test_terminal_failed_observed_fit_is_retained_and_not_scheduled_again(self):
        self.release()
        self.h.queue["jobs"] = [self.h.active_array(task_ids=(1,))]
        rows = [self.h.fixture.make_row(seed, engine) for seed in (100, 102)
                for engine in self.h.fixture.config["engines"]]
        rows[0].update(status="error", observed={})
        self.h.fixture.write_chunk(0, 2, rows=rows)
        p = self.h.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertEqual(p["runs"][0]["chunks"][0]["state"], "terminal_with_failures")
        self.assertEqual(p["runs"][0]["chunks"][0]["failed_observed"], 1)

    def test_new_receipt_for_released_member_restores_reservation_and_blocks_duplicate(self):
        self.release()
        self.h.queue["jobs"] = [self.h.active_array(task_ids=(1,))]
        new = self.h.plan()["submission_plans"][0]["ledger_intent"]
        self.h.ledger["submissions"].append(new)
        p = self.h.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertEqual(p["user_job_budget"]["used_before_plan"], 2)
        self.assertEqual(p["cpu_budget"]["unresolved_submission_reserved"], 16)

    def test_terminal_failed_permutation_remains_in_attempted_denominator(self):
        definition = {"stratify": "diagnosis", "require_max": True}
        self.h.fixture.config.update(bperm=2, alpha=.05, rule="le", full_p_first_n=4,
                                     schemes={"diagnosis": definition})
        self.h.fixture.save_config()
        self.h.seal()
        self.h.ledger["submissions"] = []
        self.old = self.h.plan()["submission_plans"][0]["ledger_intent"] | {"state": "submitted", "job_id": "1234"}
        self.h.ledger["submissions"] = [self.old]
        self.release()
        self.h.queue["jobs"] = [self.h.active_array(task_ids=(1,))]
        rows = [self.h.fixture.make_row(seed, engine) for seed in (100, 102)
                for engine in self.h.fixture.config["engines"]]
        for row in rows:
            row["schemes"] = {"diagnosis": {"definition": definition, "requested_nperm": 2,
                "nperm": 2, "attempted": 2, "complete": True, "failed": 0,
                "gt": [2]*3, "eq": [0]*3, "maxgt": 2, "maxeq": 0}}
        rows[0]["schemes"]["diagnosis"].update(nperm=1, failed=1, complete=False, gt=[1]*3, maxgt=1)
        self.h.fixture.write_chunk(0, 2, rows=rows)
        p = self.h.plan()
        self.assertEqual(p["submission_plans"], [])
        chunk = p["runs"][0]["chunks"][0]
        self.assertEqual(chunk["state"], "terminal_with_failures")
        self.assertEqual(chunk["failed_permutations"], 1)
        saved = self.h.fixture.root/"cell_0/chunk_0/rows.jsonl"
        self.assertEqual(json.loads(saved.read_text().splitlines()[0])["schemes"]["diagnosis"]["attempted"], 2)


class PartialDispatchTests(unittest.TestCase):
    def setUp(self):
        self.h = dispatch_fixtures.DispatchTests()
        self.h.setUp()
        self.addCleanup(self.h.doCleanups)
        self.old = self.h.plan["submission_plans"][0]["ledger_intent"] | {"state": "submitted", "job_id": "1234"}
        self.old[A.FIELD] = evidence(self.old, self.h.now)
        self.h.ledger["submissions"].append(self.old)
        self.h.queue["jobs"].append({"job_id": "1234", "state": "RUNNING", "project": "ebmcal_v2",
            "run_id": "fixture", "snapshot_sha256": self.old["snapshot_sha256"], "array_task_ids": [1],
            "cpus_per_task": 16, "max_concurrent": 1, "running_tasks": 1})
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 32, "max_user_submitted_jobs": 2,
            "runs": [{"run_id": "fixture", "snapshot": str(self.h.fixture.frozen),
                      "snapshot_sha256": self.old["snapshot_sha256"], "nchunks": 2,
                      "output_root": str(self.h.fixture.root)}]}
        self.h.plan = M.make_plan(self.campaign, self.h.queue, self.h.ledger, self.h.now, weekly_used=12)
        self.assertEqual(self.h.plan["submission_plans"][0]["array_task_ids"], [0])
        self.h.save_inputs()

    def test_locked_dispatch_uses_same_release_set_and_preserves_prior_receipt(self):
        before = copy.deepcopy(self.old)
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")
        self.assertEqual(self.h.dispatch(apply=True)["status"], "submitted")
        saved = json.loads(self.h.ledger_path.read_text())["submissions"]
        self.assertEqual(saved[0], before)
        self.assertEqual(saved[1]["array_task_ids"], [0])
        self.assertEqual(saved[1]["dispatch_generation"], 2)
        with self.assertRaises(ValueError):
            self.h.dispatch(apply=True)
        self.assertEqual(len(self.h.calls), 1)

    def test_accounting_replaced_with_failure_before_dispatch_blocks_retry(self):
        self.old[A.FIELD].update(query_status="failed", allocations=[], error="sacct unavailable")
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "duplicated, active or awaiting"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])

    def test_member_becomes_active_before_dispatch_even_in_other_parent_blocks_retry(self):
        new = copy.deepcopy(self.h.queue["jobs"][0])
        new.update(job_id="9999", array_task_ids=[0])
        self.h.queue["jobs"].append(new)
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "duplicated, active or awaiting"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])


if __name__ == "__main__":
    unittest.main()
