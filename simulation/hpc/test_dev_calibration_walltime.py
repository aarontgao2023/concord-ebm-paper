"""Exact old development snapshot policy, tested with static source and toy plans."""
import copy
import hashlib
import json
from pathlib import Path
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M
import test_dispatch_campaign as fixtures

PROJECT = Path(__file__).resolve().parents[2]
RUN = "dev_calibration_a"


class DevCalibrationWalltimeTests(unittest.TestCase):
    def setUp(self):
        self.h = fixtures.DispatchTests()
        self.h.setUp()
        self.addCleanup(self.h.doCleanups)
        f = self.h.fixture
        source = PROJECT/"runs/v2/snapshots"/RUN
        f.config_path.unlink()
        f.config_path = f.frozen/"configs/v2"/(RUN+".json")
        f.config = json.loads((source/"configs/v2"/(RUN+".json")).read_text())
        f.save_config()
        f.root = self.h.root/"runs"/RUN
        self.script = f.frozen/"hpc/v2/run_array.sbatch"
        self.script.write_text((source/"hpc/v2/run_array.sbatch").read_text().replace(
            "root=/N/slate/tg11/ebmcal_v2", f"root={self.h.root}"))
        self.runner = f.frozen/"scripts/v2/run_v2.py"
        self.runner.write_bytes((source/"scripts/v2/run_v2.py").read_bytes())
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 32,
            "runs": [{"run_id": RUN, "snapshot": str(f.frozen), "nchunks": 8,
                      "output_root": str(f.root), "walltime_s": 3600,
                      "execution_walltime_policy": M.DEV_CALIBRATION_1H_POLICY}]}
        self.seal()
        patch = mock.patch.object(M, "DEV_CALIBRATION_1H_PIN", self.pin())
        patch.start(); self.addCleanup(patch.stop)
        patch = mock.patch.dict(M.DEV_CALIBRATION_1H_SOURCE_PINS, self.source_pins(), clear=True)
        patch.start(); self.addCleanup(patch.stop)

    def seal(self):
        root = self.h.fixture.frozen
        paths = sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        body = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root)}\n" for p in paths)
        (root/"SHA256SUMS").write_text(body)
        self.campaign["runs"][0]["snapshot_sha256"] = hashlib.sha256(body.encode()).hexdigest()

    def pin(self):
        return (self.campaign["runs"][0]["snapshot_sha256"],
                hashlib.sha256(self.h.fixture.config_path.read_bytes()).hexdigest())

    def source_pins(self):
        return {"hpc/v2/run_array.sbatch": hashlib.sha256(self.script.read_bytes()).hexdigest(),
                "scripts/v2/run_v2.py": hashlib.sha256(self.runner.read_bytes()).hexdigest()}

    def plan(self):
        return M.make_plan(self.campaign, self.h.queue, self.h.ledger, self.h.now, weekly_used=12)

    def install(self):
        self.h.plan = self.plan()
        self.assertEqual(len(self.h.plan["submission_plans"]), 1, self.h.plan)
        self.h.save_inputs()
        return self.h.plan["submission_plans"][0]

    def test_policy_retains_cell_major24_tasks_with_one_engine_per_job(self):
        s = self.install()
        self.assertEqual(s["array_task_ids"], list(range(24)))
        self.assertEqual((s["mapping"], s["nchunks"]), ("cell_major", 8))
        self.assertEqual(s["argv"][4], "--time=01:00:00")
        evidence = s["execution_evidence"]
        self.assertEqual(evidence["configured_engines"], ["original", "repaired"])
        self.assertEqual(evidence["engines_per_fit_job"], 1)
        self.assertEqual((evidence["batch_max_jobs"], evidence["nominal_batch_drain_s"],
                          evidence["nominal_checkpoint_reserve_s"]), (32, 600, 900))
        chunk = next(v for v in self.h.plan["runs"][0]["chunks"] if v["array_index"] == 13)
        self.assertEqual((chunk["cell_index"], chunk["chunk"], chunk["seeds"]), (1, 5, [31000205]))
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")
        self.assertEqual(self.h.calls, [])
        self.assertEqual(self.h.dispatch(apply=True)["status"], "submitted")
        saved = json.loads(self.h.ledger_path.read_text())["submissions"][0]
        self.assertEqual(saved["execution_evidence"], evidence)

    def test_four_hour_default_and_existing_receipt_are_unchanged(self):
        run = self.campaign["runs"][0]
        run.pop("execution_walltime_policy")
        self.assertEqual(self.plan()["submission_plans"], [])
        run["walltime_s"] = 14400
        old = self.install()["ledger_intent"] | {"state": "submitted", "job_id": "7777"}
        before = copy.deepcopy(old)
        self.h.ledger["submissions"].append(old)
        run.update(walltime_s=3600, execution_walltime_policy=M.DEV_CALIBRATION_1H_POLICY)
        self.assertEqual(self.plan()["submission_plans"], [])
        self.assertEqual(old, before)
        self.assertNotIn("execution_evidence", old)

    def test_confirmation_policy_cannot_be_used_for_the_development_snapshot(self):
        self.campaign["runs"][0]["execution_walltime_policy"] = M.SINGLE_ENGINE_1H_POLICY
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertNotIn(RUN, M.SINGLE_ENGINE_1H_PINS)

    def test_changed_scientific_config_or_new_manifest_is_not_allowlisted(self):
        self.h.fixture.config["bperm"] = 600
        self.h.fixture.save_config()
        self.seal()
        self.assertEqual(self.plan()["submission_plans"], [])
        # Even a separately pinned test config cannot bypass the fixed B599 invariant.
        with mock.patch.object(M, "DEV_CALIBRATION_1H_PIN", self.pin()):
            self.assertEqual(self.plan()["submission_plans"], [])

    def test_wrong_workers_walltime_layout_and_paired_fit_are_rejected(self):
        for key, value in (("cpus_per_task", 8), ("walltime_s", 7200), ("nchunks", 4)):
            with self.subTest(key=key):
                old = copy.deepcopy(self.campaign["runs"][0])
                self.campaign["runs"][0][key] = value
                self.assertEqual(self.plan()["submission_plans"], [])
                self.campaign["runs"][0] = old
        self.h.fixture.config["paired_standard"] = True
        self.h.fixture.save_config(); self.seal()
        with mock.patch.object(M, "DEV_CALIBRATION_1H_PIN", self.pin()):
            self.assertEqual(self.plan()["submission_plans"], [])

    def test_stop_signal_exec_and_job_bound_markers_are_required_even_with_toy_pins(self):
        for path, before, after in ((self.script, 'exec "$HOME/', '"$HOME/'),
                                   (self.script, "B:USR1@1500", "B:USR1@900"),
                                   (self.runner, "while not stopping():", "while True:"),
                                   (self.runner, "jobs[:max(args.workers*2, 8)]", "jobs[:999]"),
                                   (self.runner, "imap_unordered(fit_job, jobs", "imap_unordered(fit_paired_job, jobs")):
            with self.subTest(marker=before):
                old = path.read_text()
                path.write_text(old.replace(before, after)); self.seal()
                with mock.patch.object(M, "DEV_CALIBRATION_1H_PIN", self.pin()), \
                        mock.patch.dict(M.DEV_CALIBRATION_1H_SOURCE_PINS, self.source_pins(), clear=True):
                    self.assertEqual(self.plan()["submission_plans"], [])
                path.write_text(old); self.seal()

    def test_dispatch_recomputes_development_timing_evidence_after_plan_rehash(self):
        s = self.install()
        s["execution_evidence"]["nominal_batch_drain_s"] = 10
        identity = D.submission_identity(s)
        s["plan_id"] = M.I.digest(identity)
        s["argv"][5] = "--comment=ebmcal_v2:"+s["plan_id"]
        s["ledger_intent"] = identity | {"plan_id": s["plan_id"], "state": "intent", "created_utc": self.h.now.isoformat()}
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "Execution evidence differs"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])


class ProductionDevPinTests(unittest.TestCase):
    def test_original_development_bytes_accept_only_the_separate_policy(self):
        f = M.snapshot_info(PROJECT/"runs/v2/snapshots"/RUN, RUN, M.DEV_CALIBRATION_1H_PIN[0])
        cfg = M.read_json(f["config_path"])
        e = M.validate_execution_walltime(f, cfg, 16, 3600, M.DEV_CALIBRATION_1H_POLICY, nchunks=8)
        self.assertEqual(e["nominal_checkpoint_reserve_s"], 900)
        self.assertIsNone(M.validate_execution_walltime(f, cfg, 16, 14400))
        with self.assertRaises(ValueError):
            M.validate_execution_walltime(f, cfg, 16, 3600, M.SINGLE_ENGINE_1H_POLICY, nchunks=8)
        for other in ("confirm_core_a", "confirm_mechanism_b", "confirm_stage_a", "dev_mechanism_a"):
            with self.subTest(run=other), self.assertRaisesRegex(ValueError, "exact dev_calibration_a"):
                M.validate_execution_walltime(f, cfg | {"run_id": other}, 16, 3600,
                                              M.DEV_CALIBRATION_1H_POLICY, nchunks=8)


if __name__ == "__main__":
    unittest.main()
