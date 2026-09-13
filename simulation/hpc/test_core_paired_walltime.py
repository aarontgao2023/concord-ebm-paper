"""Exact core policy: inert source copies, toy/mock dispatch, no HPC or fitting."""
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
RUN = "confirm_core_a"


class CorePairedWalltimeTests(unittest.TestCase):
    def setUp(self):
        self.h = fixtures.DispatchTests()
        self.h.setUp()
        self.addCleanup(self.h.doCleanups)
        f = self.h.fixture
        source = PROJECT / "runs/v2/snapshots" / RUN
        f.config_path.unlink()
        f.config_path = f.frozen / "configs/v2" / (RUN + ".json")
        f.config = json.loads((source / "configs/v2" / (RUN + ".json")).read_text())
        f.config["stop_file"] = str(self.h.root / "STOP")
        f.save_config()
        f.root = self.h.root / "runs" / RUN
        self.script = f.frozen / "hpc/v2/run_array.sbatch"
        self.script.write_text((source / "hpc/v2/run_array.sbatch").read_text().replace(
            "root=/N/slate/tg11/ebmcal_v2", f"root={self.h.root}"))
        self.runner = f.frozen / "scripts/v2/run_v2.py"
        self.runner.write_bytes((source / "scripts/v2/run_v2.py").read_bytes())
        self.paired = f.frozen / "scripts/v2/paired_engine_v2.py"
        self.paired.write_bytes((source / "scripts/v2/paired_engine_v2.py").read_bytes())
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 32, "max_tasks_per_submission": 2,
            "max_user_submitted_jobs": 10, "runs": [{"run_id": RUN, "snapshot": str(f.frozen), "nchunks": 500,
                "output_root": str(f.root), "walltime_s": 3600, "max_parallel_tasks": 1,
                "execution_walltime_policy": M.CORE_PAIRED_1H_POLICY}]}
        self.seal()
        p = mock.patch.object(M, "CORE_PAIRED_1H_PIN", self.pin())
        p.start(); self.addCleanup(p.stop)
        p = mock.patch.dict(M.CORE_PAIRED_1H_SOURCE_PINS, self.source_pins(), clear=True)
        p.start(); self.addCleanup(p.stop)

    def seal(self):
        root = self.h.fixture.frozen
        paths = sorted(p for p in root.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        body = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(root)}\n" for p in paths)
        (root / "SHA256SUMS").write_text(body)
        self.campaign["runs"][0]["snapshot_sha256"] = hashlib.sha256(body.encode()).hexdigest()

    def pin(self):
        return (self.campaign["runs"][0]["snapshot_sha256"], hashlib.sha256(self.h.fixture.config_path.read_bytes()).hexdigest())

    def source_pins(self):
        return {str(path.relative_to(self.h.fixture.frozen)): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (self.script, self.runner, self.paired)}

    def plan(self):
        return M.make_plan(self.campaign, self.h.queue, self.h.ledger, self.h.now, weekly_used=12)

    def install(self):
        self.h.plan = self.plan()
        self.assertEqual(len(self.h.plan["submission_plans"]), 1, self.h.plan)
        self.h.save_inputs()
        return self.h.plan["submission_plans"][0]

    def reidentify(self, s):
        identity = D.submission_identity(s)
        s["plan_id"] = M.I.digest(identity)
        s["argv"][5] = "--comment=ebmcal_v2:" + s["plan_id"]
        s["ledger_intent"] = identity | {"plan_id": s["plan_id"], "state": "intent", "created_utc": self.h.now.isoformat()}
        self.h.save_inputs()

    def test_exact_policy_requests_13argv_binds_evidence_and_fake_receipt(self):
        s = self.install()
        self.assertEqual(len(s["argv"]), 13)
        self.assertEqual((s["argv"][4], s["argv"][8]), ("--time=01:00:00", "--signal=B:USR1@2100"))
        self.assertEqual(s["array_task_ids"], [0, 1])
        self.assertEqual((s["nchunks"], s["mapping"], s["max_concurrent"]), (500, "chunk_major", 1))
        e = s["execution_evidence"]
        self.assertEqual((e["paired_total_timeout_s"], e["batch_max_jobs"], e["nominal_batch_waves"],
                          e["nominal_batch_drain_s"], e["nominal_checkpoint_reserve_s"]), (600, 32, 2, 1200, 900))
        self.assertEqual(e["earliest_requested_signal_elapsed_s"], 1440)
        self.assertIs(e["signal_effective_verified"], False)
        self.assertEqual(e["signal_verification_status"], "unverified_until_first_execution")
        self.assertEqual(s["plan_id"], M.I.digest(D.submission_identity(s)))
        before = self.h.ledger_path.read_bytes()
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")
        self.assertEqual(self.h.ledger_path.read_bytes(), before)
        self.assertEqual(self.h.calls, [])
        self.assertEqual(self.h.dispatch(apply=True)["status"], "submitted")
        saved = json.loads(self.h.ledger_path.read_text())["submissions"][0]
        self.assertEqual(saved["execution_evidence"], e)
        self.assertEqual(saved["argv"], s["argv"])

    def test_old_two_hour_identity_argv_and_receipt_are_unchanged(self):
        run = self.campaign["runs"][0]
        run.pop("execution_walltime_policy")
        self.assertEqual(self.plan()["submission_plans"], [])
        run["walltime_s"] = 7200
        s = self.install()
        self.assertEqual(len(s["argv"]), 12)
        self.assertNotIn("execution_evidence", s)
        self.assertEqual(s["plan_id"], M.I.digest({k: s[k] for k in D.IDENTITY_FIELDS}))
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")
        old = copy.deepcopy(s["ledger_intent"] | {"state": "submitted", "job_id": "8888"})
        self.h.ledger["submissions"].append(old)
        run.update(walltime_s=3600, execution_walltime_policy=M.CORE_PAIRED_1H_POLICY)
        plan = self.plan()
        submitted = set(old["array_task_ids"])
        self.assertTrue(all(not submitted.intersection(p["array_task_ids"]) for p in plan["submission_plans"]))
        self.assertEqual(self.h.ledger["submissions"][0], old)

    def test_signal_spelling_order_duplicate_and_extra_args_rejected(self):
        original = copy.deepcopy(self.install())
        bad_argv = []
        for option in ("--signal=B:USR1@1500", "--signal=USR1@2100", "--signal=B:USR1@2099", "--signal=B:TERM@2100"):
            argv = original["argv"][:]
            argv[8] = option
            bad_argv.append(argv)
        bad_argv += [original["argv"][:8] + original["argv"][9:],
                     original["argv"][:8] + [original["argv"][8]] + original["argv"][8:],
                     original["argv"][:8] + original["argv"][9:] + [original["argv"][8]],
                     original["argv"] + ["--wrap=echo unsafe"]]
        for argv in bad_argv:
            with self.subTest(argv=argv):
                self.h.plan["submission_plans"][0] = copy.deepcopy(original)
                self.h.plan["submission_plans"][0]["argv"] = argv
                self.h.save_inputs()
                with self.assertRaises(ValueError):
                    self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])

    def test_ordinary_policy_cannot_add_signal_and_core_cannot_omit_policy(self):
        s = self.install()
        s.pop("execution_walltime_policy")
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "requires both"):
            self.h.dispatch()
        run = self.campaign["runs"][0]
        run.pop("execution_walltime_policy")
        run["walltime_s"] = 7200
        s = self.install()
        s["argv"].insert(8, M.CORE_PAIRED_1H_SIGNAL_OPTION)
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "12-argument"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])

    def test_rehashed_false_signal_verification_or_timing_evidence_rejected(self):
        original = copy.deepcopy(self.install())
        for key, value in (("signal_effective_verified", True), ("requested_signal", "B:USR1@1500"),
                           ("nominal_batch_drain_s", 600), ("nominal_checkpoint_reserve_s", 1500)):
            with self.subTest(key=key):
                s = copy.deepcopy(original)
                self.h.plan["submission_plans"][0] = s
                s["execution_evidence"][key] = value
                self.reidentify(s)
                with self.assertRaisesRegex(ValueError, "Execution evidence differs"):
                    self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])

    def test_wrong_wall_workers_layout_and_policy_rejected(self):
        run = self.campaign["runs"][0]
        original = copy.deepcopy(run)
        for key, value in (("walltime_s", 7200), ("walltime_s", 3599), ("cpus_per_task", 8),
                           ("nchunks", 250), ("execution_walltime_policy", M.SINGLE_ENGINE_1H_POLICY),
                           ("execution_walltime_policy", M.DEV_CALIBRATION_1H_POLICY)):
            with self.subTest(key=key, value=value):
                self.campaign["runs"][0] = original | {key: value}
                self.assertEqual(self.plan()["submission_plans"], [])

    def test_changed_scientific_config_rejected_even_with_toy_manifest_pin(self):
        f = self.h.fixture
        original = copy.deepcopy(f.config)
        for key, value in (("fit_timeout_s", 600), ("paired_standard", False), ("datasets", 500),
                           ("engines", ["repaired"]), ("bperm", 600), ("full_p_first_n", 99)):
            with self.subTest(key=key):
                f.config = original | {key: value}
                f.save_config(); self.seal()
                with mock.patch.object(M, "CORE_PAIRED_1H_PIN", self.pin()):
                    self.assertEqual(self.plan()["submission_plans"], [])
        f.config = original

    def test_source_pins_independent_of_snapshot_pin(self):
        self.runner.write_text(self.runner.read_text() + "\n# changed inert source\n")
        self.seal()
        with mock.patch.object(M, "CORE_PAIRED_1H_PIN", self.pin()):
            p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertIn("source hash differs", p["runs"][0]["issues"][0])

    def test_exec_handler_paired_timeout_and_saved_half_markers_are_required(self):
        for path, before, after in (
                (self.script, 'exec "$HOME/', '"$HOME/'),
                (self.script, "B:USR1@1500", "B:USR1@2100"),
                (self.runner, "total_budget = 2*config.get", "total_budget = 3*config.get"),
                (self.runner, "signal.signal(sig, request_stop)", "pass"),
                (self.runner, "[:max(args.workers*2, 8)]", "[:999]"),
                (self.runner, "completed = getattr(exc, 'completed_results', {})", "completed = {}"),
                (self.runner, "replay_audits[candidate_key] = replay_consistent(records[candidate_key], candidate)", "pass"),
                (self.paired, "abort.cause.completed_results = results", "pass")):
            with self.subTest(marker=before):
                old = path.read_text()
                path.write_text(old.replace(before, after)); self.seal()
                with mock.patch.object(M, "CORE_PAIRED_1H_PIN", self.pin()), \
                        mock.patch.dict(M.CORE_PAIRED_1H_SOURCE_PINS, self.source_pins(), clear=True):
                    self.assertEqual(self.plan()["submission_plans"], [])
                path.write_text(old); self.seal()

    def test_stop_guard_remains_in_force(self):
        self.install()
        (self.h.root / "STOP").write_text("toy stop")
        with self.assertRaisesRegex(ValueError, "STOP"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])


class ProductionCorePairedPinTests(unittest.TestCase):
    def test_actual_frozen_hashes_and_readonly_marker_evidence(self):
        frozen = M.snapshot_info(PROJECT / "runs/v2/snapshots" / RUN, RUN, M.CORE_PAIRED_1H_PIN[0])
        cfg = M.read_json(frozen["config_path"])
        e = M.validate_execution_walltime(frozen, cfg, 16, 3600, M.CORE_PAIRED_1H_POLICY, nchunks=500)
        self.assertEqual(e["frozen_default_signal"], "B:USR1@1500")
        self.assertEqual(e["requested_signal"], "B:USR1@2100")
        self.assertEqual(e["paired_engine_evidence"]["completed_results_on_abort"]["line"], 242)
        self.assertIsNone(M.validate_execution_walltime(frozen, cfg, 16, 7200))
        with self.assertRaises(ValueError):
            M.validate_execution_walltime(frozen, cfg, 16, 3600)
        for name in ("confirm_stage_a", "confirm_mechanism_b", "dev_calibration_a"):
            with self.subTest(run=name), self.assertRaisesRegex(ValueError, "exact confirm_core_a"):
                M.validate_execution_walltime(frozen, cfg | {"run_id": name}, 16, 3600, M.CORE_PAIRED_1H_POLICY, nchunks=500)


if __name__ == "__main__":
    unittest.main()
