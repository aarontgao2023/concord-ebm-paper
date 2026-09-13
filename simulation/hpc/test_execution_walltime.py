"""Static-source and toy-plan tests; never import a simulator or call sbatch.

Integration uses temporary, explicitly monkeypatched fixture pins. Production
pins are checked read-only, without opening simulation rows or executing seeds.
"""
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


class ExecutionWalltimeTests(unittest.TestCase):
    def setUp(self):
        self.h = fixtures.DispatchTests()
        self.h.setUp()
        self.addCleanup(self.h.doCleanups)
        f = self.h.fixture
        # Reading verified source text is not importing/executing that source.
        source = PROJECT/"runs/v2/snapshots/confirm_stage_a"
        self.script = f.frozen/"hpc/v2/run_array.sbatch"
        self.script.write_text((source/"hpc/v2/run_array.sbatch").read_text().replace(
            "root=/N/slate/tg11/ebmcal_v2", f"root={self.h.root}"))
        self.runner = f.frozen/"scripts/v2/run_v2.py"
        self.runner.write_bytes((source/"scripts/v2/run_v2.py").read_bytes())
        f.config.update(phase="confirmation", base_seed=31000000, engines=["repaired"],
                        paired_standard=False, fit_timeout_s=300)
        f.save_config()
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 256,
            "runs": [{"run_id": "fixture", "snapshot": str(f.frozen), "nchunks": 2,
                      "output_root": str(f.root), "walltime_s": 3600,
                      "execution_walltime_policy": M.SINGLE_ENGINE_1H_POLICY}]}
        self.seal()
        self.pin_patch = mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": self.pin()}, clear=True)
        self.source_patch = mock.patch.dict(M.SINGLE_ENGINE_1H_SOURCE_PINS, self.source_pins(), clear=True)
        self.pin_patch.start()
        self.source_patch.start()
        self.addCleanup(self.pin_patch.stop)
        self.addCleanup(self.source_patch.stop)

    def seal(self):
        frozen = self.h.fixture.frozen
        files = sorted(p for p in frozen.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        sha = "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(frozen)}\n" for p in files)
        (frozen/"SHA256SUMS").write_text(sha)
        self.campaign["runs"][0]["snapshot_sha256"] = hashlib.sha256(sha.encode()).hexdigest()

    def pin(self):
        return (self.campaign["runs"][0]["snapshot_sha256"],
                hashlib.sha256(self.h.fixture.config_path.read_bytes()).hexdigest(), 2)

    def source_pins(self):
        return {"hpc/v2/run_array.sbatch": hashlib.sha256(self.script.read_bytes()).hexdigest(),
                "scripts/v2/run_v2.py": hashlib.sha256(self.runner.read_bytes()).hexdigest()}

    def plan(self):
        return M.make_plan(self.campaign, self.h.queue, self.h.ledger, self.h.now, weekly_used=12)

    def install_plan(self):
        self.h.plan = self.plan()
        self.assertEqual(len(self.h.plan["submission_plans"]), 1, self.h.plan)
        self.h.save_inputs()
        return self.h.plan["submission_plans"][0]

    def reidentify(self, submission):
        identity = D.submission_identity(submission)
        submission["plan_id"] = M.I.digest(identity)
        submission["argv"][5] = f"--comment=ebmcal_v2:{submission['plan_id']}"
        submission["ledger_intent"] = identity | {"plan_id": submission["plan_id"], "state": "intent",
                                                 "created_utc": self.h.now.isoformat()}
        self.h.save_inputs()

    def test_explicit_policy_binds_exact_1h_evidence_and_fake_submit_intent(self):
        s = self.install_plan()
        evidence = s["execution_evidence"]
        self.assertEqual(evidence["nominal_batch_drain_s"], 600)
        self.assertEqual(evidence["nominal_checkpoint_reserve_s"], 900)
        self.assertEqual(evidence["batch_max_jobs"], 32)
        self.assertEqual(evidence["earliest_signal_elapsed_s"], 2040)
        self.assertEqual(s["argv"][4], "--time=01:00:00")
        self.assertEqual(s["plan_id"], M.I.digest(D.submission_identity(s)))
        before = self.h.ledger_path.read_bytes()
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")
        self.assertEqual(self.h.ledger_path.read_bytes(), before)
        self.assertEqual(self.h.calls, [])
        self.assertEqual(self.h.dispatch(apply=True)["status"], "submitted")
        receipt = json.loads(self.h.ledger_path.read_text())["submissions"][0]
        self.assertEqual(receipt["execution_evidence"], evidence)
        self.assertEqual(receipt["walltime_s"], 3600)

    def test_two_hour_default_keeps_original_identity_and_one_hour_needs_opt_in(self):
        self.campaign["runs"][0].pop("execution_walltime_policy")
        self.assertEqual(self.plan()["submission_plans"], [])
        self.campaign["runs"][0]["walltime_s"] = 7200
        s = self.install_plan()
        self.assertNotIn("execution_evidence", s)
        self.assertNotIn("execution_walltime_policy", s)
        self.assertEqual(s["plan_id"], M.I.digest({k: s[k] for k in D.IDENTITY_FIELDS}))
        self.assertEqual(self.h.dispatch()["status"], "validated_dry_run")

    def test_unallowlisted_run_and_new_snapshot_are_rejected(self):
        with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {}, clear=True):
            self.assertEqual(self.plan()["submission_plans"], [])
        self.runner.write_text(self.runner.read_text()+"\n# different frozen version\n")
        self.seal()
        p = self.plan()
        self.assertEqual(p["submission_plans"], [])
        self.assertIn("exact frozen snapshot", p["runs"][0]["issues"][0])

    def test_wrong_timeout_engine_workers_wall_or_chunk_layout_rejected(self):
        for field, value in (("fit_timeout_s", 600), ("paired_standard", True),
                             ("engines", ["original", "repaired"])):
            with self.subTest(field=field):
                original = copy.deepcopy(self.h.fixture.config)
                self.h.fixture.config[field] = value
                self.h.fixture.save_config()
                self.seal()
                # Even an independently allowlisted config must meet the invariant.
                with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": self.pin()}, clear=True):
                    self.assertEqual(self.plan()["submission_plans"], [])
                self.h.fixture.config = original
                self.h.fixture.save_config()
                self.seal()
        for field, value in (("cpus_per_task", 8), ("walltime_s", 3599), ("walltime_s", 7200), ("nchunks", 1)):
            with self.subTest(field=field):
                run = self.campaign["runs"][0]
                old = copy.deepcopy(run)
                run[field] = value
                self.assertEqual(self.plan()["submission_plans"], [])
                self.campaign["runs"][0] = old

    def test_resealed_lookalike_source_does_not_bypass_exact_source_pin(self):
        self.runner.write_text(self.runner.read_text().replace("while not stopping():", "while True:"))
        self.seal()
        with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": self.pin()}, clear=True):
            p = self.plan()
        self.assertIn("source hash differs", p["runs"][0]["issues"][0])

    def test_config_pin_is_required_independently_of_manifest_pin(self):
        self.h.fixture.config["datasets"] = 6
        self.h.fixture.save_config()
        self.seal()
        old_config_hash = M.SINGLE_ENGINE_1H_PINS["fixture"][1]
        bad_pin = (self.pin()[0], old_config_hash, 2)
        with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": bad_pin}, clear=True):
            self.assertIn("configuration hash", self.plan()["runs"][0]["issues"][0])

    def test_missing_final_exec_and_changed_signal_are_rejected_even_with_toy_pin(self):
        original = self.script.read_text()
        for changed in (original.replace('exec "$HOME/', '"$HOME/'),
                        original+"\necho done\n", original.replace("B:USR1@1500", "B:USR1@900")):
            with self.subTest(changed=changed[-40:]):
                self.script.write_text(changed)
                self.seal()
                with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": self.pin()}, clear=True), \
                        mock.patch.dict(M.SINGLE_ENGINE_1H_SOURCE_PINS, self.source_pins(), clear=True):
                    self.assertEqual(self.plan()["submission_plans"], [])
        self.script.write_text(original)

    def test_missing_handler_or_loop_or_fsync_rejected_even_with_toy_pin(self):
        original = self.runner.read_text()
        for before, after in (("signal.signal(sig, request_stop)", "pass"),
                              ("while not stopping():", "while True:"),
                              ("[:max(args.workers*2, 8)]", "[:999]"),
                              ("os.fsync(stream.fileno())", "pass")):
            with self.subTest(marker=before):
                self.runner.write_text(original.replace(before, after))
                self.seal()
                with mock.patch.dict(M.SINGLE_ENGINE_1H_PINS, {"fixture": self.pin()}, clear=True), \
                        mock.patch.dict(M.SINGLE_ENGINE_1H_SOURCE_PINS, self.source_pins(), clear=True):
                    p = self.plan()
                    self.assertEqual(p["submission_plans"], [])
                    self.assertIn("evidence differs", p["runs"][0]["issues"][0])

    def test_dispatch_recomputes_evidence_even_if_tampered_plan_is_rehashed(self):
        s = self.install_plan()
        s["execution_evidence"]["nominal_checkpoint_reserve_s"] = 99999
        self.reidentify(s)
        with self.assertRaisesRegex(ValueError, "Execution evidence differs"):
            self.h.dispatch(apply=True)
        self.assertEqual(self.h.calls, [])

    def test_dispatch_rejects_missing_policy_identity_and_ledger_only_override(self):
        s = self.install_plan()
        s.pop("execution_walltime_policy")
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "requires both"):
            self.h.dispatch()
        self.campaign["runs"][0].pop("execution_walltime_policy")
        self.campaign["runs"][0]["walltime_s"] = 7200
        s = self.install_plan()
        s["ledger_intent"]["execution_walltime_policy"] = M.SINGLE_ENGINE_1H_POLICY
        self.h.save_inputs()
        with self.assertRaisesRegex(ValueError, "unbound execution"):
            self.h.dispatch()

    def test_existing_accepted_two_hour_receipt_is_unchanged_and_blocks_replacement(self):
        self.campaign["runs"][0].pop("execution_walltime_policy")
        self.campaign["runs"][0]["walltime_s"] = 7200
        old = self.plan()["submission_plans"][0]["ledger_intent"] | {"state": "submitted", "job_id": "777"}
        self.h.ledger["submissions"].append(copy.deepcopy(old))
        self.campaign["runs"][0].update(walltime_s=3600, execution_walltime_policy=M.SINGLE_ENGINE_1H_POLICY)
        self.assertEqual(self.plan()["submission_plans"], [])
        self.assertEqual(self.h.ledger["submissions"], [old])


class ProductionPinTests(unittest.TestCase):
    def test_eight_real_frozen_sources_are_verified_read_only(self):
        self.assertEqual(len(M.SINGLE_ENGINE_1H_PINS), 8)
        for name, pin in M.SINGLE_ENGINE_1H_PINS.items():
            with self.subTest(run=name):
                frozen = M.snapshot_info(PROJECT/"runs/v2/snapshots"/name, name, pin[0])
                config = M.read_json(frozen["config_path"])
                e = M.validate_execution_walltime(frozen, config, 16, 3600,
                                                 M.SINGLE_ENGINE_1H_POLICY, nchunks=pin[2])
                self.assertEqual(e["nominal_checkpoint_reserve_s"], 900)
                self.assertLess(e["runner_evidence"]["handler_installation"]["line"],
                                e["runner_evidence"]["stop_loop"]["line"])
                self.assertIsNone(M.validate_execution_walltime(frozen, config, 16, 7200))

    def test_core_mechanism_and_other_single_engine_runs_remain_excluded(self):
        for name in ("confirm_core_a", "confirm_mechanism_b", "confirm_n4_a",
                     "confirm_pair_partial_e4_a", "confirm_pair_partial_e2_a"):
            with self.subTest(run=name):
                self.assertNotIn(name, M.SINGLE_ENGINE_1H_PINS)
                with self.assertRaisesRegex(ValueError, "allowlisted run"):
                    M.validate_execution_walltime({"sha256_manifest": "a"*64},
                        {"run_id": name}, 16, 3600, M.SINGLE_ENGINE_1H_POLICY, nchunks=2)


if __name__ == "__main__":
    unittest.main()
