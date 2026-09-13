"""Independent core signal-policy regressions; no scheduler or scientific calls.

Fixtures copy inert frozen source/configuration into an isolated temporary root.
Only the root literal and corresponding in-process pins change. All submissions
are injected mocks and all ledgers live beneath that temporary root.
"""
import copy
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M

PROJECT = Path(__file__).resolve().parents[2]
RUN = "confirm_core_a"


class IndependentCorePolicyTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(prefix="independent_core_signal_")
        self.addCleanup(tmp.cleanup)
        self.root = (Path(tmp.name) / "ebmcal_v2").resolve()
        self.snapshot = self.root / "snapshots" / RUN
        shutil.copytree(PROJECT / "runs/v2/snapshots" / RUN, self.snapshot,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        script = self.snapshot / "hpc/v2/run_array.sbatch"
        body = script.read_text()
        self.assertEqual(body.count("root=/N/slate/tg11/ebmcal_v2"), 1)
        script.write_text(body.replace("root=/N/slate/tg11/ebmcal_v2", f"root={self.root}"))
        manifest = self.snapshot / "SHA256SUMS"
        manifest.write_text("".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(self.snapshot)}\n"
            for path in sorted(self.snapshot.rglob("*"))
            if path.is_file() and path != manifest))
        pin = (hashlib.sha256(manifest.read_bytes()).hexdigest(),
               hashlib.sha256((self.snapshot / "configs/v2" / f"{RUN}.json").read_bytes()).hexdigest())
        source_pins = {name: hashlib.sha256((self.snapshot / name).read_bytes()).hexdigest()
                       for name in M.CORE_PAIRED_1H_SOURCE_PINS}
        patch = mock.patch.object(M, "CORE_PAIRED_1H_PIN", pin)
        patch.start(); self.addCleanup(patch.stop)
        patch = mock.patch.dict(M.CORE_PAIRED_1H_SOURCE_PINS, source_pins, clear=True)
        patch.start(); self.addCleanup(patch.stop)
        self.now = datetime(2026, 9, 8, 12, 0, tzinfo=timezone.utc)
        self.queue = {"schema_version": "campaign_queue_v1", "captured_utc": self.now.isoformat(),
                      "scope": "all_user_jobs", "complete": True, "jobs": []}
        self.ledger = {"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": []}
        self.campaign = {"project": "ebmcal_v2", "max_project_cpus": 32,
                         "max_user_submitted_jobs": 900, "max_tasks_per_submission": 2,
                         "runs": [{"run_id": RUN, "snapshot": str(self.snapshot), "snapshot_sha256": pin[0],
                                   "output_root": str(self.root / "runs" / RUN), "nchunks": 500,
                                   "cpus_per_task": 16, "max_parallel_tasks": 1, "walltime_s": 3600,
                                   "execution_walltime_policy": M.CORE_PAIRED_1H_POLICY}]}
        self.plan = M.make_plan(self.campaign, self.queue, self.ledger, self.now, weekly_used=12)
        self.assertEqual(len(self.plan["submission_plans"]), 1, self.plan)
        self.plan_path, self.queue_path, self.ledger_path = [self.root / name for name in
                                                          ("plan.json", "queue.json", "ledger.json")]
        self.save()
        self.calls = []

    def save(self, *, ledger=True):
        pairs = [(self.plan_path, self.plan), (self.queue_path, self.queue)]
        if ledger:
            pairs.append((self.ledger_path, self.ledger))
        for path, data in pairs:
            path.write_text(json.dumps(data))

    def submit(self, argv, **kwargs):
        saved = json.loads(self.ledger_path.read_text())
        self.assertEqual(saved["submissions"][-1]["state"], "intent")
        self.assertEqual(argv, saved["submissions"][-1]["argv"])
        self.assertIs(kwargs["shell"], False)
        self.calls.append(argv[:])
        return subprocess.CompletedProcess(argv, 0, "987654\n", "")

    def dispatch(self, *, apply=True, submitter=None, **kwargs):
        return D.dispatch(self.plan_path, self.ledger_path, self.queue_path, self.root,
                          apply=apply, submitter=submitter or self.submit, now_func=lambda: self.now, **kwargs)

    def reidentify(self, submission):
        identity = D.submission_identity(submission)
        submission["plan_id"] = M.I.digest(identity)
        submission["argv"][5] = "--comment=ebmcal_v2:" + submission["plan_id"]
        submission["ledger_intent"] = identity | {"plan_id": submission["plan_id"], "state": "intent",
                                                   "created_utc": self.now.isoformat()}

    def test_exact_signal_durable_intent_and_duplicate_plan(self):
        s = self.plan["submission_plans"][0]
        self.assertEqual((s["array_task_ids"], s["nchunks"], s["max_concurrent"]), ([0, 1], 500, 1))
        self.assertEqual(len(s["argv"]), 13)
        self.assertEqual(s["argv"][8:13], ["--signal=B:USR1@2100",
            str(self.snapshot / "hpc/v2/run_array.sbatch"), str(self.snapshot), RUN, "500"])
        self.assertIs(s["execution_evidence"]["signal_effective_verified"], False)
        before = self.ledger_path.read_bytes()
        self.assertEqual(self.dispatch(apply=False)["status"], "validated_dry_run")
        self.assertEqual(self.ledger_path.read_bytes(), before)
        self.assertEqual(self.calls, [])
        self.assertEqual(self.dispatch()["status"], "submitted")
        receipt = json.loads(self.ledger_path.read_text())["submissions"][0]
        self.assertEqual(D.submission_identity(receipt), D.submission_identity(s))
        self.assertEqual(receipt["plan_id"], M.I.digest(D.submission_identity(receipt)))
        with self.assertRaisesRegex(ValueError, "already has a ledger"):
            self.dispatch()
        self.assertEqual(len(self.calls), 1)

    def test_signal_deletion_duplication_mutation_and_argument_shift_fail_before_intent(self):
        base = copy.deepcopy(self.plan["submission_plans"][0])
        variants = []
        for signal in ("--signal=B:USR1@1500", "--signal=USR1@2100", "--signal=B:USR2@2100"):
            argv = base["argv"][:]; argv[8] = signal; variants.append(argv)
        variants += [base["argv"][:8] + base["argv"][9:],
                     base["argv"][:8] + [base["argv"][8]] + base["argv"][8:],
                     base["argv"][:8] + base["argv"][9:] + [base["argv"][8]],
                     base["argv"][:8] + [base["argv"][9], base["argv"][8]] + base["argv"][10:]]
        for argv in variants:
            with self.subTest(argv=argv):
                self.plan["submission_plans"][0] = copy.deepcopy(base) | {"argv": argv}
                self.save()
                with self.assertRaises(ValueError):
                    self.dispatch()
                self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])
        self.assertEqual(self.calls, [])

    def test_resigned_forged_evidence_and_fake_policy_are_reverified(self):
        base = copy.deepcopy(self.plan["submission_plans"][0])
        forged_hashes = base["execution_evidence"]["source_sha256"] | {"scripts/v2/run_v2.py": "0" * 64}
        mutations = [("signal_effective_verified", True), ("requested_signal", "B:USR1@1500"),
                     ("nominal_batch_drain_s", 600), ("source_sha256", forged_hashes)]
        for key, value in mutations + [("__policy__", "pretend_core_paired_signal_1h_v1")]:
            with self.subTest(key=key):
                s = copy.deepcopy(base)
                if key == "__policy__":
                    s["execution_walltime_policy"] = value
                    s["argv"].pop(8)  # Even a valid-length command cannot authorize a fake policy.
                else:
                    s["execution_evidence"][key] = value
                self.reidentify(s)
                self.plan["submission_plans"][0] = s
                self.save()
                with self.assertRaises(ValueError):
                    self.dispatch()
                self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])
        self.assertEqual(self.calls, [])

    def test_unknown_submission_blocks_same_targets_under_new_generation_or_policy(self):
        def uncertain(argv, **kwargs):
            self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"][-1]["state"], "intent")
            self.calls.append(argv[:])
            raise subprocess.TimeoutExpired(argv, 30, output=b"987654\n")
        result = self.dispatch(submitter=uncertain)
        self.assertEqual(result["status"], "submission_unknown_reconciliation_required")
        old_bytes = self.ledger_path.read_bytes()
        s = self.plan["submission_plans"][0]
        s["dispatch_generation"] = 2
        # Switching back to the historical two-hour format is not a retry escape.
        s.pop("execution_walltime_policy"); s.pop("execution_evidence")
        s["walltime_s"] = 7200; s["argv"][4] = "--time=02:00:00"; s["argv"].pop(8)
        self.reidentify(s); self.save(ledger=False)
        with self.assertRaises(ValueError):
            self.dispatch()
        self.assertEqual(self.ledger_path.read_bytes(), old_bytes)
        self.assertEqual(len(self.calls), 1)

    def test_ordinary_identity_remains_nine_fields_and_thirteenth_option_rejected(self):
        run = self.campaign["runs"][0]
        run.pop("execution_walltime_policy"); run["walltime_s"] = 7200
        self.plan = M.make_plan(self.campaign, self.queue, self.ledger, self.now, weekly_used=12)
        s = self.plan["submission_plans"][0]
        self.assertEqual(len(s["argv"]), 12)
        self.assertEqual(set(D.submission_identity(s)), set(D.IDENTITY_FIELDS))
        self.assertEqual(s["plan_id"], M.I.digest({key: s[key] for key in D.IDENTITY_FIELDS}))
        self.save()
        self.assertEqual(self.dispatch(apply=False)["status"], "validated_dry_run")
        s["argv"].insert(8, "--signal=B:USR1@2100"); self.save()
        with self.assertRaisesRegex(ValueError, "12-argument"):
            self.dispatch()
        self.assertEqual(self.calls, [])

    def test_latest_guard_runs_after_work_scan_and_before_intent(self):
        recheck = D.recheck_work
        count = 0
        def slow_second_scan(*args, **kwargs):
            nonlocal count
            result = recheck(*args, **kwargs)
            count += 1
            if count == 2:
                self.now += timedelta(seconds=61)
            return result
        with mock.patch.object(D, "recheck_work", side_effect=slow_second_scan):
            result = self.dispatch()
        self.assertEqual(result["status"], "stopped_before_next_submission")
        self.assertEqual(result["outcomes"][0]["state"], "not_dispatched")
        self.assertIn("stale", result["outcomes"][0]["reason"])
        self.assertEqual(count, 2)
        self.assertEqual(self.calls, [])
        self.assertEqual(json.loads(self.ledger_path.read_text())["submissions"], [])


class IndependentProductionMetadataTest(unittest.TestCase):
    def test_source_and_config_pins_are_real_bytes_and_nominal_evidence_is_unverified(self):
        snapshot = PROJECT / "runs/v2/snapshots" / RUN
        cfg_path = snapshot / "configs/v2" / f"{RUN}.json"
        self.assertEqual(M.CORE_PAIRED_1H_PIN, (
            "e4701a864230fbc6faf64709defb4dec5fa5a3d07ac3625b95cc36cb06847464",
            "11b716e46bb5055f8b24a3ee5b756683508f66165fe3c361b359b477edf30a52"))
        self.assertEqual(hashlib.sha256((snapshot / "SHA256SUMS").read_bytes()).hexdigest(), M.CORE_PAIRED_1H_PIN[0])
        self.assertEqual(hashlib.sha256(cfg_path.read_bytes()).hexdigest(), M.CORE_PAIRED_1H_PIN[1])
        for name, expected in M.CORE_PAIRED_1H_SOURCE_PINS.items():
            self.assertEqual(hashlib.sha256((snapshot / name).read_bytes()).hexdigest(), expected)
        frozen = M.snapshot_info(snapshot, RUN, M.CORE_PAIRED_1H_PIN[0])
        evidence = M.validate_execution_walltime(frozen, json.loads(cfg_path.read_text()), 16, 3600,
                                                M.CORE_PAIRED_1H_POLICY, nchunks=500)
        self.assertEqual((evidence["nominal_batch_waves"], evidence["paired_total_timeout_s"],
                          evidence["nominal_batch_drain_s"], evidence["nominal_checkpoint_reserve_s"]),
                         (2, 600, 1200, 900))
        self.assertIs(evidence["signal_effective_verified"], False)


if __name__ == "__main__":
    unittest.main()
