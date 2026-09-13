"""Fast paired-job scheduling, durable append, and per-engine inference tests.

Run: python -m unittest discover -s scripts/v2 -p test_paired_runner_v2.py -v
The scientific imports are required, but every fit and worker pool is a fixture.
"""
from contextlib import contextmanager, redirect_stdout
from copy import deepcopy
import hashlib
import io
from itertools import product
import json
from pathlib import Path
import signal
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np

from engine_v2 import FitResult
import paired_engine_v2 as P
import run_v2 as R


def config(**overrides):
    return {"run_id": "paired_runner_fixture", "phase": "development",
            "base_seed": 502, "datasets": 1, "engines": ["original", "repaired"],
            "cells": [{"name": "REF_H0"}], "schemes": {}, "bperm": 0,
            "paired_standard": True, "fit_timeout_s": 30,
            "full_p_seeds": [502], "permutation_batch": 2} | overrides


def job(cfg=None, cell=None, seed=502, engine="__paired__", scheme=None, pid=None):
    cfg = config() if cfg is None else cfg
    return cfg, cfg["cells"][0] if cell is None else cell, seed, engine, scheme, pid


def shared_id(task):
    cfg, cell, seed, _, scheme, pid = task
    return hashlib.sha256(json.dumps([cfg["run_id"], cell, seed, scheme, pid],
                                     sort_keys=True).encode()).hexdigest()


def event(task, engine, *, status="ok", taus=None, seconds=1.):
    _, _, seed, _, scheme, pid = task
    return {"kind": "observed" if scheme is None else "permutation", "seed": seed,
            "engine": engine, "scheme": scheme, "perm_id": pid, "status": status,
            "orderings": [[0, 1], [1, 0], [0, 1]] if status == "ok" else None,
            "taus": ([.5] * 3 if taus is None else taus) if status == "ok" else None,
            "truth": {"pair_null": [True] * 3}, "diagnostics": {},
            "fit_seconds": seconds, "shared_fit_id": shared_id(task)}


class FakePaired:
    def __init__(self, result=None):
        self.calls = []
        self.result = result

    def __call__(self, task):
        self.calls.append(task)
        if self.result is not None:
            return self.result(task)
        return [event(task, engine, seconds=float(len(self.calls)))
                for engine in ("original", "repaired")]


class FakePool:
    """Exercise execute_job and append code without spawning numerical workers."""
    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def imap_unordered(self, function, jobs, chunksize):
        # Different completion order must not change identity or pairing.
        for task in reversed(jobs):
            yield function(task)


@contextmanager
def runner_fixture(cfg, output, fake):
    cfg_path = output.parent / "config.json"
    cfg_path.write_text(json.dumps(cfg))
    argv = ["run_v2.py", "--config", str(cfg_path), "--output", str(output),
            "--workers", "1"]
    old_signals = {sig: signal.getsignal(sig) for sig in
                   (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1)}
    try:
        with patch.object(R, "get_context", return_value=SimpleNamespace(
                Pool=lambda *args: FakePool())), \
                patch.object(R, "verify_pyebm", return_value={"versions": {"fixture": "1"}}), \
                patch.object(R, "source_hashes", return_value={"run_v2.py": "a" * 64}), \
                patch.object(R, "fit_paired_job", side_effect=fake), \
                patch.object(R, "STOP_REQUESTED", False), \
                patch.object(R.sys, "argv", argv), redirect_stdout(io.StringIO()):
            yield
    finally:
        for sig, handler in old_signals.items():
            signal.signal(sig, handler)


def read_jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


class PairedSchedulingTests(unittest.TestCase):
    def test_one_standard_pair_is_one_job_and_oracles_stay_separate(self):
        cfg = config(engines=["original", "repaired", "oracle_equal_original"])
        tasks = [job(cfg, engine=engine) for engine in cfg["engines"]]
        collapsed = R.collapse_paired_jobs(tasks, cfg)
        self.assertEqual([task[3] for task in collapsed], ["__paired__", "oracle_equal_original"])
        self.assertEqual(R.collapse_paired_jobs([tasks[1]], cfg)[0][3], "__paired__")

    def test_disabled_pairing_retains_jobs_and_missing_engine_is_rejected(self):
        cfg = config(paired_standard=False)
        tasks = [job(cfg, engine=engine) for engine in cfg["engines"]]
        self.assertEqual(R.collapse_paired_jobs(tasks, cfg), tasks)
        with self.assertRaisesRegex(ValueError, "requires both"):
            R.collapse_paired_jobs(tasks, config(engines=["original"]))

    def test_pairs_never_merge_different_datasets_schemes_or_permutations(self):
        cfg = config()
        identities = [job(cfg), job(cfg, seed=503),
                      job(cfg, cell={"name": "STAGE_H0"}),
                      job(cfg, cell={"name": "REF_H0", "overrides": {"n_multiplier": 4}}),
                      job(config(run_id="different_run")),
                      job(cfg, scheme="dx", pid=0), job(cfg, scheme="dx", pid=1),
                      job(cfg, scheme="unstratified", pid=0)]
        tasks = [tuple(task[:3]) + (engine,) + tuple(task[4:])
                 for task in identities for engine in ("repaired", "original")]
        collapsed = R.collapse_paired_jobs(tasks, cfg)
        self.assertEqual(len(collapsed), len(identities))
        self.assertEqual({shared_id(task) for task in collapsed},
                         {shared_id(task) for task in identities})

    def test_execute_job_returns_a_batch_for_either_path(self):
        pair, single = [object(), object()], object()
        with patch.object(R, "fit_paired_job", return_value=pair) as paired, \
                patch.object(R, "fit_job", return_value=single) as separate:
            self.assertIs(R.execute_job(job()), pair)
            paired.assert_called_once()
            separate.assert_not_called()
            self.assertEqual(R.execute_job(job(engine="oracle_equal_original")), [single])
            separate.assert_called_once()


class PairedFitRecordTests(unittest.TestCase):
    def setUp(self):
        self.truth = {"group_orderings": {"e2": [1, 0], "e33": [0, 1], "e4": [0, 1]},
                      "biomarker_names": ["a", "b"]}
        self.results = {engine: FitResult([np.array([0, 1])] * 3, "ok", {
            "optimizer": [{"success": True, "status": 0}],
            "consensus": [{"mode": engine}], "marker": engine})
            for engine in ("original", "repaired")}

    @contextmanager
    def fitted(self, result=None, error=None):
        with patch.object(R, "resolve", return_value="fixture_design"), \
                patch.object(R, "simulate", side_effect=lambda *_: ("observed_frame", deepcopy(self.truth))), \
                patch.object(R, "permute", return_value="permuted_frame") as permute, \
                patch.object(R.signal, "signal"), patch.object(R.signal, "alarm"), \
                patch.object(P, "fit_paired_orderings", return_value=self.results if result is None else result,
                             side_effect=error) as fit:
            yield fit, permute

    def test_both_observed_results_share_one_input_and_retain_truth(self):
        with self.fitted() as (fit, permute):
            records = R.fit_paired_job(job())
        fit.assert_called_once()
        permute.assert_not_called()
        self.assertEqual(fit.call_args.args, ("observed_frame",))
        self.assertEqual(fit.call_args.kwargs["abort_exceptions"], (R.FitTimeout,))
        self.assertEqual([r["engine"] for r in records], ["original", "repaired"])
        self.assertEqual(len({r["shared_fit_id"] for r in records}), 1)
        self.assertEqual(records[0]["fit_seconds"], records[1]["fit_seconds"])
        for record in records:
            self.assertEqual(record["truth"]["pair_null"], [False, False, True])
            self.assertEqual(record["distance_to_truth"], [1., 0., 0.])
            self.assertEqual(record["taus"], [0.] * 3)
            self.assertTrue(record["diagnostics"]["paired_standard"])
            json.dumps(record, allow_nan=False)

    def test_each_permutation_refits_once_with_its_own_input_and_compact_diagnostics(self):
        cfg = config(schemes={"dx": {"stratify": "diagnosis"}})
        with self.fitted() as (fit, permute):
            first = R.fit_paired_job(job(cfg, scheme="dx", pid=0))
            second = R.fit_paired_job(job(cfg, scheme="dx", pid=1))
        self.assertEqual(fit.call_count, 2)
        self.assertEqual([call.args[4] for call in permute.call_args_list], [0, 1])
        self.assertTrue(all(call.args == ("permuted_frame",) for call in fit.call_args_list))
        self.assertNotEqual(first[0]["shared_fit_id"], second[0]["shared_fit_id"])
        for record in first + second:
            self.assertNotIn("truth", record)
            self.assertNotIn("optimizer", record["diagnostics"])
            self.assertEqual(record["diagnostics"]["optimizer_calls"], 1)

    def test_timeout_retains_completed_original_and_marks_only_missing_repaired_error(self):
        exc = R.FitTimeout("fixture second-engine timeout")
        exc.completed_results = {"original": self.results["original"]}
        with self.fitted(error=exc):
            records = R.fit_paired_job(job())
        self.assertEqual([r["status"] for r in records], ["ok", "error"])
        self.assertEqual(records[0]["distance_to_truth"], [1., 0., 0.])
        self.assertIsNone(records[1]["orderings"])
        self.assertIn("FitTimeout", records[1]["diagnostics"]["error"])
        for record in records:
            self.assertEqual(record["truth"]["pair_null"], [False, False, True])
            self.assertEqual(record["diagnostics"]["paired_total_timeout_s"], 60)

    def test_failure_before_any_completed_engine_records_two_failures(self):
        with self.fitted(error=ValueError("fixture fit failure")):
            records = R.fit_paired_job(job())
        self.assertEqual([r["status"] for r in records], ["error", "error"])
        self.assertEqual(len({r["shared_fit_id"] for r in records}), 1)
        self.assertTrue(all("fixture fit failure" in r["diagnostics"]["error"] for r in records))


class PairedDurabilityTests(unittest.TestCase):
    def test_torn_second_append_replays_pair_retains_first_and_does_not_duplicate(self):
        cfg, fake = config(), FakePaired()
        with TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            with runner_fixture(cfg, output, fake):
                R.main()
                log = output / "fit_records.jsonl"
                original_bytes = log.read_bytes().splitlines(keepends=True)
                first = json.loads(original_bytes[0])
                tail = original_bytes[1][:37]
                log.write_bytes(original_bytes[0] + tail)
                R.main()
                restored = R.load_events(log)
                self.assertEqual(len(fake.calls), 2)
                self.assertEqual(len(restored), 2)
                self.assertEqual(restored[R.key(first)], first)
                self.assertEqual(len(read_jsonl(log)), 2)
                self.assertEqual(read_jsonl(output / "replay_audit.jsonl")[0]["policy"],
                                 "keep_original_record")
                quarantines = list(output.glob("fit_records.jsonl.torn-*"))
                self.assertEqual(len(quarantines), 1)
                self.assertEqual(quarantines[0].read_bytes(), tail)
                progress = json.loads((output / "progress.json").read_text())
                self.assertTrue(progress["done"])
                self.assertEqual((progress["fit_records"], progress["computational_jobs"]), (2, 1))
                saved_bytes = log.read_bytes()
                R.main()
                self.assertEqual(len(fake.calls), 2)
                self.assertEqual(log.read_bytes(), saved_bytes)

    def test_existing_repaired_half_is_verified_without_duplicate(self):
        cfg, fake = config(), FakePaired()
        with TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            output.mkdir()
            kept = event(job(cfg), "repaired", seconds=17.)
            log = output / "fit_records.jsonl"
            log.write_text(json.dumps(kept) + "\n")
            with runner_fixture(cfg, output, fake):
                R.main()
            records = R.load_events(log)
            self.assertEqual(records[R.key(kept)], kept)
            self.assertEqual(len(records), 2)
            self.assertEqual(len(fake.calls), 1)
            audit = read_jsonl(output / "replay_audit.jsonl")
            self.assertEqual(audit[0]["key"], list(R.key(kept)))

    def test_saved_failed_observed_result_is_not_replaced_by_successful_replay(self):
        cfg, fake = config(), FakePaired()
        with TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            output.mkdir()
            kept = event(job(cfg), "original", status="error")
            log = output / "fit_records.jsonl"
            log.write_text(json.dumps(kept) + "\n")
            with runner_fixture(cfg, output, fake):
                R.main()
            records = R.load_events(log)
            self.assertEqual(records[R.key(kept)], kept)
            self.assertEqual(records[(502, "repaired", None, None)]["status"], "ok")
            audit = read_jsonl(output / "replay_audit.jsonl")[0]
            self.assertEqual((audit["old_status"], audit["new_status"]), ("error", "ok"))
            rows = {r["engine"]: r for r in read_jsonl(output / "rows.jsonl")}
            self.assertEqual(rows["original"]["status"], "error")

    def test_replay_disagreement_cannot_append_new_companion_before_verification(self):
        # Either half may have survived. A rejected replay must not leave a newly
        # appended companion paired with a saved result from a different fit.
        for existing_engine, changed_field in product(("original", "repaired"), ("taus", "orderings")):
            with self.subTest(existing_engine=existing_engine, changed_field=changed_field), \
                    TemporaryDirectory() as folder:
                cfg, fake = config(), FakePaired()
                output = Path(folder) / "run"
                output.mkdir()
                kept = event(job(cfg), existing_engine)
                kept[changed_field] = [.1, .2, .3] if changed_field == "taus" else [[0, 1]] * 3
                log = output / "fit_records.jsonl"
                before = (json.dumps(kept) + "\n").encode()
                log.write_bytes(before)
                with runner_fixture(cfg, output, fake):
                    with self.assertRaisesRegex(ValueError, "Deterministic paired replay disagrees"):
                        R.main()
                self.assertEqual(log.read_bytes(), before)

    def test_full_budget_counts_failures_and_compute_units_remain_per_engine(self):
        cfg = config(schemes={"dx": {"stratify": "diagnosis"}}, bperm=3)
        statistics = {"original": [[.6, .5, .4], [.4, .6, .5], [.9, .9, .9]],
                      "repaired": [[.3, .4, .5], [.2, .6, .7], [.1, .2, .3]]}
        fake = FakePaired(lambda task: [event(task, engine, taus=statistics[engine][task[5]])
                                      for engine in ("original", "repaired")])
        with TemporaryDirectory() as folder:
            output = Path(folder) / "run"
            output.mkdir()
            initial = [event(job(cfg), "original", taus=[.5] * 3),
                       event(job(cfg), "repaired", taus=[.2, .4, .6]),
                       event(job(cfg, scheme="dx", pid=2), "original", status="error")]
            log = output / "fit_records.jsonl"
            log.write_text("".join(json.dumps(r) + "\n" for r in initial))
            with runner_fixture(cfg, output, fake):
                R.main()
                R.main()
            self.assertEqual(len(fake.calls), 3)
            records = R.load_events(log)
            self.assertEqual(records[R.key(initial[-1])], initial[-1])
            rows = {r["engine"]: r for r in read_jsonl(output / "rows.jsonl")}
            original, repaired = (rows[e]["schemes"]["dx"] for e in ("original", "repaired"))
            self.assertEqual((original["requested_nperm"], original["attempted"],
                              original["nperm"], original["failed"], original["complete"]),
                             (3, 3, 2, 1, False))
            self.assertEqual((repaired["requested_nperm"], repaired["attempted"],
                              repaired["nperm"], repaired["failed"], repaired["complete"]),
                             (3, 3, 3, 0, True))
            self.assertEqual((original["gt"], original["eq"]), ([1, 1, 0], [0, 1, 1]))
            self.assertEqual((repaired["gt"], repaired["eq"]), ([1, 1, 1], [1, 1, 0]))
            progress = json.loads((output / "progress.json").read_text())
            self.assertEqual(progress["fit_records"], 8)
            self.assertEqual(progress["computational_jobs"], 4)
            self.assertEqual(progress["statuses"], {"ok": 7, "error": 1})
            self.assertEqual(len(read_jsonl(output / "replay_audit.jsonl")), 1)


if __name__ == "__main__":
    unittest.main()
