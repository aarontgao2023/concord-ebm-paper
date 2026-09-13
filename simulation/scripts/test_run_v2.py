"""Fixed-budget inference, permutation coupling, and durable resume checks.

Run with the scientific environment and pinned pyebm wheel on PYTHONPATH:
python -m unittest discover -s scripts/v2 -p test_run_v2.py -v
Only one integration test fits DEBM; the exhaustive/resume tests use synthetic
statistics so their independent reference decisions remain inspectable.
"""

from contextlib import redirect_stdout
from fractions import Fraction
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
import pandas as pd

import design_v2 as D
from engine_v2 import FitResult
import run_v2 as R


def full_rejection(exceedances, budget, alpha, comparisons, rule):
    value = Fraction(1 + exceedances, 1 + budget)
    threshold = Fraction(str(alpha)) / comparisons
    return value < threshold if rule == "strict" else value <= threshold


def consensus(values):
    choices = set(values)
    return next(iter(choices)) if len(choices) == 1 else None


def event(seed=502, engine="original", scheme=None, pid=None, taus=None, status="ok"):
    return {"kind": "observed" if scheme is None else "permutation", "seed": seed,
            "engine": engine, "scheme": scheme, "perm_id": pid, "status": status,
            "taus": [.5, .5, .5] if taus is None else taus,
            "orderings": [[0, 1], [1, 0], [0, 1]] if status == "ok" else None,
            "truth": {"pair_null": [True, True, True]}, "diagnostics": {},
            "fit_seconds": 0.}


class BoundTests(unittest.TestCase):
    def test_bounds_equal_every_possible_completion_in_small_budgets(self):
        # Ties and strict exceedances both count; unknown fits can have either
        # outcome. This enumerates all completions, not the implementation's bounds.
        for budget in range(1, 8):
            for successful in range(budget + 1):
                for known in range(successful + 1):
                    for tied in range(known + 1):
                        gt, eq = [known - tied] * 3, [tied] * 3
                        for alpha in (.05, .3, .6):
                            answer = R.decision_bounds(gt, eq, successful, budget, alpha)
                            unknown = budget - successful
                            for rule in ("strict", "le"):
                                expected = consensus(full_rejection(known + sum(tail), budget,
                                                                     alpha, 3, rule)
                                                     for tail in product((0, 1), repeat=unknown))
                                self.assertEqual(answer[rule]["pair_reject"], [expected] * 3)
                                expected_max = consensus(full_rejection(known + sum(tail), budget,
                                                                         alpha, 1, rule)
                                                         for tail in product((0, 1), repeat=unknown))
                                self.assertIs(R.single_decision(known, successful, budget, alpha, rule),
                                              expected_max)

    def test_b599_boundary_and_b1999_early_stop(self):
        full599 = R.decision_bounds([8, 9, 10], [0] * 3, 599, 599)
        self.assertEqual(full599["strict"]["pair_reject"], [True, False, False])
        self.assertEqual(full599["le"]["pair_reject"], [True, True, False])
        before = R.decision_bounds([32] * 3, [0] * 3, 32, 1999)
        after = R.decision_bounds([20] * 3, [13] * 3, 33, 1999)
        for rule in ("strict", "le"):
            self.assertEqual(before[rule]["pair_reject"], [None] * 3)
            self.assertEqual(after[rule]["pair_reject"], [False] * 3)

    def test_failed_fit_remains_unknown_in_the_planned_denominator(self):
        records = [event()]
        records += [event(scheme="dx", pid=i, taus=[0., 0., 0.]) for i in range(97)]
        records += [event(scheme="dx", pid=i, status="error") for i in (97, 98)]
        state = R.state_for({R.key(r): r for r in records}, 502, "original", "dx", 99, .05)
        self.assertEqual((state["attempted"], state["nperm"], state["failed"]), (99, 97, 2))
        self.assertFalse(state["complete"])
        self.assertEqual(state["exact_decision"]["strict"]["pair_reject"], [None] * 3)
        # Dropping two errors would falsely certify rejection: 1/98 < .05/3.
        self.assertTrue(full_rejection(0, 97, .05, 3, "strict"))
        self.assertFalse(full_rejection(2, 99, .05, 3, "strict"))

    def test_requested_endpoints_determine_stopping(self):
        state = {"exact_decision": {"strict": {"pair_reject": [False, None, None],
                                                "max_reject": None}}}
        self.assertFalse(R.resolved_for(state, {}, "strict"))
        self.assertTrue(R.resolved_for(state, {"tested_pair_index": 0}, "strict"))
        self.assertFalse(R.resolved_for(state, {"tested_pair_index": 0, "require_max": True}, "strict"))
        state["exact_decision"]["strict"]["pair_reject"] = [False, False, False]
        self.assertTrue(R.resolved_for(state, {}, "strict"))
        self.assertFalse(R.resolved_for(state, {"require_max": True}, "strict"))


class PermutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df, cls.truth = D.simulate(D.resolve("REF_H0"), 9643)

    def test_all_three_strata_preserve_their_required_counts(self):
        for stratify in ("none", "diagnosis", "oracle_dx_stage"):
            before = self.df.copy(deep=True)
            after = R.permute(self.df, self.truth, 9643, "scheme", 17, {"stratify": stratify})
            if stratify == "none":
                strata = np.zeros(len(before), dtype=int)
            elif stratify == "diagnosis":
                strata = before.Diagnosis.to_numpy()
            else:
                # Diagnosis enters the DEBM fit. Preserving k alone is insufficient
                # when group membership still depends on diagnosis conditional on k.
                strata = np.asarray([f"{dx}:{stage}" for dx, stage in
                                     zip(before.Diagnosis, self.truth["latent_stage"])])
            for level in np.unique(strata):
                np.testing.assert_array_equal(np.bincount(before.APOE.to_numpy()[strata == level], minlength=3),
                                              np.bincount(after.APOE.to_numpy()[strata == level], minlength=3))
            pd.testing.assert_frame_equal(before.drop(columns="APOE"), after.drop(columns="APOE"), check_exact=True)
            pd.testing.assert_frame_equal(before, self.df, check_exact=True)

    def test_observed_count_stratification_preserves_diagnosis_by_count_tables(self):
        df, truth = D.simulate(D.resolve('MISS_ALL_H0'), 31000300)
        names = list(truth['biomarker_names'])
        permuted = R.permute(df, truth, 31000300, 'diagnosis_observed_count', 3,
                           {'stratify': 'diagnosis_observed_count'})
        self.assertTrue(permuted[names].equals(df[names]))
        self.assertFalse((permuted.APOE.to_numpy() == df.APOE.to_numpy()).all())
        counts = lambda f: f.assign(nobs=f[names].notna().sum(axis=1)).groupby(
            ['Diagnosis', 'nobs', 'APOE']).size().sort_index()
        self.assertTrue(counts(permuted).equals(counts(df)))

    def test_pair_only_permutation_leaves_third_group_fixed(self):
        for allowed in ([0, 1], [0, 2], [1, 2]):
            for stratify in ("none", "diagnosis", "oracle_dx_stage"):
                after = R.permute(self.df, self.truth, 9643, "pair", 9,
                                  {"stratify": stratify, "groups": allowed})
                excluded = ~self.df.APOE.isin(allowed)
                pd.testing.assert_frame_equal(after[excluded], self.df[excluded], check_exact=True)
                self.assertEqual(after.APOE.value_counts().to_dict(), self.df.APOE.value_counts().to_dict())

    def test_engine_jobs_share_data_and_permutation_indices(self):
        captured = []

        def capture(frame, config):
            captured.append((frame.copy(deep=True), config.mode))
            return FitResult([np.arange(config.expected_events) for _ in range(3)], "ok", {})

        config = {"fit_timeout_s": 30, "schemes": {"dx": {"stratify": "diagnosis"}}}
        cell = {"name": "PWR_E2_K14"}
        previous = signal.getsignal(signal.SIGALRM)
        try:
            with patch.object(R, "fit_orderings", side_effect=capture):
                left = R.fit_job((config, cell, 8241, "original", "dx", 9))
                right = R.fit_job((config, cell, 8241, "repaired", "dx", 9))
                observed = R.fit_job((config, cell, 8241, "original", None, None))
                replay = R.fit_job((config, cell, 8241, "original", "dx", 9))
        finally:
            signal.signal(signal.SIGALRM, previous)
        self.assertEqual([left["status"], right["status"], observed["status"], replay["status"]], ["ok"] * 4)
        pd.testing.assert_frame_equal(captured[0][0], captured[1][0], check_exact=True)
        pd.testing.assert_frame_equal(captured[0][0], captured[3][0], check_exact=True)
        self.assertEqual(captured[0][1:], ("original",))
        self.assertEqual(captured[1][1:], ("repaired",))
        self.assertEqual(observed["truth"]["pair_null"], [False, False, True])
        self.assertEqual(observed["truth"]["realized_inversions"]["e2"], 14)
        pd.testing.assert_frame_equal(captured[0][0].drop(columns="APOE"),
                                      captured[2][0].drop(columns="APOE"), check_exact=True)

    def test_observed_failure_preserves_known_truth(self):
        previous = signal.getsignal(signal.SIGALRM)
        try:
            with patch.object(R, "fit_orderings", side_effect=R.FitTimeout("test timeout")):
                result = R.fit_job(({"fit_timeout_s": 30}, {"name": "PWR_E4_K27"},
                                    8641, "original", None, None))
        finally:
            signal.signal(signal.SIGALRM, previous)
        self.assertEqual(result["status"], "error")
        self.assertEqual(result["truth"]["pair_null"], [True, False, False])
        self.assertEqual(result["truth"]["realized_inversions"]["e4"], 27)


class ResumeTests(unittest.TestCase):
    def test_truncated_tail_is_quarantined_and_can_be_replayed(self):
        first = event()
        good = (json.dumps(first) + "\n").encode()
        for tail in (b'{"seed":502,"engine":"orig', b'{"error":"\xe7\xbb'):
            with TemporaryDirectory() as folder:
                path = Path(folder) / "events.jsonl"
                path.write_bytes(good + tail)
                restored = R.load_events(path)
                self.assertEqual(restored, {R.key(first): first})
                self.assertEqual(path.read_bytes(), good)
                quarantines = list(path.parent.glob(path.name + ".torn-*"))
                self.assertEqual(len(quarantines), 1)
                self.assertEqual(quarantines[0].read_bytes(), tail)
                second = event(scheme="dx", pid=0)
                with path.open("a") as handle:
                    handle.write(json.dumps(second) + "\n")
                self.assertEqual(len(R.load_events(path)), 2)

    def test_complete_last_json_without_newline_is_safe_to_append(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / "events.jsonl"
            first, second = event(), event(scheme="dx", pid=0)
            path.write_text(json.dumps(first))
            self.assertEqual(len(R.load_events(path)), 1)
            self.assertTrue(path.read_bytes().endswith(b"\n"))
            with path.open("a") as handle:
                handle.write(json.dumps(second) + "\n")
            self.assertEqual(len(R.load_events(path)), 2)

    def test_corrupted_middle_and_duplicate_records_are_rejected(self):
        line = json.dumps(event()).encode() + b"\n"
        for contents in (b'{"bad":\n' + line, line + b'{"bad":\n', line + line):
            with TemporaryDirectory() as folder:
                path = Path(folder) / "events.jsonl"
                path.write_bytes(contents)
                with self.assertRaises(ValueError):
                    R.load_events(path)
                self.assertEqual(path.read_bytes(), contents)

    def test_resuming_complete_run_does_not_recompute_and_rejects_config_drift(self):
        class SerialPool:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def imap_unordered(self, function, jobs, chunksize):
                yield from map(function, jobs)

        calls = []

        def cheap_fit(job):
            _, _, seed, engine, scheme, pid = job
            calls.append((seed, engine, scheme, pid))
            return event(seed, engine, scheme, pid)

        with TemporaryDirectory() as folder:
            root = Path(folder)
            cfg_path, output = root / "config.json", root / "run"
            cfg = {"run_id": "runner_test", "phase": "development", "base_seed": 502,
                   "datasets": 1, "engines": ["original"], "cells": [{"name": "REF_H0"}],
                   "schemes": {"dx": {"stratify": "diagnosis"}}, "bperm": 3,
                   "full_p_seeds": [502], "permutation_batch": 2}
            cfg_path.write_text(json.dumps(cfg))
            args = ["run_v2.py", "--config", str(cfg_path), "--output", str(output), "--workers", "1"]
            old_signals = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1)}
            try:
                with patch.object(R, "get_context", return_value=SimpleNamespace(Pool=lambda *a, **k: SerialPool())), \
                        patch.object(R, "verify_pyebm", return_value={"versions": {"test": "1"}}) as verify_environment, \
                        patch.object(R, "fit_job", side_effect=cheap_fit), \
                        patch.object(R, "STOP_REQUESTED", False), \
                        patch.object(R.sys, "argv", args), redirect_stdout(io.StringIO()):
                    R.main()
                    self.assertEqual(len(calls), 4)
                    R.main()
                    self.assertEqual(len(calls), 4)
                    progress = json.loads((output / "progress.json").read_text())
                    self.assertTrue(progress["done"])
                    row = json.loads((output / "rows.jsonl").read_text())
                    self.assertTrue(row["diagnostics"]["full_p_selected"])
                    self.assertEqual(row["schemes"]["dx"]["nperm"], 3)
                    self.assertTrue(row["schemes"]["dx"]["complete"])
                    log = output / "fit_records.jsonl"
                    saved_lines = log.read_bytes().splitlines(keepends=True)
                    lost = json.loads(saved_lines[-1])
                    log.write_bytes(b"".join(saved_lines[:-1]) + saved_lines[-1][:20])
                    R.main()
                    self.assertEqual(len(calls), 5)
                    self.assertEqual(calls[-1], R.key(lost))
                    self.assertEqual(len(R.load_events(log)), 4)
                    verify_environment.return_value = {"versions": {"test": "2"}}
                    with self.assertRaisesRegex(ValueError, "Numerical environment changed"):
                        R.main()
                    self.assertEqual(len(calls), 5)
                    verify_environment.return_value = {"versions": {"test": "1"}}
                    cfg["bperm"] = 4
                    cfg_path.write_text(json.dumps(cfg))
                    with self.assertRaisesRegex(ValueError, "Source/config changed"):
                        R.main()
            finally:
                for sig, handler in old_signals.items():
                    signal.signal(sig, handler)


class FitIntegrationTest(unittest.TestCase):
    def test_actual_debm_job_uses_v2_truth_and_biomarker_mapping(self):
        previous = signal.getsignal(signal.SIGALRM)
        cell = {"name": "PWR_E2_K6", "overrides": {"missing": False,
                "biomarker_names": ["ABETA", "PTAU", "TAU", "NG"]}}
        try:
            result = R.fit_job(({"fit_timeout_s": 60}, cell, 9702, "original", None, None))
        finally:
            signal.signal(signal.SIGALRM, previous)
        self.assertEqual(result["status"], "ok", result["diagnostics"])
        self.assertEqual(result["truth"]["pair_null"], [False, False, True])
        self.assertEqual(result["truth"]["pair_distances"], [1., 1., 0.])
        self.assertEqual(len(result["distance_to_truth"]), 3)
        for order in result["orderings"]:
            self.assertEqual(sorted(order), list(range(4)))
        json.dumps(result, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
