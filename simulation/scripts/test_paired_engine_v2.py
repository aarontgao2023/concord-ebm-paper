"""Same-input mixture replay equivalence and failure/restoration checks."""
from __future__ import annotations

from collections import namedtuple
from copy import deepcopy
import json
from time import monotonic
import unittest
from unittest.mock import patch
import warnings

import numpy as np
from scipy.special import expit

from engine_v2 import EngineConfig, fit_orderings
from paired_engine_v2 import _array_digest, _clone_mixture_output, fit_paired_orderings
from test_engine_v2 import small_dataset


def core_record(result):
    record = result.to_dict()
    record["diagnostics"] = {key: value for key, value in record["diagnostics"].items()
                             if key not in ("seconds", "paired_engine", "model_artifacts")}
    return record


def independent_reference(frame, config):
    """Capture raw output without bypassing or replacing any numerical work."""
    from pyebm import core_utilities as cu
    saved = {}
    original_mixture, original_compile = cu.do_mixturemodel, cu.compile_model_output

    def mixture(*args, **kwargs):
        output = original_mixture(*args, **kwargs)
        saved["p_yes"] = output[1].copy()
        return output

    def compile_model(*args, **kwargs):
        model = original_compile(*args, **kwargs)
        saved["event_centers"] = deepcopy(model.EventCenters)
        saved["parameters"] = [{field: deepcopy(getattr(parameters, field))
                                 for field in ("Control", "Disease", "Mixing")}
                                for parameters in model.BiomarkerParameters]
        return model

    with patch.object(cu, "do_mixturemodel", mixture), patch.object(cu, "compile_model_output", compile_model):
        result = fit_orderings(frame, config)
    return result, saved


def cheap_mixture(*args, **kwargs):
    """Only for error-path tests: skip costly GMM, retain the real DEBM pipeline."""
    data, group_values = args[3], args[5][0]
    events, groups = data.shape[1], len(np.unique(group_values))
    parameters = namedtuple("BiomarkerParams", "Control Disease Mixing")
    parameters.Control = [np.tile([0., 1.], (events, 1)) for _ in range(groups)]
    parameters.Disease = [np.tile([2., 1.], (events, 1)) for _ in range(groups)]
    parameters.Mixing = [np.full(events, .5) for _ in range(groups)]
    p_yes = expit(2 * data[:, :, 0] - 2)
    return parameters, p_yes, 1 - p_yes, p_yes.copy(), 1 - p_yes


class PairedRealEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases = {}
        standard = small_dataset()
        standard.loc[[1, 35, 59], "b1"] = np.nan
        permuted = standard.copy()
        permuted["APOE"] = np.random.default_rng(50126).permutation(permuted.APOE.to_numpy())
        for label, frame in (("standard", standard), ("permuted", permuted)):
            config = EngineConfig(expected_events=4, audit_original_neighbors=True,
                                  biomarker_names=("b3", "b2", "b1", "b0"))
            references = {}
            began = monotonic()
            for mode in ("original", "repaired"):
                options = dict(vars(config), mode=mode)
                references[mode] = independent_reference(frame, EngineConfig(**options))
            independent_seconds = monotonic() - began
            paired = fit_paired_orderings(frame, config)
            cls.cases[label] = (references, paired)
            timing = paired["original"].diagnostics["paired_engine"]
            print(f"paired timing {label}: independent={independent_seconds:.3f}s "
                  f"paired={timing['total_wall_seconds']:.3f}s mixture={timing['mixture_seconds']:.3f}s")

    def test_both_complete_fitresults_equal_two_independent_pipelines(self):
        for label, (references, paired) in self.cases.items():
            for mode in ("original", "repaired"):
                expected, _ = references[mode]
                self.assertTrue(expected.ok, (label, mode, expected.diagnostics))
                self.assertTrue(paired[mode].ok, (label, mode, paired[mode].diagnostics))
                self.assertEqual(core_record(paired[mode]), core_record(expected), (label, mode))

    def test_model_centers_parameters_and_objectives_are_exactly_equal(self):
        for label, (references, paired) in self.cases.items():
            for mode in ("original", "repaired"):
                reference, raw = references[mode]
                details = paired[mode].diagnostics["model_artifacts"]
                self.assertEqual(details["posterior_sha256"], _array_digest(raw["p_yes"]))
                np.testing.assert_array_equal(details["event_centers"], raw["event_centers"])
                for actual, expected in zip(details["biomarker_parameters"], raw["parameters"]):
                    for field in ("Control", "Disease", "Mixing"):
                        np.testing.assert_array_equal(actual[field], expected[field])
                self.assertEqual(paired[mode].diagnostics["consensus"], reference.diagnostics["consensus"])
            # Two independently fitted pipelines have identical mixture posteriors.
            np.testing.assert_array_equal(references["original"][1]["p_yes"],
                                          references["repaired"][1]["p_yes"])
            self.assertEqual(paired["original"].diagnostics["model_artifacts"]["posterior_sha256"],
                             paired["repaired"].diagnostics["model_artifacts"]["posterior_sha256"])

    def test_shared_cost_is_explicit_and_each_record_is_serializable(self):
        for _, paired in self.cases.values():
            for mode, result in paired.items():
                d = result.diagnostics["paired_engine"]
                self.assertEqual((d["mixture_compute_calls"], d["mixture_replay_calls"]), (1, 1))
                self.assertEqual(d["optimizer_records_replayed"], mode == "repaired")
                self.assertGreater(d["total_wall_seconds"], d["mixture_seconds"])
                json.dumps(result.to_dict(), allow_nan=False)


class PairedFailureTests(unittest.TestCase):
    def setUp(self):
        from pyebm import core_utilities as cu, debm
        from pyebm.central_ordering.generalized_mallows import weighted_mallows
        from pyebm.mixture_model import gaussian_mixture_model as gmm
        self.cu, self.debm, self.mallows, self.gmm = cu, debm, weighted_mallows, gmm
        self.frame = small_dataset()
        self.config = EngineConfig(expected_events=4)

    def state(self):
        return (self.cu.do_mixturemodel, self.debm.fit,
                self.mallows.__dict__["consensus"], self.gmm.opt)

    def assert_restored(self, before):
        for actual, expected in zip(self.state(), before):
            self.assertIs(actual, expected)

    def test_dynamic_parameter_class_outputs_are_really_independent_copies(self):
        parameters = namedtuple("BiomarkerParams", "Control Disease Mixing")
        parameters.Control = [np.array([[0., 1.]])]
        parameters.Disease = [np.array([[2., 1.]])]
        parameters.Mixing = [np.array([.4])]
        output = (parameters, np.array([[.2], [.8]]), np.array([[.8], [.2]]), [], [])
        a, b = _clone_mixture_output(output), _clone_mixture_output(output)
        self.assertIsNot(a[0], parameters)
        a[0].Control[0][0, 0] = 999.
        a[0].Mixing[0][0] = .99
        a[1][0, 0] = .99
        self.assertEqual(b[0].Control[0][0, 0], 0.)
        self.assertEqual(b[0].Mixing[0][0], .4)
        self.assertEqual(b[1][0, 0], .2)

    def test_each_call_and_changed_permutation_compute_their_own_mixture(self):
        calls = []
        def recorded(*args, **kwargs):
            calls.append(np.asarray(args[5][0]).copy())
            return cheap_mixture(*args, **kwargs)
        with patch.object(self.cu, "do_mixturemodel", recorded):
            before = self.state()
            a = fit_paired_orderings(self.frame, self.config)
            permuted = self.frame.copy()
            permuted["APOE"] = np.random.default_rng(776).permutation(permuted.APOE.to_numpy())
            b = fit_paired_orderings(permuted, self.config)
            self.assert_restored(before)
        self.assertEqual(len(calls), 2)
        self.assertFalse(np.array_equal(calls[0], calls[1]))
        self.assertTrue(all(result.ok for pair in (a, b) for result in pair.values()))

    def test_changed_inputs_between_passes_refuse_replay(self):
        original_fit, calls = self.debm.fit, []
        def mutating_fit(*args, **kwargs):
            output = original_fit(*args, **kwargs)
            calls.append(1)
            if len(calls) == 1:
                self.frame.loc[[0, 30], "APOE"] = self.frame.loc[[30, 0], "APOE"].to_numpy()
            return output
        with patch.object(self.cu, "do_mixturemodel", cheap_mixture), patch.object(self.debm, "fit", mutating_fit):
            before = self.state()
            results = fit_paired_orderings(self.frame, self.config)
            self.assert_restored(before)
        self.assertTrue(results["original"].ok)
        self.assertEqual(results["repaired"].status, "error")
        self.assertIn("Mixture inputs differ", results["repaired"].diagnostics["error"])
        self.assertEqual(results["repaired"].diagnostics["paired_engine"]["mixture_replay_calls"], 0)

    def test_mixture_errors_optimizer_records_and_warnings_are_replayed(self):
        calls = []
        def failed(*args, **kwargs):
            calls.append(1)
            self.gmm.opt.records.append({"success": False, "status": 9, "method": "test"})
            warnings.warn("known mixture warning", RuntimeWarning)
            raise RuntimeError("known mixture failure")
        with patch.object(self.cu, "do_mixturemodel", failed):
            before = self.state()
            results = fit_paired_orderings(self.frame, self.config)
            self.assert_restored(before)
        self.assertEqual(len(calls), 1)
        for result in results.values():
            self.assertEqual(result.status, "error")
            self.assertEqual(result.diagnostics["error"], "RuntimeError: known mixture failure")
            self.assertEqual(result.diagnostics["optimizer_calls"], 1)
            self.assertEqual(result.diagnostics["optimizer_failures"], 1)
            self.assertIn("optimizer_nonconvergence", result.diagnostics["quality_flags"])
            self.assertEqual(result.diagnostics["warnings_count"], 1)
            self.assertEqual(result.diagnostics["warnings"], ["RuntimeWarning: known mixture warning"])

    def test_original_consensus_failure_does_not_block_repaired_consensus(self):
        def failed_consensus(*args, **kwargs):
            raise RuntimeError("original consensus failed")
        with patch.object(self.cu, "do_mixturemodel", cheap_mixture), patch.object(
                self.mallows, "consensus", classmethod(failed_consensus)):
            before = self.state()
            results = fit_paired_orderings(self.frame, self.config)
            self.assert_restored(before)
        self.assertEqual(results["original"].status, "error")
        self.assertTrue(results["repaired"].ok, results["repaired"].diagnostics)

    def test_timeout_before_mixture_completion_aborts_and_restores(self):
        def timed_out(*args, **kwargs):
            raise TimeoutError("joint budget exhausted")
        with patch.object(self.cu, "do_mixturemodel", timed_out):
            before = self.state()
            with self.assertRaisesRegex(TimeoutError, "joint budget exhausted") as caught:
                fit_paired_orderings(self.frame, self.config)
            self.assert_restored(before)
        self.assertEqual(caught.exception.completed_results, {})

    def test_second_pass_custom_timeout_keeps_complete_original_and_restores(self):
        class CustomTimeout(Exception):
            pass
        stage_calls = []
        original_staging = self.cu.do_patient_staging
        def staging(*args, **kwargs):
            stage_calls.append(1)
            if len(stage_calls) == 2:
                raise CustomTimeout("second pass timeout")
            return original_staging(*args, **kwargs)
        with patch.object(self.cu, "do_mixturemodel", cheap_mixture), patch.object(
                self.cu, "do_patient_staging", staging):
            before = self.state()
            with self.assertRaisesRegex(CustomTimeout, "second pass timeout") as caught:
                fit_paired_orderings(self.frame, self.config, abort_exceptions=(CustomTimeout,))
            self.assert_restored(before)
        complete = caught.exception.completed_results
        self.assertEqual(set(complete), {"original"})
        self.assertTrue(complete["original"].ok)
        self.assertIn("paired_engine", complete["original"].diagnostics)

    def test_keyboard_interrupt_restores_all_patches(self):
        def interrupted(*args, **kwargs):
            raise KeyboardInterrupt()
        with patch.object(self.cu, "do_mixturemodel", interrupted):
            before = self.state()
            with self.assertRaises(KeyboardInterrupt):
                fit_paired_orderings(self.frame, self.config)
            self.assert_restored(before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
