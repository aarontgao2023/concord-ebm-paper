"""Exact objective/trajectory and context-restoration checks for the fast normal-density likelihood."""
from __future__ import annotations

import unittest
import warnings

import numpy as np
import scipy.stats

from benchmark_fast_likelihood_v2 import assert_exact, full_fit_equivalence
from engine_v2 import EngineConfig
from fast_likelihood_v2 import calculate_likelihood_gmm, fast_likelihood_context
from paired_engine_v2 import fit_paired_orderings
from test_engine_v2 import small_dataset


class FastLikelihoodTests(unittest.TestCase):
    def test_pdf_arrays_are_exact_for_random_missing_and_boundary_parameters(self):
        rng = np.random.default_rng(491781)
        data = np.r_[rng.normal(size=100), np.nan, -np.inf, np.inf, -1e6, 1e6]
        parameters = [(rng.normal(), np.exp(rng.uniform(-12, 5))) for _ in range(60)]
        parameters += [(0., 0.), (0., -1.), (np.nan, 1.), (0., np.nan), (0., np.inf), (np.inf, 1.)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for loc, scale in parameters:
                a = scipy.stats.norm(loc=loc, scale=scale).pdf(data)
                b = scipy.stats.norm.pdf(data, loc=loc, scale=scale)
                self.assertTrue(np.array_equal(a, b, equal_nan=True), (loc, scale))

    def test_objective_exact_for_both_branches_and_all_missing_groups(self):
        from pyebm.mixture_model import gaussian_mixture_model as gmm
        rng = np.random.default_rng(813397)
        data = rng.normal(size=90)
        groups = np.repeat([0., 1., 2.], 30)
        data[:30] = np.nan
        data[34] = np.nan
        groups[-1] = np.nan
        params = [np.array([rng.normal(), np.exp(rng.uniform(-5, 2)), rng.normal()+2,
                            np.exp(rng.uniform(-5, 2)), rng.uniform()]) for _ in range(60)]
        params += [np.array([0., s, 2., s, q]) for s in (0., 1e-8, 1.) for q in (0., 1e-6, 1-1e-6, 1.)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for values in (data, np.full(90, np.nan)):
                for mixing in ([], [0., .5, 1.]):
                    args = (values, ["APOE"], [groups], mixing)
                    for param in params:
                        assert_exact(gmm.calculate_likelihood_gmm(param, *args),
                                     calculate_likelihood_gmm(param, *args))

    def test_disabled_nested_and_exception_contexts_restore_identity(self):
        from pyebm.mixture_model import gaussian_mixture_model as gmm
        original = gmm.calculate_likelihood_gmm
        with fast_likelihood_context(False):
            self.assertIs(gmm.calculate_likelihood_gmm, original)
        with self.assertRaisesRegex(TimeoutError, "test timeout"):
            with fast_likelihood_context():
                self.assertIs(gmm.calculate_likelihood_gmm, calculate_likelihood_gmm)
                self.assertEqual(gmm.calculate_likelihood_gmm.__name__, original.__name__)
                with fast_likelihood_context():
                    self.assertIs(gmm.calculate_likelihood_gmm, calculate_likelihood_gmm)
                self.assertIs(gmm.calculate_likelihood_gmm, calculate_likelihood_gmm)
                raise TimeoutError("test timeout")
        self.assertIs(gmm.calculate_likelihood_gmm, original)

    def test_full_original_fit_matches_every_optimizer_evaluation_and_output(self):
        frame = small_dataset()
        frame.loc[[1, 35, 59], "b1"] = np.nan
        result = full_fit_equivalence(frame, EngineConfig(expected_events=4, audit_original_neighbors=True))
        self.assertTrue(result["all_objective_evaluations_exact_equal"])
        self.assertGreater(result["objective_evaluations"], 0)

    def test_fast_context_nests_with_paired_engine_without_leaking(self):
        from pyebm.mixture_model import gaussian_mixture_model as gmm
        original_objective, original_opt = gmm.calculate_likelihood_gmm, gmm.opt
        frame = small_dataset()
        config = EngineConfig(expected_events=4)
        baseline = fit_paired_orderings(frame, config)
        with fast_likelihood_context():
            fast = fit_paired_orderings(frame, config)
        self.assertIs(gmm.calculate_likelihood_gmm, original_objective)
        self.assertIs(gmm.opt, original_opt)
        for mode in baseline:
            expected, actual = baseline[mode].to_dict(), fast[mode].to_dict()
            for record in (expected, actual):
                record["diagnostics"] = {key: value for key, value in record["diagnostics"].items()
                                          if key not in ("seconds", "paired_engine")}
            assert_exact(expected, actual)


if __name__ == "__main__":
    unittest.main(verbosity=2)
