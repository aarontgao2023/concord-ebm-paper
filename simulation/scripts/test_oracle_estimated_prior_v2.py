"""Analytic likelihood, identification and no-truth-order oracle checks."""
from __future__ import annotations

from copy import deepcopy
import json
import unittest

import numpy as np
from scipy.stats import norm

from design_v2 import resolve, simulate
from engine_v2 import EngineConfig
from oracle_estimated_prior_v2 import (
    Q_LOWER, Q_UPPER, estimate_abnormal_prior,
    fit_oracle_estimated_prior_orderings, mixture_loglik_derivatives,
)
from test_oracle_engine_v2 import simple_fixture


COMPONENTS = {"normal_mean": 0., "abnormal_mean": 2., "sd": 1.}


class EstimatedPriorTests(unittest.TestCase):
    def test_derivatives_match_direct_densities_and_finite_differences(self):
        values = np.array([-1.8, -.2, .5, 1.1, 2.3, 3.9])
        f0, f1 = norm.pdf(values), norm.pdf(values, loc=2.)
        logratios = norm.logpdf(values, loc=2.) - norm.logpdf(values)
        q, h = .37, 1e-5
        ll, score, curvature = mixture_loglik_derivatives(q, logratios)
        mixture = (1 - q) * f0 + q * f1
        self.assertAlmostEqual(ll + np.log(f0).sum(), np.log(mixture).sum(), places=12)
        self.assertAlmostEqual(score, np.sum((f1 - f0) / mixture), places=12)
        self.assertAlmostEqual(curvature, -np.sum(((f1 - f0) / mixture)**2), places=12)
        before = mixture_loglik_derivatives(q - h, logratios)
        after = mixture_loglik_derivatives(q + h, logratios)
        self.assertAlmostEqual(score, (after[0] - before[0]) / (2 * h), places=7)
        self.assertAlmostEqual(curvature, (after[1] - before[1]) / (2 * h), places=6)

    def test_global_concavity_and_grid_maximum_in_several_samples(self):
        rng = np.random.default_rng(99218)
        for q_true, delta in ((.08, .7), (.4, 2.), (.85, 4.)):
            values = rng.normal(size=100) + delta * (rng.uniform(size=100) < q_true)
            component = dict(COMPONENTS, abnormal_mean=delta)
            fit = estimate_abnormal_prior(values, component)
            ratios = delta * (values - delta / 2)
            grid = np.linspace(Q_LOWER, Q_UPPER, 1501)
            triples = np.asarray([mixture_loglik_derivatives(q, ratios) for q in grid])
            self.assertTrue(np.all(triples[:, 2] < 0))
            self.assertTrue(np.all(np.diff(triples[:, 1]) <= 0))
            self.assertGreaterEqual(fit.relative_loglik_at_estimate + 1e-11, triples[:, 0].max())
            self.assertTrue(fit.identified)
            if fit.status == "interior":
                self.assertAlmostEqual(fit.score_at_estimate, 0., places=8)
                self.assertGreater(fit.score_lower, 0.)
                self.assertLess(fit.score_upper, 0.)

    def test_extreme_ratios_and_symmetric_likelihood_remain_finite(self):
        # f0/f1 underflow in direct PDF space, but the constrained mixture is finite.
        values = np.array([-1000., -1000., -1000., 1000.])
        fit = estimate_abnormal_prior(values, COMPONENTS)
        self.assertEqual(fit.status, "interior")
        self.assertAlmostEqual(fit.abnormal_prior, .25, places=12)
        self.assertTrue(np.isfinite(fit.relative_loglik_at_estimate))
        symmetric = estimate_abnormal_prior(np.array([-4., 0., 2., 6.]), COMPONENTS)
        self.assertAlmostEqual(symmetric.abnormal_prior, .5, places=12)
        _, _, curvature = mixture_loglik_derivatives(.5, np.array([-1e-9, 1e-9]))
        self.assertLess(curvature, 0.)

    def test_boundary_empty_identical_and_flat_cases_are_explicit(self):
        lower = estimate_abnormal_prior(np.array([-4., -3., -2.]), COMPONENTS)
        upper = estimate_abnormal_prior(np.array([4., 5., 6.]), COMPONENTS)
        self.assertEqual((lower.abnormal_prior, lower.status), (Q_LOWER, "lower_boundary"))
        self.assertEqual((upper.abnormal_prior, upper.status), (Q_UPPER, "upper_boundary"))
        self.assertLess(lower.score_at_estimate, 0.)
        self.assertGreater(upper.score_at_estimate, 0.)
        cases = [
            (np.array([np.nan, np.nan]), COMPONENTS, "no_observations"),
            (np.array([1., 2.]), dict(COMPONENTS, abnormal_mean=0.), "identical_components"),
            (np.array([1., 1.]), COMPONENTS, "flat_observed_likelihood"),
        ]
        for values, params, expected in cases:
            fit = estimate_abnormal_prior(values, params)
            self.assertEqual(fit.status, expected)
            self.assertEqual(fit.abnormal_prior, .5)
            self.assertFalse(fit.identified)
            json.dumps(fit.to_dict(), allow_nan=False)

    def test_missing_values_do_not_add_likelihood_or_pseudo_observations(self):
        values = np.array([-1., .2, .4, 1.2, 2., 3.4])
        a = estimate_abnormal_prior(values, COMPONENTS)
        b = estimate_abnormal_prior(np.r_[values, np.full(20, np.nan)], COMPONENTS)
        self.assertEqual(a.abnormal_prior, b.abnormal_prior)
        self.assertEqual(a.relative_loglik_at_estimate, b.relative_loglik_at_estimate)
        self.assertEqual((b.n_observed, b.n_missing), (6, 20))
        with self.assertRaises(ValueError):
            estimate_abnormal_prior([np.inf], COMPONENTS)
        with self.assertRaises(ValueError):
            estimate_abnormal_prior([0.], dict(COMPONENTS, sd=0.))
        with self.assertRaises(ValueError):
            mixture_loglik_derivatives(0., np.array([1.]))

    def test_large_known_mixture_recovers_abnormal_fraction(self):
        rng = np.random.default_rng(86152)
        for true_q in (.1, .7):
            component = dict(COMPONENTS, abnormal_mean=3.)
            values = rng.normal(size=40000) + 3. * (rng.uniform(size=40000) < true_q)
            fit = estimate_abnormal_prior(values, component)
            self.assertLess(abs(fit.abnormal_prior - true_q), .01)

    def test_real_mallows_uses_estimated_q_and_no_diagnosis_or_order_truth(self):
        frame, truth = simple_fixture()
        # The wrapper must work without all forbidden truth fields or Diagnosis.
        minimal_truth = {key: truth[key] for key in ("biomarker_names", "group_order", "component_parameters")}
        frame = frame.drop(columns="Diagnosis")
        for mode in ("original", "repaired"):
            fit = fit_oracle_estimated_prior_orderings(frame, minimal_truth,
                                                       EngineConfig(mode=mode, expected_events=4))
            self.assertTrue(fit.ok, fit.diagnostics)
            np.testing.assert_array_equal(fit.orderings, [[3, 1, 0, 2]] * 3)
            d = fit.diagnostics
            self.assertFalse(d["inference_eligible"])
            self.assertFalse(d["truth_order_used_in_prior_and_initialization"])
            self.assertFalse(d["individual_latent_stage_used"])
            self.assertFalse(d["diagnosis_used_in_prior_estimation"])
            self.assertFalse(d["truth_dx_counts_used"])
            self.assertEqual(d["prior_fit_calls"], 12)
            self.assertEqual(d["optimizer_calls"], 0)
            for code, group in enumerate(truth["group_order"]):
                priors = np.asarray(d["abnormal_priors"][group])
                np.testing.assert_allclose(d["normal_mixing_params"][group], 1 - priors)
                np.testing.assert_array_equal(d["initial_orderings"][group], np.argsort(1 - priors))
                for index, name in enumerate(truth["biomarker_names"]):
                    direct = estimate_abnormal_prior(frame.loc[frame.APOE == code, name],
                                                     truth["component_parameters"][group][name])
                    self.assertAlmostEqual(priors[index], direct.abnormal_prior, places=14)
            if mode == "repaired":
                self.assertTrue(all(x["local_optimal"] for x in d["consensus"]))
            json.dumps(fit.to_dict(), allow_nan=False)

    def test_truth_poison_and_input_name_remapping_leave_biomarker_results_unchanged(self):
        frame, truth = simple_fixture()
        baseline = fit_oracle_estimated_prior_orderings(frame, truth)
        altered = deepcopy(truth)
        for forbidden in ("group_orderings", "latent_stage", "group_dx_counts", "stage_beta_parameters"):
            altered[forbidden] = object()
        frame["Diagnosis"] = "UNUSABLE_DIAGNOSIS"
        expected_names = tuple(reversed(truth["biomarker_names"]))
        changed = fit_oracle_estimated_prior_orderings(
            frame, altered, EngineConfig(mode="repaired", biomarker_names=expected_names))
        self.assertTrue(changed.ok, changed.diagnostics)
        self.assertEqual(changed.diagnostics["abnormal_priors"], baseline.diagnostics["abnormal_priors"])
        for order in changed.orderings:
            self.assertEqual([expected_names[i] for i in order], ["d", "b", "a", "c"])

    def test_all_missing_component_remains_flagged_in_fit_result(self):
        frame, truth = simple_fixture()
        frame.loc[frame.APOE == 0, "a"] = np.nan
        fit = fit_oracle_estimated_prior_orderings(frame, truth)
        record = fit.diagnostics["prior_fits"]["e2"]["a"]
        self.assertEqual(record["status"], "no_observations")
        self.assertEqual(record["n_observed"], 0)
        self.assertIn("prior_no_observations:e2:a", fit.diagnostics["quality_flags"])
        # The finite ordering, if returned, must retain the identification flag.
        self.assertGreaterEqual(fit.diagnostics["prior_unidentified_fits"], 1)

    def test_design_v2_truth_is_compatible_on_four_biomarkers(self):
        # Use actual declared names rather than requiring a particular historical alias.
        from design_v2 import DesignConfig
        config = resolve("STAGE_H0", biomarker_names=DesignConfig().biomarker_names[:4])
        frame, truth = simulate(config, 51000821)
        result = fit_oracle_estimated_prior_orderings(frame, truth,
                                                     EngineConfig(mode="repaired", expected_events=4))
        self.assertTrue(result.ok, result.diagnostics)
        self.assertEqual(len(result.orderings), 3)
        self.assertEqual(result.diagnostics["prior_fit_calls"], 12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
