"""Cheap analytic-prior and oracle mechanism checks, independent of large fits."""
from __future__ import annotations

from copy import deepcopy
import unittest

import numpy as np
import pandas as pd

from design_v2 import resolve, simulate
from engine_v2 import EngineConfig
from oracle_engine_v2 import analytic_group_priors, fit_oracle_orderings, oracle_posteriors


def simple_fixture():
    names, order = ["a", "b", "c", "d"], [3, 1, 0, 2]
    rows = []
    for group in range(3):
        for stage in range(5):
            for _ in range(4):
                rows.append({"APOE": group, "PTID": str(len(rows)), "Diagnosis": "MCI",
                             **{name: 8.0 * (order.index(i) < stage) for i, name in enumerate(names)}})
    frame = pd.DataFrame(rows)
    frame.loc[10, "b"] = np.nan
    groups = ["e2", "e33", "e4"]
    truth = {"biomarker_names": names, "group_order": groups,
             "group_orderings": {group: order for group in groups},
             "group_dx_counts": {group: {"MCI": 20} for group in groups},
             "stage_beta_parameters": {group: {"MCI": [1.0, 1.0]} for group in groups},
             "component_parameters": {group: {name: {"normal_mean": 0.0, "abnormal_mean": 8.0, "sd": 1.0}
                                                for name in names} for group in groups}}
    return frame, truth


class OracleTests(unittest.TestCase):
    def test_uniform_stage_closed_form_and_rank_mapping(self):
        _, truth = simple_fixture()
        q = analytic_group_priors(truth)
        ranks = np.argsort(truth["group_orderings"]["e2"])
        np.testing.assert_allclose(q["e2"], 1 - (ranks + 0.5) / 4)
        np.testing.assert_array_equal(np.argsort(1 - q["e2"]), [3, 1, 0, 2])

    def test_beta_formula_matches_large_cheap_stage_simulation(self):
        _, truth = simulate(resolve("STAGE_H0"), 9617)
        q = analytic_group_priors(truth)
        rng = np.random.default_rng(272)
        draws = 150000
        events = len(truth["biomarker_names"])
        for group in truth["group_order"]:
            dxs = list(truth["group_dx_counts"][group])
            counts = np.array([truth["group_dx_counts"][group][dx] for dx in dxs])
            sampled_dx = rng.choice(len(dxs), draws, p=counts / counts.sum())
            stage = np.empty(draws, dtype=int)
            for index, diagnosis in enumerate(dxs):
                mask = sampled_dx == index
                a, b = truth["stage_beta_parameters"][group][diagnosis]
                stage[mask] = np.rint(events * rng.beta(a, b, mask.sum())).astype(int)
            ranks = np.argsort(truth["group_orderings"][group])
            estimated = (stage[:, None] > ranks[None, :]).mean(axis=0)
            np.testing.assert_allclose(estimated, q[group], atol=0.008, rtol=0)

    def test_component_posterior_uses_abnormal_q_not_normal_mixing(self):
        from pyebm.mixture_model.gaussian_mixture_model import calculate_prob_mm
        values = np.array([[-1.0], [0.0], [2.0], [4.0], [np.nan]])
        params = {"a": {"normal_mean": 0.0, "abnormal_mean": 2.0, "sd": 1.0}}
        posterior = oracle_posteriors(values, params, ["a"], np.array([0.25]))
        upstream, *_ = calculate_prob_mm(values, np.array([[0.0, 1.0, 2.0, 1.0, 0.75]]),
                                        val_invalid=np.nan)
        np.testing.assert_allclose(posterior, upstream, equal_nan=True)
        self.assertTrue(np.isnan(posterior[-1, 0]))

    def test_simple_order_recovery_local_optimum_and_no_individual_stage_use(self):
        frame, truth = simple_fixture()
        config = EngineConfig(mode="repaired", expected_events=4)
        for mode in ("equal_prior", "known_group_prior"):
            result = fit_oracle_orderings(frame, truth, mode, config)
            self.assertTrue(result.ok, result.diagnostics)
            np.testing.assert_array_equal(result.orderings, [[3, 1, 0, 2]] * 3)
            self.assertFalse(result.diagnostics["inference_eligible"])
            self.assertEqual(result.diagnostics["optimizer_calls"], 0)
            self.assertTrue(all(record["local_optimal"] for record in result.diagnostics["consensus"]))
            altered = deepcopy(truth)
            altered["latent_stage"] = [None] * len(frame)
            repeated = fit_oracle_orderings(frame, altered, mode, config)
            np.testing.assert_array_equal(repeated.orderings, result.orderings)

    def test_output_names_and_original_mode(self):
        frame, truth = simple_fixture()
        expected_names = tuple(reversed(truth["biomarker_names"]))
        config = EngineConfig(mode="original", biomarker_names=expected_names)
        result = fit_oracle_orderings(frame, truth, "known_group_prior", config)
        self.assertTrue(result.ok, result.diagnostics)
        for ordering in result.orderings:
            self.assertEqual([expected_names[i] for i in ordering], ["d", "b", "a", "c"])
        self.assertEqual([record["mode"] for record in result.diagnostics["consensus"]], ["original"] * 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
