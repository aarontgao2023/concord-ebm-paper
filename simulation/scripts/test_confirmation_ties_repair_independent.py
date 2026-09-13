"""Independent invented count fixtures; no generator, fitted data or HPC.

These fixtures do not reuse the author's test builders. The pure summary API
assumes provenance-validated rows; production-loader integration is separately
covered by the author and inspected in this review.
"""
import copy
import hashlib
import math
import unittest

import analyze_confirmation_ties_repair as T


BASE = 31008600


def plan(n=4, selected=2, explicit=()):
    return {"run_id": "confirm_core_a", "phase": "development", "base_seed": BASE,
            "datasets": n, "full_p_first_n": selected, "full_p_seeds": list(explicit),
            "engines": ["original", "repaired"], "cells": [{"name": "REF_H0"}],
            "paired_standard": True, "fast_likelihood": True, "bperm": 599,
            "alpha": .05, "rule": "le", "schemes": {
                "unrestricted": {"stratify": "none", "require_max": True},
                "diagnosis": {"stratify": "diagnosis", "require_max": True}}}


def observation(cfg, offset, engine, decisions=(False, False), ties=None):
    seed = cfg["base_seed"] + offset
    selected = set(range(cfg["base_seed"], cfg["base_seed"] + cfg["full_p_first_n"])) | set(cfg["full_p_seeds"])
    schemes = {}
    for index, (name, decision) in enumerate(zip(("unrestricted", "diagnosis"), decisions)):
        # With one failed planned permutation and nine inclusive exceedances,
        # the pair <= test is unknown. Its global max test can still reject.
        n = 598 if decision is None else 599
        exceed = 9 if decision is None else 0 if decision else 40
        eq = 0 if ties is None else ties[index]
        gt = exceed if ties is None else 0
        schemes[name] = {"requested_nperm": 599, "nperm": n, "attempted": 599,
            "failed": 599 - n, "complete": n == 599, "gt": [gt] * 3, "eq": [eq] * 3,
            "maxgt": gt, "maxeq": eq, "definition": copy.deepcopy(cfg["schemes"][name])}
    return {"run_id": cfg["run_id"], "phase": cfg["phase"], "cell": "REF_H0",
            "seed": seed, "engine": engine, "status": "ok", "schemes": schemes,
            "truth": {"pair_null": [True] * 3,
                "manifest": {"seed": seed, "data_sha256": hashlib.sha256(str(seed).encode()).hexdigest()},
                "ordering_description": "same three-group order"},
            "diagnostics": {"full_p_selected": seed in selected}}


def report(cfg, rows):
    return T.summarize_run(cfg, rows, bootstrap_repeats=25, bootstrap_seed=31337)


def tie(result, engine="original", scheme="unrestricted", metric="within_dataset_mean_pair"):
    return next(x for x in result["ties_by_scheme"] if x["engine"] == engine and x["scheme"] == scheme)["metrics"][metric]


def contrast(result, comparison="U_minus_D_gap", endpoint="bonferroni"):
    return next(x for x in result["repair_rejection_contrasts"] if x["comparison"] == comparison and x["endpoint"] == endpoint)["summary"]


class IndependentTiesRepairTests(unittest.TestCase):
    def test_production_contract_pins_absolute_cohort_and_selection(self):
        # Metadata only: these identities are validated, never generated.
        cfg = plan(1000, 100)
        cfg.update(base_seed=42100000, phase="confirmation", cells=[{"name": "IID_H0"}, {"name": "REF_H0"}])
        cells, seeds, chosen = T._validate_config(cfg, frozen=True)
        self.assertEqual(cells, ["IID_H0", "REF_H0"])
        self.assertEqual((seeds[0], seeds[-1], len(seeds)), (42100000, 42100999, 1000))
        self.assertEqual(chosen, list(range(42100000, 42100100)))
        variants = [{"datasets": 999}, {"full_p_first_n": 99}, {"base_seed": 41100000},
                    {"full_p_first_n": 0, "full_p_seeds": list(range(42100100, 42100200))},
                    {"cells": list(reversed(cfg["cells"]))}]
        for change in variants:
            with self.subTest(change=change), self.assertRaises(ValueError):
                T._validate_config(cfg | change, frozen=True)

    def test_100_selected_and_1000_planned_do_not_become_999_completers(self):
        cfg = plan(1000, 100)
        rows = [observation(cfg, i, engine, ties=(599 if i >= 100 or i % 2 else 0, 0))
                for i in range(999) for engine in T.ENGINES]
        result = report(cfg, rows)
        selected = tie(result)
        self.assertEqual((selected["planned_datasets"], selected["known_datasets"]), (100, 100))
        self.assertEqual(selected["estimate"], .5)
        self.assertAlmostEqual(selected["mcse_mean"], math.sqrt(.25 / 99))
        self.assertEqual(len(result["selected_dataset_ties"]), 400)
        self.assertEqual(result["selected_absolute_seeds"], list(range(BASE, BASE + 100)))
        whole = contrast(result)
        self.assertEqual((whole["planned_datasets"], whole["known_datasets"], whole["unknown_datasets"]), (1000, 999, 1))
        self.assertEqual(whole["all_planned_mean_bounds"], [-.002, .002])
        for key in ("estimate", "mcse_mean", "bootstrap_ci95", "hoeffding_ci95"):
            self.assertIsNone(whole[key])

    def test_hand_worked_DID_and_missing_fifth_unit(self):
        cfg = plan(4, 2)
        patterns = [(1, 0, 0, 1), (0, 1, 0, 0), (1, 1, 0, 0), (0, 0, 1, 1)]
        rows = [observation(cfg, i, engine, tuple(bool(v) for v in pattern[2*j:2*j+2]))
                for i, pattern in enumerate(patterns) for j, engine in enumerate(T.ENGINES)]
        result = report(cfg, rows)
        # Matched dataset differences are +2,-1,0,0, with sample variance19/12.
        output = contrast(result)
        self.assertEqual(output["support"], [-2, 2])
        self.assertEqual(output["estimate"], .25)
        self.assertAlmostEqual(output["mcse_mean"], math.sqrt(19 / 48))
        self.assertEqual(contrast(result, "unrestricted")["estimate"], .25)
        self.assertEqual(contrast(result, "diagnosis")["estimate"], 0.)
        self.assertEqual(result, report(cfg, list(reversed(rows))))
        incomplete = report(plan(5, 2), rows)
        self.assertEqual(contrast(incomplete)["all_planned_mean_bounds"], [-.2, .6])
        self.assertIsNone(contrast(incomplete)["estimate"])

    def test_nonprefix_explicit_fullp_set_and_missing_selected_engine(self):
        cfg = plan(6, 0, [BASE + 1, BASE + 4])
        rows = [observation(cfg, i, engine, ties=(599 if i != 4 else 0, 0))
                for i in range(6) for engine in T.ENGINES if not (i == 4 and engine == "original")]
        result = report(cfg, rows)
        self.assertEqual(result["selected_absolute_seeds"], [BASE + 1, BASE + 4])
        self.assertEqual(tie(result)["all_planned_mean_bounds"], [.5, 1.])
        self.assertIsNone(tie(result)["estimate"])
        self.assertIsNone(tie(result)["distribution"])
        self.assertEqual(tie(result, "repaired")["estimate"], .5)
        self.assertEqual(tie(result)["unknown_seeds"], [BASE + 4])

    def test_terminal_failed_engine_with_stale_positive_fields_stays_unknown(self):
        cfg = plan(1, 1)
        original = observation(cfg, 0, "original", (True, False))
        failed = observation(cfg, 0, "repaired", (True, True))
        failed.update(status="invalid", truth={})
        result = report(cfg, [failed, original])
        self.assertEqual(contrast(result)["all_planned_mean_bounds"], [0., 2.])
        self.assertEqual(contrast(result, "unrestricted")["all_planned_mean_bounds"], [0., 1.])
        self.assertEqual(contrast(result, "diagnosis")["all_planned_mean_bounds"], [-1., 0.])
        self.assertEqual(result["status_or_unavailable_reason_counts"]["REF_H0:repaired"], {"invalid": 1})
        self.assertEqual(tie(result, "repaired")["unknown_datasets"], 1)

    def test_same_data_hash_but_contradictory_truth_masks_both_engines(self):
        cfg = plan(1, 1)
        original = observation(cfg, 0, "original", (True, False))
        repaired = observation(cfg, 0, "repaired", (False, True))
        repaired["truth"]["ordering_description"] = "different order declaration"
        before = copy.deepcopy((cfg, original, repaired))
        result = report(cfg, [original, repaired])
        self.assertEqual(contrast(result)["all_planned_mean_bounds"], [-2., 2.])
        self.assertEqual(tie(result)["unknown_datasets"], 1)
        self.assertEqual(tie(result, "repaired")["unknown_datasets"], 1)
        self.assertEqual(len(result["pairing_integrity_issues"]), 1)
        self.assertEqual((cfg, original, repaired), before)

    def test_598_successes_keep_unknown_pair_decision_and_no_tie_fraction(self):
        cfg = plan(1, 1)
        rows = [observation(cfg, 0, "original", (True, False)),
                observation(cfg, 0, "repaired", (None, None))]
        result = report(cfg, rows)
        self.assertEqual(contrast(result)["all_planned_mean_bounds"], [0., 2.])
        self.assertIsNone(contrast(result)["estimate"])
        # Both repaired global tests are resolved positive despite unknown
        # pair decisions; global DID is independently 1-0-1+1=1.
        self.assertEqual(contrast(result, endpoint="global_max")["estimate"], 1.)
        self.assertEqual(tie(result, "repaired")["known_datasets"], 0)
        self.assertFalse(result["endpoint_complete"])

    def test_degenerate_1000_differences_have_nonzero_range_correct_CI(self):
        cfg = plan(1000, 100)
        rows = [observation(cfg, i, engine, ties=(599, 599)) for i in range(1000) for engine in T.ENGINES]
        result = report(cfg, rows)
        output = contrast(result)
        radius = 4 * math.sqrt(math.log(40) / 2000)
        self.assertTrue(output["endpoint_complete"])
        self.assertTrue(output["bootstrap_degenerate_sample"])
        self.assertIsNone(output["bootstrap_ci95"])
        self.assertEqual(output["mcse_mean"], 0.)
        self.assertEqual(output["all_planned_mean_bounds"], [0., 0.])
        self.assertAlmostEqual(output["hoeffding_ci95"][0], -radius)
        self.assertAlmostEqual(output["hoeffding_ci95"][1], radius)
        self.assertEqual(tie(result)["estimate"], 1.)
        self.assertLess(tie(result)["hoeffding_ci95"][0], 1.)
        self.assertIsNone(tie(result)["bootstrap_ci95"])


if __name__ == "__main__":
    unittest.main()
