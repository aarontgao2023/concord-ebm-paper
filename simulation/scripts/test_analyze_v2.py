"""Statistical boundary, missingness, and historical-result regression checks."""
import copy
import csv
import json
from pathlib import Path
import statistics
import tempfile
import unittest

import analyze_v2 as A


def scheme(n=599, budget=599, gt=None, eq=None):
    return {"nperm": n, "requested_nperm": budget, "complete": n == budget,
            "gt": gt or [0, 0, 0], "eq": eq or [0, 0, 0]}


def row(seed=1, null=None, entries=None):
    return {"run_id": "test", "engine": "test_engine", "cell": "test_cell", "seed": seed,
            "status": "ok", "truth": {"pair_null": null or [True, True, True]},
            "observed": {}, "schemes": entries or {}, "source": "fixture"}


class StatisticalBoundaryTests(unittest.TestCase):
    def test_strict_and_weak_b599_boundary(self):
        self.assertEqual(A.rejection_rank(599, rule="strict"), 9)
        self.assertEqual(A.rejection_rank(599, rule="le"), 10)
        self.assertFalse(A.rejects(9, 599, rule="strict"))
        self.assertTrue(A.rejects(9, 599, rule="le"))
        self.assertEqual(A.rejection_rank(99, rule="strict"), 1)

    def test_per_row_budget_changes_boundary(self):
        # Reusing a group-level B599 denominator would incorrectly reject B598.
        self.assertFalse(A.evaluate_scheme(scheme(598, 598, [9, 9, 9]), .05, "le")["pair_reject"][0])
        self.assertTrue(A.evaluate_scheme(scheme(599, 599, [9, 9, 9]), .05, "le")["pair_reject"][0])

    def test_early_stop_keeps_exact_decision_but_not_p(self):
        result = A.evaluate_scheme(scheme(20, 599, [10, 11, 12]), .05, "le")
        self.assertEqual(result["pair_reject"], [False] * 3)
        self.assertIsNone(result["pge"])

    def test_unfinished_promising_tail_is_unevaluable(self):
        result = A.evaluate_scheme(scheme(20, 599), .05, "le")
        self.assertEqual(result["pair_reject"], [None] * 3)
        self.assertIsNone(A.any_known(result["pair_reject"]))

    def test_missing_one_fit_can_have_certified_decisions(self):
        result = A.evaluate_scheme(scheme(598, 599, [502, 153, 106]), .05, "strict")
        self.assertEqual(result["pair_reject"], [False] * 3)
        self.assertIsNone(result["pge"])

    def test_exact_declaration_cannot_conflict_or_omit_alpha(self):
        entry = scheme(20, 599, [10, 11, 12])
        entry["exact_decision"] = {"alpha": .05, "le": {"pair_reject": [True, False, False]}}
        with self.assertRaisesRegex(ValueError, "conflicts"):
            A.evaluate_scheme(entry, .05, "le")
        entry["exact_decision"].pop("alpha")
        with self.assertRaisesRegex(ValueError, "specify their alpha"):
            A.evaluate_scheme(entry, .05, "le")

    def test_counts_do_not_round_or_truncate(self):
        self.assertEqual(A.integer("45.0", "count"), 45)
        for value in ("45.1", -1, True):
            with self.assertRaises(ValueError):
                A.integer(value, "count")

    def test_duplicate_identity_raises_instead_of_deduplicating(self):
        first = row()
        with self.assertRaisesRegex(ValueError, "duplicate dataset identity"):
            A.validate_identities([first, copy.deepcopy(first)])
        second = copy.deepcopy(first)
        second["run_id"] = "different_run"
        A.validate_identities([first, second])

    def test_partial_null_separates_global_true_and_false_rejections(self):
        sample = row(null=[True, False, False])
        result = A.endpoints(sample, {"pair_reject": [True, False, False], "max_reject": None})
        self.assertTrue(result["global_reject"])
        self.assertFalse(result["any_true_alternative_reject"])
        self.assertTrue(result["any_true_null_reject"])

    def test_failures_not_counted_as_nonrejections(self):
        good = row(entries={"stratified": scheme(gt=[50, 50, 50])})
        failed = row(2)
        failed["status"] = "failed_observed"
        summary = A.summarize_group([good, failed], bootstrap_repeats=0)
        metric = summary["schemes"]["stratified"]["metrics"]["global_reject"]
        self.assertEqual(metric["n_evaluable"], 1)
        self.assertEqual(metric["n_unevaluable"], 1)
        self.assertEqual(metric["all_attempted_rate_bounds"], [0, .5])
        self.assertEqual(summary["schemes"]["stratified"]["metrics"]["any_true_alternative_reject"]["n_inapplicable"], 2)

    def test_cluster_shift_uses_one_value_per_dataset(self):
        a = row(1, entries={"unstratified": scheme(gt=[0, 0, 0]),
                            "stratified": scheme(gt=[60, 120, 180])})
        b = row(2, entries={"unstratified": scheme(gt=[0, 0, 0]),
                            "stratified": scheme(gt=[120, 180, 240])})
        summary = A.summarize_group([a, b], bootstrap_repeats=10)
        self.assertEqual(summary["paired"]["p_shift"]["n_datasets"], 2)
        self.assertAlmostEqual(summary["paired"]["p_shift"]["estimate"], .25)

    def test_new_scheme_names_populate_paired_with_explicit_direction(self):
        sample = row(entries={"unrestricted": scheme(), "diagnosis": scheme(gt=[30, 40, 50])})
        sample["diagnostics"] = {"full_p_selected": True, "sequential_stopping": True}
        summary = A.summarize_group([sample], bootstrap_repeats=0)
        self.assertEqual(summary["primary_comparison"], "diagnosis_minus_unrestricted")
        self.assertEqual(summary["paired"]["global_reject"]["estimate"], -1)
        self.assertEqual(summary["paired"]["global_reject"]["direction"], "diagnosis minus unrestricted")
        self.assertEqual(summary["paired"]["mean_pair_p_shift"]["n_datasets"], 1)

    def test_mean_pair_shift_differs_from_legacy_median_and_uses_seed_mcse(self):
        rows = []
        for seed, counts in enumerate(([0, 0, 540], [60, 120, 180], [60, 60, 60]), 1):
            sample = row(seed, entries={"unrestricted": scheme(), "diagnosis": scheme(gt=list(counts))})
            sample["diagnostics"] = {"full_p_selected": True}
            rows.append(sample)
        summary = A.summarize_group(rows, bootstrap_repeats=50)
        self.assertAlmostEqual(summary["paired"]["p_shift"]["estimate"], .1)
        metric = summary["paired"]["mean_pair_p_shift"]
        self.assertAlmostEqual(metric["estimate"], .2)
        self.assertEqual(metric["n_datasets"], 3)
        self.assertAlmostEqual(metric["mcse_mean"], statistics.stdev([.3, .2, .1])/(3**.5))
        self.assertIsNotNone(metric["bootstrap_ci95"])
        self.assertAlmostEqual(summary["paired"]["pair_p_shift_e33_vs_e4"]["estimate"], (.9+.3+.1)/3)

    def test_new_mean_endpoint_requires_explicit_selection_and_complete_both_budgets(self):
        rows = []
        for seed in range(1, 5):
            sample = row(seed, entries={"unrestricted": scheme(), "diagnosis": scheme(gt=[60]*3)})
            if seed != 3:
                sample["diagnostics"] = {"full_p_selected": seed != 2, "sequential_stopping": True}
            if seed == 4:
                sample["schemes"]["diagnosis"] = scheme(n=598, gt=[60]*3)
            rows.append(sample)
        summary = A.summarize_group(rows, bootstrap_repeats=0)
        metric = summary["paired"]["mean_pair_p_shift"]
        self.assertEqual(metric["n_datasets"], 1)
        self.assertEqual(metric["n_selected"], 2)
        self.assertEqual(metric["n_not_selected"], 2)
        self.assertEqual(metric["n_unevaluable"], 1)
        self.assertAlmostEqual(metric["estimate"], .1)
        # The historical endpoint intentionally keeps nonsequential full rows.
        self.assertEqual(summary["paired"]["p_shift"]["n_datasets"], 2)

    def test_stage_comparison_is_separate_and_primary_can_be_selected(self):
        u, d, o = scheme(), scheme(gt=[60]*3), scheme(gt=[180]*3)
        for entry, count in ((u, 0), (d, 30), (o, 150)):
            entry.update(maxgt=count, maxeq=0)
        sample = row(entries={"unrestricted": u, "diagnosis": d, "oracle_dx_stage": o})
        sample["diagnostics"] = {"full_p_selected": True}
        summary = A.summarize_group([sample], bootstrap_repeats=0)
        self.assertEqual(summary["primary_comparison"], "diagnosis_minus_unrestricted")
        stage = summary["scheme_comparisons"]["oracle_dx_stage_minus_diagnosis"]["metrics"]
        self.assertAlmostEqual(stage["mean_pair_p_shift"]["estimate"], .2)
        self.assertAlmostEqual(stage["max_p_shift"]["estimate"], .2)
        changed = A.summarize_group([sample], bootstrap_repeats=0,
                                    primary_comparison="oracle_dx_stage_minus_diagnosis")
        self.assertEqual(changed["paired"], changed["scheme_comparisons"]["oracle_dx_stage_minus_diagnosis"]["metrics"])
        self.assertIn("diagnosis_minus_unrestricted", changed["scheme_comparisons"])

    def test_paired_inapplicable_truth_endpoint_is_not_a_failure(self):
        null = row(1, entries={"unrestricted": scheme(), "diagnosis": scheme()})
        alt = row(2, null=[False]*3, entries={"unrestricted": scheme(), "diagnosis": scheme()})
        failed = row(3, entries={"unrestricted": scheme(), "diagnosis": scheme()})
        failed["status"] = "error"
        summary = A.summarize_group([null, alt, failed], bootstrap_repeats=0)
        metric = summary["paired"]["any_true_null_reject"]
        self.assertEqual(metric["n_datasets"], 1)
        self.assertEqual(metric["n_inapplicable"], 1)
        self.assertEqual(metric["n_unevaluable"], 1)
        alternative = summary["paired"]["any_true_alternative_reject"]
        self.assertEqual(alternative["n_inapplicable"], 2)
        self.assertEqual(alternative["n_unevaluable"], 0)

    def test_comparison_csv_preserves_new_endpoints(self):
        sample = row(entries={"unrestricted": scheme(), "diagnosis": scheme(gt=[60]*3)})
        sample["diagnostics"] = {"full_p_selected": True}
        summary = A.summarize_group([sample], bootstrap_repeats=0)
        with tempfile.TemporaryDirectory() as directory:
            A.write_outputs({"groups": [summary], "mechanism_comparisons": []}, directory)
            with (Path(directory)/"scheme_comparisons.csv").open() as handle:
                values = list(csv.DictReader(handle))
            metric = next(value for value in values if value["endpoint"] == "mean_pair_p_shift")
            self.assertEqual(metric["comparison"], "diagnosis_minus_unrestricted")
            self.assertEqual(metric["n_datasets"], "1")

    def test_mechanism_only_row_does_not_invent_test(self):
        sample = row()
        sample["observed"] = {"taus": [.1, .2, .3], "distance_to_truth": [.4, .5, .6]}
        summary = A.summarize_group([sample], bootstrap_repeats=0)
        self.assertEqual(summary["schemes"], {})
        self.assertEqual(summary["observed_statistics"]["distance_to_truth"]["e4"]["mean"], .6)

    def test_sequentially_completed_rows_do_not_select_p_distribution(self):
        a = row(1, entries={"unstratified": scheme(), "stratified": scheme(gt=[60, 120, 180])})
        b = row(2, entries={"unstratified": scheme(), "stratified": scheme(gt=[120, 180, 240])})
        a["diagnostics"] = {"sequential_stopping": True, "full_p_selected": True}
        b["diagnostics"] = {"sequential_stopping": True, "full_p_selected": False}
        summary = A.summarize_group([a, b], bootstrap_repeats=0)
        self.assertEqual(summary["schemes"]["stratified"]["n_full_p_datasets"], 2)
        self.assertEqual(summary["schemes"]["stratified"]["n_p_summary_datasets"], 1)
        self.assertEqual(summary["paired"]["p_shift"]["n_datasets"], 1)
        self.assertAlmostEqual(summary["paired"]["p_shift"]["estimate"], .2)

    def test_mechanism_engine_comparison_pairs_seeds_and_labels_oracle(self):
        a = row(1)
        a["engine"] = "original"
        a["observed"] = {"distance_to_truth": [.6, .4, .2]}
        b = row(1)
        b["engine"] = "repaired"
        b["observed"] = {"distance_to_truth": [.3, .2, .1]}
        c = row(1)
        c["engine"] = "oracle"
        c["observed"] = {"distance_to_truth": [.1, .1, .1]}
        unpaired = copy.deepcopy(b)
        unpaired["seed"] = 2
        result = A.mechanism_comparisons([a, b, c, unpaired], ["original", "repaired", "oracle"], 10)
        comparison = next(r for r in result if r["engine_b"] == "repaired" and r["group"] == "group_mean")
        self.assertEqual(comparison["n_shared_seeds"], 1)
        self.assertEqual(comparison["n_seed_union"], 2)
        self.assertAlmostEqual(comparison["paired_difference"]["estimate"], -.2)
        for comparison in result:
            if "oracle" in (comparison["engine_a"], comparison["engine_b"]):
                self.assertEqual(comparison["comparison_role"], "oracle_mechanism_diagnostic")
                self.assertFalse(comparison["oracle_is_inference_candidate"])

    def test_single_pair_scheme_cannot_claim_other_pairs_or_global(self):
        entry = scheme(gt=[100, 0, 0])
        entry["maxgt"], entry["maxeq"] = 0, 0
        entry["definition"] = {"tested_pair_index": 0, "family": "pairwise_diagnosis"}
        sample = row(entries={"PDX01": entry})
        result = A.summarize_group([sample], bootstrap_repeats=0)
        detail = result["schemes"]["PDX01"]
        self.assertEqual(detail["scope"], "single_pair_only")
        self.assertEqual(detail["metrics"]["pair_e2_vs_e33"]["rate"], 0)
        for endpoint in ("global_reject", "max_global_reject", "pair_e2_vs_e4"):
            self.assertEqual(detail["metrics"][endpoint]["n_inapplicable"], 1)
            self.assertIsNone(detail["metrics"][endpoint]["rate"])
        family = result["pairwise_families"]["pairwise_diagnosis"]
        self.assertIsNone(family["metrics"]["global_reject"]["rate"])
        self.assertEqual(family["n_full_p_datasets"], 0)

    def test_pairwise_family_uses_each_member_budget_and_preserves_partial_null(self):
        entries = {}
        for name, index, budget, counts in (("PDX01", 0, 599, [100, 0, 0]),
                                             ("PDX02", 1, 598, [0, 9, 0]),
                                             ("PDX12", 2, 99, [0, 0, 0])):
            entry = scheme(budget, budget, counts)
            entry["definition"] = {"tested_pair_index": index, "family": "pairwise_diagnosis"}
            entries[name] = entry
        sample = row(null=[True, False, False], entries=entries)
        result = A.summarize_group([sample], bootstrap_repeats=0)
        family = result["pairwise_families"]["pairwise_diagnosis"]
        self.assertEqual(family["metrics"]["global_reject"]["rate"], 1)
        self.assertEqual(family["metrics"]["any_true_alternative_reject"]["rate"], 1)
        self.assertEqual(family["metrics"]["any_true_null_reject"]["rate"], 0)
        self.assertEqual(family["metrics"]["pair_e2_vs_e4"]["rate"], 0)
        self.assertEqual(family["n_full_p_datasets"], 1)
        self.assertEqual(family["full_p_summary"]["e2_vs_e4"]["mean"], 10 / 599)


class LegacyRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        directory = Path(__file__).resolve().parents[2] / "results" / "by_cell"
        if not directory.exists():
            raise unittest.SkipTest("legacy data not present")
        cls.rows = A.load_legacy(directory)

    def test_legacy_counts_and_b599_both_decision_rules(self):
        self.assertEqual(len(self.rows), 1700)
        A.validate_identities(self.rows)
        reference = [r for r in self.rows if r["run_id"] == "legacy_B599"]
        for rule, expected in (("strict", 19), ("le", 23)):
            result = A.summarize_group(reference, rule=rule, bootstrap_repeats=0)
            self.assertEqual(result["schemes"]["unstratified"]["metrics"]["global_reject"]["reject_count"], expected)
            self.assertEqual(result["schemes"]["stratified"]["metrics"]["global_reject"]["reject_count"], 5)
            self.assertEqual(result["paired"]["p_shift"]["n_datasets"], 99)
            self.assertEqual(result["legacy_audit"]["recorded_p_shift_dataset_median"]["n_datasets"], 100)

    def test_legacy_partial_null_and_failed_distance_target(self):
        data = [r for r in self.rows if r["cell"] == "PWR_0.30"]
        result = A.summarize_group(data, rule="strict", bootstrap_repeats=0)
        metrics = result["schemes"]["stratified"]["metrics"]
        self.assertEqual(metrics["global_reject"]["reject_count"], 11)
        self.assertEqual(metrics["any_true_alternative_reject"]["reject_count"], 10)
        self.assertEqual(metrics["any_true_null_reject"]["reject_count"], 2)
        data = [r for r in self.rows if r["cell"] == "PWR_0.60"]
        result = A.summarize_group(data, rule="strict", bootstrap_repeats=0)
        self.assertEqual(result["target_distance_check"], {"n_assessable": 100, "n_achieved": 89})


if __name__ == "__main__":
    unittest.main()
