"""Hand-built count fixtures only; no generated confirmation or fitted data."""
import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import analyze_v2 as A
import analyze_confirmation_primary as P

SEED = 31009700


def config(run="confirm_core_a", datasets=1, full=1):
    names = P.SUPPORTED[run]
    return {"run_id": run, "phase": "confirmation", "base_seed": SEED, "datasets": datasets,
            "engines": ["repaired"], "cells": [{"name": "REF_H0" if run == "confirm_core_a" else "STAGE_H0"}],
            "bperm": 599, "alpha": .05, "rule": "le", "full_p_first_n": full,
            "schemes": {name: {"stratify": "none" if name == "unrestricted" else name,
                               "require_max": True} for name in names}}


def row(cfg, offset=0, a=(0, 100, 100), b=(10, 100, 100), *, nperm=599, max_a=0, max_b=10):
    schemes = {}
    for name, counts, maximum in zip(P.SUPPORTED[cfg["run_id"]], (a, b), (max_a, max_b)):
        schemes[name] = {"definition": copy.deepcopy(cfg["schemes"][name]), "requested_nperm": 599,
                         "nperm": nperm, "complete": nperm == 599, "gt": list(counts), "eq": [0, 0, 0],
                         "maxgt": maximum, "maxeq": 0}
    return {"run_id": cfg["run_id"], "cell": cfg["cells"][0]["name"], "engine": "repaired",
            "seed": SEED+offset, "status": "ok", "truth": {"pair_null": [True]*3},
            "schemes": schemes,
            "diagnostics": {"sequential_stopping": True, "full_p_selected": SEED+offset in P._selected(cfg)}}


def group(cfg, rows, **kwargs):
    return P.summarize_run(cfg, rows, bootstrap_repeats=40, **kwargs)["groups"][0]


class PrimaryTests(unittest.TestCase):
    def test_core_fwer_reverses_legacy_comparison_but_p_shift_does_not(self):
        cfg = config()
        r = row(cfg)
        output = group(cfg, [r])
        frozen_summary = A.summarize_group([r], bootstrap_repeats=40)
        legacy = frozen_summary["scheme_comparisons"]["diagnosis_minus_unrestricted"]["metrics"]
        effect = output["primary_bonferroni_rejection_difference"]
        self.assertEqual(effect["direction"], "unrestricted minus diagnosis")
        self.assertEqual(effect["estimate"], 1)
        self.assertEqual(legacy["global_reject"]["estimate"], -1)
        self.assertEqual(effect["complete_case"]["bootstrap_ci95"], [-legacy["global_reject"]["bootstrap_ci95"][1],
                                                                    -legacy["global_reject"]["bootstrap_ci95"][0]])
        self.assertIsNone(effect["bootstrap_ci95"])  # one constant dataset is not precise population evidence
        self.assertLess(effect["hoeffding_ci95"][0], effect["hoeffding_ci95"][1])
        p = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertEqual(p["direction"], "diagnosis minus unrestricted")
        self.assertAlmostEqual(p["estimate"], 1/180)
        self.assertAlmostEqual(p["estimate"], legacy["mean_pair_p_shift"]["estimate"])
        self.assertEqual(legacy["p_shift"]["estimate"], 0)  # mean is not the legacy median endpoint

    def test_shuffled_rows_keep_exact_seed_pairing_and_inputs_unchanged(self):
        cfg = config(datasets=3, full=3)
        rows = [row(cfg, 0), row(cfg, 1, a=(10, 10, 10), b=(0, 0, 0)), row(cfg, 2, a=(20, 20, 20), b=(30, 30, 30))]
        before = copy.deepcopy(rows)
        first, second = group(cfg, rows), group(cfg, list(reversed(rows)))
        self.assertEqual(first, second)
        self.assertEqual(rows, before)
        self.assertEqual(first["primary_bonferroni_rejection_difference"]["paired_table"], {"10": 1, "01": 1, "00": 1})

    def test_partial_known_scheme_tightens_planned_difference_bounds(self):
        cfg = config(datasets=3, full=2)
        first, second = row(cfg, 0), row(cfg, 1)
        del second["schemes"]["diagnosis"]
        output = group(cfg, [first, second])  # third planned seed is missing entirely
        effect = output["primary_bonferroni_rejection_difference"]
        self.assertEqual((effect["known_datasets"], effect["unknown_datasets"], effect["planned_datasets"]), (1, 2, 3))
        self.assertIsNone(effect["estimate"])
        self.assertIsNone(effect["bootstrap_ci95"])
        self.assertEqual(effect["complete_case"]["estimate"], 1)
        self.assertEqual(effect["all_planned_mean_bounds"], [0, 1])  # +1, [0,1], [-1,1]
        self.assertEqual(output["status_or_unavailable_reason_counts"], {"ok": 2, "missing_row": 1})
        self.assertIsNone(output["scheme_bonferroni_rejection_rates"]["diagnosis"]["estimate"])
        self.assertEqual(output["scheme_bonferroni_rejection_rates"]["diagnosis"]["all_planned_rate_bounds"], [0, 2/3])

    def test_full_p100_stays_100_when_only_one_selected_and_one_natural_completion_exist(self):
        cfg = config(datasets=101, full=100)
        output = group(cfg, [row(cfg, 0), row(cfg, 100, a=(599, 599, 599), b=(0, 0, 0))])
        p = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertEqual((p["planned_datasets"], p["known_datasets"], p["unknown_datasets"]), (100, 1, 99))
        self.assertFalse(p["complete"])
        self.assertIsNone(p["estimate"])
        self.assertAlmostEqual(p["complete_case"]["estimate"], 1/180)
        self.assertEqual(output["preselected_p"]["absolute_seeds"], list(range(SEED, SEED+100)))

    def test_selection_flag_mismatch_is_not_silently_used(self):
        cfg = config(datasets=2, full=1)
        for offset in (0, 1):
            r = row(cfg, offset)
            r["diagnostics"]["full_p_selected"] = not r["diagnostics"]["full_p_selected"]
            with self.subTest(offset=offset), self.assertRaisesRegex(ValueError, "selection"):
                group(cfg, [r])

    def test_early_exact_decisions_count_for_rejection_not_full_p(self):
        cfg = config()
        r = row(cfg, a=(10, 20, 30), b=(11, 21, 31), nperm=40, max_a=30, max_b=31)
        output = group(cfg, [r])
        self.assertTrue(output["primary_bonferroni_rejection_difference"]["complete"])
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["estimate"], 0)
        p = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertFalse(p["complete"])
        self.assertEqual(p["known_datasets"], 0)
        self.assertIsNone(p["estimate"])
        self.assertAlmostEqual(p["all_planned_mean_bounds"][0], -558/600)
        self.assertAlmostEqual(p["all_planned_mean_bounds"][1], 560/600)

    def test_three_pairs_are_averaged_before_seed_level_mcse(self):
        cfg = config(datasets=2, full=2)
        output = group(cfg, [row(cfg, 0, a=(0, 0, 0), b=(0, 360, 360)),
                             row(cfg, 1, a=(360, 0, 0), b=(0, 0, 0))])
        p = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertEqual(p["known_datasets"], 2)
        self.assertAlmostEqual(p["estimate"], .1)  # dataset means .4 and -.2
        self.assertAlmostEqual(p["mcse_mean"], .3)

    def test_zero_discordance_does_not_report_zero_width_primary_interval(self):
        cfg = config(datasets=2, full=2)
        output = group(cfg, [row(cfg, i, a=(10, 10, 10), b=(10, 10, 10)) for i in range(2)])
        effect = output["primary_bonferroni_rejection_difference"]
        self.assertEqual(effect["estimate"], 0)
        self.assertIsNone(effect["bootstrap_ci95"])
        self.assertTrue(effect["bootstrap_degenerate_sample"])
        self.assertEqual(effect["hoeffding_ci95"], [-1, 1])

    def test_saved_decision_does_not_override_uncertified_counts(self):
        cfg = config()
        r = row(cfg, a=(0, 0, 0), b=(0, 0, 0), nperm=0, max_a=0, max_b=0)
        for entry in r["schemes"].values():
            entry["exact_decision"] = {"alpha": .05, "le": {"pair_reject": [False]*3, "max_reject": False}}
        output = group(cfg, [r])
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["known_datasets"], 0)
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["all_planned_mean_bounds"], [-1, 1])

    def test_invalid_or_failed_fit_never_supplies_stale_success_counts(self):
        cfg = config()
        for status in ("invalid", "error", "pending", "running"):
            r = row(cfg)
            r["status"] = status
            output = group(cfg, [r])
            self.assertEqual(output["primary_bonferroni_rejection_difference"]["known_datasets"], 0)
            self.assertEqual(output["preselected_p"]["metrics"]["mean_pair_p_shift"]["known_datasets"], 0)

    def test_stage_direction_and_diagnostic_p_role_are_explicit(self):
        cfg = config("confirm_stage_a")
        output = group(cfg, [row(cfg)])
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["direction"], "oracle_dx_stage minus diagnosis")
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["estimate"], -1)
        self.assertEqual(output["preselected_p"]["role"], "execution_diagnostic_only_no_population_p_shift_claim")

    def test_global_max_remains_separate_from_bonferroni(self):
        cfg = config()
        output = group(cfg, [row(cfg, max_a=30, max_b=0)])
        self.assertEqual(output["primary_bonferroni_rejection_difference"]["estimate"], 1)
        self.assertEqual(output["global_max_separate"]["difference"]["estimate"], -1)

    def test_unplanned_duplicate_and_wrong_budget_rows_rejected(self):
        cfg = config()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            group(cfg, [row(cfg), row(cfg)])
        with self.assertRaisesRegex(ValueError, "Unplanned"):
            group(cfg, [row(cfg, 2)])
        r = row(cfg)
        r["schemes"]["diagnosis"]["requested_nperm"] = 99
        with self.assertRaisesRegex(ValueError, "budget"):
            group(cfg, [r])

    def test_empty_input_has_planned_denominators_and_no_estimated_effects(self):
        cfg = config(datasets=4, full=2)
        output = P.summarize_run(cfg, [], unavailable_reasons={("REF_H0", "repaired", SEED): "run_integrity_invalid"})
        g = output["groups"][0]
        self.assertEqual((output["input_row_count"], output["planned_row_count"]), (0, 4))
        self.assertEqual(g["primary_bonferroni_rejection_difference"]["all_planned_mean_bounds"], [-1, 1])
        self.assertIsNone(g["primary_bonferroni_rejection_difference"]["estimate"])
        self.assertEqual(g["status_or_unavailable_reason_counts"], {"run_integrity_invalid": 1, "missing_row": 3})
        json.dumps(output, allow_nan=False)
        self.assertIn("unavailable", P.render_markdown({"runs": [output]}))

    def test_cli_report_adapter_preserves_loader_rejections_as_planned_unknowns(self):
        cfg = config(datasets=2, full=2)
        run = {"cfg": cfg, "names": ["REF_H0"], "seeds": {SEED, SEED+1}, "snapshot": Path("toy_snapshot"),
               "trusted": True, "environment": None, "issues": [],
               "inventory": {"ready_for_final_analysis": False}}
        loader = SimpleNamespace(__file__=__file__,
            frozen_run=lambda *args: run, load_run=lambda *args, **kwargs: run,
            get_record=lambda run, cell, engine, seed: (row(cfg), None) if seed == SEED else (None, "observed_error"),
            sha=lambda path: "0"*64)
        entry = {"run_id": cfg["run_id"], "snapshot_sha256": "1"*64, "config_sha256": "2"*64}
        with TemporaryDirectory() as directory:
            protocol = Path(directory)/"protocol.json"
            protocol.write_text(json.dumps({"runs": [entry]}))
            with patch.dict("sys.modules", {"analyze_confirmation_cohorts": loader}):
                report = P.build_report(protocol, directory, run_ids=[cfg["run_id"]], bootstrap_repeats=10)
        result = report["runs"][0]["groups"][0]
        self.assertEqual(result["status_or_unavailable_reason_counts"], {"ok": 1, "observed_error": 1})
        self.assertIsNone(result["primary_bonferroni_rejection_difference"]["estimate"])
        self.assertFalse(report["fitting_or_generation_performed"])
        self.assertEqual(report["cross_run_environment_status"], "not_applicable_single_run")


if __name__ == "__main__":
    unittest.main()
