"""Tiny invented records only; no scientific dataset generation or fitting."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

import analyze_confirmation_cohorts as C


class ConfirmationCohortTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.design = C.resolved_design(C.design_defaults(C.PROJECT / "scripts/v2/design_v2.py"), {"name": "REF_H0"})
        self.env = {"versions": {name: "1" for name in C.I.NUMERICAL_PACKAGES},
                    "source_sha256": {"x.py": "a" * 64}, "wheel_sha256": "b" * 64}
        self.sources = {name: "c" * 64 for name in C.I.RUNTIME_FILES}

    def tearDown(self):
        self.temp.cleanup()

    def row(self, seed=100, status="ok"):
        order = list(range(14))
        truth = {"config": self.design, "seed": seed, "manifest": {"seed": seed, "config_sha256": C.canonical_sha(self.design),
                 "source_sha256": self.sources["design_v2.py"], "data_sha256": "d" * 64, "numpy_version": "1", "pandas_version": "1"},
                 "biomarker_names": self.design["biomarker_names"], "group_order": list(C.GROUPS), "base_order": order,
                 "group_orderings": {g: order[:] for g in C.GROUPS}, "planned_inversions": {g: 0 for g in C.GROUPS},
                 "realized_inversions": {g: 0 for g in C.GROUPS}, "pair_null": [True] * 3, "h1_distance": 0.,
                 "pair_truth": {pair: {"null": True, "inversions": 0, "normalized_distance": 0.} for pair in C.PAIR_NAMES}}
        return {"run_id": "toy", "phase": "confirmation", "cell": "REF_H0", "engine": "repaired", "seed": seed, "status": status,
                "truth": truth, "observed": {"orderings": [order[:] for _ in range(3)], "taus": [0.] * 3, "distance_to_truth": [0.] * 3},
                "diagnostics": {"provenance": self.env, "quality_flags": ["optimizer_unsuccessful"]}, "schemes": {}}

    def make_run(self, rows=None, seeds=range(100, 102), name="toy"):
        return {"run_id": name, "cfg": {"run_id": name, "engines": ["repaired"], "bperm": 599, "alpha": .05, "rule": "le", "schemes": {}},
                "seeds": set(seeds), "full_p": set(), "names": ["REF_H0"], "sources": self.sources,
                "designs": {"REF_H0": self.design}, "trusted": True, "records": {(r["cell"], r["engine"], r["seed"]):
                    [{"row": r, "reason": None}] for r in rows or []}}

    def test_fixed_denominator_unknown_never_valid_only(self):
        result = C.summarize([(1, True, None), (2, False, None), (3, None, "missing"), (4, None, "timeout")])
        self.assertEqual((result["planned_R"], result["reject_count"], result["unknown_count"]), (4, 1, 2))
        self.assertEqual(result["rate_bounds"], [.25, .75])
        self.assertIsNone(result["rate"])
        self.assertIsNone(result["wilson_ci95"])
        self.assertIsNone(result["mcse"])

    def test_complete_endpoint_has_mc_uncertainty(self):
        result = C.summarize([(1, True, None), (2, False, None)])
        self.assertEqual(result["rate"], .5)
        self.assertAlmostEqual(result["mcse"], (.25 / 2) ** .5)
        self.assertLess(result["wilson_ci95"][0], .5)
        self.assertGreater(result["wilson_ci95"][1], .5)

    def test_three_valued_any_and_both_discovery(self):
        self.assertIs(C.any3([True, None]), True)
        self.assertIs(C.all3([False, None]), False)
        self.assertIsNone(C.all3([True, None]))
        self.assertIsNone(C.any3([False, None]))
        self.assertIs(C.all3([True, True]), True)

    def test_paired_unknown_bounds_use_all_seed_pairs(self):
        result = C.paired_summary([(1, True, None), (2, None, "missing")], [(1, False, None), (2, False, None)])
        self.assertEqual(result["difference_bounds"], [-1., -.5])
        self.assertEqual(result["planned_pairs"], 2)
        self.assertIsNone(result["difference"])
        self.assertIsNone(result["paired_mcse"])
        with self.assertRaises(ValueError):
            C.paired_summary([(1, True, None)], [(2, False, None)])

    def test_quality_flags_are_not_exclusion_criteria(self):
        row = self.row()
        C.validate_truth(row, self.design, self.sources["design_v2.py"], self.env)
        run = self.make_run([row])
        self.assertIs(C.get_record(run, "REF_H0", "repaired", 100)[0], row)

    def test_duplicate_rows_never_selected(self):
        run = self.make_run([self.row()])
        run["records"][("REF_H0", "repaired", 100)].append(deepcopy(run["records"][("REF_H0", "repaired", 100)][0]))
        self.assertEqual(C.get_record(run, "REF_H0", "repaired", 100), (None, "duplicate_identity"))

    def test_truth_inversions_and_observed_tau_are_recomputed(self):
        row = self.row()
        row["truth"]["realized_inversions"]["e4"] = 1
        with self.assertRaisesRegex(ValueError, "inversion"):
            C.validate_truth(row, self.design, self.sources["design_v2.py"], self.env)
        row = self.row()
        row["observed"]["taus"][0] = 1 / 91
        with self.assertRaisesRegex(ValueError, "inconsistency"):
            C.validate_truth(row, self.design, self.sources["design_v2.py"], self.env)

    def test_integer_counts_le_boundary_and_failed_permutations(self):
        row = self.row()
        run = self.make_run([row])
        definition = {"stratify": "diagnosis", "groups": [0, 1], "tested_pair_index": 0, "require_max": False}
        run["cfg"]["schemes"] = {"local": definition}
        entry = {"definition": definition, "requested_nperm": 599, "nperm": 599, "attempted": 599,
                 "failed": 0, "complete": True, "gt": [9, 0, 0], "eq": [0, 0, 0]}
        row["schemes"]["local"] = entry
        self.assertEqual(C.decision(run, "REF_H0", "repaired", 100, "local", 0), (True, None))
        run["cfg"]["rule"] = "strict"
        self.assertEqual(C.decision(run, "REF_H0", "repaired", 100, "local", 0), (False, None))
        run["cfg"]["rule"] = "le"
        entry.update(nperm=598, failed=1, complete=False)
        self.assertIsNone(C.decision(run, "REF_H0", "repaired", 100, "local", 0)[0])
        entry["gt"][0] = 20
        self.assertEqual(C.decision(run, "REF_H0", "repaired", 100, "local", 0), (False, None))
        with self.assertRaises(KeyError):
            C.decision(run, "REF_H0", "repaired", 100, "unregistered", 0)

    def source_pair(self, status="ok", missing_truth=False):
        a, b = self.row(), self.row()
        a["status"] = status
        if missing_truth:
            a["truth"] = {}
        source = self.make_run([a], range(100, 101), "source")
        target = self.make_run([b], range(100, 101), C.LOCAL_RUNS[0])
        target["cfg"]["data_and_fit_reuse"] = {"source_run_id": "source", "source_cells": ["REF_H0"], "engine": "repaired",
                                               "source_seed_start": 100, "source_seed_end_inclusive": 100}
        return source, target

    def test_source_failure_never_replaced_and_missing_truth_not_local_exclusion(self):
        source, target = self.source_pair("timeout", True)
        result = C.compare_sources({"source": source, C.LOCAL_RUNS[0]: target})[0]
        self.assertFalse(result["complete_crosscheck"])
        self.assertEqual(result["mismatches"], [])
        self.assertEqual(C.get_record(source, "REF_H0", "repaired", 100), (None, "observed_timeout"))
        self.assertIsNotNone(C.get_record(target, "REF_H0", "repaired", 100)[0])

    def test_successful_source_mismatch_invalidates_both_without_selecting_winner(self):
        source, target = self.source_pair()
        target["records"][("REF_H0", "repaired", 100)][0]["row"]["observed"]["taus"][0] = .1
        result = C.compare_sources({"source": source, C.LOCAL_RUNS[0]: target})[0]
        self.assertEqual(len(result["mismatches"]), 1)
        self.assertEqual(C.get_record(source, "REF_H0", "repaired", 100)[1], "source_local_crosscheck_mismatch")
        self.assertEqual(C.get_record(target, "REF_H0", "repaired", 100)[1], "source_local_crosscheck_mismatch")

    def test_partial_segments_disjoint_and_one_target_only(self):
        design = deepcopy(self.design)
        design.update(name="PWR_E2_K27", alternative_group="e2", target_inversions=27)
        first, top = self.make_run(seeds=range(100, 102), name="first"), self.make_run(seeds=range(102, 104), name=C.PARTIAL_RUNS[0])
        definition = {"stratify": "diagnosis", "groups": [1, 2], "tested_pair_index": 2, "require_max": False}
        for run in (first, top):
            run["designs"] = {"PWR_E2_K27": design}
            run["cfg"]["schemes"] = {"target": deepcopy(definition)}
        first["cfg"]["schemes"]["target"]["family"] = C.FAMILY
        cohort = {"cohort_id": "toy4", "cell": "PWR_E2_K27", "engine": "repaired", "target_scheme": "target", "target_pair_index": 2,
                  "combined_target_R": 4, "segments": [{"run_id": "first", "seed_start": 100, "seed_end_inclusive": 101, "R": 2},
                  {"run_id": C.PARTIAL_RUNS[0], "seed_start": 102, "seed_end_inclusive": 103, "R": 2}]}
        top["cfg"]["cohort_definition"] = cohort
        runs = {"first": first, C.PARTIAL_RUNS[0]: top}
        result = C.partial_cohorts(runs)[0]
        self.assertEqual(result["summary"]["planned_R"], 4)
        self.assertEqual(result["summary"]["unknown_count"], 4)
        self.assertEqual(result["target_pair_index"], 2)
        cohort["segments"][1].update(seed_start=101, seed_end_inclusive=102)
        with self.assertRaisesRegex(ValueError, "overlapping"):
            C.partial_cohorts(runs)

    def test_family_power_keeps_first500_and_both_endpoint(self):
        run = self.make_run(seeds=range(100, 600), name=C.LOCAL_RUNS[1])
        design = deepcopy(self.design)
        design.update(name="PWR_E4_K27", alternative_group="e4", target_inversions=27)
        run["names"] = ["PWR_E4_K27"]
        run["designs"] = {"PWR_E4_K27": design}
        run["cfg"]["family_definition"] = {"id": C.FAMILY, "level": .05, "pair_order": list(C.PAIR_NAMES)}
        run["cfg"]["schemes"] = {f"p{idx}": {"tested_pair_index": idx, "family": C.FAMILY, "groups": list(C.PAIRS[idx]),
                                                       "stratify": "diagnosis", "require_max": False} for idx in range(3)}
        result = C.local_families({C.LOCAL_RUNS[1]: run})[0]
        for summary in result["endpoints"].values():
            self.assertEqual(summary["planned_R"], 500)
            self.assertEqual(summary["unknown_count"], 500)
        self.assertIn("both_true_alternatives_discovery", result["endpoints"])
        self.assertEqual(result["pair_null"], [True, False, False])

    def test_real_frozen_empty_results_yields_full_planned_bounds(self):
        result = C.analyze(C.PROJECT / "runs/v2/confirmation_protocol.json", self.root / "absent")
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertFalse(result["collection_complete"])
        self.assertEqual(len(result["global_power_panel"]), 6)
        self.assertEqual([item["summary"]["planned_R"] for item in result["partial_null_cohorts"]], [1000] * 3)
        for item in result["global_power_panel"]:
            self.assertEqual(item["summary"]["unknown_count"], 500)
            self.assertEqual(item["summary"]["rate_bounds"], [0., 1.])
        self.assertEqual(sum(r["plan"]["planned_rows"] for r in result["runs"]), 16800)

    def test_zero_empirical_paired_difference_does_not_report_zero_width_ci(self):
        values = [(i, False, None) for i in range(20)]
        result = C.paired_summary(values, values)
        self.assertTrue(result["degenerate_empirical_difference"])
        self.assertIsNone(result["paired_normal_ci95"])
        self.assertLess(result["paired_hoeffding_ci95"][0], 0)
        self.assertGreater(result["paired_hoeffding_ci95"][1], 0)

    def test_local_target_cannot_relabel_other_group_distribution(self):
        definition = {"tested_pair_index": 0, "groups": [1, 2], "stratify": "diagnosis", "require_max": False}
        with self.assertRaisesRegex(ValueError, "designated"):
            C.validate_local_definition(definition, 0)

    def loaded_fixture(self, failure=False, missing=False, tamper_source=False):
        # Uses invented 100/101 rows, not a confirmation namespace dataset.
        self.sources = {name: C.sha(C.PROJECT / "scripts/v2" / name) for name in C.I.RUNTIME_FILES + C.I.OPTIONAL_RUNTIME_FILES}
        dep = C.literal_assignments(C.PROJECT / "scripts/v2/dependency_v2.py")
        self.env["versions"]["pyebm"] = dep["PYEBM_VERSION"]
        self.env.update(source_sha256=dep["SOURCE_SHA256"], wheel_sha256=dep["PYEBM_WHEEL_SHA256"])
        rows = [self.row(100), self.row(101)]
        if failure:
            rows[0].update(status="timeout", truth={}, observed={})
        if missing:
            rows = rows[:1]
        run = self.make_run()
        run["cfg"].update(phase="confirmation", base_seed=100, datasets=2, cells=[{"name": "REF_H0"}], bperm=0)
        run.update(config_path=self.root / "toy.json", snapshot=C.PROJECT, dependency=dep, issues=[])
        run["config_path"].write_text(json.dumps(run["cfg"]))
        chunk = self.root / "raw/toy/cell_0/chunk_0"
        chunk.mkdir(parents=True)
        identity = {"config": run["cfg"], "cell_index": 0, "chunk": 0, "nchunks": 1, "source_sha256": dict(self.sources)}
        if tamper_source:
            identity["source_sha256"]["run_v2.py"] = "f" * 64
        (chunk / "manifest.json").write_text(json.dumps(identity | {"signature": C.I.digest(identity), "environment": self.env, "seeds": [100, 101]}))
        (chunk / "rows.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
        (chunk / "progress.json").write_text(json.dumps({"done": True, "rows": len(rows), "expected_rows": 2}))
        return C.load_run(run, self.root / "raw")

    def test_loader_done_does_not_supply_missing_dataset(self):
        run = self.loaded_fixture(missing=True)
        self.assertFalse(run["inventory"]["collection_complete"])
        self.assertEqual(run["inventory"]["totals"]["missing_planned_rows"], 1)
        self.assertIsNotNone(C.get_record(run, "REF_H0", "repaired", 100)[0])
        self.assertEqual(C.get_record(run, "REF_H0", "repaired", 101), (None, "missing_row"))

    def test_loader_terminal_failure_missing_truth_retains_planned_failure(self):
        run = self.loaded_fixture(failure=True)
        self.assertTrue(run["inventory"]["collection_complete"])
        self.assertFalse(run["inventory"]["ready_for_final_analysis"])
        self.assertTrue(run["trusted"])
        self.assertEqual(C.get_record(run, "REF_H0", "repaired", 100), (None, "observed_timeout"))
        self.assertEqual(run["issues"][0]["code"], "failed_row_truth_unavailable_or_invalid")
        self.assertEqual(run["inventory"]["totals"]["observed_failed_or_invalid"], 1)

    def test_loader_rejects_resigned_manifest_with_wrong_frozen_source(self):
        run = self.loaded_fixture(tamper_source=True)
        self.assertFalse(run["trusted"])
        self.assertEqual(C.get_record(run, "REF_H0", "repaired", 100), (None, "run_integrity_invalid"))


if __name__ == "__main__":
    unittest.main()
