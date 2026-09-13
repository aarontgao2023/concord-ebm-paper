"""Inert count/metadata fixtures only; no simulated data or model fits."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import analyze_confirmation_operator_comparison as O

C = O.C
PROJECT = Path(__file__).resolve().parents[2]
PROTOCOL = PROJECT / "runs/v2/confirmation_protocol.json"
CONTRACT = PROJECT / "docs/review_v2/operator_comparison_contract_20260908.json"


def empty_runs():
    entries = {e["run_id"]: e for e in json.loads(PROTOCOL.read_text())["runs"]}
    runs = {name: C.frozen_run(entries[name], PROJECT) for name in O.RUN_IDS}
    for run in runs.values():
        run.update(trusted=True, environment={"fixture": "same-numerical-environment"},
                   inventory={"collection_complete": False, "ready_for_final_analysis": False, "plan": {}, "totals": {}})
    return runs


def entry(definition, counts=(9, 9, 9), *, nperm=599, attempted=None, max_count=100):
    attempted = nperm if attempted is None else attempted
    # Split the inclusive count into gt/eq so omitting ties changes the boundary.
    value = {"definition": definition, "requested_nperm": 599, "nperm": nperm,
             "attempted": attempted, "failed": attempted-nperm, "complete": nperm == 599,
             "gt": [max(k-1, 0) for k in counts], "eq": [int(k > 0) for k in counts],
             "maxgt": max_count, "maxeq": 0}
    value["exact_decision"] = {"alpha": .05}
    for rule in ("strict", "le"):
        value["exact_decision"][rule] = {
            "pair_reject": [C.I.bounds(k, nperm, 599, .05, 3, rule) for k in counts],
            "max_reject": C.I.bounds(max_count, nperm, 599, .05, 1, rule)}
    return value


def row(run, cell, seed, counts=(9, 9, 9)):
    return {"run_id": run["run_id"], "cell": cell, "engine": "repaired", "seed": seed,
            "status": "ok", "truth": {"config": run["designs"][cell],
                "biomarker_names": run["designs"][cell]["biomarker_names"],
                "manifest": {"data_sha256": "a"*64}},
            "observed": {"orderings": [[0, 1]]*3, "taus": [0., 0., 0.], "distance_to_truth": [0., 0., 0.]},
            "schemes": {name: entry(definition, counts) for name, definition in run["cfg"]["schemes"].items()},
            "diagnostics": {"quality_flags": ["legal_optimizer_unsuccessful"]}}


def add(run, cell, seed, value=None):
    value = value or row(run, cell, seed)
    run["records"].setdefault((cell, "repaired", seed), []).append({"row": value, "reason": None, "path": "inert-fixture:1"})
    return value


class OperatorComparisonTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.base = empty_runs()

    def setUp(self):
        self.runs = copy.deepcopy(self.base)

    def test_contract_and_exact_four_run_plan(self):
        O.validate_contract(json.loads(CONTRACT.read_text()), PROTOCOL.read_bytes())
        O.validate_plan(self.runs)
        self.runs["topup"] = self.runs[O.RUN_IDS[1]]
        with self.assertRaisesRegex(ValueError, "exactly the four"):
            O.validate_plan(self.runs)

    def test_fixed_B599_inclusive_boundary_is_pair_not_max(self):
        run = self.runs[O.RUN_IDS[0]]
        r = add(run, "REF_H0", 42100000)
        self.assertEqual(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0), (True, None))
        self.assertIs(r["schemes"]["diagnosis"]["exact_decision"]["le"]["max_reject"], False)
        r["schemes"]["diagnosis"] = entry(run["cfg"]["schemes"]["diagnosis"], (10, 10, 10))
        self.assertEqual(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0), (False, None))

    def test_failed_fit_is_unknown_not_attempted_as_success(self):
        run = self.runs[O.RUN_IDS[0]]
        r = add(run, "REF_H0", 42100000)
        r["schemes"]["diagnosis"] = entry(run["cfg"]["schemes"]["diagnosis"], nperm=598, attempted=599)
        self.assertEqual(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0), (None, "permutation_decision_unknown"))
        # Even an exhausted failed budget can conclusively establish nonrejection.
        r["schemes"]["diagnosis"] = entry(run["cfg"]["schemes"]["diagnosis"], (10, 10, 10), nperm=598, attempted=599)
        self.assertEqual(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0), (False, None))

    def test_nonselected_early_stop_and_stale_exact_rejected(self):
        run = self.runs[O.RUN_IDS[0]]
        r = add(run, "REF_H0", 42100999)
        r["schemes"]["diagnosis"] = entry(run["cfg"]["schemes"]["diagnosis"], (20, 20, 20), nperm=40, max_count=35)
        self.assertEqual(O.fullgroup_pair_decision(run, "REF_H0", 42100999, 1), (False, None))
        r["schemes"]["diagnosis"]["exact_decision"]["le"]["pair_reject"] = [True]*3
        value, reason = O.fullgroup_pair_decision(run, "REF_H0", 42100999, 1)
        self.assertIsNone(value)
        self.assertIn("disagree", reason)

    def test_wrong_operator_missing_definition_and_malformed_exact_unknown(self):
        for change in ("local", "definition", "exact", "nperm"):
            with self.subTest(change=change):
                run = copy.deepcopy(self.runs[O.RUN_IDS[0]])
                r = add(run, "REF_H0", 42100000)
                item = r["schemes"]["diagnosis"]
                if change == "local":
                    run["cfg"]["schemes"]["diagnosis"]["tested_pair_index"] = 0
                elif change == "definition":
                    item.pop("definition")
                elif change == "exact":
                    item["exact_decision"]["strict"] = []
                else:
                    item["nperm"] = True
                self.assertIsNone(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0)[0])

    def test_null_mask_and_both_discovery_logic(self):
        outcomes = [(True, None), (None, "unfinished"), (False, None)]
        self.assertEqual(O.endpoint_value(outcomes, (False, False, True), "any_true_null_rejection"), (False, None))
        self.assertEqual(O.endpoint_value(outcomes, (False, False, True), "any_false_null_rejection"), (True, None))
        self.assertEqual(O.endpoint_value(outcomes, (False, False, True), "all_false_null_rejection"), (None, "unfinished"))
        self.assertEqual(O.endpoint_value(outcomes, (True, False, False), "all_false_null_rejection"), (False, None))

    def test_missing_full_plan_has_seven_endpoints_and_no_primary_estimates(self):
        result = O.summarize_comparisons(self.runs)
        self.assertEqual(len(result["endpoints"]), 7)
        self.assertEqual(result["planned_dataset_cell_units"], 2000)
        self.assertFalse(result["ready_for_final_analysis"])
        for endpoint in result["endpoints"]:
            for name in ("fullgroup", "local"):
                s = endpoint[name]
                self.assertEqual(s["unknown_count"], endpoint["planned_R"])
                self.assertEqual(s["rate_bounds"], [0, 1])
                self.assertIsNone(s["rate"]); self.assertIsNone(s["wilson_ci95"])
            self.assertEqual(endpoint["paired_difference"]["difference_bounds"], [-1, 1])
            self.assertIsNone(endpoint["paired_difference"]["difference"])

    def test_mismatch_invalidates_both_without_modifying_caller(self):
        a = add(self.runs[O.RUN_IDS[0]], "REF_H0", 42100000)
        b = add(self.runs[O.RUN_IDS[1]], "REF_H0", 42100000)
        b["truth"]["manifest"]["data_sha256"] = "b"*64
        result = O.summarize_comparisons(self.runs)
        self.assertEqual(len(result["source_local_crosschecks"][0]["mismatches"]), 1)
        for name in ("fullgroup", "local"):
            self.assertEqual(result["endpoints"][0][name]["unknown_count"], 1000)
        self.assertIsNone(self.runs[O.RUN_IDS[0]]["records"][("REF_H0", "repaired", 42100000)][0]["reason"])

    def test_failed_source_null_truth_not_replaced_by_valid_local_or_flags(self):
        a = add(self.runs[O.RUN_IDS[0]], "REF_H0", 42100000)
        a.update(status="timeout", truth=None)
        add(self.runs[O.RUN_IDS[1]], "REF_H0", 42100000)
        result = O.summarize_comparisons(self.runs)
        endpoint = result["endpoints"][0]
        self.assertEqual(endpoint["fullgroup"]["unknown_count"], 1000)
        self.assertEqual(endpoint["local"]["reject_count"], 1)
        self.assertEqual(endpoint["paired_difference"]["unknown_pairs"], 1000)

    def test_duplicate_and_environment_mismatch_never_pick_winner(self):
        source = self.runs[O.RUN_IDS[0]]
        add(source, "REF_H0", 42100000); add(source, "REF_H0", 42100000)
        add(self.runs[O.RUN_IDS[1]], "REF_H0", 42100000)
        result = O.summarize_comparisons(self.runs)
        self.assertEqual(result["endpoints"][0]["fullgroup"]["unknown_count"], 1000)
        self.runs[O.RUN_IDS[1]]["environment"] = {"fixture": "other"}
        result = O.summarize_comparisons(self.runs)
        self.assertEqual(result["endpoints"][0]["local"]["unknown_count"], 1000)
        self.assertTrue(result["cross_run_issues"])

    def test_complete_paired_rate_direction_and_nondegenerate_interval(self):
        for cell, source_id, local_id, base, n, null in O.COHORTS:
            source, local = self.runs[source_id], self.runs[local_id]
            for seed in range(base, base+n):
                add(source, cell, seed, row(source, cell, seed, (10, 10, 10)))
                add(local, cell, seed, row(local, cell, seed, (9, 9, 9)))
        result = O.summarize_comparisons(self.runs)
        self.assertTrue(result["all_endpoint_decisions_complete"])
        self.assertFalse(result["ready_for_final_analysis"])  # Entire source run is not claimed done.
        for e in result["endpoints"]:
            self.assertEqual((e["fullgroup"]["rate"], e["local"]["rate"]), (0, 1))
            p = e["paired_difference"]
            self.assertEqual(p["difference"], 1)
            self.assertEqual(p["difference_orientation"], "local_minus_fullgroup")
            self.assertIsNone(p["paired_normal_ci95"])
            self.assertLess(p["paired_hoeffding_ci95"][0], 1)

    def test_no_overwrite_and_actual_input_protection(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            raw = root/"raw"; raw.mkdir()
            snapshot = root/"custom-snapshot"; snapshot.mkdir()
            original = snapshot/"config.json"; original.write_text("immutable")
            result = {"protected_snapshot_roots": [str(snapshot)], "input_sha256": {str(original): "unused"}}
            with self.assertRaisesRegex(ValueError, "existing output"):
                O.write_outputs(original.with_suffix(""), result, PROTOCOL, CONTRACT, raw, root)
            with self.assertRaisesRegex(ValueError, "scientific inputs"):
                O.write_outputs(snapshot/"new", result, PROTOCOL, CONTRACT, raw, root)
            self.assertEqual(original.read_text(), "immutable")
            alias = root/"alias.json"; alias.symlink_to(original)
            with self.assertRaisesRegex(ValueError, "existing output"):
                O.write_outputs(root/"alias", result, PROTOCOL, CONTRACT, raw, root)

    def test_changed_protocol_capture_cannot_release_an_estimate(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            protocol = root/"protocol.json"; protocol.write_bytes(PROTOCOL.read_bytes())
            contract = root/"contract.json"; contract.write_bytes(CONTRACT.read_bytes())
            with mock.patch.object(O, "load_compact_run", side_effect=lambda run, _: self.runs[run["run_id"]]):
                baseline = O.analyze(protocol, root/"raw", contract, PROJECT)
            self.assertTrue(baseline["inputs_stable"])
            def change_protocol(run, unused):
                if run["run_id"] == O.RUN_IDS[-1]:
                    protocol.write_bytes(protocol.read_bytes()+b" ")
                return self.runs[run["run_id"]]
            with mock.patch.object(O, "load_compact_run", side_effect=change_protocol):
                with self.assertRaisesRegex(ValueError, "Inputs or source changed"):
                    O.analyze(protocol, root/"raw", contract, PROJECT)


if __name__ == "__main__":
    unittest.main()
