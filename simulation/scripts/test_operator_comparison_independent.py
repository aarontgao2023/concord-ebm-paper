"""Independent count/identity fixtures; no data generator, fit or HPC calls."""
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import analyze_confirmation_operator_comparison as O

C = O.C
PROJECT = Path(__file__).resolve().parents[2]
PROTOCOL = PROJECT / "runs/v2/confirmation_protocol.json"
CONTRACT = PROJECT / "docs/review_v2/operator_comparison_contract_20260908.json"


def counts(definition, values, nperm=599, attempted=None, ties=False):
    attempted = nperm if attempted is None else attempted
    item = {"definition": copy.deepcopy(definition), "requested_nperm": 599,
            "nperm": nperm, "attempted": attempted, "failed": attempted-nperm,
            "complete": nperm == 599, "gt": [0, 0, 0] if ties else list(values),
            "eq": list(values) if ties else [0, 0, 0]}
    if definition.get("require_max"):
        item.update(maxgt=100 if nperm >= 100 else 0, maxeq=0)
    return item


def put(run, cell, seed, values=(10, 10, 10), status="ok"):
    row = {"run_id": run["run_id"], "cell": cell, "engine": "repaired", "seed": seed,
           "status": status, "truth": {"config": run["designs"][cell],
            "biomarker_names": run["designs"][cell]["biomarker_names"],
            "manifest": {"data_sha256": hashlib.sha256(f"inert:{cell}:{seed}".encode()).hexdigest()}},
           "observed": {"orderings": [list(range(14))]*3, "taus": [0., 0., 0.], "distance_to_truth": [0., 0., 0.]},
           "schemes": {name: counts(definition, values) for name, definition in run["cfg"]["schemes"].items()}}
    wrapper = {"row": row, "reason": None, "path": "independent-count-fixture"}
    run["records"].setdefault((cell, "repaired", seed), []).append(wrapper)
    return row


class IndependentOperatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        entries = {x["run_id"]: x for x in json.loads(PROTOCOL.read_text())["runs"]}
        cls.empty = {name: C.frozen_run(entries[name], PROJECT) for name in O.RUN_IDS}
        for run in cls.empty.values():
            run.update(trusted=True, environment={"inert_fixture": "common"},
                       inventory={"collection_complete": False, "ready_for_final_analysis": False, "plan": {}, "totals": {}})

    def setUp(self):
        self.runs = copy.deepcopy(self.empty)

    def test_inclusive_pair_counts_and_failed_budget_boundaries_independent_of_max(self):
        run = self.runs[O.RUN_IDS[0]]
        row = put(run, "REF_H0", 42100000)
        definition = run["cfg"]["schemes"]["diagnosis"]
        row["schemes"]["diagnosis"] = counts(definition, (10, 9, 0), ties=True)
        self.assertEqual([O.fullgroup_pair_decision(run, "REF_H0", 42100000, i)[0] for i in range(3)], [False, True, True])
        # Global-max count100 is nonrejecting while the latter two pairs reject.
        self.assertEqual(row["schemes"]["diagnosis"]["maxgt"], 100)
        for successes, expected in [(590, True), (589, None)]:
            with self.subTest(successes=successes):
                row["schemes"]["diagnosis"] = counts(definition, (0, 0, 0), nperm=successes, attempted=599)
                self.assertIs(O.fullgroup_pair_decision(run, "REF_H0", 42100000, 0)[0], expected)

    def test_asymmetric_complete_direction_then_one_missing_pair_keeps_R(self):
        source, local = self.runs[O.RUN_IDS[0]], self.runs[O.RUN_IDS[1]]
        for offset in range(1000):
            seed = 42100000+offset
            put(source, "REF_H0", seed, (0, 10, 10) if offset == 0 else (10, 10, 10))
            put(local, "REF_H0", seed, (0, 10, 10) if offset in (1, 2) else (10, 10, 10))
        report = O.summarize_comparisons(self.runs)
        endpoint = report["endpoints"][0]; paired = endpoint["paired_difference"]
        self.assertEqual(endpoint["planned_R"], 1000)
        self.assertEqual((endpoint["fullgroup"]["rate"], endpoint["local"]["rate"]), (.001, .002))
        self.assertEqual(paired["difference_orientation"], "local_minus_fullgroup")
        self.assertEqual(paired["paired_table"], {"10": 1, "01": 2, "00": 997})
        self.assertAlmostEqual(paired["difference"], .001)
        self.assertAlmostEqual(paired["paired_mcse"], math.sqrt((3-1/1000)/(999*1000)))
        self.assertEqual(report["planned_dataset_cell_units"], 2000)
        self.assertIn("not pooled as independent", report["counting_scope"])
        local["records"].pop(("REF_H0", "repaired", 42100999))
        endpoint = O.summarize_comparisons(self.runs)["endpoints"][0]
        self.assertEqual(endpoint["fullgroup"]["rate"], .001)
        self.assertIsNone(endpoint["local"]["rate"])
        self.assertEqual(endpoint["paired_difference"]["planned_pairs"], 1000)
        self.assertEqual(endpoint["paired_difference"]["unknown_pairs"], 1)
        self.assertEqual(endpoint["paired_difference"]["difference_bounds"], [.001, .002])
        for field in ("difference", "paired_mcse", "paired_normal_ci95", "paired_hoeffding_ci95"):
            self.assertIsNone(endpoint["paired_difference"][field])

    def test_observed_mismatch_invalidates_both_but_failure_never_uses_counterpart(self):
        source, local = self.runs[O.RUN_IDS[0]], self.runs[O.RUN_IDS[1]]
        first = put(source, "REF_H0", 42100000, (0, 10, 10))
        second = put(local, "REF_H0", 42100000, (0, 10, 10))
        second["observed"]["taus"] = [1., 0., 0.]
        report = O.summarize_comparisons(self.runs)
        self.assertEqual(len(report["source_local_crosschecks"][0]["mismatches"]), 1)
        for name in ("fullgroup", "local"):
            self.assertEqual(report["endpoints"][0][name]["unknown_count"], 1000)
        self.assertIsNone(source["records"][("REF_H0", "repaired", 42100000)][0]["reason"])
        first.update(status="timeout", truth=None)
        report = O.summarize_comparisons(self.runs)
        endpoint = report["endpoints"][0]
        self.assertEqual(endpoint["fullgroup"]["unknown_reason_counts"]["observed_timeout"], 1)
        self.assertEqual(endpoint["local"]["reject_count"], 1)
        self.assertEqual(endpoint["paired_difference"]["unknown_pairs"], 1000)
        self.assertIsNone(first["truth"])

    def test_new_input_file_during_capture_is_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); protocol = root/"protocol.json"; protocol.write_bytes(PROTOCOL.read_bytes())
            contract = root/"contract.json"; contract.write_bytes(CONTRACT.read_bytes())
            raw = root/"raw"
            def load(run, unused):
                if run["run_id"] == O.RUN_IDS[-1]:
                    path = raw/O.RUN_IDS[0]/"cell_1/chunk_0/rows.jsonl"
                    path.parent.mkdir(parents=True); path.write_text("\n")
                return self.runs[run["run_id"]]
            with mock.patch.object(O, "load_compact_run", side_effect=load):
                with self.assertRaisesRegex(ValueError, "Inputs or source changed"):
                    O.analyze(protocol, raw, contract, PROJECT)

    def test_output_created_at_publication_cannot_overwrite_victim(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); raw = root/"raw"; raw.mkdir()
            victim = raw/"immutable.json"; victim.write_text("immutable")
            report = O.summarize_comparisons(self.runs)
            report.update(protected_snapshot_roots=[], input_sha256={str(victim): "fixture"})
            original = os.link
            def race(source, target, *args, **kwargs):
                Path(target).symlink_to(victim)
                return original(source, target, *args, **kwargs)
            with mock.patch.object(O.os, "link", side_effect=race), self.assertRaises(FileExistsError):
                O.write_outputs(root/"report", report, PROTOCOL, CONTRACT, raw, root)
            self.assertEqual(victim.read_text(), "immutable")


if __name__ == "__main__":
    unittest.main()
