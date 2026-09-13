"""Independent primary/real-loader integration with invented counts only.

All saved toy rows use seed 31009700 onwards. The empty-results test reads only
the actual frozen configuration/source metadata; it never opens actual run data.
No generator, estimator, or scientific-data fixture is imported or executed.
"""
from copy import deepcopy
import itertools
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import analyze_confirmation_cohorts as C
import analyze_confirmation_primary as P
from test_analyze_confirmation_primary import config, row, SEED


class PrimaryLoaderIntegration(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.snapshot = self.root / "toy_snapshot"
        self.source = self.snapshot / "scripts/v2"
        self.source.mkdir(parents=True)
        for name in C.I.RUNTIME_FILES:
            (self.source / name).write_text(f"# Never executed toy source: {name}\n")
        self.sources = {name: C.sha(self.source / name) for name in C.I.RUNTIME_FILES}
        self.cfg = config(datasets=3, full=2)
        self.config_path = self.snapshot / "configs/v2/toy.json"
        self.config_path.parent.mkdir(parents=True)
        self.config_path.write_text(json.dumps(self.cfg))
        self.design = C.resolved_design(C.design_defaults(C.PROJECT / "scripts/v2/design_v2.py"), {"name": "REF_H0"})
        self.environment = {"versions": {name: "toy1" for name in C.I.NUMERICAL_PACKAGES},
                            "source_sha256": {"upstream.py": "a" * 64}, "wheel_sha256": "b" * 64}
        self.protocol = self.root / "toy_protocol.json"
        self.protocol.write_text(json.dumps({"runs": [{"run_id": self.cfg["run_id"],
            "snapshot_sha256": "c" * 64, "config_sha256": C.sha(self.config_path)}]}))
        self.results_root = self.root / "toy_results"

    def tearDown(self):
        self.temporary.cleanup()

    def frozen_toy_run(self):
        return {"run_id": self.cfg["run_id"], "cfg": deepcopy(self.cfg), "config_path": self.config_path,
                "snapshot": self.snapshot, "sources": self.sources,
                "designs": {"REF_H0": self.design}, "names": ["REF_H0"],
                "seeds": set(range(SEED, SEED + 3)), "full_p": {SEED, SEED + 1},
                "dependency": {"PYEBM_VERSION": "toy1", "SOURCE_SHA256": self.environment["source_sha256"],
                               "PYEBM_WHEEL_SHA256": self.environment["wheel_sha256"]},
                "records": {}, "issues": [], "environment": None}

    def make_row(self, offset, early=False):
        r = row(self.cfg, offset, a=(10, 20, 30) if early else (0, 100, 100),
                b=(10, 20, 30) if early else (10, 100, 100),
                nperm=40 if early else 599, max_a=30, max_b=30)
        r["phase"] = self.cfg["phase"]
        order = list(range(14))
        r["truth"] = {"config": self.design, "seed": r["seed"], "manifest": {
            "seed": r["seed"], "config_sha256": C.canonical_sha(self.design),
            "source_sha256": self.sources["design_v2.py"], "data_sha256": "d" * 64,
            "numpy_version": "toy1", "pandas_version": "toy1"},
            "biomarker_names": self.design["biomarker_names"], "group_order": list(C.GROUPS),
            "base_order": order, "group_orderings": {g: order[:] for g in C.GROUPS},
            "planned_inversions": {g: 0 for g in C.GROUPS}, "realized_inversions": {g: 0 for g in C.GROUPS},
            "pair_null": [True] * 3, "h1_distance": 0.,
            "pair_truth": {p: {"null": True, "inversions": 0, "normalized_distance": 0.} for p in C.PAIR_NAMES},
            "latent_stage": [0, 1, 2]}
        r["observed"] = {"orderings": [order[:] for _ in C.GROUPS], "taus": [0.] * 3, "distance_to_truth": [0.] * 3}
        r["diagnostics"].update(provenance=deepcopy(self.environment), optimizer={"history": ["toy unused detail"]})
        for entry in r["schemes"].values():
            entry.update(attempted=entry["nperm"], failed=0)
        return r

    def write_rows(self, rows):
        path = self.results_root / self.cfg["run_id"] / "cell_0/chunk_0"
        path.mkdir(parents=True, exist_ok=True)
        identity = {"config": self.cfg, "cell_index": 0, "chunk": 0, "nchunks": 1, "source_sha256": self.sources}
        manifest = identity | {"signature": C.I.digest(identity), "environment": self.environment,
                               "seeds": list(range(SEED, SEED + 3))}
        (path / "manifest.json").write_text(json.dumps(manifest))
        (path / "rows.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
        (path / "progress.json").write_text(json.dumps({"done": True, "stopped": False, "rows": len(rows), "expected_rows": 3}))

    def build_toy_report(self):
        # Only the frozen-registry adapter is replaced for namespace31 fixtures;
        # the real inventory, truth checks, load_run and get_record all execute.
        with patch.object(C, "frozen_run", side_effect=lambda *args: self.frozen_toy_run()):
            return P.build_report(self.protocol, self.results_root, run_ids=[self.cfg["run_id"]], bootstrap_repeats=20)

    def test_compaction_preserves_primary_statistics_through_real_loader(self):
        self.write_rows([self.make_row(i, early=i == 2) for i in range(3)])
        summaries = []
        for compact in (False, True):
            run = C.load_run(self.frozen_toy_run(), self.results_root, compact=compact)
            self.assertTrue(run["trusted"])
            self.assertTrue(run["inventory"]["ready_for_final_analysis"])
            rows = [C.get_record(run, "REF_H0", "repaired", seed)[0] for seed in sorted(run["seeds"])]
            self.assertEqual("latent_stage" in rows[0]["truth"], not compact)
            self.assertEqual("optimizer" in rows[0]["diagnostics"], not compact)
            summaries.append(P.summarize_run(self.cfg, rows, bootstrap_repeats=20))
        self.assertEqual(*summaries)
        with patch.object(C, "load_run", wraps=C.load_run) as loader:
            output = self.build_toy_report()["runs"][0]
        self.assertIs(loader.call_args.kwargs["compact"], True)
        self.assertEqual(output["groups"], summaries[0]["groups"])
        self.assertTrue(output["collection_complete"])
        self.assertTrue(output["endpoint_complete"])
        self.assertEqual(output["status"], "complete")

    def test_terminal_observed_failure_retains_denominators_and_unestimated_mean(self):
        rows = [self.make_row(i, early=i == 2) for i in range(3)]
        rows[1].update(status="error", truth={}, observed=None, schemes={})
        self.write_rows(rows)
        loaded = C.load_run(self.frozen_toy_run(), self.results_root, compact=True)
        self.assertTrue(loaded["inventory"]["collection_complete"])
        self.assertFalse(loaded["inventory"]["ready_for_final_analysis"])
        report = self.build_toy_report()
        run = report["runs"][0]
        self.assertTrue(run["collection_complete"])
        self.assertTrue(run["provenance"]["inventory_collection_complete"])
        self.assertFalse(run["endpoint_complete"])
        self.assertEqual(run["status"], "terminal_with_unknown")
        self.assertIn("Collection is finished; unknown endpoints", P.render_markdown(report))
        output = run["groups"][0]
        effect = output["primary_bonferroni_rejection_difference"]
        self.assertEqual((effect["planned_datasets"], effect["known_datasets"], effect["unknown_datasets"]), (3, 2, 1))
        self.assertEqual(effect["all_planned_mean_bounds"], [0., 2 / 3])
        self.assertIsNone(effect["estimate"])
        p = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertEqual((p["planned_datasets"], p["known_datasets"]), (2, 1))
        self.assertIsNone(p["estimate"])
        self.assertEqual(output["status_or_unavailable_reason_counts"], {"ok": 2, "observed_error": 1})

    def test_inventory_rejects_invalid_attempted_counts_before_primary_summary(self):
        rows = [self.make_row(i, early=i == 2) for i in range(3)]
        rows[0]["schemes"]["diagnosis"]["attempted"] = 598  # nperm is 599
        self.write_rows(rows)
        output = self.build_toy_report()["runs"][0]
        self.assertFalse(output["provenance"]["trusted_run"])
        self.assertFalse(output["endpoint_complete"])
        self.assertEqual(output["status"], "integrity_invalid")
        effect = output["groups"][0]["primary_bonferroni_rejection_difference"]
        self.assertEqual(effect["unknown_datasets"], 3)
        self.assertEqual(effect["all_planned_mean_bounds"], [-1., 1.])
        self.assertIsNone(effect["estimate"])

    def test_failed_permutation_keeps_fixed_budget_and_can_leave_terminal_unknown(self):
        rows = [self.make_row(i, early=i == 2) for i in range(3)]
        for name, count, maximum in (("unrestricted", 9, 29), ("diagnosis", 10, 30)):
            rows[0]["schemes"][name].update(nperm=598, attempted=599, failed=1, complete=False,
                                            gt=[count] * 3, maxgt=maximum)
        self.write_rows(rows)
        loaded = C.load_run(self.frozen_toy_run(), self.results_root, compact=True)
        self.assertTrue(loaded["trusted"])
        self.assertTrue(loaded["inventory"]["collection_complete"])
        self.assertFalse(loaded["inventory"]["ready_for_final_analysis"])
        run = self.build_toy_report()["runs"][0]
        self.assertTrue(run["collection_complete"])
        self.assertFalse(run["endpoint_complete"])
        self.assertEqual(run["status"], "terminal_with_unknown")
        output = run["groups"][0]
        effect = output["primary_bonferroni_rejection_difference"]
        self.assertEqual((effect["planned_datasets"], effect["known_datasets"]), (3, 2))
        self.assertEqual(effect["all_planned_mean_bounds"], [1 / 3, 2 / 3])
        self.assertIsNone(effect["estimate"])
        selected = output["preselected_p"]["metrics"]["mean_pair_p_shift"]
        self.assertEqual((selected["planned_datasets"], selected["known_datasets"]), (2, 1))
        self.assertIsNone(selected["estimate"])
        # Dataset0 shift bounds are [0,2/600], not recomputed with denominator599.
        self.assertAlmostEqual(selected["all_planned_mean_bounds"][0], (1 / 180) / 2)
        self.assertAlmostEqual(selected["all_planned_mean_bounds"][1], (1 / 180 + 2 / 600) / 2)

    def test_terminal_full_p_failure_does_not_hide_complete_rejection_endpoint(self):
        rows = [self.make_row(i, early=i == 2) for i in range(3)]
        for name, count in (("unrestricted", 20), ("diagnosis", 21)):
            rows[0]["schemes"][name].update(nperm=598, attempted=599, failed=1, complete=False,
                                            gt=[count] * 3, maxgt=30)
        self.write_rows(rows)
        report = self.build_toy_report()
        run = report["runs"][0]
        self.assertTrue(run["collection_complete"])
        self.assertFalse(run["endpoint_complete"])
        self.assertEqual(run["status"], "terminal_with_unknown")
        group = run["groups"][0]
        self.assertTrue(group["endpoint_completion"]["primary_bonferroni_rejection_difference"])
        self.assertFalse(group["endpoint_completion"]["preselected_mean_pair_p_shift"])
        self.assertIsNotNone(group["primary_bonferroni_rejection_difference"]["estimate"])
        self.assertIsNone(group["preselected_p"]["metrics"]["mean_pair_p_shift"]["estimate"])
        markdown = P.render_markdown(report)
        self.assertIn("Rejection endpoint complete", markdown)
        self.assertIn("p endpoint complete", markdown)

    def test_real_frozen_registry_and_absent_results_keep_1000_and_100_or_10(self):
        output = P.build_report(C.PROJECT / "runs/v2/confirmation_protocol.json", self.root / "absent", bootstrap_repeats=1)
        self.assertEqual(sum(run["planned_row_count"] for run in output["runs"]), 5000)
        for run in output["runs"]:
            self.assertFalse(run["collection_complete"])
            self.assertFalse(run["endpoint_complete"])
            self.assertEqual(run["status"], "collecting")
            for group in run["groups"]:
                effect = group["primary_bonferroni_rejection_difference"]
                self.assertEqual((effect["planned_datasets"], effect["known_datasets"]), (1000, 0))
                self.assertEqual(effect["all_planned_mean_bounds"], [-1., 1.])
                self.assertIsNone(effect["estimate"])
                selected = group["preselected_p"]["metrics"]["mean_pair_p_shift"]
                self.assertEqual(selected["planned_datasets"], 100 if run["run_id"] == "confirm_core_a" else 10)
                self.assertEqual(selected["known_datasets"], 0)
        self.assertIn("Results are still missing or unfinished.", P.render_markdown(output))
        json.dumps(output, allow_nan=False)

    def test_all_binary_pair_completion_bounds_are_sharp(self):
        # Two dataset pairs, independently missing or observed in either scheme.
        for pattern in itertools.product((False, True, None), repeat=4):
            positive, negative = pattern[:2], pattern[2:]
            result = P._difference(positive, negative, positive_name="U", negative_name="D",
                                   context="independent toy exhaustive", repeats=1, bootstrap_seed=1)
            completions = []
            for candidate in itertools.product(*[(False, True) if v is None else (v,) for v in pattern]):
                completions.append(sum(int(candidate[i]) - int(candidate[i + 2]) for i in range(2)) / 2)
            self.assertEqual(result["all_planned_mean_bounds"], [min(completions), max(completions)])
            if None in pattern:
                self.assertIsNone(result["estimate"])


if __name__ == "__main__":
    unittest.main()
