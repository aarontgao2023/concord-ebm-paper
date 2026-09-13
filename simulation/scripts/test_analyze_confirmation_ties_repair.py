"""Pure namespace31 count fixtures; never generate/fitted confirmation data."""
import copy
import itertools
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import analyze_confirmation_cohorts as C
import analyze_confirmation_ties_repair as T
from test_analyze_confirmation_primary import config as base_config, row as base_row, SEED
import test_confirmation_primary_integration as IntegrationFixture


def config(n=2, full=2):
    cfg = base_config(datasets=n, full=full)
    cfg["engines"] = ["original", "repaired"]
    return cfg


def row(cfg, offset=0, engine="original", *, u=False, d=False, nperm=599, eq_u=None, eq_d=None):
    r = base_row(cfg, offset, a=([0]*3 if u else [20]*3), b=([0]*3 if d else [20]*3),
                 max_a=0 if u else 30, max_b=0 if d else 30, nperm=nperm)
    r["engine"] = engine
    r["truth"]["manifest"] = {"data_sha256": "a"*64, "seed": r["seed"]}
    for scheme, ties in zip(T.SCHEMES, (eq_u, eq_d)):
        entry = r["schemes"][scheme]
        entry.update(attempted=nperm, failed=0)
        if ties is not None:
            entry.update(gt=[0]*3, eq=list(ties), maxgt=0, maxeq=max(ties))
    return r


def summarize(cfg, rows):
    return T.summarize_run(cfg, rows, bootstrap_repeats=20)


def metric(report, *, engine="original", scheme="unrestricted", name="e2_vs_e33"):
    return next(g for g in report["ties_by_scheme"] if g["engine"] == engine and g["scheme"] == scheme)["metrics"][name]


def contrast(report, name="U_minus_D_gap", kind="bonferroni"):
    return next(c for c in report["repair_rejection_contrasts"] if c["comparison"] == name and c["endpoint"] == kind)["summary"]


class TiesRepairTests(unittest.TestCase):
    def test_ties_use_E_over599_and_difference_is_D_minus_U(self):
        cfg = config(1, 1)
        rows = [row(cfg, engine=e, eq_u=[0, 299, 599], eq_d=[1, 300, 599]) for e in T.ENGINES]
        report = summarize(cfg, rows)
        self.assertEqual(metric(report)["estimate"], 0.)
        self.assertEqual(metric(report, name="global_max")["estimate"], 1.)
        self.assertAlmostEqual(metric(report, name="within_dataset_mean_pair")["estimate"], 898/(3*599))
        shift = report["tie_differences"][0]
        self.assertEqual(shift["direction"], "diagnosis minus unrestricted")
        self.assertAlmostEqual(shift["metrics"]["within_dataset_mean_pair"]["estimate"], 2/(3*599))
        self.assertEqual(report["planned_tie_datasets_per_cell"], 1)

    def test_dataset_unit_mean_distribution_and_MCSE(self):
        cfg = config()
        rows = [row(cfg, i, e, eq_u=[599*i]*3, eq_d=[0]*3) for i in range(2) for e in T.ENGINES]
        output = metric(summarize(cfg, rows), name="within_dataset_mean_pair")
        self.assertEqual(output["known_datasets"], 2)
        self.assertEqual(output["estimate"], .5)
        self.assertAlmostEqual(output["mcse_mean"], .5)
        self.assertEqual(output["distribution"]["quantiles"]["0.25"], .25)
        self.assertEqual(output["distribution"]["quantiles"]["0.75"], .75)

    def test_selected100_and_rejection1000_do_not_shrink_or_add_natural_completions(self):
        cfg = config(1000, 100)
        report = summarize(cfg, [row(cfg, i, e) for i in (0, 100) for e in T.ENGINES])
        ties, repair = metric(report), contrast(report)
        self.assertEqual((ties["planned_datasets"], ties["known_datasets"], ties["unknown_datasets"]), (100, 1, 99))
        self.assertEqual((repair["planned_datasets"], repair["known_datasets"]), (1000, 2))
        self.assertIsNone(ties["estimate"])
        self.assertIsNone(ties["distribution"])
        self.assertIsNone(ties["bootstrap_ci95"])
        self.assertEqual(len(report["selected_dataset_ties"]), 400)
        self.assertNotIn(SEED+100, {r["seed"] for r in report["selected_dataset_ties"]})

    def test_early_decision_is_usable_for_rejection_but_not_ties(self):
        cfg = config(1, 1)
        report = summarize(cfg, [row(cfg, engine=e, nperm=40) for e in T.ENGINES])
        self.assertEqual(contrast(report)["estimate"], 0)
        self.assertTrue(contrast(report)["endpoint_complete"])
        self.assertEqual(metric(report)["known_datasets"], 0)
        self.assertFalse(report["endpoint_complete"])
        self.assertEqual(metric(report)["all_planned_mean_bounds"], [0, 1])

    def test_terminal_failure_and_stale_fields_stay_unknown(self):
        cfg = config(1, 1)
        good, failed = row(cfg), row(cfg, engine="repaired", u=True, d=True)
        failed.update(status="error", truth={})
        report = summarize(cfg, [good, failed])
        self.assertEqual(contrast(report, "unrestricted")["all_planned_mean_bounds"], [-1, 0])
        self.assertIsNone(contrast(report, "unrestricted")["estimate"])
        self.assertEqual(metric(report, engine="repaired")["known_datasets"], 0)
        self.assertEqual(report["status_or_unavailable_reason_counts"]["REF_H0:repaired"], {"error": 1})

    def test_failed_permutation_retains_B599_even_when_decision_resolved(self):
        cfg = config(1, 1)
        rows = [row(cfg, engine=e, nperm=598) for e in T.ENGINES]
        for r in rows:
            for entry in r["schemes"].values():
                entry.update(failed=1, attempted=599)
        report = summarize(cfg, rows)
        self.assertEqual(contrast(report)["estimate"], 0)
        self.assertEqual(metric(report)["known_datasets"], 0)
        self.assertIsNone(metric(report)["distribution"])

    def test_DID_sign_support_and_separate_max_endpoint(self):
        cfg = config()
        rows = [row(cfg, 0, "original", u=True, d=False), row(cfg, 0, "repaired", u=False, d=True),
                row(cfg, 1, "original", u=False, d=True), row(cfg, 1, "repaired", u=True, d=False)]
        report = summarize(cfg, rows)
        result = contrast(report)
        self.assertEqual(result["support"], [-2, 2])
        self.assertEqual(result["estimate"], 0)
        self.assertAlmostEqual(result["mcse_mean"], 2)
        self.assertEqual(result["observed_joint_decision_counts"], {"1001": 1, "0110": 1})
        one = summarize(config(1, 1), rows[:2])
        self.assertEqual(contrast(one)["estimate"], 2)
        self.assertEqual(contrast(one, "unrestricted")["estimate"], 1)
        self.assertEqual(contrast(one, "diagnosis")["estimate"], -1)
        changed = copy.deepcopy(rows[:2])
        for r in changed:
            for entry in r["schemes"].values():
                entry["maxgt"] = 30
        self.assertEqual(contrast(summarize(config(1, 1), changed), kind="global_max")["estimate"], 0)

    def test_all_partial_DID_completion_bounds_are_sharp(self):
        for pattern in itertools.product((0, 1, None), repeat=4):
            values, bounds = T._linear([[p] for p in pattern], [1, -1, -1, 1])
            possible = [sum(c*v for c, v in zip([1, -1, -1, 1], filled))
                        for filled in itertools.product(*[(0, 1) if p is None else (p,) for p in pattern])]
            self.assertEqual(bounds[0], (min(possible), max(possible)))
            self.assertEqual(values[0], None if None in pattern else possible[0])

    def test_pairing_hash_mismatch_masks_both_without_dropping_dataset(self):
        cfg = config(1, 1)
        original, repaired = row(cfg), row(cfg, engine="repaired")
        repaired["truth"]["manifest"]["data_sha256"] = "b"*64
        output = summarize(cfg, [original, repaired])
        self.assertEqual(output["pairing_integrity_issues"][0]["code"], "paired_truth_or_data_mismatch")
        self.assertEqual(metric(output)["known_datasets"], 0)
        self.assertEqual(contrast(output)["all_planned_mean_bounds"], [-2, 2])

    def test_reordered_rows_and_inputs_unchanged_bootstrap_deterministic(self):
        cfg = config()
        rows = [row(cfg, i, e, u=i == 0, d=e == "original") for i in range(2) for e in T.ENGINES]
        before = copy.deepcopy((cfg, rows))
        self.assertEqual(summarize(cfg, rows), summarize(cfg, list(reversed(rows))))
        self.assertEqual((cfg, rows), before)

    def test_constant_sample_gets_nonzero_Hoeffding_not_zero_percentile(self):
        cfg = config()
        output = summarize(cfg, [row(cfg, i, e) for i in range(2) for e in T.ENGINES])
        for result in (metric(output), contrast(output)):
            self.assertTrue(result["bootstrap_degenerate_sample"])
            self.assertIsNone(result["bootstrap_ci95"])
            self.assertLess(*result["hoeffding_ci95"])

    def test_invalid_selection_duplicate_budget_or_counts_rejected(self):
        cfg = config(1, 1)
        r = row(cfg)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            summarize(cfg, [r, r])
        for label, mutation in [("selection", lambda x: x["diagnostics"].update(full_p_selected=False)),
                                 ("budget", lambda x: x["schemes"]["diagnosis"].update(requested_nperm=99)),
                                 ("accounting", lambda x: x["schemes"]["diagnosis"].update(attempted=598)),
                                 ("exceeds", lambda x: x["schemes"]["diagnosis"].update(eq=[600]*3))]:
            bad = copy.deepcopy(r)
            mutation(bad)
            with self.subTest(label=label), self.assertRaises(ValueError):
                summarize(cfg, [bad])
        with self.assertRaises(ValueError):
            T.summarize_run(cfg, [], bootstrap_repeats=0)

    def test_real_frozen_plan_empty_input_keeps_all_fixed_denominators(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = T.build_report(C.PROJECT/"runs/v2/confirmation_protocol.json", Path(tmp)/"absent", bootstrap_repeats=1)
        self.assertEqual(output["planned_rows"], 4000)
        self.assertEqual(output["selected_absolute_seeds"], list(range(42100000, 42100100)))
        self.assertEqual(output["status"], "collecting")
        self.assertFalse(output["collection_complete"])
        self.assertTrue(output["provenance"]["compact_rows"])
        self.assertEqual(metric(output)["planned_datasets"], 100)
        self.assertEqual(contrast(output)["planned_datasets"], 1000)
        self.assertIn("supplementary", T.render_markdown(output))
        json.dumps(output, allow_nan=False)

    def test_real_compact_loader_with_namespace31_fixture_preserves_failure_status(self):
        helper = IntegrationFixture.PrimaryLoaderIntegration()
        helper.setUp()
        try:
            helper.cfg["engines"] = list(T.ENGINES)
            helper.config_path.write_text(json.dumps(helper.cfg))
            rows = []
            for i in range(3):
                for engine in T.ENGINES:
                    r = helper.make_row(i, early=i == 2)
                    r["engine"] = engine
                    if i == 1 and engine == "repaired":
                        r.update(status="error", truth={}, observed=None, schemes={})
                    rows.append(r)
            helper.write_rows(rows)
            path = helper.results_root/helper.cfg["run_id"]/"cell_0/chunk_0/progress.json"
            progress = json.loads(path.read_text())
            progress["expected_rows"] = 6
            path.write_text(json.dumps(progress))
            original_validate = T._validate_config
            # Only frozen registry/finite plan adapter uses namespace31 toys;
            # real inventory, truth validation, compaction and loader execute.
            with patch.object(C, "frozen_run", side_effect=lambda *a: helper.frozen_toy_run()), \
                 patch.object(T, "_validate_config", side_effect=lambda cfg, **kw: original_validate(cfg)), \
                 patch.object(C, "load_run", wraps=C.load_run) as loader:
                report = T.build_report(helper.protocol, helper.results_root, bootstrap_repeats=5)
            self.assertIs(loader.call_args.kwargs["compact"], True)
            self.assertTrue(report["integrity_valid"])
            self.assertTrue(report["collection_complete"])
            self.assertEqual(report["status"], "terminal_with_unknown")
            self.assertEqual(metric(report, engine="repaired")["known_datasets"], 1)
            self.assertEqual(contrast(report)["known_datasets"], 2)
            self.assertIsNone(contrast(report)["estimate"])
        finally:
            helper.tearDown()

    def test_cli_refuses_output_hardlink_and_raw_before_reading_effects(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            raw = root/"raw"
            raw.mkdir()
            source = raw/"rows.jsonl"
            source.write_text("unchanged\n")
            alias = root/"alias.json"
            os.link(source, alias)
            protocol = root/"protocol.json"
            protocol.write_text(json.dumps({"runs": []}))
            unregistered_snapshot = root/"runs/v2/snapshots/old_development/report.json"
            relocated_snapshot = root/"relocated/old_snapshot/report.md"
            for output in (alias, source, unregistered_snapshot, relocated_snapshot):
                with patch.object(T, "build_report") as build, self.assertRaises(ValueError):
                    T.main(["--protocol", str(protocol), "--results-root", str(raw), "--project", str(root),
                            "--snapshots-root", str(root/"relocated"), "--output-json", str(output)])
                build.assert_not_called()
            self.assertEqual(source.read_text(), "unchanged\n")


if __name__ == "__main__":
    unittest.main()
