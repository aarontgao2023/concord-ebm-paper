"""Small hand-built development fixtures; no42 data, generator, fit or HPC reads."""
from copy import deepcopy
import hashlib
import json
import statistics
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import analyze_confirmation_mechanism as M
import run_record_audit as U

BASE = 31009601
ENVIRONMENT = {"versions": {name: "fixture1" for name in U.I.NUMERICAL_PACKAGES},
               "source_sha256": {"dependency": "a" * 64}, "wheel_sha256": "b" * 64}


class MechanismFixture:
    def __init__(self, root, *, cells=1, mutate=None, omit=None):
        self.root = Path(root)
        self.config = {"run_id": "mechanism_fixture", "phase": "development", "base_seed": BASE,
                       "datasets": 3, "engines": list(M.ENGINES), "schemes": {}, "bperm": 0,
                       "paired_standard": True, "cells": [{"name": "REF_H0", "id": f"FIXTURE_{i}"} for i in range(cells)]}
        self.config_path = self.root / "frozen/configs/v2/fixture.json"
        self.config_path.parent.mkdir(parents=True)
        self.config_path.write_text(json.dumps(self.config))
        source = self.root / "frozen/scripts/v2"
        source.mkdir(parents=True)
        for name in U.I.RUNTIME_FILES + U.I.OPTIONAL_RUNTIME_FILES:
            (source / name).write_text("# hand-built fixture source, never imported\n")
        (source / "design_v2.py").write_text("""DESIGN_VERSION = 'fixture1'
BIOMARKERS = (BiomarkerSpec('a', 'X', .8, 1), BiomarkerSpec('b', 'X', .8, 1), BiomarkerSpec('c', 'X', .8, 1))
class DesignConfig:
    name: str = 'REF_H0'
    diagnosis_mode: str = 'fixed'
    composition_lambda: float = 1.
    sample_scale: int = 1
    biomarker_names: tuple = ()
    fixed_base_order: object = None
    alternative_group: object = None
    target_inversions: int = 0
    stage_eta: float = 0.
    group_pdf_shift: dict = field(default_factory=dict)
raise RuntimeError('This fixture source must never be executed')
""")
        self.design_plans, _ = M._design_plans(self.config_path, self.config)
        self.hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir()}
        self.raw = self.root / "raw"
        for ci, cell in enumerate(self.config["cells"]):
            records = [self.record(U.cell_id(cell), seed, engine) for seed in range(BASE, BASE + 3) for engine in M.ENGINES]
            if mutate:
                mutate(records)
            if omit:
                records = [r for r in records if not omit(ci, r)]
            directory = self.raw / f"cell_{ci}/chunk_0"
            directory.mkdir(parents=True)
            identity = {"config": self.config, "cell_index": ci, "chunk": 0, "nchunks": 1, "source_sha256": self.hashes}
            environment = {"versions": {name: "fixture1" for name in U.I.NUMERICAL_PACKAGES},
                           "source_sha256": {"dependency": "a" * 64}, "wheel_sha256": "b" * 64}
            (directory / "manifest.json").write_text(json.dumps(identity | {"signature": U.I.digest(identity),
                "environment": environment, "seeds": list(range(BASE, BASE + 3))}))
            (directory / "fit_records.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
            rows = [{"run_id": self.config["run_id"], "phase": "development", "cell": U.cell_id(cell),
                     "seed": r["seed"], "engine": r["engine"], "status": r["status"], "truth": r["truth"],
                     "observed": {k: r.get(k) for k in ("orderings", "taus", "distance_to_truth")},
                     "diagnostics": r["diagnostics"], "schemes": {}} for r in records]
            (directory / "rows.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
            (directory / "progress.json").write_text(json.dumps({"done": True, "rows": len(rows), "expected_rows": 18}))

    def record(self, cell, seed, engine):
        truth_order = [2, 1, 0]
        order = {"original": [0, 1, 2], "repaired": [1, 0, 2] if seed % 2 else truth_order,
                 "oracle_equal_original": [0, 1, 2], "oracle_equal_repaired": [1, 0, 2],
                 "oracle_estimated_original": truth_order, "oracle_estimated_repaired": [2, 0, 1]}[engine]
        design = self.design_plans[cell]
        truth = {"biomarker_names": ["a", "b", "c"], "group_orderings": {g: truth_order for g in U.GROUPS},
                 "group_order": list(U.GROUPS), "seed": seed, "design_version": design["design_version"],
                 "config": design, "base_order": truth_order, "common_order_null": True,
                 "planned_inversions": dict.fromkeys(U.GROUPS, 0), "realized_inversions": dict.fromkeys(U.GROUPS, 0),
                 "pair_null": [True] * 3, "pair_distances": [0.] * 3, "h1_distance": 0.,
                 "pair_truth": {name: {"null": True, "inversions": 0, "normalized_distance": 0.} for name in M.PAIR_NAMES},
                 "manifest": {"data_sha256": U.I.digest([cell, seed]), "seed": seed,
                     "config_sha256": hashlib.sha256(json.dumps(design, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                     "source_sha256": self.hashes["design_v2.py"], "numpy_version": "fixture1", "pandas_version": "fixture1"}}

        optimizer = {"method": "SLSQP", "success": True, "status": 0, "nit": 3, "message": "success",
                     "finite_x": True, "finite_objective": True}
        d = {"quality_flags": [], "event_centers_finite": True, "provenance": deepcopy(ENVIRONMENT),
             "consensus": [{"mode": engine.rsplit("_", 1)[-1], "local_optimal": engine != "original",
                            "final_score": 1., "best_neighbor_score": .9 if engine == "original" else 1.,
                            "residual_improvement": .1 if engine == "original" else 0.} for _ in U.GROUPS]}
        record = {"kind": "observed", "seed": seed, "engine": engine, "scheme": None, "perm_id": None,
                  "status": "ok", "orderings": [order] * 3, "taus": [0.] * 3,
                  "distance_to_truth": [U.kendall(order, truth_order)] * 3, "truth": truth, "diagnostics": d}
        if engine in ("original", "repaired"):
            d.update(parameters_finite=True, optimizer=[deepcopy(optimizer)] * 2, optimizer_calls=2, optimizer_failures=0,
                     model_artifacts={"posterior_sha256": "c" * 64, "model_biomarkers": ["a", "b", "c"],
                                      "biomarker_parameters": [{"Control": [0], "Disease": [1], "Mixing": [.5]}]})
            record["shared_fit_id"] = U.I.digest([cell, seed, "paired"])
        else:
            estimated = engine.startswith("oracle_estimated")
            q = [.1, .5, .9] if estimated else [.5] * 3
            initial = truth_order if estimated else [0, 1, 2]
            d.update({flag: False for flag in ("inference_eligible", "truth_order_used_in_prior_and_initialization",
                                              "individual_latent_stage_used", "diagnosis_used_in_prior_estimation", "truth_dx_counts_used")})
            d.update(posterior_biomarkers=["a", "b", "c"], abnormal_priors={g: q for g in U.GROUPS},
                     normal_mixing_params={g: [1 - v for v in q] for g in U.GROUPS},
                     initial_orderings={g: initial for g in U.GROUPS},
                     oracle_prior_mode="estimated_group_prior" if estimated else "equal_prior")
            if estimated:
                d.update(prior_fit_calls=9, prior_boundary_fits=0, prior_unidentified_fits=0,
                         prior_fits={g: {name: {"abnormal_prior": value, "status": "interior", "identified": True,
                                                "n_observed": 10, "n_missing": 2}
                                         for name, value in zip(["a", "b", "c"], q)} for g in U.GROUPS})
        return record

    def analyze(self):
        bundle = U.load_run(self.config_path, self.raw)
        return M._summarize_bundle(bundle, repeats=40, bootstrap_seed=72)


def effect(result, first="original", second="repaired", cell="FIXTURE_0"):
    return next(e for e in result["effects"] if e["cell"] == cell and e["first"] == first and e["second"] == second
                and e["field"] == "distance_to_truth" and e["group_or_pair"] == "dataset_mean")


class MechanismTests(unittest.TestCase):
    @staticmethod
    def grid_values(counts, denominator=91):
        values = [count / denominator for count in counts]
        return values + [statistics.mean(values)]

    def test_integer_direction_grid_cancels_mathematically_zero_dataset_difference(self):
        first = self.grid_values([0, 0, 4])
        second = self.grid_values([0, 1, 3])
        self.assertNotEqual(first[3], second[3])  # equal total K, different rounding
        self.assertEqual(M._integer_grid_directions([first], [second], 14, 3), [0, 1, 0])
        self.assertEqual(M._integer_grid_directions([first], [second], 14, 1), [0, 0, 1])

    def test_integer_direction_grid_keeps_smallest_single_and_dataset_steps(self):
        zero, step = self.grid_values([0, 0, 0]), self.grid_values([1, 0, 0])
        for index in (0, 3):
            with self.subTest(index=index):
                self.assertEqual(M._integer_grid_directions([zero, step], [step, zero], 14, index), [1, 0, 1])

    def test_integer_direction_grid_preserves_missing_pairs_without_imputing_ties(self):
        zero, one = self.grid_values([0, 0, 0]), self.grid_values([91, 91, 91])
        self.assertEqual(M._integer_grid_directions([None, zero, one, None], [one, None, zero, None], 14, 3), [1, 0, 0])
        self.assertEqual(M._integer_grid_directions([None], [None], 14, 0), [0, 0, 0])

    def test_integer_direction_grid_rejects_off_grid_or_inconsistent_mean(self):
        good = self.grid_values([0, 1, 3])
        bad_grid = [.1, good[1], good[2], statistics.mean([.1, good[1], good[2]])]
        bad_mean = good[:3] + [good[3] + .001]
        for bad in (bad_grid, bad_mean, [float("nan")]*4, [False, 0., 0., 0.]):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                M._integer_grid_directions([bad], [None], 14, 3)

    def test_complete_mean_bounds_do_not_invert_from_count_cancellation(self):
        values = [0.24393772893772894] * 200
        lower, upper = M._bounded_mean(values)
        self.assertEqual(lower, upper)
        self.assertAlmostEqual(lower, values[0])
        self.assertEqual(M._bounded_mean([0.25, None]), [0.125, 0.625])

    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def fixture(self, **kwargs):
        return MechanismFixture(self.temp.name, **kwargs)

    def test_complete_fixture_uses_config_denominators_and_four_directional_contrasts(self):
        r = self.fixture().analyze()
        self.assertTrue(r["ready"], r["issues"])
        self.assertEqual((r["planned_cell_datasets"], r["planned_fit_records"]), (3, 18))
        self.assertEqual(len(r["effects"]), 4 * 8)
        self.assertTrue(all(e["estimable"] and e["planned_datasets"] == 3 for e in r["effects"]))
        self.assertTrue(all(e["direction_count_grid_denominator"] == (9 if e["group_or_pair"] == "dataset_mean" else 3) for e in r["effects"]))
        self.assertLess(effect(r)["estimate"], 0)
        self.assertAlmostEqual(effect(r, "oracle_estimated_original", "oracle_estimated_repaired")["estimate"], 1 / 3)
        self.assertTrue(all(x["all_six_hashes_present_and_equal"] for x in r["dataset_hash_checks"]))
        self.assertTrue(all(all(x["checks"].values()) for x in r["nuisance_pair_checks"]))
        self.assertFalse(r["fitting_or_generation_performed"])
        json.dumps(r, allow_nan=False)

    def test_incomplete_cell_blocks_all_bootstrap_but_other_complete_cell_remains_eligible(self):
        r = self.fixture(cells=2, omit=lambda ci, row: ci == 0 and row["engine"] == "oracle_equal_original" and row["seed"] == BASE).analyze()
        self.assertFalse(r["collection_complete"])
        self.assertIsNone(effect(r)["estimate"])
        self.assertEqual(effect(r)["known_paired_datasets"], 3)  # complete pair insufficient when cell unfinished
        self.assertIsNotNone(effect(r, cell="FIXTURE_1")["bootstrap_ci95"])
        self.assertEqual(r["cohorts"][0]["planned_records"], 18)

    def test_terminal_failure_preserves_denominator_and_marks_only_affected_endpoints_unestimable(self):
        def mutate(rows):
            row = next(r for r in rows if r["engine"] == "oracle_estimated_repaired" and r["seed"] == BASE)
            row.update(status="error", orderings=None, taus=None, distance_to_truth=None)
            row["diagnostics"] = {"error": "fixed timeout", "quality_flags": []}
        r = self.fixture(mutate=mutate).analyze()
        self.assertTrue(r["collection_complete"])
        self.assertTrue(effect(r)["estimable"])
        bad = effect(r, "repaired", "oracle_estimated_repaired")
        self.assertFalse(bad["estimable"])
        self.assertEqual((bad["planned_datasets"], bad["known_paired_datasets"], bad["unknown_paired_datasets"]), (3, 2, 1))
        self.assertIsNone(bad["bootstrap_ci95"])
        self.assertLess(bad["all_planned_difference_identification_bounds"][0], bad["all_planned_difference_identification_bounds"][1])

    def test_optimizer_quality_flags_do_not_select_successful_ordering_subset(self):
        def mutate(rows):
            for row in rows:
                if row["engine"] in ("original", "repaired"):
                    row["diagnostics"]["optimizer"][0]["success"] = False
                    row["diagnostics"]["optimizer_failures"] = 2  # fixture shares same optimizer object twice
                    row["diagnostics"]["quality_flags"] = ["optimizer_nonconvergence"]
        r = self.fixture(mutate=mutate).analyze()
        self.assertTrue(effect(r)["estimable"], r["issues"])
        summary = next(x for x in r["engine_summaries"] if x["engine"] == "original")
        self.assertEqual(summary["has_quality_flags"]["true"], 3)
        self.assertEqual(summary["optimizer_calls_recorded"], 6)
        self.assertEqual(summary["any_remaining_improving_neighbor"]["true"], 3)

    def test_data_hash_mismatch_blocks_effect_without_erasing_planned_records(self):
        def mutate(rows):
            rows[1]["truth"]["manifest"]["data_sha256"] = "f" * 64
        r = self.fixture(mutate=mutate).analyze()
        self.assertEqual(r["cohorts"][0]["present_records"], 18)
        self.assertFalse(effect(r)["estimable"])
        self.assertTrue(any(i["code"] == "within_dataset_data_hash_mismatch" for i in r["issues"]))

    def test_optimizer_history_mismatch_is_an_integrity_issue(self):
        def mutate(rows):
            rows[1]["diagnostics"]["optimizer"][0]["nit"] += 1
        r = self.fixture(mutate=mutate).analyze()
        self.assertFalse(effect(r)["estimable"])
        self.assertTrue(any(i["code"] == "consensus_pair_nuisance_audit_mismatch" for i in r["issues"]))

    def test_initial_match_truth_and_q_counts_use_group_units_without_inflating_R(self):
        r = self.fixture().analyze()
        oracle = [x for x in r["oracle_diagnostics"] if x["engine"] == "oracle_estimated_original"]
        self.assertEqual(len(oracle), 3)
        self.assertTrue(all(all(x["initial_matches_truth"].values()) for x in oracle))
        self.assertEqual(sum(x["n_prior_fits"] for x in oracle), 27)
        self.assertEqual(effect(r)["planned_datasets"], 3)

    def test_q_boundary_values_and_provenance_declarations_are_checked(self):
        def mutate(rows):
            row = next(x for x in rows if x["engine"] == "oracle_estimated_original")
            row["diagnostics"]["truth_order_used_in_prior_and_initialization"] = True
            row["diagnostics"]["prior_fits"]["e2"]["a"]["status"] = "lower_boundary"
        r = self.fixture(mutate=mutate).analyze()
        self.assertTrue(any(i["code"] == "oracle_diagnostic_integrity" for i in r["issues"]))
        self.assertFalse(effect(r)["estimable"])

    def test_saved_distance_must_match_legal_order_arithmetic(self):
        def mutate(rows):
            rows[0]["distance_to_truth"] = [.2] * 3
        r = self.fixture(mutate=mutate).analyze()
        self.assertTrue(any(i["code"] == "invalid_observed_endpoint" for i in r["issues"]))
        self.assertFalse(effect(r)["estimable"])
        self.assertEqual(effect(r)["unknown_paired_datasets"], 1)

    def test_derived_rows_cannot_disagree_with_raw_records(self):
        f = self.fixture()
        path = f.raw / "cell_0/chunk_0/rows.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        rows[0]["observed"]["distance_to_truth"] = [.5] * 3
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        r = f.analyze()
        self.assertTrue(any(i["code"] == "derived_row_raw_record_mismatch" for i in r["issues"]))
        self.assertFalse(effect(r)["estimable"])

    def test_uncommitted_raw_tail_is_not_a_complete_inference_input(self):
        f = self.fixture()
        path = f.raw / "cell_0/chunk_0/fit_records.jsonl"
        path.write_bytes(path.read_bytes().rstrip(b"\n"))
        before = path.read_bytes()
        r = f.analyze()
        self.assertTrue(any(i["code"] == "uncommitted_raw_log_tail" for i in r["issues"]))
        self.assertFalse(effect(r)["estimable"])
        self.assertEqual(before, path.read_bytes())

    def test_valid_q_boundaries_are_reported_without_filtering(self):
        def mutate(rows):
            for row in rows:
                if row["engine"].startswith("oracle_estimated"):
                    d = row["diagnostics"]
                    for g in U.GROUPS:
                        d["abnormal_priors"][g][0] = 1e-6
                        d["normal_mixing_params"][g][0] = 1 - 1e-6
                        d["prior_fits"][g]["a"].update(abnormal_prior=1e-6, status="lower_boundary")
                    d["prior_boundary_fits"] = 3
        r = self.fixture(mutate=mutate).analyze()
        self.assertTrue(r["ready"], r["issues"])
        summary = next(x for x in r["oracle_summaries"] if x["engine"] == "oracle_estimated_original")
        self.assertEqual(summary["q_boundary_fits_recorded"], 9)
        self.assertTrue(effect(r, "repaired", "oracle_estimated_repaired")["estimable"])

    def test_nonfinite_consensus_audit_is_reportable_json_and_blocks_CI(self):
        def mutate(rows):
            rows[0]["diagnostics"]["consensus"][0]["final_score"] = float("nan")
        r = self.fixture(mutate=mutate).analyze()
        self.assertFalse(effect(r)["estimable"])
        self.assertTrue(any(i["code"] == "incomplete_consensus_optimality_audit" for i in r["issues"]))
        json.dumps(r, allow_nan=False)

    def test_shared_dataset_bootstrap_is_deterministic_and_preserves_joint_columns(self):
        a = M._bootstrap_columns({"x": [0., 1., 4.], "2x": [0., 2., 8.]}, 60, 72, "cell")
        b = M._bootstrap_columns({"2x": [0., 2., 8.], "x": [0., 1., 4.]}, 60, 72, "cell")
        self.assertEqual(a, b)
        self.assertEqual(a["2x"], [2 * x for x in a["x"]])

    def test_six_matching_but_wrong_truth_records_do_not_pass_frozen_design_checks(self):
        def alter_config(row):
            row["truth"]["config"]["sample_scale"] = 4
        mutations = {
            "config": alter_config,
            "config_hash": lambda row: row["truth"]["manifest"].update(config_sha256="f" * 64),
            "source_hash": lambda row: row["truth"]["manifest"].update(source_sha256="f" * 64),
            "truth_seed": lambda row: row["truth"].update(seed=row["seed"] + 1),
            "generator_environment": lambda row: row["truth"]["manifest"].update(numpy_version="wrong"),
            "true_null_plan": lambda row: row["truth"].update(pair_null=[False] * 3),
        }
        for label, change in mutations.items():
            with self.subTest(label=label):
                f = MechanismFixture(Path(self.temp.name) / label, mutate=lambda rows: [change(row) for row in rows])
                r = f.analyze()
                self.assertFalse(effect(r)["estimable"])
                self.assertTrue(any(i["code"] == "truth_frozen_design_or_provenance_mismatch" for i in r["issues"]))
                self.assertEqual(effect(r)["planned_datasets"], 3)

    def test_output_guard_preserves_raw_protocol_and_all_frozen_snapshots(self):
        root = Path(self.temp.name)
        protocol = root / "confirmation_mechanism_analysis.json"
        protocol.write_text(json.dumps({"runs": [{"snapshot": "archive/frozen_mechanism"}]}))
        raw = root / "raw"
        for directory in (raw / "nested", root, root / "runs/v2/snapshots/unregistered",
                          root / "archive/frozen_mechanism/notes"):
            with self.subTest(directory=directory), self.assertRaises(ValueError):
                M._output_directory(directory, protocol, raw, root)
        self.assertEqual(M._output_directory(root / "analysis", protocol, raw, root), (root / "analysis").resolve())

    def test_public_protocol_reader_confirms_R200_four_cells_six_arms_without_raw_reads(self):
        root = Path(__file__).resolve().parents[2]
        path, _ = M._protocol_config(root / "runs/v2/confirmation_protocol.json", root)
        cfg = json.loads(path.read_text())
        self.assertEqual(cfg["datasets"] * len(cfg["cells"]) * len(cfg["engines"]), 4800)
        self.assertEqual(tuple(cfg["engines"]), M.ENGINES)
        plans, generator_hash = M._design_plans(path, cfg)
        self.assertEqual(plans["REF_H0_N4"]["sample_scale"], 4)
        self.assertEqual(plans["IID_H0"]["diagnosis_mode"], "iid")
        self.assertEqual(plans["STAGE_H0"]["stage_eta"], 1.)
        self.assertTrue(all(len(plan["biomarker_names"]) == 14 for plan in plans.values()))
        self.assertTrue(U.I.is_hash(generator_hash))


if __name__ == "__main__":
    unittest.main()
