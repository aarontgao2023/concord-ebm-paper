"""Small complete/incomplete fixtures for the post-HPC development analyses."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import run_record_audit as U


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, U.ROOT / relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


E = module("estimated_audit_fixture", "results/v2_dev_estimated_prior_a/analyze_estimated_prior.py")
B = module("bridge_audit_fixture", "results/v2_dev_paired_bridge_a/compare_bridge.py")


class DevelopmentAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def cfg(self, run, bridge=False):
        return {"run_id": run, "phase": "development", "base_seed": 100, "datasets": 2,
                "engines": ["original", "repaired"] if bridge else
                    ["repaired", "oracle_estimated_original", "oracle_estimated_repaired"],
                "cells": [{"name": "REF"}], "schemes": {"diagnosis": {"stratify": "diagnosis"}} if bridge else {},
                "bperm": 3 if bridge else 0, "full_p_first_n": 2, "paired_standard": bridge}

    def fit(self, cfg, seed, engine, scheme=None, pid=None):
        ordering = [1, 0] if engine == "oracle_estimated_repaired" else [0, 1]
        optimizer = {"method": "SLSQP", "success": True, "status": 0, "nit": 3,
                     "message": "success", "finite_x": True, "finite_objective": True}
        diagnostics = {"optimizer_calls": 1, "consensus": [{"mode": engine, "final_score": 1., "seconds": 2.}]}
        if scheme is None:
            diagnostics.update(optimizer=[optimizer], optimizer_failures=0)
        else:
            diagnostics.update(optimizer_unsuccessful=0, optimizer_status_counts={"0": 1})
        if engine.startswith("oracle_estimated"):
            prior = {"status": "interior", "identified": True, "abnormal_prior": .5,
                     "score_at_estimate": 0., "n_observed": 20, "n_missing": 0}
            diagnostics = {field: False for field in ("inference_eligible", "truth_order_used_in_prior_and_initialization",
                "individual_latent_stage_used", "diagnosis_used_in_prior_estimation", "truth_dx_counts_used")}
            diagnostics.update(prior_fits={g: {name: deepcopy(prior) for name in ("a", "b")} for g in U.GROUPS},
                               initial_orderings={g: [0, 1] for g in U.GROUPS})
        record = {"kind": "observed" if scheme is None else "permutation", "seed": seed,
                "engine": engine, "scheme": scheme, "perm_id": pid, "status": "ok",
                "orderings": [ordering] * 3, "taus": [0.] * 3,
                "distance_to_truth": [float(ordering != [1, 0])] * 3,
                "truth": {"manifest": {"data_sha256": str(seed)}, "biomarker_names": ["a", "b"],
                          "group_orderings": {g: [1, 0] for g in U.GROUPS}},
                "diagnostics": diagnostics, "fit_seconds": 2.}
        if cfg.get("paired_standard"):
            record["shared_fit_id"] = U.I.digest([seed, scheme, pid])
        return record

    def write_run(self, cfg, bridge=False, missing=None, mutate=None):
        frozen = self.base / cfg["run_id"] / "frozen"
        config_path = frozen / "configs/v2" / (cfg["run_id"] + ".json")
        source = frozen / "scripts/v2"
        config_path.parent.mkdir(parents=True)
        source.mkdir(parents=True)
        config_path.write_text(json.dumps(cfg))
        for name in U.I.RUNTIME_FILES:
            (source / name).write_text("def permute(x):\n    return x\n" if name == "run_v2.py" else "# fixture\n")
        hashes = {name: hashlib.sha256((source / name).read_bytes()).hexdigest() for name in U.I.RUNTIME_FILES}
        raw = self.base / cfg["run_id"] / "raw"
        chunk = raw / "cell_0/chunk_0"
        chunk.mkdir(parents=True)
        seeds = list(range(100, 102))
        records = [self.fit(cfg, seed, engine) for seed in seeds for engine in cfg["engines"]]
        if bridge:
            records += [self.fit(cfg, seed, engine, scheme, pid) for seed in seeds for engine in cfg["engines"]
                        for scheme in cfg["schemes"] for pid in range(3)]
        if missing:
            records = [r for r in records if not missing(r)]
        if mutate:
            mutate(records)
        (chunk / "fit_records.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
        identity = {"config": cfg, "cell_index": 0, "chunk": 0, "nchunks": 1, "source_sha256": hashes}
        environment = {"versions": {name: "1" for name in U.I.NUMERICAL_PACKAGES},
                       "source_sha256": {"dependency": "a" * 64}, "wheel_sha256": "b" * 64}
        (chunk / "manifest.json").write_text(json.dumps(identity | {
            "signature": U.I.digest(identity), "environment": environment, "seeds": seeds}))
        rows = []
        for r in records:
            if r["scheme"] is not None:
                continue
            rows.append({"run_id": cfg["run_id"], "phase": "development", "cell": "REF", "engine": r["engine"],
                         "seed": r["seed"], "status": "ok", "truth": r["truth"],
                         "observed": {name: r[name] for name in ("orderings", "taus", "distance_to_truth", "fit_seconds")},
                         "diagnostics": r["diagnostics"] | {"full_p_selected": True},
                         "schemes": {scheme: {"requested_nperm": cfg["bperm"], "nperm": 3, "attempted": 3,
                             "failed": 0, "complete": cfg["bperm"] == 3, "gt": [0] * 3, "eq": [3] * 3,
                             "definition": definition} for scheme, definition in cfg["schemes"].items()}})
        (chunk / "rows.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
        (chunk / "progress.json").write_text(json.dumps({"done": True, "rows": len(rows), "expected_rows": 2 * len(cfg["engines"]),
            "computational_jobs": len({r.get("shared_fit_id", str(U.identity("REF", r))) for r in records})}))
        return config_path, raw

    def test_estimated_complete_fixture_has_dataset_effects_and_exact_reference(self):
        config, raw = self.write_run(self.cfg("new"))
        reference, old = self.write_run(self.cfg("old") | {"engines": ["repaired"]})
        result = E.analyze(config, raw, reference, old, 20)
        self.assertTrue(result["ready"], result["issues"])
        effect = next(e for e in result["effects"] if e["first"] == "repaired" and e["second"] == "oracle_estimated_repaired"
                      and e["field"] == "distance_to_truth" and e["group_or_pair"] == "dataset_mean")
        self.assertEqual(effect["effect"]["difference_second_minus_first"]["estimate"], -1.)
        self.assertEqual(effect["effect"]["first"]["n"], 2)
        self.assertEqual(len(result["prior_fits"]), 24)

    def test_incomplete_estimated_cell_never_outputs_partial_cohort_effects(self):
        config, raw = self.write_run(self.cfg("new"), missing=lambda r: r["seed"] == 101)
        reference, old = self.write_run(self.cfg("old") | {"engines": ["repaired"]})
        result = E.analyze(config, raw, reference, old, 20)
        self.assertFalse(result["ready"])
        self.assertEqual(result["effects"], [])
        self.assertEqual(result["cohorts"][0]["planned_rows"], 6)

    def test_duplicate_fit_is_an_error_and_is_not_deduplicated_into_complete_pair(self):
        config, raw = self.write_run(self.cfg("new"), mutate=lambda rows: rows.append(deepcopy(rows[0])))
        bundle = U.load_run(config, raw)
        self.assertTrue(any(i["code"] == "duplicate_fit_identity" for i in bundle["issues"]))
        self.assertNotIn(("REF", 100, "repaired", None, None), bundle["records"])

    def test_bridge_compares_fixed_pid_subset_and_counts_shared_work_once(self):
        config, raw = self.write_run(self.cfg("new", True), bridge=True)
        reference, old = self.write_run(self.cfg("old", True) | {"bperm": 599, "paired_standard": False}, bridge=True)
        result = B.analyze(config, raw, reference, old)
        self.assertTrue(result["ready"], result["issues"])
        self.assertEqual(result["comparison_counts"], {"exact": 16})
        self.assertEqual(result["planned_computational_jobs"], 8)
        self.assertEqual(result["timing_summary"]["total_new_seconds"], 16.)
        self.assertEqual(result["timing_summary"]["total_reference_seconds"], 32.)
        self.assertFalse(result["is_calibration_experiment"])

    def test_bridge_missing_reference_does_not_become_exact(self):
        config, raw = self.write_run(self.cfg("new", True), bridge=True)
        reference, old = self.write_run(self.cfg("old", True) | {"bperm": 599, "paired_standard": False}, bridge=True,
                                        missing=lambda r: r["seed"] == 101)
        result = B.analyze(config, raw, reference, old)
        self.assertFalse(result["ready"])
        self.assertEqual(result["comparison_counts"], {"exact": 8, "missing": 8})

    def test_bridge_detects_ordering_difference_even_with_identical_taus(self):
        config, raw = self.write_run(self.cfg("new", True), bridge=True)
        def change(records):
            records[-1]["orderings"] = [[1, 0]] * 3
        reference, old = self.write_run(self.cfg("old", True) | {"bperm": 599, "paired_standard": False}, bridge=True, mutate=change)
        result = B.analyze(config, raw, reference, old)
        self.assertFalse(result["ready"])
        self.assertEqual(result["comparison_counts"]["different"], 1)

    def test_fast_bridge_counts_each_paired_side_once_and_ignores_only_seconds(self):
        def faster(records):
            for record in records:
                record["fit_seconds"] = 1.
                record["diagnostics"]["seconds"] = 999.
                record["diagnostics"]["fast_extra_metadata"] = True
                record["diagnostics"]["consensus"][0]["seconds"] = .1
        config, raw = self.write_run(self.cfg("fast", True) | {"fast_likelihood": True}, bridge=True, mutate=faster)
        reference, old = self.write_run(self.cfg("paired", True), bridge=True)
        result = B.analyze(config, raw, reference, old)
        self.assertTrue(result["ready"], result["issues"])
        self.assertTrue(result["optimizer_diagnostics_required"])
        self.assertEqual(result["matched_scope_computational_job_counts"], {"new": 8, "reference": 8})
        self.assertEqual(result["timing_summary"]["total_new_seconds"], 8.)
        self.assertEqual(result["timing_summary"]["total_reference_seconds"], 16.)
        self.assertEqual(result["reference_run_id"], "paired")
        self.assertEqual(result["reference_bperm"], 3)
        self.assertTrue(all("B599" not in line for line in result["limitations"]))
        self.assertTrue(all(c["optimizer_diagnostics"]["exact"] for c in result["comparisons"]))

    def test_fast_bridge_checks_full_optimizer_records_compact_counts_and_consensus(self):
        for field in ("optimizer_record", "optimizer_calls", "optimizer_unsuccessful", "optimizer_status_counts", "consensus"):
            with self.subTest(field=field):
                def alter(records):
                    record = records[0] if field == "optimizer_record" else records[-1]
                    d = record["diagnostics"]
                    if field == "optimizer_record":
                        d["optimizer"][0]["nit"] += 1
                    elif field == "optimizer_status_counts":
                        d[field] = {"1": 1}
                    elif field == "consensus":
                        d[field][0]["final_score"] += .01
                    else:
                        d[field] += 1
                config, raw = self.write_run(self.cfg("fast_" + field, True) | {"fast_likelihood": True}, bridge=True, mutate=alter)
                reference, old = self.write_run(self.cfg("reference_" + field, True), bridge=True)
                result = B.analyze(config, raw, reference, old)
                self.assertFalse(result["ready"])
                self.assertEqual(result["comparison_counts"]["different"], 1)
                self.assertTrue(all(c["core_successful_exact"] for c in result["comparisons"]))


if __name__ == "__main__":
    unittest.main()
