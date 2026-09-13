"""Small fixtures for denominator, provenance, and exact-completion invariants."""
import copy
from contextlib import redirect_stderr, redirect_stdout
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import inventory_run as I


class InventoryFixtures(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.base = Path(self.temporary.name)
        self.frozen = self.base / "frozen"
        self.config_path = self.frozen / "configs/v2/fixture.json"
        self.config_path.parent.mkdir(parents=True)
        self.source_root = self.frozen / "scripts/v2"
        self.source_root.mkdir(parents=True)
        for name in I.RUNTIME_FILES:
            (self.source_root / name).write_text(f"# frozen fixture {name}\n")
        self.hashes, _ = I.source_expectations(self.config_path, None)
        self.root = self.base / "run"
        self.config = {"run_id": "fixture", "phase": "development", "base_seed": 100,
                       "datasets": 4, "engines": ["original", "repaired"],
                       "cells": [{"name": "REF"}], "schemes": {}}
        self.save_config()

    def tearDown(self):
        self.temporary.cleanup()

    def save_config(self):
        self.config_path.write_text(json.dumps(self.config))

    def environment(self):
        return {"versions": {name: "1.0" for name in I.NUMERICAL_PACKAGES},
                "source_sha256": {"upstream.py": "a" * 64}, "wheel_sha256": "b" * 64,
                "pyebm_path": "/machine/environment/pyebm"}

    def make_row(self, seed, engine):
        return {"run_id": self.config["run_id"], "phase": self.config["phase"],
                "cell": "REF", "seed": seed, "engine": engine, "status": "ok",
                "truth": {"manifest": {"data_sha256": "c" * 64}},
                "observed": {"orderings": [[0, 1, 2]] * 3, "taus": [0, 0, 0]},
                "diagnostics": {"full_p_selected": seed in I.validate_config(self.config)[2]},
                "schemes": {}}

    def write_chunk(self, chunk=0, nchunks=1, rows=None, environment=None):
        directory = self.root / "cell_0" / f"chunk_{chunk}"
        directory.mkdir(parents=True, exist_ok=True)
        seeds = list(range(self.config["base_seed"], self.config["base_seed"] + self.config["datasets"]))[chunk::nchunks]
        identity = {"config": copy.deepcopy(self.config), "cell_index": 0, "chunk": chunk,
                    "nchunks": nchunks, "source_sha256": self.hashes}
        manifest = identity | {"signature": I.digest(identity), "environment": environment or self.environment(), "seeds": seeds}
        (directory / "manifest.json").write_text(json.dumps(manifest))
        if rows is None:
            rows = [self.make_row(seed, engine) for seed in seeds for engine in self.config["engines"]]
        (directory / "rows.jsonl").write_text("".join(json.dumps(row) + "\n" for row in rows))
        (directory / "progress.json").write_text(json.dumps({"done": True, "stopped": False,
            "rows": len(rows), "expected_rows": len(seeds) * len(self.config["engines"])}))
        return directory

    def inspect(self):
        return I.inventory(self.config_path, self.root)

    def configure_permutations(self, full=False, definition=None):
        self.config.update(datasets=1, engines=["original"], bperm=599, alpha=.05, rule="le",
                           schemes={"diagnosis": definition or {"stratify": "diagnosis", "require_max": True}},
                           full_p_first_n=int(full))
        self.save_config()

    def permutation_row(self, n=40, attempted=None, gt=None, maxgt=30):
        row = self.make_row(100, "original")
        row["schemes"] = {"diagnosis": {"requested_nperm": 599, "nperm": n,
             "attempted": n if attempted is None else attempted, "complete": n == 599,
             "failed": 0 if attempted is None else attempted - n,
             "gt": gt or [20, 20, 20], "eq": [0, 0, 0], "maxgt": maxgt, "maxeq": 0,
             "definition": self.config["schemes"]["diagnosis"]}}
        return row

    def test_complete_mechanism_run_uses_whole_plan(self):
        self.write_chunk(0, 2)
        self.write_chunk(1, 2)
        result = self.inspect()
        self.assertTrue(result["complete"])
        self.assertTrue(result["ready_for_final_analysis"])
        self.assertEqual(result["plan"]["planned_rows"], 8)
        self.assertEqual(result["totals"]["observed_ok"], 8)
        self.assertEqual(result["totals"]["missing_planned_rows"], 0)
        self.assertEqual(result["manifest_coverage"][0]["manifest_seed_union_count"], 4)

    def test_done_chunk_cannot_hide_unstarted_seeds(self):
        self.write_chunk(0, 2)
        result = self.inspect()
        self.assertFalse(result["complete"])
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertEqual(result["plan"]["planned_rows"], 8)
        self.assertEqual(result["totals"]["missing_planned_rows"], 4)
        self.assertEqual(result["manifest_coverage"][0]["missing_manifest_seeds"], [101, 103])
        self.assertEqual(result["manifest_coverage"][0]["missing_chunk_indices"], [1])

    def test_done_flag_cannot_hide_missing_engine_rows(self):
        rows = [self.make_row(seed, "original") for seed in range(100, 104)]
        self.write_chunk(rows=rows)
        result = self.inspect()
        self.assertFalse(result["complete"])
        repaired = next(r for r in result["cell_engines"] if r["engine"] == "repaired")
        self.assertEqual(repaired["missing_seeds"], [100, 101, 102, 103])

    def test_duplicate_identity_is_not_silently_deduplicated(self):
        rows = [self.make_row(seed, engine) for seed in range(100, 104) for engine in self.config["engines"]]
        rows.append(copy.deepcopy(rows[0]))
        self.write_chunk(rows=rows)
        result = self.inspect()
        self.assertFalse(result["complete"])
        self.assertEqual(result["totals"]["duplicate_planned_identities"], 1)
        self.assertEqual(result["totals"]["observed_ok"], 7)

    def test_signature_is_recomputed(self):
        path = self.write_chunk() / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest["signature"] = "0" * 64
        path.write_text(json.dumps(manifest))
        result = self.inspect()
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertTrue(any("signature fails" in issue["message"] for issue in result["issues"]))

    def test_manifest_seed_list_is_verified_independently_of_signature(self):
        path = self.write_chunk(0, 2) / "manifest.json"
        self.write_chunk(1, 2)
        manifest = json.loads(path.read_text())
        # Seeds are outside the runner's signed identity; verify the partition.
        manifest["seeds"] = [100, 101, 102, 103]
        path.write_text(json.dumps(manifest))
        result = self.inspect()
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertTrue(any("strided chunk allocation" in issue["message"] for issue in result["issues"]))

    def test_numerical_environment_difference_blocks_readiness(self):
        self.write_chunk(0, 2)
        environment = self.environment()
        environment["versions"]["numpy"] = "different"
        self.write_chunk(1, 2, environment=environment)
        result = self.inspect()
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertTrue(any(i["code"] == "environment_inconsistency" for i in result["issues"]))

    def test_installation_path_difference_is_not_numerical_difference(self):
        self.write_chunk(0, 2)
        environment = self.environment()
        environment["pyebm_path"] = "/another/machine/pyebm"
        self.write_chunk(1, 2, environment=environment)
        self.assertTrue(self.inspect()["ready_for_final_analysis"])

    def test_frozen_source_hash_mismatch_blocks_readiness(self):
        path = self.write_chunk() / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest["source_sha256"]["run_v2.py"] = "d" * 64
        identity = {name: manifest[name] for name in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
        manifest["signature"] = I.digest(identity)
        path.write_text(json.dumps(manifest))
        result = self.inspect()
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertTrue(any(i["code"] == "frozen_source_mismatch" for i in result["issues"]))

    def test_optional_runtime_sources_are_checked_when_frozen_files_exist(self):
        for name in I.OPTIONAL_RUNTIME_FILES:
            (self.source_root / name).write_text(f"# frozen optional module {name}\n")
        self.hashes, _ = I.source_expectations(self.config_path, None)
        self.assertEqual(set(self.hashes), set(I.RUNTIME_FILES + I.OPTIONAL_RUNTIME_FILES))
        self.write_chunk()
        self.assertTrue(self.inspect()["ready_for_final_analysis"])
        for name in I.OPTIONAL_RUNTIME_FILES:
            with self.subTest(name=name):
                original = (self.source_root / name).read_text()
                (self.source_root / name).write_text("# different frozen source\n")
                result = self.inspect()
                self.assertFalse(result["ready_for_final_analysis"])
                self.assertTrue(any(i["code"] == "frozen_source_mismatch" for i in result["issues"]))
                (self.source_root / name).write_text(original)

    def test_manifest_cannot_omit_a_present_optional_frozen_module(self):
        for name in I.OPTIONAL_RUNTIME_FILES:
            with self.subTest(name=name):
                # The manifest remains a correctly signed historical five-file
                # profile, but the supplied snapshot includes a sixth source.
                (self.source_root / name).write_text("# newly required optional source\n")
                self.write_chunk()
                result = self.inspect()
                self.assertFalse(result["ready_for_final_analysis"])
                self.assertTrue(any(i["code"] == "frozen_source_mismatch" for i in result["issues"]))
                (self.source_root / name).unlink()
        self.assertTrue(self.inspect()["ready_for_final_analysis"])

    def test_quiet_output_saves_json_and_preserves_readiness_exit_status(self):
        self.write_chunk(0, 2)  # Missing chunk: quiet must not conceal exit 3.
        output = self.base / "inventory.json"
        args = ["--config", str(self.config_path), "--root", str(self.root),
                "--output", str(output), "--require-ready"]
        default_stdout = io.StringIO()
        with redirect_stdout(default_stdout):
            self.assertEqual(I.main(args), 3)
        self.assertEqual(json.loads(default_stdout.getvalue()), json.loads(output.read_text()))
        quiet_stdout = io.StringIO()
        with redirect_stdout(quiet_stdout):
            self.assertEqual(I.main(args + ["--quiet"]), 3)
        self.assertEqual(quiet_stdout.getvalue(), "")
        self.assertFalse(json.loads(output.read_text())["ready_for_final_analysis"])

    def test_quiet_without_output_fails_instead_of_discarding_the_inventory(self):
        stderr = io.StringIO()
        with redirect_stderr(stderr), self.assertRaises(SystemExit) as caught:
            I.main(["--config", str(self.config_path), "--root", str(self.root), "--quiet"])
        self.assertEqual(caught.exception.code, 2)
        self.assertIn("--quiet requires --output", stderr.getvalue())

    def test_early_stopped_certified_decisions_are_complete(self):
        self.configure_permutations()
        self.write_chunk(rows=[self.permutation_row()])
        result = self.inspect()
        self.assertTrue(result["ready_for_final_analysis"])
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["states"], {"early_stopped_resolved": 1})

    def test_strict_and_weak_integer_boundaries(self):
        self.assertFalse(I.bounds(9, 599, 599, .05, 3, "strict"))
        self.assertTrue(I.bounds(9, 599, 599, .05, 3, "le"))
        self.assertIsNone(I.bounds(9, 598, 599, .05, 3, "le"))
        self.assertFalse(I.bounds(9, 598, 599, .05, 3, "strict"))

    def test_a_forged_exact_flag_does_not_mark_pending_work_complete(self):
        self.configure_permutations()
        row = self.permutation_row(n=8, gt=[0, 0, 0], maxgt=0)
        row["schemes"]["diagnosis"]["exact_decision"] = {
            "alpha": .05, "le": {"pair_reject": [False, False, False], "max_reject": False}}
        self.write_chunk(rows=[row])
        result = self.inspect()
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertTrue(any(issue["code"] == "invalid_scheme" for issue in result["issues"]))

    def test_done_cannot_replace_unresolved_permutation_decision(self):
        self.configure_permutations()
        self.write_chunk(rows=[self.permutation_row(n=8, gt=[0, 0, 0], maxgt=0)])
        result = self.inspect()
        self.assertFalse(result["collection_complete"])
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["states"], {"pending": 1})

    def test_forced_full_p_obligation_blocks_early_stop_completion(self):
        self.configure_permutations(full=True)
        self.write_chunk(rows=[self.permutation_row()])
        result = self.inspect()
        self.assertFalse(result["complete"])
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["states"], {"pending_full_p": 1})

    def test_failed_budget_exhaustion_is_terminal_but_not_resolved(self):
        self.configure_permutations()
        self.write_chunk(rows=[self.permutation_row(n=598, attempted=599, gt=[9, 99, 99], maxgt=99)])
        result = self.inspect()
        self.assertTrue(result["collection_complete"])
        self.assertTrue(result["terminal_with_failures"])
        self.assertFalse(result["ready_for_final_analysis"])
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["failed_permutations"], 1)

    def test_resolved_failed_permutation_is_retained_and_counted(self):
        self.configure_permutations()
        self.write_chunk(rows=[self.permutation_row(n=598, attempted=599, gt=[99, 99, 99], maxgt=99)])
        result = self.inspect()
        self.assertTrue(result["ready_for_final_analysis"])
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["states"], {"budget_exhausted_resolved": 1})
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["failed_permutations"], 1)

    def test_single_pair_scope_does_not_require_untested_pairs(self):
        self.configure_permutations(definition={"stratify": "diagnosis", "tested_pair_index": 0})
        row = self.permutation_row(n=20, gt=[10, 0, 0], maxgt=0)
        self.write_chunk(rows=[row])
        self.assertTrue(self.inspect()["ready_for_final_analysis"])

    def test_observed_failure_uses_planned_denominator(self):
        self.configure_permutations()
        row = self.make_row(100, "original")
        row.update(status="error", observed={})
        self.write_chunk(rows=[row])
        result = self.inspect()
        self.assertTrue(result["collection_complete"])
        self.assertFalse(result["complete"])
        self.assertEqual(result["totals"]["observed_failed_or_invalid"], 1)
        self.assertEqual(result["totals"]["schemes"]["diagnosis"]["states"], {"observed_failed_or_invalid": 1})

    def test_fit_records_are_never_opened(self):
        directory = self.write_chunk()
        (directory / "fit_records.jsonl").write_text("POISON: DO NOT READ\n")
        original = Path.open
        def guarded(path, *args, **kwargs):
            if path.name.startswith("fit_records"):
                raise AssertionError("Inventory tried to open a permutation event log")
            return original(path, *args, **kwargs)
        with mock.patch.object(Path, "open", guarded):
            self.assertTrue(self.inspect()["ready_for_final_analysis"])


if __name__ == "__main__":
    unittest.main()
