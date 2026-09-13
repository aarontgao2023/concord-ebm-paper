"""Synthetic small fixtures only. No confirmation data generation or model fitting."""
import ast
import copy
import hashlib
import importlib.util
from itertools import permutations
import json
import math
import marshal
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import struct
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

import analyze_stage_exchangeability as A

SEED = 31009501  # development-only fixture identity, never a confirmation seed
DX = ["CN"] * 5 + ["AD"] * 2
STAGE = [0, 0, 0, 1, 1, 0, 0]
LABEL = [0, 0, 1, 2, 2, 1, 2]


class MetricsTests(unittest.TestCase):
    def test_hand_counted_block_and_subject_denominators(self):
        m = A.block_metrics(DX, STAGE, LABEL, max_stage=1)
        o = m["oracle_dx_stage"]
        self.assertEqual((o["n"], o["possible_blocks"], o["occupied_blocks"], o["empty_blocks"]), (7, 6, 3, 3))
        self.assertEqual(o["single_group_blocks"], 1)
        self.assertAlmostEqual(o["single_group_block_fraction_among_occupied"], 1 / 3)
        self.assertAlmostEqual(o["single_group_sample_fraction"], 2 / 7)
        self.assertEqual(o["movable_samples"], 5)
        self.assertAlmostEqual(o["movable_sample_fraction"], 5 / 7)
        self.assertAlmostEqual(o["expected_label_changes_one_uniform_permutation"], 7 / 3)
        self.assertAlmostEqual(o["log_distinct_label_assignments"], math.log(6))
        self.assertAlmostEqual(m["diagnosis"]["log_distinct_label_assignments"], math.log(60))
        self.assertAlmostEqual(o["pair_overlap"]["e2_vs_e33"]["pair_sample_overlap_fraction"], 3 / 4)
        self.assertEqual(o["pair_overlap"]["e2_vs_e4"]["pair_sample_overlap_fraction"], 0)
        self.assertEqual(o["by_group"]["e4"]["movable_samples"], 1)

    def test_distinct_labels_and_expected_changes_match_complete_enumeration(self):
        labels = [0, 0, 1, 2]
        assignments = set(permutations(labels))
        expected_changes = sum(sum(a != b for a, b in zip(labels, p)) for p in assignments) / len(assignments)
        o = A.block_metrics(["CN"] * 4, [0] * 4, labels, max_stage=0)["oracle_dx_stage"]
        self.assertEqual(len(assignments), 12)  # 4!/(2!1!1!), not24 index permutations
        self.assertAlmostEqual(o["log_distinct_label_assignments"], math.log(len(assignments)))
        self.assertAlmostEqual(o["expected_label_changes_one_uniform_permutation"], expected_changes)
        self.assertFalse(o["randomization_space_degenerate"])

    def test_pure_blocks_are_degenerate_even_with_all_groups_in_dataset(self):
        m = A.block_metrics(["CN"] * 6, [0, 0, 1, 1, 2, 2], [0, 0, 1, 1, 2, 2], max_stage=2)
        o = m["oracle_dx_stage"]
        self.assertTrue(o["randomization_space_degenerate"])
        self.assertEqual(o["log_distinct_label_assignments"], 0)
        self.assertEqual(o["movable_samples"], 0)
        self.assertEqual(o["single_group_block_fraction_among_occupied"], 1)
        self.assertAlmostEqual(m["diagnosis"]["log_distinct_label_assignments"], math.log(90))
        self.assertIn("cannot establish", m["interpretation"])

    def test_invalid_or_misaligned_arrays_rejected(self):
        for dx, stages, labels in [([], [], []), (["CN"], [], [0]), (["CN"], [0.2], [0]),
                                  (["CN"], [2], [0]), (["unknown"], [0], [0]), (["CN"], [0], [3])]:
            with self.subTest(dx=dx, stages=stages, labels=labels), self.assertRaises(ValueError):
                A.block_metrics(dx, stages, labels, max_stage=1)

    def test_frozen_runner_grouping_matches_fixture_blocks(self):
        # Extract only the pure permutation function; importing the runner would
        # import fitting engines unnecessarily. No simulate/fit functions execute.
        runner = Path(__file__).resolve().parents[2] / "runs/v2/snapshots/confirm_stage_a/scripts/v2/run_v2.py"
        tree = ast.parse(runner.read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "permute")
        namespace = {"np": np, "hashlib": hashlib}
        exec(compile(ast.Module(body=[node], type_ignores=[]), str(runner), "exec"), namespace)
        df = pd.DataFrame({"Diagnosis": DX, "APOE": LABEL})
        before = A.block_metrics(DX, STAGE, LABEL, max_stage=1)["oracle_dx_stage"]["blocks"]
        for pid in range(12):
            out = namespace["permute"](df, {"latent_stage": STAGE}, SEED, "oracle_dx_stage", pid,
                                       {"stratify": "oracle_dx_stage", "require_max": True})
            after = A.block_metrics(out.Diagnosis, STAGE, out.APOE, max_stage=1)["oracle_dx_stage"]["blocks"]
            self.assertEqual(before, after)
            self.assertEqual(out.Diagnosis.tolist(), DX)
            self.assertEqual(out.APOE.iloc[3:5].tolist(), [2, 2])


class EvidenceFixture:
    def __init__(self, root):
        self.root = Path(root)
        self.snapshot = self.root / "snapshot"
        self.chunk = self.root / "chunk"
        self.chunk.mkdir()
        self.config = {"run_id": "confirm_stage_a", "base_seed": SEED, "datasets": 2,
                       "engines": ["repaired"], "cells": [{"name": "FIXTURE_STAGE"}],
                       "schemes": {"oracle_dx_stage": {"stratify": "oracle_dx_stage", "require_max": True}}}
        cfg_rel = "configs/v2/confirm_stage_a.json"
        files = {cfg_rel: json.dumps(self.config), "scripts/v2/design_v2.py": "# nonexecuted synthetic fixture source\n",
                 "scripts/v2/run_v2.py": "# nonexecuted synthetic fixture runner\n"}
        self.sources = {}
        lines = []
        for name, text in files.items():
            path = self.snapshot / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
            digest = A._sha(path)
            lines.append(f"{digest}  {name}")
            if name.startswith("scripts/"):
                self.sources[Path(name).name] = digest
        (self.snapshot / "SHA256SUMS").write_text("\n".join(lines) + "\n")
        protocol = {"runs": [{"run_id": "confirm_stage_a", "snapshot": "snapshot", "config": cfg_rel,
                               "snapshot_sha256": A._sha(self.snapshot / "SHA256SUMS"),
                               "config_sha256": A._sha(self.snapshot / cfg_rel)}]}
        self.protocol = self.root / "protocol.json"
        self.protocol.write_text(json.dumps(protocol))
        identity = {"config": self.config, "cell_index": 0, "chunk": 0, "nchunks": 1,
                    "source_sha256": self.sources}
        manifest = identity | {"signature": hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
                               "environment": {"versions": {"numpy": np.__version__, "pandas": pd.__version__}},
                               "seeds": [SEED, SEED + 1]}
        (self.chunk / "manifest.json").write_text(json.dumps(manifest))
        self.df = pd.DataFrame({"PTID": [f"fixture{i}" for i in range(7)], "Diagnosis": DX,
                                "a": [float(i) for i in range(7)], "b": [0., 1., np.nan, 3., 4., 5., 6.], "APOE": LABEL})
        self.resolved = {"fixture_only": True}
        config_sha = hashlib.sha256(json.dumps(self.resolved, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.truth = {"seed": SEED, "config": self.resolved, "group_order": list(A.GROUPS),
                      "biomarker_names": ["a", "b"], "latent_stage": STAGE,
                      "latent_row_ptid": self.df.PTID.tolist(),
                      "manifest": {"seed": SEED, "source_sha256": self.sources["design_v2.py"],
                                   "config_sha256": config_sha,
                                   "numpy_version": np.__version__, "pandas_version": pd.__version__,
                                   "data_sha256": hashlib.sha256(self.df.to_csv(index=False, float_format="%.17g").encode()).hexdigest()}}
        self.record = {"kind": "observed", "seed": SEED, "engine": "repaired", "scheme": None,
                       "perm_id": None, "status": "ok", "truth": self.truth}
        self.write_records([self.record])
        self.generator = Mock(return_value=(self.df.copy(), copy.deepcopy(self.truth)))
        self.module = SimpleNamespace(resolve=lambda *a, **k: SimpleNamespace(to_dict=lambda: self.resolved),
                                      simulate=self.generator)

    def write_records(self, records, tail=b""):
        (self.chunk / "fit_records.jsonl").write_bytes(b"".join(json.dumps(r).encode() + b"\n" for r in records) + tail)

    def analyze(self, seeds=None):
        return A.analyze_chunk(self.protocol, self.chunk, seeds or [SEED], workspace=self.root)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.f = EvidenceFixture(self.temp.name)

    def test_missing_any_selected_observed_record_blocks_all_reconstruction(self):
        with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "generation forbidden"):
            self.f.analyze([SEED, SEED + 1])
        load.assert_not_called()
        self.f.generator.assert_not_called()

    def test_matching_evidence_reconstructs_only_explicit_seed_and_keeps_failed_fit(self):
        self.f.record["status"] = "error"
        self.f.write_records([self.f.record])
        with patch.object(A, "_load_frozen_design", return_value=self.f.module):
            result = self.f.analyze()
        self.f.generator.assert_called_once()
        self.assertEqual(self.f.generator.call_args.args[1], SEED)
        self.assertEqual(result["datasets"][0]["observed_status"], "error")
        self.assertEqual(result["verified_datasets"], 1)
        self.assertFalse(result["fitting_performed"])
        json.dumps(result, allow_nan=False)
        self.assertIn("1/1", A.render_markdown(result))

    def test_error_without_generation_hash_retained_without_generator_import(self):
        self.f.record.update(status="error", truth={})
        self.f.write_records([self.f.record])
        with patch.object(A, "_load_frozen_design") as load:
            result = self.f.analyze()
        load.assert_not_called()
        self.assertEqual(result["verified_datasets"], 0)
        self.assertFalse(result["datasets"][0]["reconstructed"])

    def test_data_hash_and_stage_alignment_mismatches_are_rejected(self):
        for field in ("data", "stage", "ptid", "truth_ptid", "source", "config"):
            with self.subTest(field=field):
                df, truth = self.f.df.copy(), copy.deepcopy(self.f.truth)
                if field == "data":
                    df.loc[0, "a"] = 999.
                elif field == "stage":
                    truth["latent_stage"] = [1] * 7
                elif field == "ptid":
                    df.loc[0, "PTID"] = "unexpected_row"
                elif field == "truth_ptid":
                    truth["latent_row_ptid"] = list(reversed(truth["latent_row_ptid"]))
                elif field == "source":
                    truth["manifest"]["source_sha256"] = "0"*64
                else:
                    truth["manifest"]["config_sha256"] = "0"*64
                self.f.generator.return_value = df, truth
                with patch.object(A, "_load_frozen_design", return_value=self.f.module), self.assertRaises(ValueError):
                    self.f.analyze()

    def test_snapshot_or_chunk_source_tampering_blocks_import(self):
        path = self.f.snapshot / "scripts/v2/design_v2.py"
        path.write_text("# changed\n")
        with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "source bytes"):
            self.f.analyze()
        load.assert_not_called()

    def test_uncommitted_tail_never_authorizes_seed_and_is_not_repaired(self):
        second = self.f.record | {"seed": SEED + 1}
        self.f.write_records([self.f.record], json.dumps(second).encode())
        path = self.f.chunk / "fit_records.jsonl"
        before = path.read_bytes()
        with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "generation forbidden"):
            self.f.analyze([SEED + 1])
        load.assert_not_called()
        self.assertEqual(path.read_bytes(), before)
        with patch.object(A, "_load_frozen_design", return_value=self.f.module):
            result = self.f.analyze()
        self.assertIsNotNone(result["evidence"]["ignored_uncommitted_tail_sha256"])

    def test_duplicate_and_corrupt_committed_records_rejected(self):
        for records, tail in [([self.f.record, self.f.record], b""), ([self.f.record], b"{bad}\n")]:
            with self.subTest(tail=tail):
                self.f.write_records(records, tail)
                with patch.object(A, "_load_frozen_design") as load, self.assertRaises(ValueError):
                    self.f.analyze()
                load.assert_not_called()

    def test_wrong_environment_rejected_before_loading_frozen_module(self):
        with self.assertRaisesRegex(ValueError, "original NumPy/Pandas"):
            A._load_frozen_design(self.f.snapshot / "scripts/v2/design_v2.py",
                                  {"numpy_version": "wrong", "pandas_version": "wrong"})

    def test_unchecked_bytecode_is_ignored_and_module_registry_restored(self):
        path = Path(self.temp.name)/"toy_source.py"
        path.write_text('marker = "verified_source"\n')
        cache = Path(importlib.util.cache_from_source(str(path)))
        cache.parent.mkdir()
        code = compile('marker = "unchecked_bytecode"\n', str(path), "exec")
        cache.write_bytes(importlib.util.MAGIC_NUMBER +
                          struct.pack("<III", 0, int(path.stat().st_mtime), path.stat().st_size) + marshal.dumps(code))
        cache_before = cache.read_bytes()
        name = "_stage_diagnostic_frozen_design_" + A._sha(path)
        previous = SimpleNamespace(marker="prior_registry_entry")
        with patch.dict(sys.modules, {name: previous}):
            module = A._load_frozen_design(path, {"numpy_version": np.__version__, "pandas_version": pd.__version__}, A._sha(path))
            self.assertEqual(module.marker, "verified_source")
            self.assertIs(sys.modules[name], previous)
        self.assertEqual(cache.read_bytes(), cache_before)

    def test_source_change_between_snapshot_check_and_import_rejected(self):
        path = self.f.snapshot / "scripts/v2/design_v2.py"
        expected = A._sha(path)
        path.write_text('raise AssertionError("must not execute")\n')
        with self.assertRaisesRegex(ValueError, "before source compilation"):
            A._load_frozen_design(path, {"numpy_version": np.__version__, "pandas_version": pd.__version__}, expected)

    def test_nonterminal_record_cannot_authorize_generation(self):
        for status in (None, "pending", "running", "unknown"):
            with self.subTest(status=status):
                self.f.record["status"] = status
                self.f.write_records([self.f.record])
                with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "terminal"):
                    self.f.analyze()
                load.assert_not_called()

    def test_success_without_data_hash_is_not_silently_a_failed_diagnostic(self):
        self.f.record["truth"] = {}
        self.f.write_records([self.f.record])
        with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "Successful"):
            self.f.analyze()
        load.assert_not_called()

    def test_saved_config_hash_or_chunk_environment_mismatch_blocks_import(self):
        original = copy.deepcopy(self.f.truth)
        for field in ("config_sha256", "numpy_version", "pandas_version"):
            with self.subTest(field=field):
                self.f.record["truth"] = copy.deepcopy(original)
                self.f.record["truth"]["manifest"][field] = "incorrect"
                self.f.write_records([self.f.record])
                with patch.object(A, "_load_frozen_design") as load, self.assertRaises(ValueError):
                    self.f.analyze()
                load.assert_not_called()

    def test_missing_runtime_source_is_rejected_even_with_recomputed_manifest_signature(self):
        path = self.f.chunk/"manifest.json"
        manifest = json.loads(path.read_text())
        del manifest["source_sha256"]["run_v2.py"]
        identity = {k: manifest[k] for k in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
        manifest["signature"] = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        path.write_text(json.dumps(manifest))
        with patch.object(A, "_load_frozen_design") as load, self.assertRaisesRegex(ValueError, "closure"):
            self.f.analyze()
        load.assert_not_called()

    def test_all_selected_configs_checked_before_first_reconstruction(self):
        second = copy.deepcopy(self.f.record)
        second["seed"] = second["truth"]["seed"] = second["truth"]["manifest"]["seed"] = SEED+1
        second["truth"]["config"] = {"different_design": True}
        second["truth"]["manifest"]["config_sha256"] = hashlib.sha256(json.dumps(
            second["truth"]["config"], sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        self.f.write_records([self.f.record, second])
        with patch.object(A, "_load_frozen_design", return_value=self.f.module), self.assertRaisesRegex(ValueError, "Resolved"):
            self.f.analyze([SEED, SEED+1])
        self.f.generator.assert_not_called()

    def test_requested_failures_are_counted_including_unavailable_generation(self):
        failed = {"kind": "observed", "seed": SEED+1, "engine": "repaired", "scheme": None,
                  "perm_id": None, "status": "error", "truth": {}}
        self.f.write_records([self.f.record, failed])
        with patch.object(A, "_load_frozen_design", return_value=self.f.module):
            result = self.f.analyze([SEED, SEED+1])
        self.assertEqual(result["observed_status_counts"], {"ok": 1, "error": 1})
        self.assertEqual(len(result["datasets"]), 2)
        self.assertEqual(result["verified_datasets"], 1)
        self.f.generator.assert_called_once()

    def test_output_paths_cannot_overwrite_inputs_or_alias_each_other(self):
        forbidden = [self.f.protocol, self.f.chunk/"fit_records.jsonl",
                     self.f.chunk/"manifest.json", self.f.snapshot/"SHA256SUMS",
                     self.f.snapshot/"new_output.json", self.f.root/"configs/v2/confirm_stage_a.json"]
        for target in forbidden:
            with self.subTest(target=target), self.assertRaisesRegex(ValueError, "overwrite"):
                A.validate_output_paths(self.f.protocol, self.f.chunk, [target], workspace=self.f.root)
        same = self.f.root/"report.json"
        with self.assertRaisesRegex(ValueError, "different"):
            A.validate_output_paths(self.f.protocol, self.f.chunk, [same, same], workspace=self.f.root)
        A.validate_output_paths(self.f.protocol, self.f.chunk, [same, self.f.root/"report.md"], workspace=self.f.root)

    def test_per_seed_cli_rejects_hardlinked_raw_output_before_reconstruction(self):
        log = self.f.chunk / "fit_records.jsonl"
        before = log.read_bytes()
        output = self.f.root / "aliased_report.json"
        os.link(log, output)
        args = ["analyze_stage_exchangeability", "--protocol", str(self.f.protocol), "--chunk-dir", str(self.f.chunk),
                "--seed", str(SEED), "--workspace", str(self.f.root), "--output-json", str(output)]
        with patch("sys.argv", args), patch.object(A, "analyze_chunk") as analyze, self.assertRaisesRegex(ValueError, "hard links"):
            A.main()
        analyze.assert_not_called()
        self.assertEqual(log.read_bytes(), before)

    def test_output_alias_created_during_analysis_is_rejected_at_publication(self):
        log = self.f.chunk / "fit_records.jsonl"
        before = log.read_bytes()
        output = self.f.root / "report.json"
        A.validate_output_paths(self.f.protocol, self.f.chunk, [output], workspace=self.f.root)
        os.link(log, output)
        with self.assertRaisesRegex(ValueError, "hard links"):
            A.write_outputs(self.f.protocol, self.f.chunk, [(output, "new content")], workspace=self.f.root)
        self.assertEqual(log.read_bytes(), before)

    def test_atomic_report_publication_does_not_truncate_late_hardlink_alias(self):
        log = self.f.chunk / "fit_records.jsonl"
        before = log.read_bytes()
        output = self.f.root / "report.json"
        replace = os.replace
        def raced(source, destination):
            # Alias appears after the final destination check. Atomic replace
            # detaches this directory entry instead of truncating the raw inode.
            os.link(log, destination)
            replace(source, destination)
        with patch.object(A.os, "replace", side_effect=raced):
            A.write_outputs(self.f.protocol, self.f.chunk, [(output, "new content")], workspace=self.f.root)
        self.assertEqual(output.read_text(), "new content")
        self.assertEqual(log.read_bytes(), before)
        self.assertFalse(output.samefile(log))
        self.assertEqual(list(self.f.root.glob(".report.json.*.tmp")), [])

    def test_atomic_publication_rechecks_alias_after_temporary_file_creation(self):
        log = self.f.chunk / "fit_records.jsonl"
        before = log.read_bytes()
        output = self.f.root / "report.json"
        validate = A.validate_output_paths
        calls = 0
        def raced(*args, **kwargs):
            nonlocal calls
            calls += 1
            if calls == 2:
                os.link(log, output)
            return validate(*args, **kwargs)
        with patch.object(A, "validate_output_paths", side_effect=raced), self.assertRaisesRegex(ValueError, "hard links"):
            A.write_outputs(self.f.protocol, self.f.chunk, [(output, "new content")], workspace=self.f.root)
        self.assertEqual(log.read_bytes(), before)
        self.assertEqual(list(self.f.root.glob(".report.json.*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
