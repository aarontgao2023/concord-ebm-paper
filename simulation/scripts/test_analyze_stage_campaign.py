"""Pure hand-built stage campaign fixtures; no real confirmation data or fits.

Evidence uses namespace31 and a fake generator returning seven invented rows.
The complete-R1000 aggregation fixture is only repeated hand-counted metrics.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import analyze_stage_campaign as W
import analyze_stage_exchangeability as S
from test_analyze_stage_exchangeability import EvidenceFixture, SEED, DX, STAGE, LABEL


class StageCampaignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.f = EvidenceFixture(self.temp.name)
        self.workspace = Path(self.temp.name)
        self.raw = self.workspace / "raw"
        self.run = self.raw / W.RUN_ID
        self.f.config.update(phase="confirmation", datasets=1000, bperm=599, alpha=.05, rule="le",
                             full_p_first_n=10, cells=[{"name": "STAGE_H0", "overrides": {"stage_eta": 1.0}}],
                             schemes={"diagnosis": {"stratify": "diagnosis", "require_max": True},
                                      "oracle_dx_stage": {"stratify": "oracle_dx_stage", "require_max": True}})
        self.reseal()

    def reseal(self):
        cfg_rel = "configs/v2/confirm_stage_a.json"
        (self.f.snapshot / cfg_rel).write_text(json.dumps(self.f.config))
        files = sorted(p for p in self.f.snapshot.rglob("*") if p.is_file() and p.name != "SHA256SUMS")
        (self.f.snapshot / "SHA256SUMS").write_text("".join(f"{S._sha(p)}  {p.relative_to(self.f.snapshot)}\n" for p in files))
        cfg = self.f.config
        entry = {"run_id": W.RUN_ID, "snapshot": "snapshot", "config": cfg_rel,
                 "snapshot_sha256": S._sha(self.f.snapshot / "SHA256SUMS"), "config_sha256": S._sha(self.f.snapshot / cfg_rel),
                 "base_seed": cfg["base_seed"], "datasets_per_cell": cfg["datasets"], "cells": cfg["cells"],
                 "engines": cfg["engines"], "bperm": cfg["bperm"], "full_p_first_n": cfg["full_p_first_n"],
                 "full_p_seeds": cfg.get("full_p_seeds")}
        self.f.protocol.write_text(json.dumps({"runs": [entry]}))

    def record(self, offset=0, status="ok", no_hash=False):
        record = copy.deepcopy(self.f.record)
        seed = SEED + offset
        record.update(seed=seed, status=status)
        record["truth"]["seed"] = record["truth"]["manifest"]["seed"] = seed
        if no_hash:
            record["truth"] = {}
        return record

    def chunk(self, index=0, nchunks=1, records=(), tail=b""):
        path = self.run / "cell_0" / f"chunk_{index}"
        path.mkdir(parents=True, exist_ok=True)
        identity = {"config": self.f.config, "cell_index": 0, "chunk": index,
                    "nchunks": nchunks, "source_sha256": self.f.sources}
        manifest = identity | {"signature": hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest(),
            "environment": {"versions": {"numpy": self.f.truth["manifest"]["numpy_version"],
                                           "pandas": self.f.truth["manifest"]["pandas_version"]}},
            "seeds": list(range(SEED, SEED + 1000))[index::nchunks]}
        (path / "manifest.json").write_text(json.dumps(manifest))
        (path / "fit_records.jsonl").write_bytes(b"".join(json.dumps(r).encode() + b"\n" for r in records) + tail)
        return path

    def analyze(self):
        def toy_generator(cfg, seed):
            return self.f.df.copy(), self.record(seed - SEED)["truth"]
        self.f.generator.side_effect = toy_generator
        with patch.object(S, "_load_frozen_design", return_value=self.f.module):
            return W.analyze(self.f.protocol, self.raw, workspace=self.workspace)

    def test_empty_inputs_never_call_chunk_or_generator_and_retain_1000(self):
        with patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.f.generator.assert_not_called()
        self.assertEqual((len(result["datasets"]), result["planned_datasets"], result["pending_datasets"]), (1000, 1000, 1000))
        self.assertIsNone(result["cohort_summary"]["statistics"])
        self.assertFalse(result["cohort_summary"]["available"])
        self.assertIn("Cohort mean/quantiles: **unavailable**", W.render_markdown(result))

    def test_only_committed_terminal_observed_seed_is_reconstructed(self):
        path = self.chunk(records=[self.record(status="error"), self.record(2, status="running")],
                          tail=json.dumps(self.record(1)).encode())
        before = (path / "fit_records.jsonl").read_bytes()
        result = self.analyze()
        self.assertEqual([call.args[1] for call in self.f.generator.call_args_list], [SEED])
        self.assertEqual((result["verified_datasets"], result["committed_terminal_observed_count"], result["pending_datasets"]), (1, 1, 999))
        self.assertEqual(result["datasets"][0]["observed_status"], "error")
        self.assertEqual(result["datasets"][1]["diagnostic_status"], "pending_no_committed_observed")
        self.assertEqual(result["datasets"][2]["diagnostic_status"], "pending_nonterminal_observed")
        self.assertIsNone(result["cohort_summary"]["statistics"])
        self.assertEqual((path / "fit_records.jsonl").read_bytes(), before)
        json.dumps(result, allow_nan=False)

    def test_terminal_failure_without_hash_is_unavailable_without_import(self):
        self.chunk(records=[self.record(status="error", no_hash=True)])
        with patch.object(S, "_load_frozen_design") as load:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        load.assert_not_called()
        self.assertEqual((result["unavailable_datasets"], result["verified_datasets"], result["pending_datasets"]), (1, 0, 999))
        self.assertFalse(result["datasets"][0]["reconstructed"])
        self.assertIsNone(result["cohort_summary"]["statistics"])

    def test_corrupt_chunk_keeps_its_planned_seeds_and_continues_next_chunk(self):
        bad = self.chunk(index=0, nchunks=2)
        (bad / "fit_records.jsonl").write_bytes(b"corrupt committed line\n")
        self.chunk(index=1, nchunks=2, records=[self.record(1)])
        result = self.analyze()
        self.assertEqual([call.args[1] for call in self.f.generator.call_args_list], [SEED + 1])
        self.assertEqual((result["error_datasets"], result["verified_datasets"], result["pending_datasets"]), (500, 1, 499))
        self.assertEqual(len(result["datasets"]), 1000)
        self.assertEqual(result["errors"][0]["affected_planned_seeds"], list(range(SEED, SEED + 1000, 2)))
        self.assertIsNone(result["cohort_summary"]["statistics"])

    def test_missing_fitlog_and_manifest_stay_pending_without_generation(self):
        path = self.chunk(index=0, nchunks=2, records=[self.record()])
        (path / "fit_records.jsonl").unlink()
        path2 = self.chunk(index=1, nchunks=2, records=[self.record(1)])
        (path2 / "manifest.json").unlink()
        with patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.assertEqual(result["pending_datasets"], 1000)
        self.assertEqual(result["verified_datasets"], 0)

    def test_inconsistent_manifests_block_all_reconstruction_before_first_chunk(self):
        self.chunk(index=0, nchunks=2, records=[self.record()])
        self.chunk(index=1, nchunks=3, records=[self.record(1)])
        with patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.assertEqual(result["errors"][0]["stage"], "inconsistent_layout")
        self.assertEqual(result["verified_datasets"], 0)
        self.assertEqual(len(result["datasets"]), 1000)

    def test_source_or_seed_layout_mismatch_never_calls_per_chunk_reconstruction(self):
        path = self.chunk(records=[self.record()])
        original = json.loads((path / "manifest.json").read_text())
        for field in ("source", "seeds", "directory"):
            value = copy.deepcopy(original)
            if field == "source":
                value["source_sha256"]["design_v2.py"] = "0" * 64
            elif field == "seeds":
                value["seeds"] = [SEED]
            else:
                value["chunk"] = 1
            (path / "manifest.json").write_text(json.dumps(value))
            with self.subTest(field=field), patch.object(S, "analyze_chunk") as chunk:
                result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
            chunk.assert_not_called()
            self.assertEqual(result["errors"][0]["stage"], "manifest")
            self.assertEqual(len(result["datasets"]), 1000)
            self.assertIsNone(result["cohort_summary"]["statistics"])

    def test_successful_fit_missing_hash_errors_without_import(self):
        self.chunk(records=[self.record(no_hash=True)])
        with patch.object(S, "_load_frozen_design") as load:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        load.assert_not_called()
        self.assertEqual(result["error_datasets"], 1)
        self.assertEqual(result["datasets"][0]["diagnostic_status"], "error_chunk_diagnostic")
        self.assertIsNone(result["datasets"][0]["reconstructed"])
        self.assertEqual(result["pending_datasets"], 999)

    def test_only_complete_fixed1000_metrics_get_cohort_descriptive_statistics(self):
        metrics = S.block_metrics(DX, STAGE, LABEL, max_stage=2)
        rows = [{"seed": SEED + i, "diagnostic_status": "verified", "metrics": metrics} for i in range(1000)]
        complete = W.cohort_summary(rows, [])
        self.assertTrue(complete["available"])
        key = "oracle_dx_stage.movable_sample_fraction"
        self.assertAlmostEqual(complete["statistics"][key]["mean"], 5 / 7)
        self.assertEqual(complete["statistics"][key]["n_datasets"], 1000)
        self.assertEqual(set(complete["statistics"][key]["quantiles"]), {"0.05", "0.25", "0.5", "0.75", "0.95"})
        self.assertFalse(W.cohort_summary(rows[:999], [])["available"])
        self.assertFalse(W.cohort_summary(rows, [{"message": "unresolved source error"}])["available"])
        rows[0] = {"seed": SEED, "diagnostic_status": "unavailable_saved_generation_hash_missing"}
        self.assertIsNone(W.cohort_summary(rows, [])["statistics"])

    def test_changed_frozen_denominator_is_rejected_before_any_evidence_discovery(self):
        self.f.config["datasets"] = 2
        self.reseal()
        with patch.object(S, "analyze_chunk") as chunk, self.assertRaisesRegex(ValueError, "R1000"):
            W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()

    def test_cli_rejects_output_into_raw_before_analysis(self):
        args = ["analyze_stage_campaign", "--protocol", str(self.f.protocol), "--raw-root", str(self.raw),
                "--workspace", str(self.workspace), "--output-json", str(self.raw / "overwrite.json")]
        with patch("sys.argv", args), patch.object(W, "analyze") as analyze, self.assertRaisesRegex(ValueError, "overwrite"):
            W.main()
        analyze.assert_not_called()

    def test_symlinked_fit_log_outside_raw_is_not_reconstructed(self):
        path = self.chunk(records=[self.record()])
        outside = self.workspace / "outside_fit_records.jsonl"
        (path / "fit_records.jsonl").rename(outside)
        (path / "fit_records.jsonl").symlink_to(outside)
        with patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.assertEqual(result["error_datasets"], 1000)
        self.assertEqual(result["errors"][0]["stage"], "discovery")

    def test_protocol_changed_after_discovery_never_authorizes_reconstruction(self):
        self.chunk(records=[self.record()])
        discover = W._discover
        def changed(path, planned):
            value = discover(path, planned)
            self.f.protocol.write_text(self.f.protocol.read_text() + "\n")
            return value
        with patch.object(W, "_discover", side_effect=changed), patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.assertEqual(result["error_datasets"], 1000)
        self.assertIsNone(result["cohort_summary"]["statistics"])
        self.assertTrue(any(error["stage"] == "protocol_changed" for error in result["errors"]))

    def test_resigned_manifest_layout_changed_after_discovery_blocks_reconstruction(self):
        path = self.chunk(records=[self.record()])
        discover = W._discover
        def changed(log, planned):
            value = discover(log, planned)
            manifest = json.loads((path / "manifest.json").read_text())
            manifest["nchunks"] = 2
            manifest["seeds"] = manifest["seeds"][::2]
            identity = {k: manifest[k] for k in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
            manifest["signature"] = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
            (path / "manifest.json").write_text(json.dumps(manifest))
            return value
        with patch.object(W, "_discover", side_effect=changed), patch.object(S, "analyze_chunk") as chunk:
            result = W.analyze(self.f.protocol, self.raw, workspace=self.workspace)
        chunk.assert_not_called()
        self.assertEqual((result["verified_datasets"], result["error_datasets"]), (0, 1000))
        self.assertIsNone(result["cohort_summary"]["statistics"])

    def test_manifest_changed_during_per_chunk_call_never_accepts_verified(self):
        path = self.chunk(records=[self.record()])
        analyze_chunk = S.analyze_chunk
        def changed(*args, **kwargs):
            result = analyze_chunk(*args, **kwargs)
            manifest = path / "manifest.json"
            manifest.write_text(manifest.read_text() + "\n")
            return result
        with patch.object(S, "analyze_chunk", side_effect=changed):
            result = self.analyze()
        self.assertEqual([call.args[1] for call in self.f.generator.call_args_list], [SEED])
        self.assertEqual(result["verified_datasets"], 0)
        self.assertNotIn("metrics", result["datasets"][0])
        self.assertEqual(result["error_datasets"], 1000)

    def test_later_chunk_cannot_leave_earlier_changed_manifest_verified(self):
        earlier = self.chunk(index=0, nchunks=2, records=[self.record()])
        self.chunk(index=1, nchunks=2, records=[self.record(1)])
        analyze_chunk = S.analyze_chunk
        def changed(protocol, path, seeds, **kwargs):
            result = analyze_chunk(protocol, path, seeds, **kwargs)
            if path.name == "chunk_1":
                manifest = earlier / "manifest.json"
                manifest.write_text(manifest.read_text() + "\n")
            return result
        with patch.object(S, "analyze_chunk", side_effect=changed):
            result = self.analyze()
        self.assertEqual((result["verified_datasets"], result["error_datasets"], result["pending_datasets"]), (1, 500, 499))
        self.assertEqual(result["datasets"][0]["diagnostic_status"], "error_manifest_changed")
        self.assertNotIn("metrics", result["datasets"][0])
        self.assertEqual(result["datasets"][1]["diagnostic_status"], "verified")

    def test_per_chunk_return_must_preserve_planned_chunk_denominator(self):
        self.chunk(records=[self.record()])
        analyze_chunk = S.analyze_chunk
        def changed(*args, **kwargs):
            result = analyze_chunk(*args, **kwargs)
            result["planned_chunk_datasets"] = 999
            return result
        with patch.object(S, "analyze_chunk", side_effect=changed):
            result = self.analyze()
        self.assertEqual((result["verified_datasets"], result["error_datasets"]), (0, 1))
        self.assertNotIn("metrics", result["datasets"][0])
        self.assertEqual(len(result["datasets"]), 1000)

    def test_campaign_cli_rejects_hardlinked_raw_output_before_analysis(self):
        path = self.chunk(records=[self.record()])
        log = path / "fit_records.jsonl"
        before = log.read_bytes()
        output = self.workspace / "aliased_report.json"
        os.link(log, output)
        args = ["analyze_stage_campaign", "--protocol", str(self.f.protocol), "--raw-root", str(self.raw),
                "--workspace", str(self.workspace), "--output-json", str(output)]
        with patch("sys.argv", args), patch.object(W, "analyze") as analyze, self.assertRaisesRegex(ValueError, "hard links"):
            W.main()
        analyze.assert_not_called()
        self.assertEqual(log.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
