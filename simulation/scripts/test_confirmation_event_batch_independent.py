"""Bounded independent wrapper checks; raw audit arithmetic is not retested.

Reuse the already-reviewed inert signed-log builder. All changes made by these
tests affect temporary invented evidence only, never actual snapshots or data.
"""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit_confirmation_event_batch as B
from test_audit_confirmation_events import FrozenToy


class IndependentBatchBindingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)

    def whole_chunk(self):
        t = FrozenToy(self.temp.name, R=2)
        t.records = [t.event(pid, seed=seed) for seed in sorted(t.run["seeds"]) for pid in (None, 0, 1)]
        t.write()
        return t

    def split_chunks(self, n=2):
        t = FrozenToy(self.temp.name, R=n)
        cell = t.chunk.parent
        t.chunks = []
        for index, seed in enumerate(sorted(t.run["seeds"])):
            t.chunk = cell / f"chunk_{index}"
            t.chunk.mkdir(exist_ok=True)
            signed = {"config": t.cfg, "cell_index": 0, "chunk": index, "nchunks": n,
                      "source_sha256": t.run["sources"]}
            manifest = signed | {"signature": B.A.I.digest(signed), "seeds": [seed], "environment": t.env}
            (t.chunk / "manifest.json").write_text(json.dumps(manifest))
            t.records = [t.event(pid, seed=seed) for pid in (None, 0, 1)]
            t.write()
            self.progress(t.chunk, expected_rows=1)
            t.chunks.append(t.chunk)
        return t

    def progress(self, chunk, **changes):
        path = chunk / "progress.json"
        value = json.loads(path.read_text())
        value.update(changes)
        path.write_text(json.dumps(value))

    def reshard_whole_chunk(self, t):
        """Legally re-sign n1/two seeds into n2/one seed; snapshot unchanged."""
        seed = min(t.run["seeds"])
        manifest = json.loads((t.chunk / "manifest.json").read_text())
        manifest.update(nchunks=2, seeds=[seed])
        signed = {k: manifest[k] for k in ("config", "cell_index", "chunk", "nchunks", "source_sha256")}
        manifest["signature"] = B.A.I.digest(signed)
        (t.chunk / "manifest.json").write_text(json.dumps(manifest))
        t.records = [t.event(pid, seed=seed) for pid in (None, 0, 1)]
        t.write()
        self.progress(t.chunk, expected_rows=1)

    def run_batch(self, t, **options):
        return B.run_batch(t.protocol, ["toy"], t.root / "raw", t.root / "out", project=t.root, **options)

    def test_legally_resharded_candidate_is_not_scanned_as_original_selection(self):
        t = self.whole_chunk()
        real = B.exclusive_write
        def publish(path, value):
            real(path, value)
            if path.name == "selection.json":
                self.reshard_whole_chunk(t)
        with patch.object(B, "exclusive_write", side_effect=publish), \
                patch.object(B.A, "analyze", side_effect=AssertionError("Different manifest must not be scanned")):
            result = self.run_batch(t)
        selection = json.loads((t.root / "out/selection.json").read_text())
        self.assertEqual(len(selection["items"][0]["planned_seeds"]), 2)
        self.assertEqual(result["audits_started"], 0)
        self.assertEqual(result["batch_checkpoint_certified_count"], 0)
        self.assertEqual(result["final_category_counts"], {"deferred": 1})

    def test_manifest_changed_after_recheck_cannot_certify_replacement_audit(self):
        t = self.whole_chunk()
        real = B.A.analyze
        def audit(*args, **kwargs):
            self.reshard_whole_chunk(t)
            return real(*args, **kwargs)
        with patch.object(B.A, "analyze", side_effect=audit):
            result = self.run_batch(t)
        self.assertEqual(result["audits_started"], 1)
        self.assertEqual(result["batch_checkpoint_certified_count"], 0)
        artifact = json.loads((t.root / "out/toy__cell_0__chunk_0.json").read_text())
        self.assertFalse(artifact["checkpoint_certified"])
        self.assertIn("audit_manifest_differs_from_fixed_selection", [x["code"] for x in artifact["deferred_reasons"]])

    def test_later_chunk_can_revoke_earlier_certificate_only_in_final_summary(self):
        t = self.split_chunks()
        real = B.A.analyze
        original_first_artifact = []
        def audit(protocol, run_id, chunk, *args, **kwargs):
            result = real(protocol, run_id, chunk, *args, **kwargs)
            if Path(chunk).name == "chunk_1":
                original_first_artifact.append((t.root / "out/toy__cell_0__chunk_0.json").read_bytes())
                path = t.chunks[0] / "manifest.json"
                path.write_text(path.read_text() + "\n")  # Even same semantic identity has different selected bytes.
            return result
        with patch.object(B.A, "analyze", side_effect=audit):
            result = self.run_batch(t)
        self.assertTrue(result["source_stable"])
        self.assertFalse(result["selected_all_checkpoint_certified"])
        self.assertEqual(result["batch_checkpoint_certified_count"], 1)
        self.assertEqual(result["final_category_counts"], {"certified_done": 1, "deferred": 1})
        self.assertFalse(result["outcomes"][0]["checkpoint_certified"])
        self.assertEqual((t.root / "out/toy__cell_0__chunk_0.json").read_bytes(), original_first_artifact[0])

    def test_offset_selection_stays_fixed_when_live_candidates_change(self):
        t = self.split_chunks(4)
        self.progress(t.chunks[0], done=False, stopped=False)
        real = B.exclusive_write
        def publish(path, value):
            real(path, value)
            if path.name == "selection.json":
                self.progress(t.chunks[0], done=True, stopped=False)
                self.progress(t.chunks[2], done=False, stopped=False)
        with patch.object(B, "exclusive_write", side_effect=publish), \
                patch.object(B.A, "analyze", side_effect=AssertionError("Do not refill a now-live selected slot")):
            result = self.run_batch(t, offset=1, max_chunks=1)
        selection = json.loads((t.root / "out/selection.json").read_text())
        self.assertEqual(selection["selected_order"], [str(t.chunks[2].resolve())])
        self.assertEqual((result["candidate_chunks"], result["selected_chunks"], result["audits_started"]), (3, 1, 0))
        self.assertEqual(result["batch_checkpoint_certified_count"], 0)
        self.assertEqual(result["final_category_counts"], {"live": 2, "not_selected_max_chunks": 1, "not_selected_offset": 1})


if __name__ == "__main__":
    unittest.main()
