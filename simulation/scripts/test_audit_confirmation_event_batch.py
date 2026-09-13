"""Small invented log fixtures only; no real run data, fitting or HPC."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit_confirmation_event_batch as B
from test_audit_confirmation_events import FrozenToy


class BatchAuditTests(unittest.TestCase):
    def toy(self, chunks=3):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        t = FrozenToy(temp.name, R=chunks)
        root = t.chunk.parent
        t.chunks = []
        for i, seed in enumerate(sorted(t.run["seeds"])):
            t.chunk = root / f"chunk_{i}"
            t.chunk.mkdir(exist_ok=True)
            signed = {"config": t.cfg, "cell_index": 0, "chunk": i, "nchunks": chunks, "source_sha256": t.run["sources"]}
            manifest = signed | {"signature": B.A.I.digest(signed), "seeds": [seed], "environment": t.env}
            (t.chunk / "manifest.json").write_text(json.dumps(manifest))
            t.records = [t.event(seed=seed)] + [t.event(pid, seed=seed) for pid in range(t.cfg["bperm"])]
            t.write()
            progress = json.loads((t.chunk / "progress.json").read_text())
            progress["expected_rows"] = 1
            (t.chunk / "progress.json").write_text(json.dumps(progress))
            t.chunks.append(t.chunk)
        return t

    def run_batch(self, t, **kwargs):
        return B.run_batch(t.protocol, ["toy"], t.root / "raw", t.root / "outputs", project=t.root, **kwargs)

    def test_offset_and_limit_leave_unselected_uncertified_and_immutable(self):
        t = self.toy()
        inputs = {p: p.read_bytes() for p in (t.root / "raw").rglob("*") if p.is_file()}
        r = self.run_batch(t, offset=1, max_chunks=1)
        self.assertTrue(r["selected_all_checkpoint_certified"], r)
        self.assertEqual(r["audits_started"], 1)
        self.assertEqual(r["batch_done_certified_count"], 1)
        self.assertEqual(r["final_category_counts"], {"certified_done": 1, "not_selected_max_chunks": 1, "not_selected_offset": 1})
        selection = json.loads((t.root / "outputs/selection.json").read_text())
        self.assertEqual(selection["selected_order"], [str(t.chunks[1].resolve())])
        self.assertEqual(len(selection["items"]), 3)
        self.assertFalse(r["whole_run_completion_assessed"])
        self.assertEqual({p: p.read_bytes() for p in inputs}, inputs)
        saved = {p: p.read_bytes() for p in (t.root / "outputs").iterdir()}
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.run_batch(t)
        self.assertEqual({p: p.read_bytes() for p in saved}, saved)

    def test_live_missing_metadata_and_invalid_manifest_separately_accounted(self):
        t = self.toy(chunks=4)
        p = t.chunks[1] / "progress.json"
        progress = json.loads(p.read_text())
        progress.update(done=False, stopped=False)
        p.write_text(json.dumps(progress))
        (t.chunks[2] / "progress.json").unlink()
        p = t.chunks[3] / "manifest.json"
        manifest = json.loads(p.read_text())
        manifest["signature"] = "a"*64
        p.write_text(json.dumps(manifest))
        r = self.run_batch(t)
        self.assertEqual(r["final_category_counts"], {"certified_done": 1, "deferred": 1, "error": 1, "live": 1})
        self.assertEqual((r["discovered_chunks"], r["candidate_chunks"], r["audits_started"]), (4, 1, 1))

    def test_writer_restarts_after_selection_is_live_and_never_scanned(self):
        t = self.toy(chunks=1)
        real = B.exclusive_write
        def publish(path, value):
            real(path, value)
            if path.name == "selection.json":
                p = t.chunks[0] / "progress.json"
                progress = json.loads(p.read_text())
                progress.update(done=False, stopped=False)
                p.write_text(json.dumps(progress))
        with patch.object(B, "exclusive_write", side_effect=publish), patch.object(B.A, "analyze", side_effect=AssertionError("Must not scan live chunk")):
            r = self.run_batch(t)
        self.assertEqual(r["final_category_counts"], {"live": 1})
        self.assertEqual(r["audits_started"], 0)
        self.assertFalse(r["selected_all_checkpoint_certified"])

    def test_walltime_budget_stops_between_chunks_not_counted_as_pass(self):
        t = self.toy()
        now = [0.]
        real = B.A.analyze
        def audit(*args, **kwargs):
            result = real(*args, **kwargs)
            now[0] = 11.
            return result
        with patch.object(B.time, "monotonic", side_effect=lambda: now[0]), patch.object(B.A, "analyze", side_effect=audit):
            r = self.run_batch(t, walltime_budget_s=10)
        self.assertEqual(r["final_category_counts"], {"certified_done": 1, "unprocessed_walltime_budget": 2})
        self.assertEqual(r["audits_started"], 1)
        self.assertFalse(r["selected_all_checkpoint_certified"])

    def test_source_change_revokes_batch_and_leaves_remaining_unprocessed(self):
        t = self.toy()
        real = B.A.analyze
        def audit(*args, **kwargs):
            result = real(*args, **kwargs)
            t.protocol.write_text('{"replacement": true}')
            return result
        with patch.object(B.A, "analyze", side_effect=audit):
            r = self.run_batch(t)
        self.assertFalse(r["source_stable"])
        self.assertEqual(r["batch_checkpoint_certified_count"], 0)
        self.assertEqual(r["final_category_counts"], {"deferred": 1, "unprocessed_source_changed": 2})
        self.assertIn(str(t.protocol.resolve()), r["changed_source_paths"])

    def test_budget_expiring_during_metadata_recheck_does_not_start_scan(self):
        t = self.toy(chunks=1)
        now, calls = [0.], [0]
        real = B.inspect_selection
        def inspect(*args):
            result = real(*args)
            calls[0] += 1
            if calls[0] == 2:
                now[0] = 11.
            return result
        with patch.object(B.time, "monotonic", side_effect=lambda: now[0]), patch.object(B, "inspect_selection", side_effect=inspect), \
                patch.object(B.A, "analyze", side_effect=AssertionError("Budget already expired")):
            r = self.run_batch(t, walltime_budget_s=10)
        self.assertEqual(r["audits_started"], 0)
        self.assertEqual(r["final_category_counts"], {"unprocessed_walltime_budget": 1})

    def test_audit_errors_and_deferred_outputs_are_individual_artifacts(self):
        t = self.toy()
        real = B.A.analyze
        def audit(protocol, run, chunk, *args):
            if str(chunk).endswith("chunk_0"):
                raise ValueError("Toy raw identity mismatch")
            result = real(protocol, run, chunk, *args)
            if str(chunk).endswith("chunk_1"):
                result.update(status="deferred", checkpoint_certified=False, certified_done_chunk=False)
                result["deferred_reasons"].append({"code": "toy_live"})
            return result
        with patch.object(B.A, "analyze", side_effect=audit):
            r = self.run_batch(t)
        self.assertEqual(r["final_category_counts"], {"certified_done": 1, "deferred": 1, "error": 1})
        self.assertEqual(len(list((t.root / "outputs").glob("toy__*.json"))), 3)
        self.assertEqual(len(list((t.root / "outputs").glob("toy__*.md"))), 3)

    def test_guard_snapshot_output_and_exclusive_publish(self):
        t = self.toy(chunks=1)
        with self.assertRaisesRegex(ValueError, "raw/frozen"):
            B.run_batch(t.protocol, ["toy"], t.root / "raw", t.snapshot / "new-output", project=t.root)
        target = t.root / "immutable.json"
        target.symlink_to(t.protocol)
        before = t.protocol.read_bytes()
        with self.assertRaises(FileExistsError):
            B.exclusive_write(target, {"replacement": True})
        self.assertEqual(t.protocol.read_bytes(), before)

    def test_only_explicit_registered_runs_and_existing_manifest_scope(self):
        t = self.toy(chunks=1)
        with self.assertRaisesRegex(ValueError, "distinct"):
            B.run_batch(t.protocol, ["toy", "toy"], t.root / "raw", t.root / "outputs", project=t.root)
        with self.assertRaisesRegex(ValueError, "registered"):
            B.run_batch(t.protocol, ["unregistered"], t.root / "raw", t.root / "outputs", project=t.root)
        (t.chunks[0] / "manifest.json").unlink()
        r = self.run_batch(t)
        self.assertEqual(r["discovered_chunks"], 0)
        self.assertFalse(r["selected_all_checkpoint_certified"])
        self.assertEqual(r["scopes"][0]["planned_rows"], 1)


if __name__ == "__main__":
    unittest.main()
