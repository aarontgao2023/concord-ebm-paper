"""Independent pure-toy raw-auditor checks; no simulator or real run data.

The small direct audit_chunk fixture intentionally bypasses frozen_run. It
exercises log/restart semantics against invented three-event observations;
the author's integration suite covers full signed snapshot loading.
"""
from collections import Counter
from contextlib import redirect_stderr, redirect_stdout
import copy
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit_confirmation_events as A


class Toy:
    def __init__(self, root, *, full=True, budget=5):
        self.root = Path(root).resolve()
        self.chunk = self.root / "raw/toy/cell_0/chunk_0"
        self.chunk.mkdir(parents=True)
        self.seed = 42999991  # Invented identity; no data generation.
        self.cfg = {"run_id": "toy", "phase": "confirmation", "base_seed": self.seed,
                    "datasets": 1, "engines": ["repaired"], "cells": [{"name": "REF_H0"}],
                    "schemes": {"diagnosis": {"stratify": "diagnosis", "require_max": True}},
                    "bperm": budget, "alpha": .05, "rule": "le", "fast_likelihood": True}
        self.env = {"versions": {name: "toy" for name in A.I.NUMERICAL_PACKAGES},
                    "source_sha256": {"toy.py": "b" * 64}, "wheel_sha256": "c" * 64}
        self.env["versions"]["pyebm"] = "2.0.3"
        self.design = {"biomarker_names": ["a", "b", "c"], "target_inversions": 0,
                       "alternative_group": None}
        snap = self.root / "snapshot"
        (snap / "scripts/v2").mkdir(parents=True)
        (snap / "scripts/v2/fast_likelihood_v2.py").write_text("FAST_LIKELIHOOD_VERSION = 'toy-fast'\n")
        (snap / "scripts/v2/design_v2.py").write_text("# Inert toy evidence; never imported.\n")
        config = snap / "config.json"
        config.write_text(json.dumps(self.cfg))
        sources = {p.name: A.C.sha(p) for p in (snap / "scripts/v2").iterdir()}
        (snap / "SHA256SUMS").write_text("".join(f"{A.C.sha(p)}  {p.relative_to(snap)}\n"
            for p in sorted(snap.rglob("*")) if p.is_file() and p.name != "SHA256SUMS"))
        self.run = {"cfg": self.cfg, "names": ["REF_H0"], "designs": {"REF_H0": self.design},
                    "seeds": {self.seed}, "full_p": {self.seed} if full else set(),
                    "sources": sources, "snapshot": snap,
                    "config_path": config, "dependency": {"PYEBM_VERSION": "2.0.3",
                        "SOURCE_SHA256": self.env["source_sha256"], "PYEBM_WHEEL_SHA256": self.env["wheel_sha256"]}}
        signed = {"config": self.cfg, "cell_index": 0, "chunk": 0, "nchunks": 1,
                  "source_sha256": self.run["sources"]}
        manifest = signed | {"signature": A.I.digest(signed), "seeds": [self.seed], "environment": self.env}
        (self.chunk / "manifest.json").write_text(json.dumps(manifest))
        identity = [0, 1, 2]
        self.truth = {"config": self.design, "seed": self.seed,
                      "manifest": {"seed": self.seed, "config_sha256": A.C.canonical_sha(self.design),
                          "source_sha256": sources["design_v2.py"], "data_sha256": "d" * 64,
                          "numpy_version": "toy", "pandas_version": "toy"},
                      "biomarker_names": ["a", "b", "c"], "group_order": list(A.C.GROUPS),
                      "base_order": identity, "group_orderings": {g: identity for g in A.C.GROUPS},
                      "planned_inversions": {g: 0 for g in A.C.GROUPS},
                      "realized_inversions": {g: 0 for g in A.C.GROUPS},
                      "pair_null": [True] * 3, "h1_distance": 0.,
                      "pair_truth": {name: {"inversions": 0, "null": True, "normalized_distance": 0.}
                                     for name in A.C.PAIR_NAMES}}
        self.records = [self.record(None)]

    def record(self, pid, *, status="ok", orders=None):
        orders = orders or [[0, 1, 2]] * 3
        record = {"kind": "observed" if pid is None else "permutation", "seed": self.seed,
                  "engine": "repaired", "scheme": None if pid is None else "diagnosis",
                  "perm_id": pid, "status": status, "likelihood_implementation": "toy-fast",
                  "orderings": orders if status == "ok" else None,
                  "taus": [v / 3 for v in A.pair_counts(orders, 3)] if status == "ok" else None,
                  "fit_seconds": 1., "diagnostics": {"provenance": self.env, "mode": "repaired"}}
        if pid is None:
            record.update(truth=copy.deepcopy(self.truth), distance_to_truth=[A.C.inversions(x, [0, 1, 2]) / 3 for x in orders] if status == "ok" else None)
        return record

    def row(self, records):
        observed = records[0]
        good = [r for r in records[1:] if r["status"] == "ok"]
        n, B = len(good), self.cfg["bperm"]
        reference = observed["taus"]
        schemes = {}
        if observed["status"] == "ok":
            gt = [sum(r["taus"][i] > reference[i] + 1e-9 for r in good) for i in range(3)]
            eq = [sum(abs(r["taus"][i] - reference[i]) <= 1e-9 for r in good) for i in range(3)]
            maxgt = sum(max(r["taus"]) > max(reference) + 1e-9 for r in good)
            maxeq = sum(abs(max(r["taus"]) - max(reference)) <= 1e-9 for r in good)
            exact = {"alpha": .05}
            for rule in ("strict", "le"):
                exact[rule] = {"pair_reject": [A.I.bounds(g + e, n, B, .05, 3, rule) for g, e in zip(gt, eq)],
                               "max_reject": A.I.bounds(maxgt + maxeq, n, B, .05, 1, rule)}
            schemes["diagnosis"] = {"requested_nperm": B, "nperm": n, "attempted": len(records) - 1,
                "failed": len(records) - 1 - n, "complete": n == B, "gt": gt, "eq": eq,
                "maxgt": maxgt, "maxeq": maxeq, "exact_decision": exact, "definition": self.cfg["schemes"]["diagnosis"]}
        return {"run_id": "toy", "phase": "confirmation", "cell": "REF_H0", "seed": self.seed,
                "engine": "repaired", "status": observed["status"], "truth": observed.get("truth", {}),
                "observed": {k: observed.get(k) for k in ("orderings", "taus", "distance_to_truth", "fit_seconds")},
                "diagnostics": observed["diagnostics"] | {"sequential_stopping": True, "full_p_selected": self.seed in self.run["full_p"]},
                "schemes": schemes}

    def write(self, *, derived_records=None, done=False, stopped=True):
        (self.chunk / "fit_records.jsonl").write_text("".join(json.dumps(r) + "\n" for r in self.records))
        row = self.row(self.records if derived_records is None else derived_records)
        (self.chunk / "rows.jsonl").write_text(json.dumps(row) + "\n")
        progress = {"done": done, "stopped": stopped, "fit_records": len(self.records), "rows": 1,
                    "expected_rows": 1, "computational_jobs": len(self.records),
                    "statuses": dict(Counter(r["status"] for r in self.records))}
        (self.chunk / "progress.json").write_text(json.dumps(progress))


class IndependentEventAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.toy = Toy(self.temp.name)

    def test_nonprefix_success_and_failure_preserve_exact_missing_ids(self):
        t = self.toy
        t.records += [t.record(4), t.record(1, status="error"), t.record(3)]
        t.write()
        result = A.audit_chunk(t.run, t.chunk)
        self.assertEqual(result["status"], "certified_stopped_checkpoint", result)
        item = result["continuation"][0]
        self.assertEqual(item["required_resume_id_ranges"], [[0, 0], [2, 2]])
        self.assertEqual(item["failed_permutation_id_ranges"], [[1, 1]])
        self.assertEqual((item["attempted"], item["successful"], item["failed"]), (3, 2, 1))

    def test_pair_and_max_tie_counts_match_hand_worked_permutations(self):
        t = self.toy
        reference = [[0, 1, 2], [0, 2, 1], [1, 0, 2]]  # Counts 1,1,2.
        t.records = [t.record(None, orders=reference),
                     t.record(4),  # Counts 0,0,0.
                     t.record(1, orders=reference),  # Counts 1,1,2.
                     t.record(3, orders=[[0, 1, 2], [2, 1, 0], [0, 1, 2]]),  # 3,0,3.
                     t.record(0, orders=[[0, 1, 2], [1, 2, 0], [2, 0, 1]]),  # 2,2,2.
                     t.record(2, status="error")]
        t.write(done=True, stopped=False)
        result = A.audit_chunk(t.run, t.chunk)
        self.assertTrue(result["certified_done_chunk"], result)
        item = result["continuation"][0]
        raw = item["raw_recomputed"]
        self.assertEqual(raw["gt"], [2, 1, 1])
        self.assertEqual(raw["eq"], [1, 1, 2])
        self.assertEqual((raw["maxgt"], raw["maxeq"]), (1, 2))
        self.assertEqual(item["action"], "budget_exhausted_failures_preserved")
        self.assertFalse(item["full_p_available"])
        self.assertEqual(item["failed_permutation_id_ranges"], [[2, 2]])
        self.assertEqual(item["required_resume_id_ranges"], [])

    def test_actual_599_budget_equality_and_failed_unknown_completion(self):
        # 9 inclusive exceedances give p=(1+9)/600=.05/3 exactly.
        # One failed planned permutation must leave the <= decision unknown.
        for failed in (False, True):
            with self.subTest(failed=failed):
                t = Toy(Path(self.temp.name) / str(failed), budget=599)
                reference = [[0, 1, 2], [0, 2, 1], [1, 0, 2]]
                t.records = [t.record(None, orders=reference)]
                t.records += [t.record(pid, orders=reference) for pid in range(9)]
                t.records += [t.record(pid) for pid in range(9, 598)]
                t.records += [t.record(598, status="error" if failed else "ok")]
                t.write(done=True, stopped=False)
                result = A.audit_chunk(t.run, t.chunk)
                self.assertTrue(result["certified_done_chunk"], result)
                item = result["continuation"][0]
                self.assertEqual(item["raw_recomputed"]["exact_decision"]["strict"]["pair_reject"], [False] * 3)
                self.assertEqual(item["raw_recomputed"]["exact_decision"]["le"]["pair_reject"], [None if failed else True] * 3)
                self.assertEqual(item["full_p_available"], not failed)
                self.assertEqual(item["required_resume_id_ranges"], [])

    def test_old_rows_with_new_raw_and_progress_is_deferred_sync_skew(self):
        t = self.toy
        old_records = copy.deepcopy(t.records)
        t.records.append(t.record(4))
        t.write(derived_records=old_records)
        result = A.audit_chunk(t.run, t.chunk)
        self.assertEqual(result["status"], "deferred", result)
        self.assertFalse(result["checkpoint_certified"])
        self.assertEqual(result["findings_classification"], "requires_consistent_snapshot_before_corruption_claim")

    def test_stable_impossible_gt_count_remains_invalid(self):
        t = self.toy
        t.records.append(t.record(4))
        t.write()
        path = t.chunk / "rows.jsonl"
        row = json.loads(path.read_text())
        row["schemes"]["diagnosis"]["gt"][0] = 99
        path.write_text(json.dumps(row) + "\n")
        result = A.audit_chunk(t.run, t.chunk)
        self.assertEqual(result["status"], "invalid", result)

    def test_full_p_unfinished_done_is_not_certified(self):
        t = self.toy
        t.records.append(t.record(4))
        t.write(done=True, stopped=False)
        result = A.audit_chunk(t.run, t.chunk)
        self.assertFalse(result["certified_done_chunk"])
        self.assertEqual(result["continuation"][0]["action"], "resume_prespecified_full_p")

    def test_observed_mode_cannot_contradict_engine_identity(self):
        t = self.toy
        t.records[0]["diagnostics"]["mode"] = "original"
        t.write()
        result = A.audit_chunk(t.run, t.chunk)
        self.assertFalse(result["checkpoint_certified"], result)

    def test_nonfull_exact_earlystop_needs_no_missing_ids(self):
        t = self.toy
        t.run["full_p"] = set()
        t.write(done=True, stopped=False)
        result = A.audit_chunk(t.run, t.chunk)
        self.assertTrue(result["certified_done_chunk"], result)
        self.assertEqual(result["continuation"][0]["action"], "exact_early_stop_no_more_permutations_required")
        self.assertEqual(result["continuation"][0]["required_resume_id_ranges"], [])

    def test_null_truth_is_reported_without_losing_failed_identity(self):
        t = self.toy
        t.records = [t.record(None, status="error")]
        t.records[0]["truth"] = None
        t.write(done=True, stopped=False)
        result = A.audit_chunk(t.run, t.chunk)
        self.assertEqual(len(result["failed_observed"]), 1)
        self.assertEqual(result["failed_observed"][0]["status"], "error")
        self.assertFalse(result["missing_observed"])

    def test_predictable_temp_alias_cannot_write_through_raw(self):
        t = self.toy
        t.write()
        raw = t.chunk / "fit_records.jsonl"
        original = raw.read_bytes()
        protocol = t.root / "protocol.json"
        protocol.write_text("{}")
        prefix = t.root / "audit"
        trap = prefix.with_name(f"audit.json.tmp-{os.getpid()}")
        trap.symlink_to(raw)
        result = A.audit_chunk(t.run, t.chunk)
        with patch.object(A, "analyze", return_value=result), redirect_stdout(io.StringIO()):
            try:
                A.main(["--protocol", str(protocol), "--run-id", "toy", "--chunk", str(t.chunk),
                        "--project-root", str(t.root), "--output-prefix", str(prefix)])
            except SystemExit:
                pass  # Rejecting unsafe output also preserves the input.
        self.assertEqual(raw.read_bytes(), original)

    def test_registered_snapshot_outside_default_directory_is_read_only(self):
        # Reuse only the author's inert signed-file builder to exercise the
        # full CLI loader with a legitimate nondefault registered snapshot.
        from test_audit_confirmation_events import FrozenToy
        t = FrozenToy(Path(self.temp.name) / "registered")
        t.records = [t.event(), t.event(0), t.event(1)]
        t.write()
        before = t.config.read_bytes()
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            try:
                A.main(["--protocol", str(t.protocol), "--run-id", "toy", "--chunk", str(t.chunk),
                        "--project-root", str(t.root), "--output-prefix", str(t.config.with_suffix(""))])
            except SystemExit:
                pass
        self.assertEqual(t.config.read_bytes(), before)
        self.assertFalse(t.config.with_suffix(".md").exists())

    def test_protocol_change_during_audit_cannot_certify_new_hash(self):
        t = self.toy
        t.write()
        protocol = t.root / "protocol.json"
        protocol.write_text(json.dumps({"configuration_frozen_before_confirmation_data": True,
            "confirmation_namespace": "42xxxxxx", "runs": [{"run_id": "toy",
                "snapshot_sha256": A.C.sha(t.run["snapshot"] / "SHA256SUMS"),
                "config_sha256": A.C.sha(t.run["config_path"])}]}))
        original_protocol_hash = A.C.sha(protocol)
        real_audit = A.audit_chunk
        def changed(run, chunk):
            result = real_audit(run, chunk)
            protocol.write_text('{"replacement": true}')
            return result
        with patch.object(A.C, "frozen_run", return_value=t.run):
            baseline = A.analyze(protocol, "toy", t.chunk, t.root)
            self.assertTrue(baseline["checkpoint_certified"], baseline)
            with patch.object(A, "audit_chunk", side_effect=changed):
                try:
                    result = A.analyze(protocol, "toy", t.chunk, t.root)
                except ValueError:
                    return  # A fail-closed refusal also binds certification safely.
        self.assertFalse(result["checkpoint_certified"], result)
        self.assertEqual(result["protocol"]["sha256"], original_protocol_hash)
        changed_paths = [p for reason in result["deferred_reasons"] for p in reason.get("paths", [])]
        self.assertIn(str(protocol), changed_paths)


if __name__ == "__main__":
    unittest.main()
