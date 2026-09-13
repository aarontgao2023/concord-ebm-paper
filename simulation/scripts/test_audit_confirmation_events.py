"""Signed temporary snapshots and invented 3-event logs only; no scientific calls."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import audit_confirmation_events as A


def exact(exceed, n, B, comparisons, rule):
    threshold = Fraction(1, 20 * comparisons)
    values = [Fraction(1 + exceed + extra, B + 1) for extra in (0, B - n)]
    decisions = [value < threshold if rule == "strict" else value <= threshold for value in values]
    return decisions[0] if decisions[0] == decisions[1] else None


class FrozenToy:
    def __init__(self, root, B=2, *, full=True, paired=False, R=1, local=False):
        self.root = Path(root)
        self.snapshot = self.root / "snapshots/toy"
        source = self.snapshot / "scripts/v2"
        source.mkdir(parents=True)
        for name in A.I.RUNTIME_FILES + A.I.OPTIONAL_RUNTIME_FILES:
            shutil.copyfile(A.C.PROJECT / "scripts/v2" / name, source / name)
        defaults = A.C.design_defaults(source / "design_v2.py")
        self.cell = {"name": "REF_H0", "overrides": {"biomarker_names": defaults["biomarker_names"][:3]}}
        self.cfg = {"run_id": "toy", "phase": "confirmation", "base_seed": 42999990,
                    "datasets": R, "cells": [self.cell], "engines": ["original", "repaired"] if paired else ["repaired"],
                    "bperm": B, "alpha": .05, "rule": "le", "full_p_first_n": R if full else 0,
                    "paired_standard": paired, "fast_likelihood": True,
                    "schemes": {"D": {"stratify": "diagnosis", "require_max": not local}}}
        if local:
            self.cfg["schemes"]["D"].update(groups=[0, 1], tested_pair_index=0)
        self.config = self.snapshot / "config.json"
        self.config.write_text(json.dumps(self.cfg))
        sums = self.snapshot / "SHA256SUMS"
        files = sorted(p for p in self.snapshot.rglob("*") if p.is_file())
        sums.write_text("".join(f"{A.C.sha(p)}  {p.relative_to(self.snapshot)}\n" for p in files))
        self.entry = {"run_id": "toy", "config": "config.json", "config_sha256": A.C.sha(self.config),
                      "snapshot": "snapshots/toy", "snapshot_sha256": A.C.sha(sums), "base_seed": self.cfg["base_seed"],
                      "datasets_per_cell": R, "cells": self.cfg["cells"], "engines": self.cfg["engines"], "bperm": B,
                      "full_p_first_n": self.cfg["full_p_first_n"], "full_p_seeds": None,
                      "paired_standard": paired, "fast_likelihood": True}
        self.protocol = self.root / "protocol.json"
        self.protocol.write_text(json.dumps({"configuration_frozen_before_confirmation_data": True,
            "confirmation_namespace": "42xxxxxx", "runs": [self.entry]}))
        self.run = A.C.frozen_run(self.entry, self.root)
        self.env = {"versions": {name: "toy" for name in A.I.NUMERICAL_PACKAGES},
                    "source_sha256": self.run["dependency"]["SOURCE_SHA256"],
                    "wheel_sha256": self.run["dependency"]["PYEBM_WHEEL_SHA256"]}
        self.env["versions"]["pyebm"] = self.run["dependency"]["PYEBM_VERSION"]
        self.fast = A.C.literal_assignments(source / "fast_likelihood_v2.py")["FAST_LIKELIHOOD_VERSION"]
        self.chunk = self.root / "raw/toy/cell_0/chunk_0"
        self.chunk.mkdir(parents=True)
        signed = {"config": self.cfg, "cell_index": 0, "chunk": 0, "nchunks": 1, "source_sha256": self.run["sources"]}
        self.manifest = signed | {"signature": A.I.digest(signed), "seeds": sorted(self.run["seeds"]), "environment": self.env}
        (self.chunk / "manifest.json").write_text(json.dumps(self.manifest))
        self.records = []

    def event(self, pid=None, *, seed=None, engine="repaired", status="ok", orderings=None):
        seed = self.cfg["base_seed"] if seed is None else seed
        orders = orderings or [[0, 1, 2] for _ in range(3)]
        record = {"kind": "observed" if pid is None else "permutation", "seed": seed, "engine": engine,
                  "scheme": None if pid is None else "D", "perm_id": pid, "status": status,
                  "orderings": orders if status == "ok" else None,
                  "taus": [A.C.inversions(orders[a], orders[b]) / 3 for a, b in A.C.PAIRS] if status == "ok" else None,
                  "fit_seconds": 1., "likelihood_implementation": self.fast,
                  "diagnostics": {"provenance": self.env, "mode": engine}}
        if self.cfg["paired_standard"]:
            shared = A.shared_id(self.cfg, self.cell, seed, record["scheme"], pid)
            record["shared_fit_id"] = shared
            record["diagnostics"]["shared_fit_id"] = shared
        if pid is None:
            design = self.run["designs"]["REF_H0"]
            record["distance_to_truth"] = [A.C.inversions(x, [0, 1, 2]) / 3 for x in orders] if status == "ok" else None
            record["truth"] = {"config": design, "seed": seed,
                "manifest": {"seed": seed, "config_sha256": A.C.canonical_sha(design), "source_sha256": self.run["sources"]["design_v2.py"],
                             "data_sha256": hashlib.sha256(str(seed).encode()).hexdigest(), "numpy_version": "toy", "pandas_version": "toy"},
                "biomarker_names": design["biomarker_names"], "group_order": list(A.C.GROUPS), "base_order": [0, 1, 2],
                "group_orderings": {g: [0, 1, 2] for g in A.C.GROUPS}, "planned_inversions": {g: 0 for g in A.C.GROUPS},
                "realized_inversions": {g: 0 for g in A.C.GROUPS}, "pair_null": [True] * 3, "h1_distance": 0.,
                "pair_truth": {name: {"inversions": 0, "null": True, "normalized_distance": 0.} for name in A.C.PAIR_NAMES}}
        return record

    def rows(self):
        result = []
        for obs in [r for r in self.records if r["scheme"] is None]:
            schemes = {}
            if obs["status"] == "ok":
                perms = [r for r in self.records if r["scheme"] == "D" and r["seed"] == obs["seed"] and r["engine"] == obs["engine"]]
                good = [r for r in perms if r["status"] == "ok"]
                n, B = len(good), self.cfg["bperm"]
                gt = [sum(r["taus"][i] > obs["taus"][i] + 1e-9 for r in good) for i in range(3)]
                eq = [sum(abs(r["taus"][i] - obs["taus"][i]) <= 1e-9 for r in good) for i in range(3)]
                mg = sum(max(r["taus"]) > max(obs["taus"]) + 1e-9 for r in good)
                me = sum(abs(max(r["taus"]) - max(obs["taus"])) <= 1e-9 for r in good)
                decisions = {"alpha": .05}
                for rule in ("strict", "le"):
                    decisions[rule] = {"pair_reject": [exact(g+e, n, B, 3, rule) for g, e in zip(gt, eq)],
                                       "max_reject": exact(mg+me, n, B, 1, rule)}
                schemes["D"] = {"requested_nperm": B, "nperm": n, "attempted": len(perms), "failed": len(perms)-n,
                    "complete": n == B, "gt": gt, "eq": eq, "maxgt": mg, "maxeq": me, "exact_decision": decisions,
                    "definition": self.cfg["schemes"]["D"]}
            result.append({"run_id": "toy", "phase": "confirmation", "cell": "REF_H0", "seed": obs["seed"],
                "engine": obs["engine"], "status": obs["status"], "truth": obs.get("truth", {}), "schemes": schemes,
                "observed": {k: obs.get(k) for k in ("orderings", "taus", "distance_to_truth", "fit_seconds")},
                "diagnostics": obs["diagnostics"] | {"sequential_stopping": True, "full_p_selected": obs["seed"] in self.run["full_p"]}})
        return result

    def write(self, *, done=True, stopped=False):
        for name, records in (("fit_records", self.records), ("rows", self.rows())):
            (self.chunk / (name + ".jsonl")).write_text("".join(json.dumps(r) + "\n" for r in records))
        keys = {r.get("shared_fit_id", (r["seed"], r["engine"], r["scheme"], r["perm_id"])) for r in self.records}
        progress = {"done": done, "stopped": stopped, "fit_records": len(self.records), "rows": len(self.rows()),
                    "expected_rows": len(self.run["seeds"])*len(self.cfg["engines"]), "computational_jobs": len(keys),
                    "statuses": dict(Counter(r["status"] for r in self.records))}
        (self.chunk / "progress.json").write_text(json.dumps(progress))

    def audit(self):
        return A.analyze(self.protocol, "toy", self.chunk, self.root)


class FrozenEventAuditTests(unittest.TestCase):
    def toy(self, **kwargs):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        return FrozenToy(temp.name, **kwargs)

    def complete(self, **kwargs):
        t = self.toy(**kwargs)
        for seed in sorted(t.run["seeds"]):
            for engine in t.cfg["engines"]:
                t.records += [t.event(seed=seed, engine=engine)]
                t.records += [t.event(pid, seed=seed, engine=engine) for pid in range(t.cfg["bperm"])]
        t.write()
        return t

    def test_signed_complete_stream_is_read_only_and_small_memory(self):
        t = self.complete()
        paths = list(t.chunk.iterdir())
        before = {p: p.read_bytes() for p in paths}
        r = t.audit()
        self.assertTrue(r["certified_done_chunk"], r)
        self.assertEqual(r["raw_recomputed_progress"]["fit_records"], 3)
        self.assertEqual(r["continuation"][0]["raw_recomputed"]["eq"], [2, 2, 2])
        self.assertEqual(r["memory_model"]["full_permutation_records_retained"], 0)
        self.assertEqual(r["memory_model"]["permutation_id_bytes"], 2)
        self.assertEqual({p: p.read_bytes() for p in paths}, before)
        self.assertEqual(r["protocol"]["sha256"], A.C.sha(t.protocol))

    def test_permutations_before_observed_hand_counts(self):
        t = self.toy()
        t.records = [t.event(0, orderings=[[1, 0, 2], [0, 1, 2], [0, 1, 2]]),
                     t.event(1, orderings=[[0, 1, 2], [2, 1, 0], [0, 1, 2]]), t.event()]
        t.write()
        r = t.audit()
        self.assertTrue(r["certified_done_chunk"], r)
        s = r["continuation"][0]["raw_recomputed"]
        self.assertEqual((s["gt"], s["eq"], s["maxgt"], s["maxeq"]), ([2, 1, 1], [0, 1, 1], 2, 0))

    def test_paired_jobs_once_statistical_fits_per_engine(self):
        t = self.complete(paired=True, R=2)
        r = t.audit()
        self.assertTrue(r["certified_done_chunk"], r)
        self.assertEqual((r["raw_recomputed_progress"]["fit_records"], r["raw_recomputed_progress"]["computational_jobs"]), (12, 6))
        self.assertEqual([s["successful"] for s in r["continuation"]], [2]*4)

    def test_wrong_shared_identity_rejected(self):
        t = self.complete(paired=True)
        t.records[1]["shared_fit_id"] = "a"*64
        t.write()
        self.assertIn("invalid_paired_shared_fit_identity", t.audit()["finding_counts"])

    def test_duplicate_failed_record_is_never_replaced(self):
        t = self.toy()
        t.records = [t.event(), t.event(0, status="error"), t.event(1)]
        t.write()
        with (t.chunk / "fit_records.jsonl").open("a") as out:
            out.write(json.dumps(t.event(0)) + "\n")
        r = t.audit()
        self.assertEqual(r["status"], "invalid", r)
        self.assertEqual(r["continuation"][0]["failed_permutation_id_ranges"], [[0, 0]])
        self.assertEqual(r["continuation"][0]["successful"], 1)

    def test_failed_observed_missing_truth_is_preserved(self):
        t = self.toy()
        t.records = [t.event(status="error")]
        t.records[0]["truth"] = {}
        t.write()
        r = t.audit()
        self.assertTrue(r["certified_done_chunk"], r)
        self.assertEqual(r["continuation"][0]["action"], "observed_failure_preserved_no_retry")
        self.assertEqual(r["continuation"][0]["required_resume_id_ranges"], [])
        self.assertEqual(len(r["failed_observed"]), 1)

    def test_nonprefix_ids_compress_missing_ranges_and_never_retry_failure(self):
        t = self.toy(B=10)
        t.records = [t.event(), t.event(0), t.event(3, status="error"), t.event(7)]
        t.write(done=False, stopped=True)
        r = t.audit()
        self.assertTrue(r["checkpoint_certified"], r)
        s = r["continuation"][0]
        self.assertEqual(s["required_resume_id_ranges"], [[1, 2], [4, 6], [8, 9]])
        self.assertEqual(s["failed_permutation_id_ranges"], [[3, 3]])

    def test_local_exact_stop_does_not_wait_for_untested_max(self):
        t = self.toy(B=599, full=False, local=True)
        t.records = [t.event()] + [t.event(i) for i in range(10)]
        t.write()
        r = t.audit()
        self.assertTrue(r["certified_done_chunk"], r)
        s = r["continuation"][0]
        self.assertEqual(s["action"], "exact_early_stop_no_more_permutations_required")
        self.assertEqual(s["missing_permutation_id_ranges"], [[10, 598]])
        self.assertIsNone(s["raw_recomputed"]["exact_decision"]["le"]["max_reject"])

    def test_same_resolved_decision_still_requires_preselected_full_budget(self):
        t = self.toy(B=599, full=True, local=True)
        t.records = [t.event()] + [t.event(i) for i in range(10)]
        t.write(done=False, stopped=True)
        r = t.audit()
        self.assertTrue(r["checkpoint_certified"], r)
        self.assertEqual(r["continuation"][0]["required_resume_id_ranges"], [[10, 598]])

    def test_success_with_forged_tau_or_illegal_ordering_invalid(self):
        for mode in ("tau", "duplicate", "boolean"):
            with self.subTest(mode=mode):
                t = self.complete()
                if mode == "tau":
                    t.records[1]["taus"][0] = .5
                else:
                    t.records[1]["orderings"][0] = [0, 0, 2] if mode == "duplicate" else [False, 1, 2]
                (t.chunk / "fit_records.jsonl").write_text("".join(json.dumps(r)+"\n" for r in t.records))
                r = t.audit()
                self.assertFalse(r["checkpoint_certified"])
                self.assertIn("invalid_successful_ordering_or_tau", r["finding_counts"])

    def test_foreign_seed_and_out_of_budget_pid_never_accepted(self):
        for changes in ({"seed": 42999989}, {"perm_id": 2}, {"engine": "foreign"}, {"scheme": "foreign"}, {"kind": "observed"}):
            with self.subTest(changes=changes):
                t = self.complete()
                t.records[1].update(changes)
                (t.chunk / "fit_records.jsonl").write_text("".join(json.dumps(r)+"\n" for r in t.records))
                self.assertIn("invalid_raw_record", t.audit()["finding_counts"])

    def test_unsigned_snapshot_or_wrong_manifest_rejected(self):
        t = self.complete()
        with (t.snapshot / "scripts/v2/run_v2.py").open("a") as f:
            f.write("\n# changed\n")
        with self.assertRaisesRegex(ValueError, "frozen file mismatch"):
            t.audit()
        t = self.complete()
        t.manifest["seeds"] = [42999989]
        (t.chunk / "manifest.json").write_text(json.dumps(t.manifest))
        with self.assertRaisesRegex(ValueError, "allocation"):
            t.audit()

    def test_uncommitted_tail_is_deferred_and_never_repaired(self):
        t = self.complete()
        p = t.chunk / "fit_records.jsonl"
        p.write_bytes(p.read_bytes().rstrip(b"\n"))
        before = p.read_bytes()
        r = t.audit()
        self.assertEqual(r["status"], "deferred")
        self.assertEqual(p.read_bytes(), before)
        self.assertFalse(r["inputs"][p.name]["final_newline"])

    def test_live_lock_or_flags_prevent_certification(self):
        t = self.complete()
        with patch.object(A, "lock_probe", return_value="held_by_writer"):
            self.assertEqual(t.audit()["status"], "deferred")
        t.write(done=False, stopped=False)
        self.assertEqual(t.audit()["status"], "deferred")

    def test_missing_raw_or_metadata_is_deferred(self):
        for name in ("rows.jsonl", "progress.json", "fit_records.jsonl", "manifest.json"):
            with self.subTest(name=name):
                t = self.complete()
                (t.chunk / name).unlink()
                self.assertEqual(t.audit()["status"], "deferred")

    def test_progress_mismatch_is_sync_deferral(self):
        t = self.complete()
        p = t.chunk / "progress.json"
        progress = json.loads(p.read_text())
        progress["fit_records"] -= 1
        p.write_text(json.dumps(progress))
        r = t.audit()
        self.assertEqual(r["status"], "deferred")
        self.assertFalse(r["continuation_manifest_available"])

    def test_protocol_changed_mid_audit_retains_original_hash_and_defers(self):
        t = self.complete()
        original = A.C.sha(t.protocol)
        real = A.audit_chunk
        def change_protocol(run, chunk):
            result = real(run, chunk)
            t.protocol.write_text('{"replacement":true}')
            return result
        with patch.object(A, "audit_chunk", side_effect=change_protocol):
            r = t.audit()
        self.assertEqual(r["status"], "deferred")
        self.assertFalse(r["checkpoint_certified"])
        self.assertEqual(r["protocol"]["sha256"], original)

    def test_malformed_scheme_containers_are_findings_not_tracebacks(self):
        for value in ([], {"exact_decision": {"strict": []}}):
            with self.subTest(value=value):
                t = self.complete()
                p = t.chunk / "rows.jsonl"
                row = json.loads(p.read_text())
                if isinstance(value, list):
                    row["schemes"]["D"] = value
                else:
                    row["schemes"]["D"].update(value)
                p.write_text(json.dumps(row) + "\n")
                r = t.audit()
                self.assertEqual(r["status"], "invalid", r)
                self.assertIn("invalid_derived_scheme_types_or_decisions", r["finding_counts"])

    def test_nonhash_data_identity_is_a_finding_not_traceback(self):
        t = self.complete()
        t.records[0]["truth"]["manifest"]["data_sha256"] = ["malformed"]
        t.write()
        r = t.audit()
        self.assertEqual(r["status"], "invalid")
        self.assertIn("invalid_observed_truth_or_provenance", r["finding_counts"])

    def test_progress_boolean_cannot_impersonate_integer_status_count(self):
        t = self.toy()
        t.records = [t.event(status="error")]
        t.write()
        p = t.chunk / "progress.json"
        progress = json.loads(p.read_text())
        progress["statuses"]["error"] = True
        p.write_text(json.dumps(progress))
        self.assertEqual(t.audit()["status"], "invalid")


if __name__ == "__main__":
    unittest.main()
