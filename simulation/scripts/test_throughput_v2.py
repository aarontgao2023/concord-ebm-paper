"""Throughput denominator, paired-cost, censoring and Slurm-margin regressions."""
import json
import unittest

import analyze_throughput_v2 as T
import test_development_record_audits as F


class ThroughputTests(unittest.TestCase):
    def setUp(self):
        self.fixture = F.DevelopmentAuditTests()
        self.fixture.setUp()

    def tearDown(self):
        self.fixture.tearDown()

    def run_fixture(self, missing=None):
        cfg = self.fixture.cfg("throughput", True) | {"bperm": 99, "full_p_first_n": 0, "paired_standard": True}
        config, raw = self.fixture.write_run(cfg, bridge=True, missing=missing)
        # This seed plan intentionally has no forced full-p cohort.
        path = raw / "cell_0/chunk_0/rows.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        for row in rows:
            row["diagnostics"]["full_p_selected"] = False
        path.write_text("".join(json.dumps(row) + "\n" for row in rows))
        return config, raw

    def test_complete_cohort_cost_uses_shared_seconds_once_and_all_planned_datasets(self):
        config, raw = self.run_fixture()
        result = T.analyze(config, raw, forecast_datasets=10, forecast_workers=2)
        self.assertTrue(result["complete_for_throughput_estimation"], result["issues"])
        self.assertEqual(result["unambiguous_statistical_fit_records"], 16)
        self.assertEqual(result["computation"]["unique_computational_jobs"], 8)
        self.assertEqual(result["computation"]["recorded_unique_job_seconds"], 16.)
        self.assertEqual(result["complete_cohort_savings"]["statistical_fit_fraction_avoided"], .96)
        self.assertEqual(result["future_run_planning_scenario"]["ideal_wall_seconds"], 40.)
        scheme = next(g for g in result["cell_engine_schemes"] if g["scheme"] == "diagnosis")
        self.assertEqual((scheme["planned_datasets"], scheme["recorded_started_datasets"],
                          scheme["decision_resolved_datasets"], scheme["finished_datasets"], scheme["full_budget_datasets"]),
                         (2, 2, 2, 2, 0))

    def test_incomplete_run_has_no_savings_or_eta_from_completed_subset(self):
        config, raw = self.run_fixture(missing=lambda r: r["seed"] == 101)
        result = T.analyze(config, raw, forecast_datasets=10, forecast_workers=2)
        self.assertTrue(result["right_censored_or_incomplete"])
        self.assertIsNone(result["complete_cohort_savings"])
        self.assertIsNone(result["future_run_planning_scenario"])
        self.assertTrue(all(g["complete_cohort_mean_statistical_fit_fraction_avoided"] is None for g in result["cell_engine_schemes"]))

    def test_timeout_durations_are_counted_and_row_log_failure_mismatch_blocks_forecast(self):
        config, raw = self.run_fixture()
        path = raw / "cell_0/chunk_0/fit_records.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        for row in rows:
            if row["seed"] == 100 and row["perm_id"] == 2:
                row.update(status="error", fit_seconds=300., diagnostics={"error": "FitTimeout: limit"})
        path.write_text("".join(json.dumps(r) + "\n" for r in rows))
        result = T.analyze(config, raw, forecast_datasets=10, forecast_workers=2)
        self.assertEqual(result["failed_statistical_fit_records"], 2)
        scheme = next(g for g in result["cell_engine_schemes"] if g["scheme"] == "diagnosis")
        self.assertEqual(scheme["fit_seconds_all_recorded_attempts"]["max"], 300.)
        self.assertEqual(scheme["fit_seconds_all_recorded_attempts"]["n"], 6)
        self.assertIsNone(result["future_run_planning_scenario"])
        self.assertTrue(any(i["code"] == "scheme_row_log_count_mismatch" for i in result["issues"]))

    def test_replayed_half_pair_cost_is_not_a_complete_cpu_consumption_estimate(self):
        config, raw = self.run_fixture()
        replay = raw / "cell_0/chunk_0/replay_audit.jsonl"
        replay.write_text('{"policy":"keep_original_record"}\n')
        result = T.analyze(config, raw, forecast_datasets=10, forecast_workers=2)
        self.assertEqual(result["paired_replay_audit_entries"], 1)
        self.assertEqual(result["computation"]["unique_computational_jobs"], 8)
        self.assertIsNone(result["future_run_planning_scenario"])

    def test_exhausted_timeout_budget_stays_unknown_instead_of_using_successful_denominator(self):
        config, raw = self.run_fixture()
        cfg = json.loads(config.read_text())
        chunk = raw / "cell_0/chunk_0"
        observed = [self.fixture.fit(cfg, seed, engine) for seed in (100, 101) for engine in cfg["engines"]]
        for record in observed:
            record["taus"] = [.5] * 3
        permutations = [self.fixture.fit(cfg, seed, engine, "diagnosis", pid)
                        for seed in (100, 101) for engine in cfg["engines"] for pid in range(99)]
        for record in permutations:
            if record["seed"] == 100 and record["perm_id"] == 98:
                record.update(status="error", fit_seconds=300., diagnostics={"error": "FitTimeout"})
        (chunk / "fit_records.jsonl").write_text("".join(json.dumps(r) + "\n" for r in observed + permutations))
        rows = [json.loads(line) for line in (chunk / "rows.jsonl").read_text().splitlines()]
        for row in rows:
            row["observed"]["taus"] = [.5] * 3
            failed = row["seed"] == 100
            row["schemes"]["diagnosis"].update(nperm=98 if failed else 99, attempted=99, failed=int(failed),
                                                complete=not failed, gt=[0] * 3, eq=[0] * 3)
        (chunk / "rows.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
        result = T.analyze(config, raw)
        self.assertEqual(result["issues"], [])
        self.assertTrue(result["inventory"]["collection_complete"])
        self.assertFalse(result["complete_for_throughput_estimation"])
        self.assertEqual(result["failed_statistical_fit_records"], 2)
        scheme = next(g for g in result["cell_engine_schemes"] if g["scheme"] == "diagnosis")
        self.assertEqual(scheme["states"], {"terminal_unresolved": 1, "full": 1})
        self.assertEqual(scheme["decision_resolved_datasets"], 1)
        self.assertEqual(scheme["statistical_fit_records"], 198)
        self.assertEqual(result["computation"]["recorded_unique_job_seconds"], 698.)
        self.assertIsNone(result["complete_cohort_savings"])

    def test_slurm_drain_uses_earliest_trigger_and_exposes_large_paired_timeout(self):
        independent = T.drain_budget(300)
        paired = T.drain_budget(300, True)
        mechanism = T.drain_budget(600, True)
        self.assertEqual(independent["nominal_drain_s"], 600)
        self.assertEqual(paired["signal_only_remaining_after_drain_s"], 300)
        self.assertEqual(mechanism["nominal_drain_s"], 2400)
        self.assertEqual(mechanism["runtime_only_remaining_after_drain_s"], 200)
        self.assertEqual(mechanism["signal_only_remaining_after_drain_s"], -900)
        self.assertEqual(T.drain_budget(600, True, startup=250)["remaining_after_nominal_drain_s"], -50)

    def test_sacct_does_not_double_count_allocation_and_batch_step(self):
        path = self.fixture.base / "sacct.txt"
        path.write_text("JobIDRaw|State|CPUTimeRAW|TotalCPU|\n1|COMPLETED|100|00:01:20|\n1.batch|COMPLETED|100|00:01:20|\n")
        summary = T.sacct_summary(path)
        self.assertEqual(summary["allocated_cpu_seconds"], 100.)
        self.assertEqual(summary["measured_cpu_utilization"], .8)
        self.assertEqual(summary["excluded_step_ids"], ["1.batch"])
        self.assertIsNone(T.sacct_summary(None))


if __name__ == "__main__":
    unittest.main()
