"""Reuse the reviewed helper's fake-Slurm fixtures with exact target substitutions."""
import hashlib
import json
from pathlib import Path
import types
import unittest

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
VERSION = PROJECT / "runs/v2/campaign_versions/confirmation_20260908_final_queue_fill"


class FixedDerivationTests(unittest.TestCase):
    def test_only_recorded_constant_and_explanation_changes(self):
        preparation = json.loads((VERSION / "preparation.json").read_text())
        base = (HERE / "reserve_core_pilot_slot_20260908.sh").read_bytes()
        self.assertEqual(hashlib.sha256(base).hexdigest(), preparation["base_helper_sha256"])
        source = base.decode()
        for change in preparation["helper_exact_text_replacements"]:
            self.assertEqual(source.count(change["old"]), 1)
            source = source.replace(change["old"], change["new"], 1)
        derived = (HERE / "reserve_final_queue_slot_20260908.sh").read_bytes()
        self.assertEqual(source.encode(), derived)
        self.assertEqual(hashlib.sha256(derived).hexdigest(), preparation["helper_sha256"])

    def test_stage_only_proposal_preserves_caps_and_existing_policy(self):
        proposal = json.loads((VERSION / "campaign.json").read_text())
        self.assertEqual(proposal["max_tasks_per_submission"], 900)
        self.assertEqual([r["run_id"] for r in proposal["runs"]], ["confirm_stage_a"])
        stage = proposal["runs"][0]
        self.assertEqual((stage["max_parallel_tasks"], stage["cpus_per_task"], stage["nchunks"], stage["walltime_s"]), (1, 16, 1000, 3600))
        self.assertEqual(stage["execution_walltime_policy"], "frozen_single_engine_signal_1h_v1")
        self.assertEqual([proposal[k] for k in ("max_project_cpus", "max_user_submitted_jobs", "max_array_size", "queue_max_age_s")], [1024, 900, 1000, 60])


def load_tests(loader, tests, pattern):
    # The original nine tests already cover identity, two pending checks,
    # stale-after-fsync, readback failure, unknown/no-retry and one-shot guards.
    # Only their fixture's exact target values change; no Slurm is invoked.
    source = (HERE / "test_reserve_core_pilot_slot.py").read_text()
    changes = [
        ("0661eee6714148115c1ae9789e84ec63590e425416276a0a0302017ab83005a6", "c997a2e98a610c02aabd45acfb81f8396cdcd7507ee61ff4b2219925b1e380bc"),
        ("10265062", "10264903"),
        ("MEMBERS = [25] + list(range(430,445))", "MEMBERS = list(range(414,430))"),
        ("confirmation_20260908_core_paired_1h_pilot/reserve_core_pilot_slot.json", "confirmation_20260908_final_queue_fill/reserve_final_queue_slot.json"),
        ('"dispatch_generation":5', '"dispatch_generation":4'),
        ("reserve_core_pilot_slot_20260908.sh", "reserve_final_queue_slot_20260908.sh"),
        ("ArrayTaskId=25,430-444%", "ArrayTaskId=414-429%"),
        ("_25", "_414"),
        ('["25"]', '["414"]'),
        ('PARENT+"_[25]"', 'PARENT+"_[414]"'),
    ]
    for old, new in changes:
        if old not in source:
            raise AssertionError("Base fixture changed; review adaptation: " + old)
        source = source.replace(old, new)
    fixture = types.ModuleType("final_queue_fake_slurm_fixtures")
    fixture.__file__ = str(HERE / "test_reserve_core_pilot_slot.py")
    exec(compile(source, fixture.__file__, "exec"), fixture.__dict__)
    tests.addTests(loader.loadTestsFromModule(fixture))
    return tests


if __name__ == "__main__":
    unittest.main()
