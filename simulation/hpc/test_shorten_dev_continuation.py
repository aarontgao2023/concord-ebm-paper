"""Run only the helper's Python body against a temporary root and fake Slurm."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M

PLAN = "6d9d64a72f26dad6f28d1ae32deb257a2771b96969384d6405762deb774ba495"
MEMBER = "10264821_13"


class ShortenDevelopmentTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = (Path(self.tmp.name)/"ebmcal_v2").resolve()
        self.ledger = self.root/"control/submission_ledger.json"
        self.evidence = self.root/"control/campaign_versions/confirmation_20260908_dev_continue_1h/pending_dev_walltime_amendment.json"
        self.evidence.parent.mkdir(parents=True)
        self.snapshot = self.root/"snapshots/dev_calibration_a"
        self.snapshot.mkdir(parents=True)
        self.cfg = self.snapshot/"config.json"
        self.cfg.write_text('{}')
        self.receipt = {"run_id": "dev_calibration_a", "snapshot_sha256": M.DEV_CALIBRATION_1H_PIN[0],
            "array_task_ids": [13], "nchunks": 8, "mapping": "cell_major", "cpus_per_task": 16,
            "max_concurrent": 1, "walltime_s": 14400, "dispatch_generation": 2,
            "state": "submitted", "plan_id": PLAN, "job_id": "10264821"}
        self.assertEqual(M.I.digest(D.submission_identity(self.receipt)), PLAN)
        self.data = {"schema_version": "campaign_ledger_v1", "project": "ebmcal_v2", "submissions": [self.receipt]}
        self.ledger.write_text(json.dumps(self.data))
        self.binding = self.root/"campaign.dispatch.binding.json"
        self.binding.write_text(json.dumps({"project": "ebmcal_v2", "ledger_path": str(self.ledger)}))
        source = (Path(__file__).parent/"shorten_dev_continuation_20260908.sh").read_text()
        body = source.split("<<'PY'\n", 1)[1].rsplit("\nPY", 1)[0]
        body = body.replace("Path('/N/slate/tg11/ebmcal_v2').resolve()", f"Path({str(self.root)!r}).resolve()", 1)
        # Safety checks must survive exactly the mode that removes assert.
        self.code = compile(body, "temporary_shorten_helper", "exec", optimize=2)
        self.before = f"{MEMBER}|PENDING|4:00:00|ebmcal_v2:{PLAN}\n999999_4|PENDING|4:00:00|phase-test:unrelated\n"
        self.after = f"{MEMBER}|RUNNING|1:00:00|ebmcal_v2:{PLAN}\n999999_4|PENDING|4:00:00|phase-test:unrelated\n"
        self.calls = []
        self.queries = 0
        self.now = 0.0
        self.timeout_update = False
        self.update_returncode = 0
        self.slow_intent = False
        self.real_durable = D.durable_json

    def fake_subprocess(self, argv, **kwargs):
        self.calls.append(argv)
        if argv[0] == "squeue":
            self.queries += 1
            self.assertEqual(argv, ['squeue','-r','-u','tg11','-h','-o','%i|%T|%l|%k'])
            return subprocess.CompletedProcess(argv, 0, self.before if self.queries == 1 else self.after, '')
        self.assertEqual(argv, ['scontrol','update','JobId='+MEMBER,'TimeLimit=01:00:00'])
        saved = json.loads(self.evidence.read_text())
        self.assertEqual(saved['status'], 'update_unknown')
        if self.timeout_update:
            raise subprocess.TimeoutExpired(argv, 30)
        return subprocess.CompletedProcess(argv, self.update_returncode, '', 'mock response')

    def durable(self, path, value):
        self.real_durable(path, value)
        if self.slow_intent and path == self.evidence and value['status'] == 'update_unknown':
            self.now = 61

    def execute(self):
        old_path = sys.path[:]
        try:
            with mock.patch.object(M, 'snapshot_info', return_value={'config_path': self.cfg}), \
                    mock.patch.object(M, 'validate_execution_walltime', return_value={'toy_validated_policy': True}), \
                    mock.patch.object(D, 'durable_json', side_effect=self.durable), \
                    mock.patch('subprocess.run', side_effect=self.fake_subprocess), \
                    mock.patch('time.monotonic', side_effect=lambda: self.now), mock.patch('builtins.print'):
                exec(self.code, {'__name__': '__main__'})
        finally:
            sys.path[:] = old_path

    def test_exact_pending_member_only_and_original_identity_preserved(self):
        before = copy.deepcopy(self.receipt)
        self.execute()
        updates = [argv for argv in self.calls if argv[0] == 'scontrol']
        self.assertEqual(len(updates), 1)
        receipt = json.loads(self.ledger.read_text())['submissions'][0]
        amendment = receipt.pop('execution_adjustments')
        self.assertEqual(receipt, before)
        self.assertEqual(amendment[0]['state'], 'verified')
        self.assertEqual(json.loads(self.evidence.read_text())['query_age_s'], 0)

    def test_running_or_absent_member_is_skipped_without_ledger_or_scheduler_change(self):
        for before in (self.before.replace('|PENDING|4:00:00|ebmcal_v2:', '|RUNNING|4:00:00|ebmcal_v2:'), ''):
            with self.subTest(before=before):
                self.before = before
                original = self.ledger.read_bytes()
                with self.assertRaises(SystemExit) as exc:
                    self.execute()
                self.assertEqual(exc.exception.code, 0)
                self.assertEqual(self.ledger.read_bytes(), original)
                self.assertEqual(json.loads(self.evidence.read_text())['status'], 'skipped_no_change')
                self.assertFalse(any(argv[0] == 'scontrol' for argv in self.calls))
                self.evidence.unlink(); self.queries = 0; self.calls.clear()

    def test_duplicate_evidence_still_blocks_with_python_optimization(self):
        self.evidence.write_text('{}')
        with self.assertRaisesRegex(ValueError, 'One-shot evidence exists'):
            self.execute()
        self.assertEqual(self.calls, [])

    def test_previous_adjustment_blocks_even_if_standalone_evidence_was_removed(self):
        self.receipt['execution_adjustments'] = [{'kind': 'pending_member_walltime_reduction', 'job_ids': [MEMBER]}]
        self.ledger.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'Prior member amendment'):
            self.execute()
        self.assertEqual(self.calls, [])

    def test_identity_or_project_tampering_rejected_before_queue(self):
        self.receipt['mapping'] = 'chunk_major'
        self.ledger.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'identity hash differs'):
            self.execute()
        self.assertEqual(self.calls, [])
        self.binding.write_text(json.dumps({'project': 'phase-test', 'ledger_path': str(self.ledger)}))
        with self.assertRaisesRegex(ValueError, 'binding differs'):
            self.execute()

    def test_malformed_empty_adjustment_history_rejected_before_scheduler_change(self):
        self.receipt['execution_adjustments'] = {}
        self.ledger.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError, 'Malformed adjustment history'):
            self.execute()
        self.assertEqual(self.calls, [])

    def test_pending_evidence_expiring_during_fsync_cannot_shorten_running_job(self):
        self.slow_intent = True
        original = self.ledger.read_bytes()
        with self.assertRaises(SystemExit):
            self.execute()
        value = json.loads(self.evidence.read_text())
        self.assertEqual(value['status'], 'skipped_stale_queue_no_change')
        self.assertIs(value['subprocess_invoked'], False)
        self.assertEqual(self.ledger.read_bytes(), original)
        self.assertFalse(any(argv[0] == 'scontrol' for argv in self.calls))

    def test_update_timeout_persists_unknown_and_appends_audit_without_retry(self):
        self.timeout_update = True
        with self.assertRaisesRegex(ValueError, 'no automatic retry'):
            self.execute()
        self.assertEqual(json.loads(self.evidence.read_text())['status'], 'update_unknown')
        receipt = json.loads(self.ledger.read_text())['submissions'][0]
        self.assertEqual(receipt['execution_adjustments'][0]['state'], 'update_unknown')
        self.assertEqual(receipt['walltime_s'], 14400)
        with self.assertRaisesRegex(ValueError, 'One-shot evidence exists'):
            self.execute()
        self.assertEqual(sum(argv[0] == 'scontrol' for argv in self.calls), 1)

    def test_nonzero_update_response_requires_reconciliation(self):
        self.update_returncode = 1
        with self.assertRaisesRegex(ValueError, 'no automatic retry'):
            self.execute()
        self.assertEqual(json.loads(self.evidence.read_text())['status'], 'needs_reconciliation')
        self.assertEqual(sum(argv[0] == 'scontrol' for argv in self.calls), 1)


if __name__ == '__main__':
    unittest.main()
