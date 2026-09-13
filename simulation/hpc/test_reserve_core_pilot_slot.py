"""Execute only the extracted helper body with a temporary root and fake Slurm."""
import copy
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import dispatch_campaign as D
import manage_campaign as M

PLAN = "0661eee6714148115c1ae9789e84ec63590e425416276a0a0302017ab83005a6"
PARENT = "10265062"
MEMBERS = [25] + list(range(430,445))


class ReserveCorePilotSlotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = (Path(self.tmp.name)/"ebmcal_v2").resolve()
        self.ledger = self.root/"control/submission_ledger.json"
        self.evidence = self.root/"control/campaign_versions/confirmation_20260908_core_paired_1h_pilot/reserve_core_pilot_slot.json"
        self.ledger.parent.mkdir(parents=True)
        self.snapshot = self.root/"snapshots/confirm_core_a"
        self.snapshot.mkdir(parents=True)
        self.cfg = self.snapshot/"config.json"
        project = Path(__file__).resolve().parents[2]
        self.cfg.write_bytes((project/"runs/v2/snapshots/confirm_core_a/configs/v2/confirm_core_a.json").read_bytes())
        (self.snapshot/"SHA256SUMS").write_bytes((project/"runs/v2/snapshots/confirm_core_a/SHA256SUMS").read_bytes())
        self.receipt = {"run_id":"confirm_core_a","snapshot_sha256":M.CORE_PAIRED_1H_PIN[0],
            "array_task_ids":MEMBERS,"nchunks":500,"mapping":"chunk_major","cpus_per_task":16,
            "max_concurrent":3,"walltime_s":7200,"dispatch_generation":5,
            "state":"submitted","plan_id":PLAN,"job_id":PARENT}
        self.assertEqual(M.I.digest(D.submission_identity(self.receipt)),PLAN)
        self.data = {"schema_version":"campaign_ledger_v1","project":"ebmcal_v2","submissions":[self.receipt]}
        self.ledger.write_text(json.dumps(self.data))
        self.binding = self.root/"campaign.dispatch.binding.json"
        self.binding.write_text(json.dumps({"project":"ebmcal_v2","ledger_path":str(self.ledger)}))
        self.helper = self.root/"reserve.sh"
        source = (Path(__file__).parent/"reserve_core_pilot_slot_20260908.sh").read_text()
        self.helper.write_text(source)
        body = source.split("<<'PY'\n",1)[1].rsplit("\nPY",1)[0]
        body = body.replace("Path('/N/slate/tg11/ebmcal_v2').resolve()",f"Path({str(self.root)!r}).resolve()",1)
        self.code = compile(body,"temporary_core_reservation_helper","exec",optimize=2)
        self.before = self.queue()
        self.recheck = self.queue()
        self.after = self.queue().replace(f"{PARENT}_25|PENDING",f"{PARENT}_25|RUNNING")
        self.control_before = self.control(3)
        self.control_after = self.control(2)
        self.calls, self.queries, self.controls = [],0,0
        self.now, self.slow_intent, self.timeout_update, self.update_returncode = 0.,False,False,0
        self.real_fsync = os.fsync

    def queue(self):
        return "".join(f"{PARENT}_{i}|PENDING|2:00:00|ebmcal_v2:{PLAN}\n" for i in MEMBERS)

    def control(self, throttle):
        return f"JobId={PARENT} ArrayJobId={PARENT} ArrayTaskId=25,430-444%{throttle} ArrayTaskThrottle={throttle} JobState=PENDING TimeLimit=02:00:00 NumCPUs=16 Comment=ebmcal_v2:{PLAN}\n"

    def fake_subprocess(self, argv, **kwargs):
        self.calls.append(argv)
        if argv[0] == "squeue":
            self.assertEqual(argv,["squeue","-r","-j",PARENT,"-h","-o","%i|%T|%l|%k"])
            values = [self.before,self.recheck,self.after]
            value = values[min(self.queries,2)]
            self.queries += 1
            return subprocess.CompletedProcess(argv,0,value,"")
        if argv[:3] == ["scontrol","show","job"]:
            self.assertEqual(argv,["scontrol","show","job",PARENT])
            value = self.control_before if self.controls == 0 else self.control_after
            self.controls += 1
            return subprocess.CompletedProcess(argv,0,value,"")
        self.assertEqual(argv,["scontrol","update","JobId=10265062","ArrayTaskThrottle=2"])
        evidence = json.loads(self.evidence.read_text())
        self.assertEqual(evidence["status"],"update_unknown")
        self.assertIsNone(evidence["subprocess_invoked"])
        if self.timeout_update:
            raise subprocess.TimeoutExpired(argv,30)
        return subprocess.CompletedProcess(argv,self.update_returncode,"","mock response")

    def fsync(self, fd):
        self.real_fsync(fd)
        if self.slow_intent and self.evidence.exists() and json.loads(self.evidence.read_text())["status"] == "update_unknown":
            self.now = 61.

    def execute(self):
        old_path = sys.path[:]
        try:
            with mock.patch.object(M,"snapshot_info",return_value={"config_path":self.cfg}), \
                    mock.patch.object(M,"validate_execution_walltime",return_value=None), \
                    mock.patch("sys.argv",["-",str(self.helper)]), \
                    mock.patch("subprocess.run",side_effect=self.fake_subprocess), \
                    mock.patch("time.monotonic",side_effect=lambda:self.now), \
                    mock.patch("os.fsync",side_effect=self.fsync), mock.patch("builtins.print"):
                exec(self.code,{"__name__":"__main__"})
        finally:
            sys.path[:] = old_path

    def updates(self):
        return [a for a in self.calls if a[:2] == ["scontrol","update"]]

    def test_exact_update_once_and_original_receipt_identity_remains_unchanged(self):
        original = copy.deepcopy(self.receipt)
        self.execute()
        self.assertEqual(len(self.updates()),1)
        saved = json.loads(self.ledger.read_text())["submissions"][0]
        adjustment = saved.pop("execution_adjustments")
        self.assertEqual(saved,original)
        self.assertEqual(adjustment[0]["state"],"verified")
        self.assertEqual((adjustment[0]["old_array_task_throttle"],adjustment[0]["new_array_task_throttle"]),(3,2))
        self.assertEqual(M.I.digest(D.submission_identity(saved)),PLAN)
        value = json.loads(self.evidence.read_text())
        self.assertEqual(value["parsed_queue_after"]["25"],"RUNNING")
        self.assertEqual(value["status"],"verified")
        self.assertIs(value["walltime_changes"],False)

    def test_one_shot_evidence_and_history_block_repeated_mutation(self):
        self.execute()
        with self.assertRaisesRegex(ValueError,"One-shot evidence exists"):
            self.execute()
        self.assertEqual(len(self.updates()),1)
        self.evidence.unlink()
        with self.assertRaisesRegex(ValueError,"Prior throttle amendment"):
            self.execute()
        self.assertEqual(len(self.updates()),1)

    def test_changed_identity_or_bound_ledger_rejected_before_queries(self):
        self.receipt["max_concurrent"] = 2
        self.ledger.write_text(json.dumps(self.data))
        with self.assertRaisesRegex(ValueError,"identity hash differs"):
            self.execute()
        self.assertEqual(self.calls,[])
        self.binding.write_text(json.dumps({"project":"phase-test","ledger_path":str(self.ledger)}))
        with self.assertRaisesRegex(ValueError,"binding differs"):
            self.execute()

    def test_member_starting_on_second_queue_or_controller_is_safe_skip(self):
        for where in ("queue","controller"):
            with self.subTest(where=where):
                if where == "queue":
                    self.recheck = self.before.replace("|PENDING|","|RUNNING|",1)
                else:
                    self.recheck = self.before
                    self.control_before = self.control(3).replace("JobState=PENDING","JobState=RUNNING")
                original = self.ledger.read_bytes()
                with self.assertRaises(SystemExit) as exc:
                    self.execute()
                self.assertEqual(exc.exception.code,0)
                self.assertEqual(json.loads(self.evidence.read_text())["status"],"skipped_no_change")
                self.assertEqual(self.ledger.read_bytes(),original)
                self.assertEqual(self.updates(),[])
                self.evidence.unlink(); self.calls.clear(); self.queries = self.controls = 0

    def test_missing_or_compressed_or_wrong_comment_members_skip(self):
        variants = [self.before.split("\n",1)[1],self.before+self.before.splitlines()[0]+"\n",
                    self.before.replace(PARENT+"_25",PARENT+"_[25]",1),self.before.replace(PLAN,"a"*64,1)]
        for value in variants:
            with self.subTest(value=value[:80]):
                self.before = value
                with self.assertRaises(SystemExit):
                    self.execute()
                self.assertEqual(self.updates(),[])
                self.assertEqual(json.loads(self.evidence.read_text())["status"],"skipped_no_change")
                self.evidence.unlink(); self.calls.clear(); self.queries = self.controls = 0

    def test_stale_after_intent_fsync_skips_and_does_not_modify_ledger(self):
        original = self.ledger.read_bytes()
        self.slow_intent = True
        with self.assertRaises(SystemExit):
            self.execute()
        value = json.loads(self.evidence.read_text())
        self.assertEqual(value["status"],"skipped_stale_or_changed_no_change")
        self.assertIs(value["subprocess_invoked"],False)
        self.assertEqual(self.ledger.read_bytes(),original)
        self.assertEqual(self.updates(),[])

    def test_update_timeout_is_unknown_preserved_and_never_retried(self):
        self.timeout_update = True
        with self.assertRaisesRegex(ValueError,"no automatic retry"):
            self.execute()
        self.assertEqual(json.loads(self.evidence.read_text())["status"],"update_unknown")
        receipt = json.loads(self.ledger.read_text())["submissions"][0]
        self.assertEqual(receipt["execution_adjustments"][0]["state"],"update_unknown")
        self.assertEqual(receipt["max_concurrent"],3)
        with self.assertRaisesRegex(ValueError,"One-shot evidence"):
            self.execute()
        self.assertEqual(len(self.updates()),1)

    def test_nonzero_update_needs_reconciliation(self):
        self.update_returncode = 1
        with self.assertRaisesRegex(ValueError,"no automatic retry"):
            self.execute()
        self.assertEqual(json.loads(self.evidence.read_text())["status"],"needs_reconciliation")
        self.assertEqual(len(self.updates()),1)

    def test_wrong_final_throttle_preserves_uncertainty_without_retry(self):
        self.control_after = self.control(3)
        with self.assertRaisesRegex(ValueError,"no automatic retry"):
            self.execute()
        value = json.loads(self.evidence.read_text())
        self.assertEqual(value["status"],"needs_reconciliation")
        self.assertEqual(value["returncode"],0)
        self.assertEqual(len(self.updates()),1)


if __name__ == "__main__":
    unittest.main()
