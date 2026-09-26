#!/usr/bin/env python3
"""AMENDMENT1: isolate first-group stopping, without changing inner SLSQP."""
from pathlib import Path
import argparse
import numpy as np
import audit_engine
OriginalState=audit_engine.AuditState
class AllGroupState(OriginalState):
 def iteration(self,*args,**kwargs):
  saved=self.mode;self.mode='original'
  super().iteration(*args,**kwargs)
  self.mode=saved;self.termination=None
  allpass=max(self.iterations[-1]['mixing_mean_abs_change_by_group'])<.01
  if allpass:self.termination='all_groups_at_upstream_tolerance';return True
  if len(self.iterations)>=100:self.termination='outer_cap';return True
  return False
if __name__=='__main__':
 audit_engine.AuditState=AllGroupState
 from run_numerical_audit import dataset_manifest,run
 ap=argparse.ArgumentParser();ap.add_argument('--index',type=int);ap.add_argument('--id');ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 tasks=dataset_manifest();task=tasks[a.index] if a.index is not None else next(t for t in tasks if t['id']==a.id)
 run(task,a.out,['allgroup01'])
