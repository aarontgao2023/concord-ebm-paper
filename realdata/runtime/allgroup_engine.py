"""Original-tolerance all-group stopping; caps fail rather than certify convergence."""
from contextlib import contextmanager
from pathlib import Path
import sys,inspect,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'vendor/w3'))
import audit_engine
from run_stop_isolation import AllGroupState

class LightState:
 def __init__(self):
  self.iterations=[];self.termination=None;self.last_parameters=None;self.function_source_sha256=None
 def iteration(self,params,mixes,old,Groups,GroupValues,Data_all):
  dx=np.asarray(mixes,float)-np.asarray(old,float)
  if dx.ndim==1:dx=dx[None,:]
  change=np.mean(np.abs(dx),axis=1)
  p=np.asarray(params,float)
  if p.ndim==3:p=p[None,...]
  self.last_parameters=p[...,0].copy()
  self.iterations.append(change.tolist())
  if not np.isfinite(change).all() or not np.isfinite(self.last_parameters).all():raise FloatingPointError('Nonfinite all-group stopping diagnostics')
  # Same arithmetic/strict inequality as W3 Amendment1; no inner optimizer change.
  if max(change.tolist())<.01:self.termination='all_groups_at_upstream_tolerance';return True
  if len(self.iterations)>=100:self.termination='outer_cap';return True
  return False
 def summary(self,include_trajectory=False):
  result={'mode':'allgroup01-light','termination':self.termination,'outer_iterations':len(self.iterations),'function_source_sha256':self.function_source_sha256,'final_mixing_mean_abs_change_by_group':self.iterations[-1] if self.iterations else None,'max_final_mixing_residual':max(self.iterations[-1]) if self.iterations else None,'final_parameter_sha256':hashlib.sha256(np.asarray(self.last_parameters,dtype='<f8').tobytes()).hexdigest() if self.last_parameters is not None else None}
  if include_trajectory:result['mixing_trajectory']=self.iterations
  return result

@contextmanager
def allgroup_stopping(detailed=True):
 if detailed:
  original_state=audit_engine.AuditState;audit_engine.AuditState=AllGroupState
  try:
   with audit_engine.instrument_mixture('allgroup01') as state:
    yield state
    if state.termination=='outer_cap':raise RuntimeError('All-group0.01 rule reached100-iteration safeguard')
  finally:audit_engine.AuditState=original_state
 else:
  import pyebm.core_utilities as cu
  original=cu.do_mixturemodel;source=inspect.getsource(original);state=LightState();state.function_source_sha256=hashlib.sha256(source.encode()).hexdigest()
  for old in ('if np.mean(np.abs(mixes-mixes_old))<10**-2:','if np.mean(np.abs(mixes[0]-mixes_old[0]))<10**-2:'):
   assert source.count(old)==1,old
   source=source.replace(old,'if _allgroup_iteration(params_opt, mixes, mixes_old, Groups, GroupValues, Data_all):')
  ns=dict(original.__globals__);ns['_allgroup_iteration']=state.iteration
  exec(compile(source,'<allgroup01-stop-only>','exec'),ns);cu.do_mixturemodel=ns['do_mixturemodel']
  try:
   yield state
   if state.termination=='outer_cap':raise RuntimeError('All-group0.01 rule reached100-iteration safeguard')
  finally:cu.do_mixturemodel=original
