"""Temporary, auditable mixture instrumentation; never edits installed pyebm.

Original mode preserves upstream calls/arithmetic and stops. Strict mode changes
only SLSQP tolerance and all-group outer convergence, as frozen in PROTOCOL.json.
"""
from __future__ import annotations
from contextlib import contextmanager
import inspect,hashlib
import numpy as np


class AuditState:
    def __init__(self,mode):
        self.mode=mode;self.iterations=[];self.optimizers=[];self.previous=None;self.streak=0
        self.termination=None;self.function_source_sha256=None
    def iteration(self,params,mixes,old,Groups,GroupValues,Data_all):
        import pyebm.mixture_model.gaussian_mixture_model as gmm
        p=np.asarray(params,dtype=float)
        if p.ndim==3:p=p[None,...]
        p=p[...,0].copy()
        groups=np.unique(GroupValues[0]) if len(Groups) else [0]
        dx=np.asarray(mixes,dtype=float)-np.asarray(old,dtype=float)
        if dx.ndim==1:dx=dx[None,:]
        mixing_change=np.mean(np.abs(dx),axis=1)
        nll=[]
        for idx,g in enumerate(groups):
            data=Data_all[np.asarray(GroupValues[0])==g] if len(Groups) else Data_all
            nll.append(float(sum(gmm.calculate_likelihood_gmm(p[idx,k],data[:,k,:],[],[],[]) for k in range(p.shape[1]))))
        if self.previous is None:parchange=objchange=None
        else:
            pp,pnll=self.previous
            parchange=np.max(np.abs(p-pp)/np.maximum(1.,np.abs(pp)),axis=(1,2)).tolist()
            objchange=(np.abs(np.array(nll)-pnll)/np.maximum(1.,np.abs(pnll))).tolist()
        row={'iteration':len(self.iterations)+1,'mixing_mean_abs_change_by_group':mixing_change.tolist(),
             'parameters':p.tolist(),'nll_by_group':nll,'max_relative_parameter_change_by_group':parchange,
             'relative_nll_change_by_group':objchange}
        self.iterations.append(row);self.previous=(p,np.array(nll))
        if self.mode=='original':
            stop=bool(mixing_change[0]<.01)
            if stop:self.termination='upstream_mixing_stop'
            return stop
        stable=(parchange is not None and max(mixing_change)<1e-4 and max(parchange)<1e-4 and max(objchange)<1e-6)
        self.streak=self.streak+1 if stable else 0
        if self.streak>=2:self.termination='all_groups_strict_converged';return True
        if len(self.iterations)>=100:self.termination='outer_cap';return True
        return False
    def summary(self):
        return {'mode':self.mode,'termination':self.termination,'outer_iterations':len(self.iterations),
                'function_source_sha256':self.function_source_sha256,'trajectory':self.iterations,
                'optimizer_calls':self.optimizers,'strict_converged':self.termination=='all_groups_strict_converged'}


class OptimizerProxy:
    def __init__(self,original,state):self.original=original;self.state=state
    def __getattr__(self,name):return getattr(self.original,name)
    def minimize(self,*args,**kwargs):
        if self.state.mode=='strict':kwargs['options']={**kwargs.get('options',{}),'ftol':1e-9,'maxiter':600}
        r=self.original.minimize(*args,**kwargs)
        x=np.asarray(r.x,float);bounds=np.asarray(kwargs.get('bounds',[]),float)
        near=np.maximum(1.,np.abs(x))*1e-7
        active_low=(np.abs(x-bounds[:,0])<=near).tolist() if len(bounds)==len(x) else []
        active_high=(np.abs(x-bounds[:,1])<=near).tolist() if len(bounds)==len(x) else []
        self.state.optimizers.append({'success':bool(r.success),'status':int(r.status),'nit':int(getattr(r,'nit',-1)),
            'message':str(r.message),'fun':float(r.fun),'x':x.tolist(),'bounds':bounds.tolist(),
            'active_lower':active_low,'active_upper':active_high,'ftol':kwargs.get('options',{}).get('ftol','scipy-default')})
        return r


@contextmanager
def instrument_mixture(mode):
    import pyebm.core_utilities as cu
    import pyebm.mixture_model.gaussian_mixture_model as gmm
    state=AuditState(mode);original=cu.do_mixturemodel;opt=gmm.opt
    source=inspect.getsource(original);state.function_source_sha256=hashlib.sha256(source.encode()).hexdigest()
    replacements=[('if np.mean(np.abs(mixes-mixes_old))<10**-2:',
                   'if _audit_iteration(params_opt, mixes, mixes_old, Groups, GroupValues, Data_all):'),
                  ('if np.mean(np.abs(mixes[0]-mixes_old[0]))<10**-2:',
                   'if _audit_iteration(params_opt, mixes, mixes_old, Groups, GroupValues, Data_all):')]
    for old,new in replacements:
        assert source.count(old)==1,old
        source=source.replace(old,new)
    ns=dict(original.__globals__);ns['_audit_iteration']=state.iteration
    exec(compile(source,'<W3-instrumented-pyebm-mixture>','exec'),ns)
    cu.do_mixturemodel=ns['do_mixturemodel'];gmm.opt=OptimizerProxy(opt,state)
    try:yield state
    finally:cu.do_mixturemodel=original;gmm.opt=opt
