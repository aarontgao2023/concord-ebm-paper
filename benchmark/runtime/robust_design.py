"""Two common-noise stresses preserving diagnosis-conditional exchangeability."""
from pathlib import Path
import sys,hashlib,json
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import t
HERE=Path(__file__).resolve().parent
# design_v2 is the simulation generator in simulation/scripts of this repository.
sys.path.insert(0,str(HERE.parents[1]/'simulation/scripts'))
import design_v2 as base
T_DF=5
T_SCALE=np.sqrt((T_DF-2)/T_DF)

def canonical_sha(obj):
 return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

@lru_cache(None)
def t_auc(delta):
 # Integrating with the unscaled t variable is more accurate than tail truncation.
 value,error=quad(lambda x:t.pdf(x,T_DF)*t.cdf(x+delta/T_SCALE,T_DF),-np.inf,np.inf,epsabs=2e-12,epsrel=2e-12,limit=200)
 return float(value),float(error)

@lru_cache(None)
def calibrated_shift(auc):
 shift=brentq(lambda d:t_auc(d)[0]-auc,0.,20.,xtol=1e-12,rtol=1e-12)
 value,error=t_auc(shift)
 return {'target_component_auc':float(auc),'shift':float(shift),'calculated_component_auc':value,'quadrature_error_bound':error,'absolute_calibration_error':abs(value-auc),'noise_df':T_DF,'noise_scale':float(T_SCALE),'noise_variance':1.}

def calibration_table():
 return {b.name:calibrated_shift(b.component_auc) for b in base.BIOMARKERS}

def simulate_robust(mechanism,cell,seed):
 if mechanism not in ('CORR','T5'):raise ValueError(mechanism)
 specs={'H0':('REF_H0',{}),'SINGLE6_G1':('PWR_E2_K6',{'geometry':'single_displacement'}),'SINGLE13_G3':('PWR_E4_K13',{'geometry':'single_displacement'})}
 name,overrides=specs[cell]
 cfg=base.resolve(name,**overrides,fixed_base_order=tuple(range(14)),correlated=mechanism=='CORR')
 df,truth=base.simulate(cfg,seed)
 baseline_data_sha=truth['manifest']['data_sha256']
 if mechanism=='T5':
  rng=np.random.default_rng(np.random.SeedSequence(seed).spawn(len(base.STREAM_NAMES))[base.STREAM_NAMES.index('measurement_noise')])
  names=truth['biomarker_names'];complete=T_SCALE*rng.standard_t(T_DF,size=(len(df),len(names)))
  stage=np.array(truth['latent_stage']);groups=df.APOE.to_numpy();delta=np.array([calibrated_shift(base._BIOMARKER_BY_NAME[n].component_auc)['shift'] for n in names])
  for gi,g in enumerate(base.GROUP_ORDER):
   mask=groups==gi;ranks=np.argsort(truth['group_orderings'][g]);complete[mask]+=((ranks[None,:]<stage[mask,None])*delta)
   truth['component_parameters'][g]={n:calibrated_shift(base._BIOMARKER_BY_NAME[n].component_auc)|{'normal_mean':0.,'abnormal_mean':float(delta[j]),'sd':1.} for j,n in enumerate(names)}
  observed=complete.copy();observed[~np.isfinite(df[names].to_numpy())]=np.nan
  df.loc[:,names]=observed
  for gi,g in enumerate(base.GROUP_ORDER):
   cn=(groups==gi)&(df.Diagnosis.to_numpy()=='CN');ad=(groups==gi)&(df.Diagnosis.to_numpy()=='AD')
   truth['realized_diagnosis_auc'][g]={n:{'complete':base._auc(complete[cn,j],complete[ad,j]),'observed':base._auc(observed[cn,j],observed[ad,j])} for j,n in enumerate(names)}
 truth['robustness']={'version':'w5-common-noise-v1','mechanism':mechanism,'cell':cell,'one_fixed_latent_layout':True,'baseline_data_sha256':baseline_data_sha,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'component_calibration':calibration_table() if mechanism=='T5' else None,'conditional_exchangeability_rationale':'Under H0: same diagnosis-specific stage, common measurement and observation laws; no group-specific perturbations.'}
 truth['manifest']['data_sha256']=hashlib.sha256(df.to_csv(index=False,float_format='%.17g').encode()).hexdigest()
 truth['manifest']['robust_config_sha256']=canonical_sha({'base_config':truth['config'],'mechanism':mechanism,'cell':cell,'noise_df':T_DF if mechanism=='T5' else None,'calibration':truth['robustness']['component_calibration']})
 assert truth['conditional_exchangeability_by_design']==(cell=='H0')
 return df,truth
