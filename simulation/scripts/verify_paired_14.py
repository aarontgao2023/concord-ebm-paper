"""One development case: exact full-pipeline equivalence in the execution environment."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse,json
from dataclasses import replace
from pathlib import Path
from time import monotonic
import numpy as np
from design_v2 import resolve,simulate
from engine_v2 import EngineConfig
from paired_engine_v2 import fit_paired_orderings,_array_digest
from test_paired_engine_v2 import independent_reference,core_record
from run_v2 import permute
from dependency_v2 import verify_pyebm
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
df,truth=simulate(resolve('REF_H0'),31009402)
df=permute(df,truth,31009402,'diagnosis',0,{'stratify':'diagnosis'})
cfg=EngineConfig(expected_events=14,biomarker_names=tuple(truth['biomarker_names']),audit_original_neighbors=True)
started=monotonic()
refs={mode:independent_reference(df,replace(cfg,mode=mode)) for mode in ('original','repaired')}
independent_seconds=monotonic()-started
paired=fit_paired_orderings(df,cfg)
for mode in refs:
    ref,raw=refs[mode]; actual=paired[mode]
    assert ref.ok and actual.ok,(mode,ref.diagnostics,actual.diagnostics)
    assert core_record(ref)==core_record(actual),mode
    details=actual.diagnostics['model_artifacts']
    assert details['posterior_sha256']==_array_digest(raw['p_yes'])
    np.testing.assert_array_equal(details['event_centers'],raw['event_centers'])
    for a,b in zip(details['biomarker_parameters'],raw['parameters']):
        for key in ('Control','Disease','Mixing'): np.testing.assert_array_equal(a[key],b[key])
report={'status':'all_exact_equal','seed':31009402,'phase':'development','input':'REF_H0 diagnosis permutation id0','environment':verify_pyebm(),'independent_two_fit_seconds':independent_seconds,'paired_diagnostics':paired['original'].diagnostics['paired_engine'],'orderings':{m:paired[m].to_dict()['orderings'] for m in paired}}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps(report,allow_nan=False),flush=True)
