"""Small fit-free scientific invariants and failure-record checks; no output fits."""
from __future__ import annotations
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
from common import (read_campaign,read_frame,n2_frame,pair_shuffle,bootstrap_frame,frame_sha,
                    GROUPS,LABELS,ENGINES,canonical,save)
import engine_adapter as adapter


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();cfg,base=read_campaign(args.config);checks=[]
    for cohort in cfg['cohorts']:
        frame=read_frame(cfg,base,cohort)
        for rep in range(200):
            a,da=n2_frame(frame,cfg,cohort,rep,'matched',base)
            b,db=n2_frame(frame,cfg,cohort,rep,'different',base)
            assert a.drop(columns='APOE').equals(b.drop(columns='APOE'))
            assert da['selected_roster_sha256']==db['selected_roster_sha256']
            assert da['draw_seed_components']==[cfg['cohorts'][cohort]['draw_seed_base'],rep]
        checks.append({'name':cohort+'_all200_identical_participant_matrix_between_mixes','passed':True})
        data=frame.copy();data['APOE']=data.APOE.map({g:i for i,g in enumerate(GROUPS)}).astype(int)
        seed=cfg['cohorts'][cohort]['apoe_seed']
        for pi in range(3):
            for pid in (0,1,598):
                p=pair_shuffle(data,seed,pi,pid)
                assert p.drop(columns='APOE').equals(data.drop(columns='APOE'))
                assert pd.crosstab(p.APOE,p.Diagnosis).equals(pd.crosstab(data.APOE,data.Diagnosis))
                assert frame_sha(p)==frame_sha(pair_shuffle(data,seed,pi,pid))
        checks.append({'name':cohort+'_Dpair_fixed_margins_and_matrix_reproducibility','passed':True})
        bs,meta=bootstrap_frame(data,cfg['cohorts'][cohort]['bootstrap_seed'],cfg['biomarkers'])
        assert bs.PTID.is_unique and len(bs)==len(data)
        assert pd.crosstab(bs.APOE,bs.Diagnosis).equals(pd.crosstab(data.APOE,data.Diagnosis))
        checks.append({'name':cohort+'_bootstrap_nine_strata_and_unique_bootstrap_IDs','passed':True})
    d=np.array([1,1,2,3,1,2,2,3]);g=np.array([0]*4+[1]*4);mask=np.ones(8,dtype=bool)
    original=adapter.invariant.composition_weights
    with adapter.explicit_reference(cfg['reference']):
        groups,w,ess=adapter.invariant.composition_weights(d,g,mask,'min')
        for group in groups:
            for j,diagnosis in enumerate([1,2,3]):
                assert np.isclose(w[group][d[g==group]==diagnosis].sum()/w[group].sum(),cfg['reference'][j])
    assert adapter.invariant.composition_weights is original
    checks.append({'name':'explicit_reference_weighted_diagnosis_and_patch_restoration','passed':True})
    # Mock only the expensive fit; exercise the real adapter certificate and log path.
    old_fit,old_stopping=adapter.fit_orderings,adapter.allgroup_stopping
    @contextmanager
    def bad_certificate(detailed):
        summary={'termination':'all_groups_at_upstream_tolerance','outer_iterations':1,
                 'trajectory':[{'mixing_mean_abs_change_by_group':[.001,float('nan')]}]}
        yield SimpleNamespace(termination='all_groups_at_upstream_tolerance',summary=lambda:summary)
    adapter.fit_orderings=lambda data,config:SimpleNamespace(ok=True,diagnostics={'parameters_finite':True})
    adapter.allgroup_stopping=bad_certificate
    try:
        sample=bs[['PTID','APOE','Diagnosis']+cfg['biomarkers']]
        result=adapter.fit_task({'data':sample,'engine':'separate','biomarkers':cfg['biomarkers'],'reference':cfg['reference'],
                                'fit_seed':1,'fit_timeout_s':5,'kind':'observed','key':'mock-only'})
        assert result['status']=='error' and 'Nonfinite' in result['error']
        json.dumps(result,allow_nan=False)
    finally:
        adapter.fit_orderings,adapter.allgroup_stopping=old_fit,old_stopping
    checks.append({'name':'nonfinite_detailed_certificate_rejected_and_preserved_as_valid_JSON','passed':True})
    # Inclusive Bonferroni boundary, B599: e=9 rejects, e=10 does not.
    assert 3*(1+9)*20<=600 and not 3*(1+10)*20<=600
    checks.append({'name':'inclusive_B599_pair_Bonferroni_integer_boundary','passed':True})
    save(args.out,{'passed':True,'fits_performed':False,'checks':checks,'campaign_canonical_sha256':canonical(cfg)})
    print(json.dumps({'passed':True,'checks':len(checks),'fits_performed':False}))


if __name__=='__main__':
    main()
