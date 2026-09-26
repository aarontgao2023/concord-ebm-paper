#!/usr/bin/env python3
"""ADNI six-event panel, step 2: verify the anchors and write the de-identified input.

Rebuilds every anchor of adni_audit.py with an independent per-participant loop, checks that
all values agree, keeps participants with the six events, age, sex and education, and writes
<out>/raw_deid.csv (fresh sequential pseudonyms; no RID, PTID or dates). The linkage to source
identifiers and dates goes to <out>/local_only only. The asserted counts are those of Table 2
and Supplementary Table 3 (1,231 anchors; 1,203 participants).

  python adni_export_core6.py --data-dir $CONCORD_ADNI_DIR --out $CONCORD_WORK_DIR/audit/adni

Run adni_audit.py with the same arguments first. The next step is ../prepare_campaign.py.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import adni_audit

EVENTS=['ABETA','TAU','PTAU','MEM','EXF','LAN']
GROUPS=['e2','e33','e4']; DX=['CN','MCI','AD']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def freeze(path, payload, private=False):
    if path.exists():
        assert path.read_bytes()==payload, 'Changed frozen output: '+path.name
    else:
        path.write_bytes(payload)
        path.chmod(0o600 if private else 0o444)


def main(argv=None):
    adni_audit.arguments(__doc__, argv)
    DATA, OUT = adni_audit.DATA, adni_audit.OUT
    audit=json.loads((OUT/'coverage_audit.json').read_text())
    frame=pd.read_csv(OUT/'local_only/candidate_frame.csv')
    c=pd.read_csv(DATA/adni_audit.FILES['csf'])
    d=pd.read_csv(DATA/adni_audit.FILES['dx'])
    g=pd.read_csv(DATA/adni_audit.FILES['apoe'])
    q=pd.read_csv(DATA/adni_audit.FILES['cog'])
    c['date']=pd.to_datetime(c.DRAWDATE,errors='coerce')
    d['date']=pd.to_datetime(d.EXAMDATE,errors='coerce')
    q['date']=pd.to_datetime(q.EXAMDATE,errors='coerce')
    for t in [c,d,q]:
        t['csv_row']=np.arange(len(t))+2
    gt=g.set_index('RID').GENOTYPE.to_dict()
    gmap={'2/2':'e2','2/3':'e2','3/3':'e33','3/4':'e4','4/4':'e4'}
    dmap={1:'CN',2:'MCI',3:'AD'}
    csfcols=['PHC_AB42','PHC_Tau','PHC_pTau181']
    cogcols=['PHC_MEM','PHC_EXF','PHC_LAN']
    db={rid:t for rid,t in d.dropna(subset=['date','PHC_Diagnosis']).groupby('RID')}
    qb={rid:t for rid,t in q.dropna(subset=['date']).dropna(subset=cogcols,how='all').groupby('RID')}
    expected={}
    for rid,visits in c.dropna(subset=['date']).dropna(subset=csfcols,how='all').groupby('RID'):
        if gt.get(rid) not in gmap or rid not in db:
            continue
        for _,v in visits.sort_values(['date','csv_row']).iterrows():
            dt=db[rid].copy()
            dt['distance']=(dt.date-v.date).dt.days.abs()
            dt=dt[dt.distance<=180].sort_values(['distance','date','csv_row'])
            if dt.empty or dt.iloc[0].PHC_Diagnosis not in dmap:
                continue
            dv=dt.iloc[0]
            qt=qb.get(rid,pd.DataFrame()).copy()
            if not qt.empty:
                qt['distance']=(qt.date-v.date).dt.days.abs()
                qt=qt[qt.distance<=180].sort_values(['distance','date','csv_row'])
            qv=qt.iloc[0] if not qt.empty else None
            values=v[csfcols].tolist()+(qv[cogcols].tolist() if qv is not None else [np.nan]*3)
            expected[int(rid)]={'anchor_id':int(v.csv_row),'Diagnosis':dmap[dv.PHC_Diagnosis],
                                'APOE':gmap[gt[rid]],'values':values,
                                'cognition_source_row':int(qv.csv_row) if qv is not None else None,
                                'diagnosis_source_row':int(dv.csv_row)}
            break
    assert len(expected)==len(frame)==1231
    assert set(expected)==set(frame.RID.astype(int))
    for _,r in frame.iterrows():
        e=expected[int(r.RID)]
        for name in ['anchor_id','Diagnosis','APOE','cognition_source_row','diagnosis_source_row']:
            assert e[name]==r[name], name
        assert np.allclose(e['values'],r[EVENTS].to_numpy(dtype=float),equal_nan=True,atol=1e-14)
    mask=np.isfinite(frame[EVENTS+['Age','Education']]).all(axis=1)&frame.Sex.isin(['Male','Female'])
    final=frame[mask].sort_values('RID').reset_index(drop=True)
    assert len(final)==1203
    counts=pd.crosstab(final.APOE,final.Diagnosis).reindex(index=GROUPS,columns=DX,fill_value=0)
    assert counts.values.tolist()==[[50,39,8],[222,273,65],[98,299,149]]
    raw=final[['APOE','Diagnosis','Age','Sex','Education']+EVENTS].copy()
    raw.insert(0,'PTID',[f'adni_csfcommon_{i:06d}' for i in range(1,len(raw)+1)])
    assert raw.PTID.is_unique and raw.notna().all().all()
    assert not any(x in raw.columns for x in ['RID','SUBJID','DRAWDATE','EXAMDATE','anchor_date'])
    csv=raw.to_csv(index=False,float_format='%.17g',lineterminator='\n').encode()
    freeze(OUT/'raw_deid.csv',csv)
    link=final[['RID','PTID','anchor_id','anchor_date','diagnosis_date','diagnosis_source_row','diagnosis_offset_days',
                'cognition_date','cognition_source_row','cognition_offset_days','Platform']].copy().rename(columns={'PTID':'source_PTID'})
    link.insert(0,'PTID',raw.PTID)
    freeze(OUT/'local_only/subject_linkage.csv',link.to_csv(index=False,float_format='%.17g',lineterminator='\n').encode(),True)
    freeze(OUT/'group_diagnosis_counts.csv',counts.rename_axis('APOE').reset_index().to_csv(index=False).encode())
    manifest={
        'cohort':'ADNI','n':len(raw),'input_file':'raw_deid.csv','input_sha256':sha(OUT/'raw_deid.csv'),
        'input_columns':raw.columns.tolist(),'events':EVENTS,
        'sensitivity_core5':{'events':[x for x in EVENTS if x!='PTAU'],'same_people_as_core6':True,'n':1203},
        'cell_counts':{g:{d:int(counts.loc[g,d]) for d in DX} for g in GROUPS},
        'diagnosis_order':DX,'group_order':GROUPS,
        'selection':audit['selection'],'matching':audit['matching'],
        'completeness':'Six raw PHC event scores plus age at CSF draw, sex and education required after anchor selection; no imputation, model fitting, reference choice or outcome-based exclusions.',
        'age_source':'PHC_Age_Biomarker at CSF anchor; finite, 18-110 years',
        'sex_source':'Nearest contemporaneous PHC DEMODX PHC_Sex; 1 Male, 2 Female',
        'education_source':'Same DEMODX record PHC_Education; finite 0-40 years',
        'event_source_mapping':dict(zip(EVENTS,csfcols+cogcols)),
        'CSF_platform':{'xMAP':1203},'CSF_definition':audit['CSF_units'],
        'cognition_definition':audit['cognition_units'],
        'preprocessing':'Unadjusted PHC scores; no covariate regression, no sign flip, no restandardization.',
        'source_sha256':audit['source_sha256'],
        'source_audit_sha256':sha(OUT/'coverage_audit.json'),
        'builder_sha256':audit['script_sha256'], 'exporter_sha256':sha(Path(__file__)),
        'privacy':{'upload_allowlist':['raw_deid.csv','preparation_manifest.json','group_diagnosis_counts.csv','validation.json'],
                   'local_only':'Source RID/PTID and acquisition dates in local_only; never upload.'},
        'limitations':['e2-AD cell has 8 participants and is retained without changing diagnosis or selecting later visits.',
                       'Index is earliest eligible CSF draw, not necessarily original study enrollment.',
                       'ADNI date verification cannot validate NACC CSF draw timing, absent from its return file.'],
    }
    validation={'independent_per_person_anchor_reconstruction':True,'nearest_diagnosis_and_cognition_rows_match':True,
                'all_raw_values_match_source':True,'complete_core6_and_covariates':True,'fresh_unique_ids_no_source_dates':True,
                'n':1203,'minimum_cell':8,'new_fits':0,'new_permutations':0}
    freeze(OUT/'preparation_manifest.json',(json.dumps(manifest,indent=2)+'\n').encode())
    freeze(OUT/'validation.json',(json.dumps(validation,indent=2)+'\n').encode())
    print(json.dumps({'n':len(raw),'cell_counts':manifest['cell_counts'],'input_sha256':manifest['input_sha256'],'independent_verification_passed':True},indent=2))


if __name__=='__main__':
    main()
