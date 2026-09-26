#!/usr/bin/env python3
"""ADNI six-event panel, step 1: anchor visits and coverage audit (no model fits, no tests).

Selects each participant's anchor visit (earliest dated ADSP-PHC CSF measurement with an
eligible APOE genotype and a CN, MCI or AD diagnosis within 180 days) and records aggregate
coverage of the candidate event panels (Supplementary Table 3, Supplementary Methods 2).

  python adni_audit.py --data-dir $CONCORD_ADNI_DIR --out $CONCORD_WORK_DIR/audit/adni

--data-dir holds the ADNI files listed in FILES (ADSP-PHC tables and the ADNI APOE table).
--out must lie outside the repository. Identifying rows and acquisition dates are written only
to <out>/local_only; all other outputs are aggregate counts, source definitions and hashes.
The next step is adni_export_core6.py with the same --data-dir and --out.
"""
import argparse
import os
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
DATA = None      # set by configure()
OUT = None
PRIVATE = None
GROUPS = ['e2', 'e33', 'e4']
DX = ['CN', 'MCI', 'AD']
CSF = {'ABETA': 'PHC_AB42', 'TAU': 'PHC_Tau', 'PTAU': 'PHC_pTau181'}
COG = {'MEM': 'PHC_MEM', 'EXF': 'PHC_EXF', 'LAN': 'PHC_LAN'}
MRI = {
    'Hippocampus': ['Left.Hippocampus_combat', 'Right.Hippocampus_combat'],
    'Ventricles': ['Left.Lateral.Ventricle_combat', 'Right.Lateral.Ventricle_combat'],
    'Entorhinal': ['lh_entorhinal_volume_combat', 'rh_entorhinal_volume_combat'],
    'Fusiform': ['lh_fusiform_volume_combat', 'rh_fusiform_volume_combat'],
    'MiddleTemporal': ['lh_middletemporal_volume_combat', 'rh_middletemporal_volume_combat'],
    'Precuneus': ['lh_precuneus_volume_combat', 'rh_precuneus_volume_combat'],
}
FILES = {
    'csf': 'ADSP_PHC_CSF_14Sep2026.csv',
    'cog': 'ADSP_PHC_COGN_14Sep2026.csv',
    'dx': 'ADSP_PHC_DEMODX_14Sep2026.csv',
    'apoe': 'APOERES_14Sep2026.csv',
    'mri': 'ADSP_PHC_T1_FS_18Oct2025.csv',
    'pet': 'ADSP_PHC_PET_Amyloid_Simple_22Aug2025.csv',
    'plasma': 'ADSP_PHC_PLASMA_14Sep2026.csv',
    'dictionary': 'ADSP_PHC_MERGED_DATADIC_20260521_14Sep2026.csv',
}


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def outside_repository(path):
    path = Path(path).resolve()
    if path == REPO or REPO in path.parents:
        raise SystemExit(f'The output folder must lie outside the repository: {path}')
    return path


def configure(data_dir, out_dir):
    global DATA, OUT, PRIVATE
    DATA = Path(data_dir).resolve()
    OUT = outside_repository(out_dir)
    PRIVATE = OUT / 'local_only'


def arguments(description, argv=None):
    ap = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--data-dir', type=Path, default=os.environ.get('CONCORD_ADNI_DIR'),
                    help='folder with the ADNI files (default: $CONCORD_ADNI_DIR)')
    ap.add_argument('--out', type=Path, default=None,
                    help='output folder outside the repository (default: $CONCORD_WORK_DIR/audit/adni)')
    args = ap.parse_args(argv)
    if args.out is None and os.environ.get('CONCORD_WORK_DIR'):
        args.out = Path(os.environ['CONCORD_WORK_DIR']) / 'audit/adni'
    if args.data_dir is None or args.out is None:
        ap.error('--data-dir and --out are required (or set CONCORD_ADNI_DIR and CONCORD_WORK_DIR)')
    configure(args.data_dir, args.out)
    return args


def load(key):
    d = pd.read_csv(DATA / FILES[key], low_memory=False)
    d['_source_row'] = np.arange(len(d)) + 2
    return d


def counts(d):
    t = pd.crosstab(d.APOE, d.Diagnosis).reindex(index=GROUPS, columns=DX, fill_value=0)
    return {g: {s: int(t.loc[g,s]) for s in DX} for g in GROUPS}


def genotype(s):
    return {'2/2':'e2', '2/3':'e2', '3/3':'e33', '3/4':'e4', '4/4':'e4'}.get(s, None)


def nearest(anchors, table, date, values, prefix, window=180):
    """One record per anchor; nearest date, earlier date, then source CSV row.

    At least one requested measurement must be available; completeness is not
    used to pick a later visit. All requested values come from that one record.
    """
    d = table[['RID', date, '_source_row'] + values].copy()
    d[date] = pd.to_datetime(d[date], errors='coerce')
    d = d.dropna(subset=[date]).dropna(subset=values, how='all')
    d = d.merge(anchors[['RID','anchor_id','anchor_date']], on='RID')
    d['_offset'] = (d[date] - d.anchor_date).dt.days
    d = d[d._offset.abs().le(window)].copy()
    d['_abs'] = d._offset.abs()
    d = d.sort_values(['anchor_id','_abs',date,'_source_row'], kind='stable').drop_duplicates('anchor_id')
    return d[['anchor_id',date,'_source_row','_offset']+values].rename(columns={
        date: prefix+'_date', '_source_row':prefix+'_source_row', '_offset':prefix+'_offset_days'})


def panel_stats(d, events, covariates, label):
    complete = np.isfinite(d[events]).all(axis=1)
    covok = d[covariates].notna().all(axis=1)
    use = d.loc[covok].copy()
    pairrows = []
    for a,b in itertools.combinations(events,2):
        n = counts(use.loc[np.isfinite(use[[a,b]]).all(axis=1)])
        for g in GROUPS:
            for s in DX:
                pairrows.append({'panel':label,'event1':a,'event2':b,'APOE':g,'Diagnosis':s,'n':n[g][s]})
    return {
        'events':events, 'covariates':covariates, 'anchor_n':len(d),
        'covariate_complete_n':int(covok.sum()), 'covariate_complete_cells':counts(use),
        'complete_events_n':int(complete.sum()),
        'complete_events_and_covariates_n':int((complete&covok).sum()),
        'complete_cells':counts(d.loc[complete&covok]),
        'minimum_event_pair_cell_support':min(x['n'] for x in pairrows),
        'zero_event_pair_cells':sum(x['n']==0 for x in pairrows),
        'event_availability':{e:int(np.isfinite(use[e]).sum()) for e in events},
    }, pairrows


def main(argv=None):
    arguments(__doc__, argv)
    OUT.mkdir(parents=True, exist_ok=True)
    PRIVATE.mkdir(parents=True, exist_ok=True)
    PRIVATE.chmod(0o700)
    csf, dx, cog, mri, apoe, pet, plasma = [load(k) for k in ['csf','dx','cog','mri','apoe','pet','plasma']]
    source_csf_rows, source_csf_people = len(csf), int(csf.RID.nunique())
    for columns,d in [(list(CSF.values()),csf),(list(COG.values()),cog)]:
        for col in columns:
            d[col] = pd.to_numeric(d[col],errors='coerce')
    # Stable genotype is validated rather than choosing among conflicting records.
    assert apoe.groupby('RID').GENOTYPE.nunique().max()==1
    gt=apoe.drop_duplicates('RID')[['RID','GENOTYPE']].copy()
    gt['APOE']=gt.GENOTYPE.map(genotype)
    csf['anchor_date']=pd.to_datetime(csf.DRAWDATE,errors='coerce')
    csf=csf.dropna(subset=['anchor_date']).dropna(subset=list(CSF.values()),how='all')
    assert not csf.duplicated(['RID','anchor_date','Platform']).any()
    csf['anchor_id']=csf._source_row.astype(int)
    csf=csf.merge(gt,on='RID',validate='many_to_one')
    csf=csf[csf.APOE.isin(GROUPS)].copy()
    # Match nearest known clinical diagnosis, including non-AD dementia before
    # excluding it. This avoids skipping a nearer non-AD visit for a desired label.
    dx=dx[dx.PHC_Diagnosis.notna()].copy()
    linked=nearest(csf,dx,'EXAMDATE',['PHC_Diagnosis','PHC_Sex','PHC_Education'],'diagnosis')
    csf=csf.merge(linked,on='anchor_id',how='left',validate='one_to_one')
    csf['Diagnosis']=csf.PHC_Diagnosis.map({1:'CN',2:'MCI',3:'AD'})
    eligible=csf[csf.Diagnosis.isin(DX)].copy()
    # Freeze anchor BEFORE cognition completeness, MRI, PET or plasma availability.
    frame=eligible.sort_values(['RID','anchor_date','_source_row'],kind='stable').drop_duplicates('RID').copy()
    frame=frame.rename(columns={v:k for k,v in CSF.items()})
    frame['Age']=pd.to_numeric(frame.PHC_Age_Biomarker,errors='coerce')
    frame['Sex']=frame.PHC_Sex.map({1:'Male',2:'Female'})
    frame['Education']=pd.to_numeric(frame.PHC_Education,errors='coerce')
    frame.loc[~frame.Education.between(0,40),'Education']=np.nan
    frame.loc[~frame.Age.between(18,110),'Age']=np.nan
    c=nearest(frame,cog,'EXAMDATE',list(COG.values()),'cognition').rename(columns={v:k for k,v in COG.items()})
    frame=frame.merge(c,on='anchor_id',how='left',validate='one_to_one')
    mfields=[v for vs in MRI.values() for v in vs]+['EstimatedTotalIntraCranialVol_combat']
    for col in mfields:
        mri[col]=pd.to_numeric(mri[col],errors='coerce')
    # Use one MRI record, without preferring its completeness.
    m=nearest(frame,mri,'PHC_SCANDATE',mfields,'mri')
    for event,fields in MRI.items():
        m[event]=m[fields].mean(axis=1,skipna=False)
        m.loc[~(m[fields]>0).all(axis=1),event]=np.nan
    m['ICV']=m.EstimatedTotalIntraCranialVol_combat.where(m.EstimatedTotalIntraCranialVol_combat.gt(0))
    frame=frame.merge(m[['anchor_id','mri_date','mri_source_row','mri_offset_days','ICV']+list(MRI)],on='anchor_id',how='left',validate='one_to_one')
    p=pet[(pet.PHC_QC_IMAGE==1)&(pet.PHC_QC_TIMING==1)].copy()
    p.PHC_CENTILOIDS=pd.to_numeric(p.PHC_CENTILOIDS,errors='coerce')
    p=nearest(frame,p,'PHC_SCANDATE',['PHC_CENTILOIDS'],'pet').rename(columns={'PHC_CENTILOIDS':'CENTILOID'})
    frame=frame.merge(p,on='anchor_id',how='left',validate='one_to_one')
    plasma.PHC_pTau217=pd.to_numeric(plasma.PHC_pTau217,errors='coerce')
    p=nearest(frame,plasma,'DRAWDATE',['PHC_pTau217'],'plasma').rename(columns={'PHC_pTau217':'PLASMA_PTAU217'})
    frame=frame.merge(p,on='anchor_id',how='left',validate='one_to_one')
    core=list(CSF)+list(COG)
    panels={'CSFcog6':core,'CSFcog5_noPTAU':[e for e in core if e!='PTAU'],
            'CSFcogMRI8':core+['Hippocampus','Ventricles'],
            'CSFcogMRI12':core+list(MRI),'CSFcogPET7':core+['CENTILOID'],
            'CSFcogPlasma7':core+['PLASMA_PTAU217'],
            'all14':core+list(MRI)+['CENTILOID','PLASMA_PTAU217']}
    summaries={};pairs=[]
    for name,events in panels.items():
        covars=['Age','Sex','Education']+(['ICV'] if any(e in MRI for e in events) else [])
        summary,pair=panel_stats(frame,events,covars,name);summaries[name]=summary;pairs+=pair
    coremask=np.isfinite(frame[core]).all(axis=1)&frame[['Age','Sex','Education']].notna().all(axis=1)
    final=frame.loc[coremask].copy()
    privatepath=PRIVATE/'candidate_frame.csv'
    frame.to_csv(privatepath,index=False,float_format='%.17g');privatepath.chmod(0o600)
    pd.DataFrame(pairs).to_csv(OUT/'event_pair_support.csv',index=False)
    dictionary=load('dictionary')
    definitions=dictionary[dictionary.VARNAME.isin(list(CSF.values())+list(COG.values())+['PHC_Diagnosis','PHC_Sex','PHC_Education','Platform','PHC_CENTILOIDS','PHC_pTau217'])].drop(columns=['_source_row','update_stamp'])
    definitions.to_csv(OUT/'source_definitions.csv',index=False)
    offsets={}
    for prefix in ['diagnosis','cognition','mri','pet','plasma']:
        x=frame[prefix+'_offset_days'].dropna().abs()
        y=final[prefix+'_offset_days'].dropna().abs()
        offsets[prefix]={'all_anchors_n':len(x),'all_anchors_max_abs_days':float(x.max()) if len(x) else None,
                         'complete_core_n':len(y),'complete_core_median_abs_days':float(y.median()) if len(y) else None,
                         'complete_core_max_abs_days':float(y.max()) if len(y) else None}
    report={
        'status':'COVERAGE_ONLY_NOT_ANALYSIS_INPUT', 'new_model_fits':0,
        'selection':'Earliest dated PHC CSF record with >=1 finite CSF z score, valid e2/e33/e4, and nearest known PHC clinical diagnosis within 180d classified CN/MCI/AD. Anchor selected before cognition/MRI/PET/plasma completeness. No converter or future diagnosis selection.',
        'matching':'Nearest single dated record within +/-180 days; ties earlier date then source CSV row. All events in a modality from one record. No replacement with more complete or later visit.',
        'diagnosis_codes':{'1':'CN','2':'prevalent MCI','3':'AD dementia'},
        'genotype_map':{'2/2':'e2','2/3':'e2','3/3':'e33','3/4':'e4','4/4':'e4','2/4':'excluded'},
        'source_CSF_rows':source_csf_rows,'source_CSF_people':source_csf_people,
        'genotype_eligible_dated_CSF_rows':len(csf),
        'with_eligible_nearby_diagnosis_CSF_rows':len(eligible),
        'anchor_people':len(frame),'anchor_cells':counts(frame),'panels':summaries,
        'platforms':frame.Platform.value_counts().to_dict(),
        'anchor_phases':frame.PHASE.value_counts().to_dict(),
        'complete_core_phases':final.PHASE.value_counts().to_dict(),
        'offsets':offsets,
        'CSF_units':'Unitless PHC z scores of log10 concentrations, raw concentrations NOT used; xMAP. Explicit pTau181 source field.',
        'cognition_units':'PHC harmonized cognitive composite scores; no PreciseFilter restriction or manual sign flip in this coverage audit.',
        'MRI_units':'PHC ComBat volume mm^3, bilateral arithmetic mean requiring both hemispheres; common ICV field. No covariate correction yet.',
        'PET_QC':'Only image and timing QC pass (both1); PHC_CENTILOIDS, no PHC_CL_FAIL substitute.',
        'PLASMA_platforms':plasma.Platform.value_counts().to_dict(),
        'source_sha256':{FILES[k]:sha(DATA/FILES[k]) for k in FILES},
        'script_sha256':sha(Path(__file__)),
        'private_candidate_sha256':sha(privatepath),
        'privacy':'Identifying RID/PTID/source records/acquisition dates only in local_only; no raw participant rows printed or cluster export produced.'}
    (OUT/'coverage_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'anchor_n':len(frame),'anchor_cells':counts(frame),'complete_core_n':len(final),'complete_core_cells':counts(final),'panel_complete_n':{k:v['complete_events_and_covariates_n'] for k,v in summaries.items()}},indent=2))


if __name__=='__main__':
    main()
