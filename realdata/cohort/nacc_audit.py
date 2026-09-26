#!/usr/bin/env python3
"""NACC six-event panel, step 1: anchor visits and coverage audit (no model fits, no tests).

Links each CSF measurement of the NACC ADSP-PHC release to its Uniform Data Set (UDS) visit
(NACCID + NACCVNUM), keeps the earliest eligible linked visit per participant (APOE e2, e3/e3
or e4; CN, prevalent MCI or AD dementia with AD aetiology) and records aggregate coverage of
the candidate event panels (Supplementary Table 3, Supplementary Methods 2). The anchor is
fixed before completeness is checked; later visits never replace an incomplete profile.

  python nacc_audit.py --phc-dir $CONCORD_NACC_DIR --uds-file $CONCORD_NACC_UDS \\
      --out $CONCORD_WORK_DIR/audit/nacc

--phc-dir is the NACC ADSP-PHC release folder (subfolders Biomarker/, Cognition/, Imaging_T1/,
Imaging_PET/). --uds-file is a CSV export of the NACC UDS investigator file with the columns
read below. --out must lie outside the repository; participant frames go to <out>/local_only
and every other output is aggregate. The next step is nacc_prepare_core6.py.
"""
import argparse
from pathlib import Path
import hashlib
import itertools
import json
import os
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = None       # set by configure()
NEW = None
KEY = ['NACCID', 'NACCVNUM']
GROUPS = ['e2', 'e33', 'e4']
DX = ['CN', 'MCI', 'AD']
CSF = {'ABETA': 'PHC_AB42', 'TAU': 'PHC_Tau', 'PTAU': 'PHC_pTau'}
COG = {'MEM': 'PHC_MEM', 'EXF': 'PHC_EXF', 'LAN': 'PHC_LAN'}
MRI = {
 'Hippocampus':['Left.Hippocampus_combat','Right.Hippocampus_combat'],
 'Ventricles':['Left.Lateral.Ventricle_combat','Right.Lateral.Ventricle_combat'],
 'Entorhinal':['lh_entorhinal_volume_combat','rh_entorhinal_volume_combat'],
 'Fusiform':['lh_fusiform_volume_combat','rh_fusiform_volume_combat'],
 'MiddleTemporal':['lh_middletemporal_volume_combat','rh_middletemporal_volume_combat'],
 'Precuneus':['lh_precuneus_volume_combat','rh_precuneus_volume_combat'],
}
ICV = 'EstimatedTotalIntraCranialVol_combat'
EVENTS = list(CSF) + list(COG) + list(MRI) + ['CENTILOID']
PANELS = {
 'CSFcog6':list(CSF)+list(COG),
 'CSFcog5_noPTAU':['ABETA','TAU']+list(COG),
 'CSFcog6_MRI2':list(CSF)+list(COG)+['Hippocampus','Ventricles'],
 'CSFcog6_MRI6':list(CSF)+list(COG)+list(MRI),
 'CSFcog6_PET':list(CSF)+list(COG)+['CENTILOID'],
 'MRIcog9':list(COG)+list(MRI),
 'PETcog4':['CENTILOID']+list(COG),
}
RELEASE_FILES = {
 'csf':'Biomarker/NACC_ADSP_PHC_Biomarker_2024.csv',
 'cog':'Cognition/NACC_ADSP_PHC_Cognition_2024.csv',
 'fs':'Imaging_T1/NACC_ADSP_PHC_T1_Freesurfer_2024.csv',
 'pet':'Imaging_PET/NACC_ADSP_PHC_Amyloid_Simple_2024.csv',
}
FILES = {}       # key -> path, set by configure()


def configure(phc_dir, uds_file, out_dir):
 global OUT, NEW
 NEW = Path(phc_dir).resolve()
 OUT = Path(out_dir).resolve()
 if OUT == REPO or REPO in OUT.parents:
  raise SystemExit(f'The output folder must lie outside the repository: {OUT}')
 FILES.clear()
 FILES.update({k:NEW/v for k,v in RELEASE_FILES.items()})
 FILES['uds'] = Path(uds_file).resolve()


def provenance_name(key):
 return RELEASE_FILES.get(key, FILES[key].name)


def arguments(description, argv=None):
 ap = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
 ap.add_argument('--phc-dir', type=Path, default=os.environ.get('CONCORD_NACC_DIR'),
                 help='NACC ADSP-PHC release folder (default: $CONCORD_NACC_DIR)')
 ap.add_argument('--uds-file', type=Path, default=os.environ.get('CONCORD_NACC_UDS'),
                 help='CSV export of the NACC UDS file (default: $CONCORD_NACC_UDS)')
 ap.add_argument('--out', type=Path, default=None,
                 help='output folder outside the repository (default: $CONCORD_WORK_DIR/audit/nacc)')
 args = ap.parse_args(argv)
 if args.out is None and os.environ.get('CONCORD_WORK_DIR'):
  args.out = Path(os.environ['CONCORD_WORK_DIR']) / 'audit/nacc'
 if args.phc_dir is None or args.uds_file is None or args.out is None:
  ap.error('--phc-dir, --uds-file and --out are required (or set CONCORD_NACC_DIR, CONCORD_NACC_UDS and CONCORD_WORK_DIR)')
 configure(args.phc_dir, args.uds_file, args.out)
 return args

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def finite(d,cols):return np.isfinite(d[cols].to_numpy(float)).all(axis=1)
def counts(d):
 c=pd.crosstab(d.APOE,d.Diagnosis).reindex(index=GROUPS,columns=DX,fill_value=0)
 return {g:{s:int(c.loc[g,s]) for s in DX} for g in GROUPS}
def read(k,cols):
 d=pd.read_csv(FILES[k],usecols=cols,low_memory=False)
 d.NACCID=d.NACCID.astype(str)
 if 'NACCVNUM' in d:
  v=pd.to_numeric(d.NACCVNUM,errors='raise')
  assert (v%1==0).all() and v.notna().all()
  d.NACCVNUM=v.astype(int)
 return d

def main(argv=None):
 arguments(__doc__, argv)
 os.umask(0o077); OUT.mkdir(parents=True,exist_ok=True)
 private=OUT/'local_only';private.mkdir(exist_ok=True);os.chmod(private,0o700)
 csf=read('csf',KEY+list(CSF.values())+['Platform'])
 cog=read('cog',KEY+list(COG.values()))
 fs=read('fs',KEY+[c for p in MRI.values() for c in p]+[ICV])
 uds=read('uds',KEY+['VISITYR','VISITMO','VISITDAY','NACCAPOE','NACCUDSD','NACCALZD','NACCAGE','SEX','EDUC'])
 pet=read('pet',['NACCID','PHC_SCANDATE','PHC_CENTILOIDS','PHC_QC_IMAGE','PHC_QC_TIMING','PHC_CL_FAIL'])
 assert csf.NACCID.is_unique and not cog.duplicated(KEY).any() and not uds.duplicated(KEY).any()
 badfs=fs.duplicated(KEY,keep=False);fs=fs.loc[~badfs].copy()
 for e,c in CSF.items():csf[e]=pd.to_numeric(csf[c],errors='coerce')
 for e,c in COG.items():cog[e]=pd.to_numeric(cog[c],errors='coerce')
 for e,p in MRI.items():fs[e]=(pd.to_numeric(fs[p[0]],errors='coerce')+pd.to_numeric(fs[p[1]],errors='coerce'))/2
 fs['ICV']=pd.to_numeric(fs[ICV],errors='coerce')
 uds['UDS_date']=pd.to_datetime(dict(year=uds.VISITYR,month=uds.VISITMO,day=uds.VISITDAY),errors='coerce')
 uds['APOE']=uds.NACCAPOE.map({1:'e33',2:'e4',3:'e2',4:'e4',5:'e2e4',6:'e2'})
 uds['Diagnosis']=np.select([uds.NACCUDSD.eq(1),uds.NACCUDSD.eq(3),uds.NACCUDSD.eq(4)&uds.NACCALZD.eq(1)],DX,default='excluded')
 uds['Age']=pd.to_numeric(uds.NACCAGE,errors='coerce').where(lambda x:x.between(18,120))
 uds['Sex']=uds.SEX.map({1:'Male',2:'Female'})
 uds['Education']=pd.to_numeric(uds.EDUC,errors='coerce').where(lambda x:x.between(0,98))
 uds=uds.loc[uds.APOE.isin(GROUPS)&uds.Diagnosis.isin(DX)].copy()
 frames={}
 for anchor,source in [('CSF',csf),('MRI',fs)]:
  a=source[KEY].merge(uds,on=KEY,how='inner',validate='one_to_one')
  assert a.UDS_date.notna().all()
  a=a.sort_values(['NACCID','UDS_date','NACCVNUM'],kind='mergesort').drop_duplicates('NACCID')
  a=a.merge(csf[KEY+list(CSF)+['Platform']],on=KEY,how='left',validate='one_to_one')
  a=a.merge(cog[KEY+list(COG)],on=KEY,how='left',validate='one_to_one')
  a=a.merge(fs[KEY+list(MRI)+['ICV']],on=KEY,how='left',validate='one_to_one').reset_index(drop=True)
  # Objective coverage bound: nearest finite PET on actual scan date within 180 days.
  # No availability-dependent change to index UDS visit; QC statuses are reported.
  pp=pet.copy();pp['PET_date']=pd.to_datetime(pp.PHC_SCANDATE,errors='coerce')
  pp['CENTILOID']=pd.to_numeric(pp.PHC_CENTILOIDS,errors='coerce');pp['source_row']=range(len(pp))
  pp=pp.loc[np.isfinite(pp.CENTILOID)&pp.PET_date.notna()]
  link=a[['NACCID','UDS_date']].merge(pp,on='NACCID',how='left')
  link['PET_gap_days']=(link.PET_date-link.UDS_date).dt.total_seconds().abs()/86400
  link=link.loc[link.PET_gap_days.le(180)].sort_values(['NACCID','PET_gap_days','PET_date','source_row'],kind='mergesort').drop_duplicates('NACCID')
  a=a.merge(link[['NACCID','CENTILOID','PET_date','PET_gap_days','PHC_QC_IMAGE','PHC_QC_TIMING','PHC_CL_FAIL']],on='NACCID',how='left',validate='one_to_one')
  a['core_covariates']=finite(a,['Age','Education'])&a.Sex.notna()
  a['mri_covariates']=a.core_covariates&np.isfinite(a.ICV)&a.ICV.gt(0)
  frames[anchor]=a
 # Save local input frame solely for deterministic later preparation. Never transfer.
 for name,a in frames.items():a.to_csv(private/f'{name.lower()}_anchor_LOCAL_ONLY.csv',index=False,float_format='%.17g',date_format='%Y-%m-%d')
 rows=[];pairrows=[];evrows=[];patrows=[];head={}
 for anchor,a in frames.items():
  head[anchor]={'anchor_n':len(a),'anchor_counts':counts(a),'core_covariate_n':int(a.core_covariates.sum()),'all_event_availability':{e:int(np.isfinite(a[e]).sum()) for e in EVENTS}}
  for panel,events in PANELS.items():
   cov=a.mri_covariates if any(e in MRI for e in events) else a.core_covariates
   masks={'anchor':np.ones(len(a),bool),'covariates':cov,'complete':cov&finite(a,events),'masked_at_least2':cov&(np.isfinite(a[events]).sum(axis=1)>=2)}
   for selection,mask in masks.items():
    sub=a.loc[mask];ct=counts(sub)
    rows.extend({'anchor':anchor,'panel':panel,'selection':selection,'APOE':g,'Diagnosis':s,'n':ct[g][s]} for g,s in itertools.product(GROUPS,DX))
   complete=a.loc[masks['complete']]
   head[anchor][panel]={'events':events,'complete_n':len(complete),'complete_counts':counts(complete),'masked_at_least2_n':int(masks['masked_at_least2'].sum())}
   masked=a.loc[masks['masked_at_least2']]
   supports=[]
   for e1,e2 in itertools.combinations(events,2):
    for g,s in itertools.product(GROUPS,DX):
     q=masked.APOE.eq(g)&masked.Diagnosis.eq(s)
     n=int((q&np.isfinite(masked[e1])&np.isfinite(masked[e2])).sum());supports.append(n)
     pairrows.append({'anchor':anchor,'panel':panel,'selection':'masked_at_least2','event1':e1,'event2':e2,'APOE':g,'Diagnosis':s,'n':n})
   head[anchor][panel]['minimum_pair_cell_n']=min(supports)
   head[anchor][panel]['zero_pair_cells']=supports.count(0)
  for e in EVENTS:
   for g,s in itertools.product(GROUPS,DX):
    q=a.APOE.eq(g)&a.Diagnosis.eq(s)
    evrows.append({'anchor':anchor,'event':e,'APOE':g,'Diagnosis':s,'n':int((q&np.isfinite(a[e])).sum()),'denominator':int(q.sum())})
  patt=np.isfinite(a[list(CSF)+list(COG)]).astype(int).astype(str).agg(''.join,axis=1)
  for pattern,n in patt.value_counts().items():patrows.append({'anchor':anchor,'event_order':';'.join(list(CSF)+list(COG)),'observed_pattern':pattern,'n':int(n)})
 pd.DataFrame(rows).to_csv(OUT/'nine_cell_coverage.csv',index=False)
 pd.DataFrame(pairrows).to_csv(OUT/'event_pair_support.csv',index=False)
 pd.DataFrame(evrows).to_csv(OUT/'event_coverage.csv',index=False)
 pd.DataFrame(patrows).to_csv(OUT/'missingness_patterns.csv',index=False)
 csfa=frames['CSF'];complete=csfa.loc[csfa.core_covariates&finite(csfa,PANELS['CSFcog6'])]
 platform={'source':{str(k):int(v) for k,v in csf.Platform.fillna('missing').value_counts().items()},'CSFcog6_complete':{str(k):int(v) for k,v in complete.Platform.fillna('missing').value_counts().items()}}
 head['source']={'csf_rows':len(csf),'csf_participants':csf.NACCID.nunique(),'excluded_ambiguous_FS_rows':int(badfs.sum()),'platform':platform,'source_date_fields':{'CSF_draw_date':False,'MRI_scan_date':False,'PET_scan_date':True}}
 dump(OUT/'coverage_summary.json',head)
 dump(OUT/'provenance.json',{'inputs':{k:{'path':provenance_name(k),'sha256':sha(v)} for k,v in FILES.items()},'script_sha256':sha(Path(__file__)),'software':{'numpy':np.__version__,'pandas':pd.__version__},'no_fitting_or_inference':True,'participant_exports':'local_only only; no HPC file created','selection':'earliest eligible modality-linked UDS before cognition/covariate completeness; no later rescue','cognition_link':'exact NACCID+NACCVNUM','MRI_link':'exact NACCID+NACCVNUM; bilateral means require both sides','PET_link':'nearest finite dated value within180days to immutable anchor; coverage only, QC recorded not filtered','clinical_rules':'CN=NACCUDSD1; MCI=NACCUDSD3 prevalent; AD=NACCUDSD4+NACCALZD1; e2 codes3/6,e33code1,e4codes2/4;exclude2/4code5','actual_CSF_and_MRI_to_cognition_intervals_verified':False})
 print(json.dumps({'CSF':head['CSF'],'MRI':{k:head['MRI'][k] for k in ['anchor_n','MRIcog9','CSFcog6']},'platform':platform},indent=2))

if __name__=='__main__':main()
