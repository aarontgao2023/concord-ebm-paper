#!/usr/bin/env python3
"""NACC six-event panel, step 2: write the unadjusted de-identified input (no model fits).

Keeps the anchors of nacc_audit.py with all six events, age, sex and education and writes
<out>/prepared/nacc_csfcog6_raw_deid.csv (fresh sequential pseudonyms; no NACCID, visit
numbers or dates). The linkage goes to <out>/local_only only. The six-event panel was fixed on
23 September 2026 from coverage alone. The asserted counts are those of Table 2 (1,584
participants).

  python nacc_prepare_core6.py --phc-dir $CONCORD_NACC_DIR --uds-file $CONCORD_NACC_UDS \\
      --out $CONCORD_WORK_DIR/audit/nacc

Use the same arguments as for nacc_audit.py. The next step is ../prepare_campaign.py.
"""
from pathlib import Path
import json
import os
import numpy as np
import pandas as pd
import nacc_audit
from nacc_audit import FILES, CSF, COG, GROUPS, DX, sha, dump, counts, finite

def main(argv=None):
 nacc_audit.arguments(__doc__, argv)
 OUT = nacc_audit.OUT
 os.umask(0o077)
 provenance=json.loads((OUT/'provenance.json').read_text())
 for name,path in FILES.items():
  assert sha(path)==provenance['inputs'][name]['sha256'],f'Input changed: {name}'
 source=OUT/'local_only/csf_anchor_LOCAL_ONLY.csv'
 a=pd.read_csv(source,float_precision='round_trip')
 events=list(CSF)+list(COG)
 eligible=finite(a,events+['Age','Education'])&a.Sex.isin(['Male','Female'])
 final=a.loc[eligible].sort_values('NACCID',kind='mergesort').reset_index(drop=True).copy()
 assert final.NACCID.is_unique and len(final)==1584
 expected={'e2':{'CN':106,'MCI':15,'AD':22},'e33':{'CN':480,'MCI':104,'AD':175},'e4':{'CN':281,'MCI':94,'AD':307}}
 assert counts(final)==expected
 final.insert(0,'PTID',[f'NCM26N{i+1:05d}' for i in range(len(final))])
 columns=['PTID','APOE','Diagnosis','Age','Sex','Education']+events
 data=final[columns].copy()
 out=OUT/'prepared';out.mkdir(exist_ok=True)
 file=out/'nacc_csfcog6_raw_deid.csv'
 data.to_csv(file,index=False,float_format='%.17g')
 reload=pd.read_csv(file,float_precision='round_trip')
 pd.testing.assert_frame_equal(data,reload,check_dtype=False,check_exact=True)
 final[['NACCID','NACCVNUM','UDS_date','Platform']+columns].to_csv(OUT/'local_only/core6_linkage_LOCAL_ONLY.csv',index=False,float_format='%.17g')
 platform=[]
 for g in GROUPS:
  for d in DX:
   x=final.loc[final.APOE.eq(g)&final.Diagnosis.eq(d),'Platform'].fillna('missing')
   for p in ['missing','Innogenetics']:
    platform.append({'APOE':g,'Diagnosis':d,'Platform':p,'n':int(x.eq(p).sum()),'cell_n':len(x)})
 pd.DataFrame(platform).to_csv(OUT/'platform_by_nine_cells.csv',index=False)
 pd.DataFrame([{'APOE':g,'Diagnosis':d,'n':expected[g][d]} for g in GROUPS for d in DX]).to_csv(out/'counts.csv',index=False)
 dictionary=[]
 for e,c in dict(CSF,**COG).items():
  dictionary.append({'event':e,'source_field':c,'scale':'PHC harmonized log10 concentration z-score' if e in CSF else 'PHC harmonized cognition score',
   'higher_abnormal':e in ['TAU','PTAU'],'manual_sign_flip':False,'covariates':['Age','SexMale'] if e in CSF else ['Age','SexMale','Education'],
   'definition_note':('CSF p-tau; NACC instrument CSFPTAU documents p-tau181P, but return lacks explicit CSFPTAU-to-PHC_pTau crosswalk; not same-assay claim' if e=='PTAU' else '')})
 manifest={
  'status':'FROZEN_CORE6_SOURCE_INPUT','freeze_date':'2026-09-23','panel_fixed_from_coverage_only':True,'cohort':'NACC CSF-linked cognition, prevalent MCI',
  'n':len(data),'counts':expected,'events':events,'nested_sensitivity':{'CSFcog5_noPTAU_same_people':[e for e in events if e!='PTAU'],'participants_and_adjusted_shared_columns_unchanged':True},
  'groups':GROUPS,'diagnosis_labels':DX,'export_schema':columns,'contains_original_identifiers_or_dates':False,
  'inputs':provenance['inputs'],'anchor_private_frame_sha256':sha(source),'input_sha256':sha(file),
  'input_path':str(file.relative_to(OUT)),'preparer_sha256':sha(Path(__file__)),
  'event_dictionary':dictionary,'preprocessing_performed':False,
  'preprocessing_plan':'prepare_campaign.py: within-cohort pooled-CN OLS, CSF Age/Sex; cognition Age/Sex/Education, apply frozen adjustment before group permutations; no imputation or extra log transform',
  'timing':'Exact supplied CSF-linked NACCID/NACCVNUM to UDS and PHC cognition; actual CSF draw date unavailable, actual acquisition interval not verified',
  'selection':'earliest eligible CSF-linked UDS before completeness, then complete core6 and Age/Sex/Education; no later rescue; no model-derived selection',
  'platform_inference':'Missing platform remains unknown; retained and audited. Shared scores fit separately within ADNI and NACC, no participant pooling across cohorts.',
  'pTau_mapping_evidence':{'source_instrument':'https://files.alz.washington.edu/documentation/biomarker-ee2-csf-ded.pdf','page':6,'definition':'CSFPTAU: P-tau181P reported value/concentration (pg/mL)','PHC_direct_crosswalk_available':False},
  'no_fitting_or_permutation':True,'CSV_roundtrip_exact':True,
 }
 dump(out/'manifest.json',manifest)
 print(json.dumps({'path':str(file),'n':len(data),'sha256':sha(file),'counts':expected,'status':manifest['status']},indent=2))

if __name__=='__main__':main()
