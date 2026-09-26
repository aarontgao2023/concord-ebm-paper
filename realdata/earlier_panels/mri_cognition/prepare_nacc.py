#!/usr/bin/env python3
"""NACC MRI-cognition population (1,886 participants; Supplementary Methods 2): write the input.

Earlier analysis population behind the NACC MRI-cognition rows of Fig. 4e. Links each ADSP-PHC
FreeSurfer (ComBat) record of the NACC release to its UDS visit (NACCID + NACCVNUM), keeps the
earliest eligible linked visit per participant (APOE e2, e3/e3 or e4; CN, prevalent MCI or AD
dementia with AD aetiology) and then participants with all eight events (MMSE, ADSP-PHC memory
and six MRI volumes) plus age, sex, education and ICV. No fitting or inference; the runner
(run_shared_panel.correct) adjusts all eight events once before the 4-event subset is taken.

  python prepare_nacc.py --phc-dir $CONCORD_NACC_DIR --uds-file $CONCORD_NACC_UDS \\
      --pet-cohort <derived>/nacc_cohort_full.csv --out $CONCORD_WORK_DIR/mri_cognition/nacc

--pet-cohort (from ../nacc/build_cohort_nacc.py) is used only to count the overlap with the
amyloid PET population. --out must lie outside the repository; it receives a pseudonymous raw
input (nacc_common_raw_deid.csv) and a local-only linkage file. The script asserts the population
size of Supplementary Methods 2 (1,886); other counts can be checked with --expected-counts, a
local JSON file {"anchors": n, "feature_complete": n, "cells": [[..3 x 3..]]} (rows e2, e3/e3,
e4; columns CN, MCI, AD). The next step is freeze_campaign.py.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

KEY = ['NACCID', 'NACCVNUM']
GROUPS = ['e2', 'e33', 'e4']
LABELS = ['CN', 'MCI', 'AD']
EVENTS = ['MMSE', 'MEM', 'Hippocampus', 'Entorhinal', 'Fusiform',
          'MiddleTemporal', 'Precuneus', 'Ventricles']
K4 = EVENTS[:4]
MRI_FIELDS = {
    'Hippocampus': ['Left.Hippocampus_combat', 'Right.Hippocampus_combat'],
    'Entorhinal': ['lh_entorhinal_volume_combat', 'rh_entorhinal_volume_combat'],
    'Fusiform': ['lh_fusiform_volume_combat', 'rh_fusiform_volume_combat'],
    'MiddleTemporal': ['lh_middletemporal_volume_combat', 'rh_middletemporal_volume_combat'],
    'Precuneus': ['lh_precuneus_volume_combat', 'rh_precuneus_volume_combat'],
    'Ventricles': ['Left.Lateral.Ventricle_combat', 'Right.Lateral.Ventricle_combat'],
}
ICV_FIELD = 'EstimatedTotalIntraCranialVol_combat'
REPO = Path(__file__).resolve().parents[3]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def normalize_keys(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if df[KEY].isna().any().any():
        raise ValueError('Missing participant/visit key; no values printed.')
    df.NACCID = df.NACCID.astype(str)
    visit = pd.to_numeric(df.NACCVNUM, errors='raise')
    if not np.equal(visit, np.floor(visit)).all():
        raise ValueError('Noninteger visit key; no values printed.')
    df.NACCVNUM = visit.astype(int)
    return df


def finite_matrix(df: pd.DataFrame, columns: list[str]) -> np.ndarray:
    return np.isfinite(df[columns].to_numpy(dtype=float)).all(axis=1)


def dump_json(path: Path, content: dict) -> None:
    path.write_text(json.dumps(content, indent=2, allow_nan=False) + '\n')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--phc-dir', type=Path, default=os.environ.get('CONCORD_NACC_DIR'),
                        help='NACC ADSP-PHC release folder (default: $CONCORD_NACC_DIR)')
    parser.add_argument('--uds-file', type=Path, default=os.environ.get('CONCORD_NACC_UDS'),
                        help='CSV export of the NACC UDS file (default: $CONCORD_NACC_UDS)')
    parser.add_argument('--pet-cohort', type=Path, required=True,
                        help='nacc_cohort_full.csv written by ../nacc/build_cohort_nacc.py')
    parser.add_argument('--out', type=Path, default=None,
                        help='output folder outside the repository (default: $CONCORD_WORK_DIR/mri_cognition/nacc)')
    parser.add_argument('--expected-counts', type=Path, help='optional local JSON with expected counts')
    args = parser.parse_args(argv)
    if args.out is None and os.environ.get('CONCORD_WORK_DIR'):
        args.out = Path(os.environ['CONCORD_WORK_DIR']) / 'mri_cognition/nacc'
    if args.phc_dir is None or args.uds_file is None or args.out is None:
        parser.error('--phc-dir, --uds-file and --out are required (or set CONCORD_NACC_DIR, CONCORD_NACC_UDS and CONCORD_WORK_DIR)')
    out = args.out.resolve()
    if out == REPO or REPO in out.parents:
        raise SystemExit(f'The output folder must lie outside the repository: {out}')
    expected = json.loads(args.expected_counts.read_text()) if args.expected_counts else {}
    out.mkdir(parents=True, exist_ok=True)
    # These outputs are private research inputs. Never print individual records.
    os.chmod(out, 0o700)
    os.umask(0o077)
    new = args.phc_dir.resolve()
    paths = {
        'freesurfer': new / 'Imaging_T1/NACC_ADSP_PHC_T1_Freesurfer_2024.csv',
        'cognition': new / 'Cognition/NACC_ADSP_PHC_Cognition_2024.csv',
        'uds': args.uds_file.resolve(),
        'existing_cohort': args.pet_cohort.resolve(),
    }
    input_provenance = {k: {'path': str(p.relative_to(new)) if p.is_relative_to(new) else p.name, 'sha256': sha256(p)}
                        for k, p in paths.items()}
    if (out / 'manifest.json').exists():
        prior_inputs = json.loads((out / 'manifest.json').read_text()).get('inputs')
        if prior_inputs != input_provenance:
            raise ValueError('Frozen input source hashes changed; existing freeze will not be overwritten.')
    fs_columns = KEY + [c for pair in MRI_FIELDS.values() for c in pair] + [ICV_FIELD]
    # Duplicate status is evaluated on keys from the full FS row inventory.
    fs = normalize_keys(pd.read_csv(paths['freesurfer'], usecols=fs_columns, low_memory=False))
    cog = normalize_keys(pd.read_csv(paths['cognition'], usecols=KEY + ['PHC_MEM'], low_memory=False))
    uds_columns = KEY + ['VISITYR', 'VISITMO', 'VISITDAY', 'NACCAPOE', 'NACCUDSD',
                         'NACCALZD', 'NACCAGE', 'NACCMMSE', 'SEX', 'EDUC']
    uds = normalize_keys(pd.read_csv(paths['uds'], usecols=uds_columns, low_memory=False))
    old = normalize_keys(pd.read_csv(paths['existing_cohort'], usecols=KEY, low_memory=False))
    if cog.duplicated(KEY).any() or uds.duplicated(KEY).any():
        raise ValueError('Cognition or UDS visit keys are ambiguous; no keys printed.')
    flow = []

    def stage(name: str, d: pd.DataFrame) -> None:
        flow.append({'stage': name, 'rows': int(len(d)), 'participants': int(d.NACCID.nunique())})

    stage('FreeSurfer source', fs)
    ambiguous = fs.duplicated(KEY, keep=False)
    duplicate_key_count = int(fs.loc[ambiguous, KEY].drop_duplicates().shape[0])
    duplicate_row_count = int(ambiguous.sum())
    fs = fs.loc[~ambiguous].copy()
    stage('Exclude all ambiguous FreeSurfer ID-visit keys', fs)
    m = fs.merge(uds, on=KEY, how='inner', validate='one_to_one')
    stage('Exact UDS participant-visit match', m)
    m['APOE'] = m.NACCAPOE.map({1: 'e33', 2: 'e4', 3: 'e2', 4: 'e4', 5: 'e2e4', 6: 'e2'}).fillna('NA')
    m['Diagnosis'] = np.select([m.NACCUDSD.eq(1), m.NACCUDSD.eq(3),
                               m.NACCUDSD.eq(4) & m.NACCALZD.eq(1)],
                              ['CN', 'MCI', 'AD'], default='excluded')
    m = m.loc[m.APOE.isin(GROUPS) & m.Diagnosis.isin(LABELS)].copy()
    stage('Known eligible APOE and existing CN-prevalentMCI-AD rule', m)
    m['UDS_date'] = pd.to_datetime(dict(year=m.VISITYR, month=m.VISITMO, day=m.VISITDAY), errors='coerce')
    if m.UDS_date.isna().any():
        raise ValueError('An eligible MRI-linked UDS visit lacks a valid date; freeze rule requires resolution.')
    # Freeze the anchor before looking for complete cognition/covariates. Do not
    # search later visits to rescue incomplete profiles.
    m = m.sort_values(['NACCID', 'UDS_date', 'NACCVNUM'], kind='mergesort').drop_duplicates('NACCID').copy()
    stage('First eligible MRI-linked UDS visit before feature completeness', m)
    if 'anchors' in expected and len(m) != expected['anchors']:
        raise AssertionError(f"MRI anchor count changed: {len(m)}; expected {expected['anchors']}.")
    m = m.merge(cog, on=KEY, how='left', validate='one_to_one')
    m['MMSE'] = pd.to_numeric(m.NACCMMSE, errors='coerce').where(lambda s: s.between(0, 30))
    m['MEM'] = pd.to_numeric(m.PHC_MEM, errors='coerce')
    for event, pair in MRI_FIELDS.items():
        # Addition propagates a missing side; DataFrame.mean would not.
        left = pd.to_numeric(m[pair[0]], errors='coerce')
        right = pd.to_numeric(m[pair[1]], errors='coerce')
        m[event] = (left + right) / 2.0
    m['Age'] = pd.to_numeric(m.NACCAGE, errors='coerce')
    m['Sex'] = m.SEX.map({1: 'Male', 2: 'Female'})
    m['Education'] = pd.to_numeric(m.EDUC, errors='coerce').where(lambda s: s.between(0, 98))
    m['ICV'] = pd.to_numeric(m[ICV_FIELD], errors='coerce')
    complete = m.loc[finite_matrix(m, EVENTS)].copy()
    stage('All eight events finite; same selected visit', complete)
    if 'feature_complete' in expected and len(complete) != expected['feature_complete']:
        raise AssertionError(f"Feature-complete count changed: {len(complete)}; expected {expected['feature_complete']}.")
    covariate_missing = {'Age': int((~np.isfinite(complete.Age)).sum()),
                         'Sex': int(complete.Sex.isna().sum()),
                         'Education': int((~np.isfinite(complete.Education)).sum()),
                         'ICV': int((~np.isfinite(complete.ICV)).sum())}
    final = complete.loc[finite_matrix(complete, ['Age', 'Education', 'ICV']) & complete.Sex.notna()].copy()
    if not final.ICV.gt(0).all():
        raise AssertionError('Nonpositive ICV in frozen candidate; no values printed.')
    final = final.sort_values('NACCID', kind='mergesort').reset_index(drop=True)
    stage('Complete age-sex-education-ICV; frozen common K4-K8 participants', final)
    counts = pd.crosstab(final.APOE, final.Diagnosis).reindex(index=GROUPS, columns=LABELS, fill_value=0)
    if len(final) != 1886 or ('cells' in expected and not np.array_equal(counts.to_numpy(), np.array(expected['cells']))):
        raise AssertionError('Frozen count/cell validation failed; inspect aggregate flow.')
    final.insert(0, 'PTID', [f'NCSP{i:05d}' for i in range(len(final))])
    schema = ['PTID', 'APOE', 'Diagnosis', 'Age', 'Sex', 'Education', 'ICV'] + EVENTS
    remote = final[schema].copy()
    linked = final[['NACCID', 'NACCVNUM', 'UDS_date'] + schema].copy()
    remote_path = out / 'nacc_common_raw_deid.csv'
    linked_path = out / 'nacc_common_linked_LOCAL_ONLY.csv'
    remote.to_csv(remote_path, index=False, float_format='%.17g')
    linked.to_csv(linked_path, index=False, float_format='%.17g', date_format='%Y-%m-%d')
    # Read back the exact HPC artifact to validate its content and field contract.
    reloaded = pd.read_csv(remote_path, float_precision='round_trip')
    pd.testing.assert_frame_equal(remote, reloaded, check_dtype=False, check_exact=True)
    counts.reset_index().to_csv(out / 'counts.csv', index=False)
    pd.DataFrame(flow).to_csv(out / 'flow.csv', index=False)
    dictionary = []
    for event in EVENTS:
        fields = MRI_FIELDS.get(event, ['NACCMMSE'] if event == 'MMSE' else ['PHC_MEM'])
        dictionary.append({'event': event, 'source_fields': ';'.join(fields),
            'construction': '(left+right)/2, both finite' if event in MRI_FIELDS else ('valid range 0-30' if event == 'MMSE' else 'supplied PHC score'),
            'source_scale': 'source volume units, PHC ComBat harmonized' if event in MRI_FIELDS else ('MMSE points' if event == 'MMSE' else 'PHC harmonized memory score'),
            'conventional_abnormal_direction': 'higher' if event == 'Ventricles' else 'lower',
            'manual_sign_flip': False, 'in_K4': event in K4, 'in_K8': True})
    pd.DataFrame(dictionary).to_csv(out / 'event_dictionary.csv', index=False)
    table = {g: {d: int(counts.loc[g, d]) for d in LABELS} for g in GROUPS}
    manifest = {
        'freeze_date': '2026-09-20', 'cohort': 'NACC common MRI-cognition, no PET selection',
        'n': len(final), 'groups': GROUPS, 'diagnosis_labels': LABELS, 'counts': table,
        'group_n': {g: int(counts.loc[g].sum()) for g in GROUPS},
        'panels': {'K4': K4, 'K8': EVENTS}, 'export_schema': schema,
        'inputs': input_provenance,
        'software': {'numpy': np.__version__, 'pandas': pd.__version__},
        'script': {'path': str(Path(__file__).resolve().relative_to(REPO)), 'sha256': sha256(Path(__file__))},
        'outputs': {
            'hpc_input': {'path': remote_path.name, 'sha256': sha256(remote_path), 'contains_original_identifiers': False},
            'local_linkage': {'path': linked_path.name, 'sha256': sha256(linked_path), 'contains_original_identifiers': True, 'transfer_to_HPC': False},
        },
        'source_version': 'ADSP-PHC December 2024 NACC return-to-cohort files',
        'ambiguous_FS_keys_excluded': duplicate_key_count, 'ambiguous_FS_rows_excluded': duplicate_row_count,
        'selection_order': ['exclude every duplicated FS ID-visit key', 'exact UDS join',
            'existing eligible APOE and diagnosis rule', 'earliest eligible linked UDS date, then visit number',
            'exact PHC cognition join', 'all K8 values finite', 'complete Age/Sex/Education/ICV'],
        'select_later_visit_to_rescue_incomplete_profile': False,
        'APOE_rule': 'codes 1=e33; 2/4=e4; 3/6=e2; exclude 5=e2e4 and unknown',
        'diagnosis_rule': 'CN=NACCUDSD1; prevalent MCI=NACCUDSD3; AD=NACCUDSD4 and NACCALZD1',
        'MRI_fields': MRI_FIELDS, 'ICV_field': ICV_FIELD,
        'feature_complete_n_before_covariates': len(complete),
        'feature_complete_missing_covariates': covariate_missing,
        'prior_PET_cohort_overlap_n': int(final.NACCID.isin(set(old.NACCID)).sum()),
        'prior_PET_cohort_exact_selected_visit_overlap_n': int(sum(k in set(map(tuple, old[KEY].values)) for k in map(tuple, final[KEY].values))),
        'timing': {'link': 'exact NACCID + NACCVNUM across MRI, UDS, cognition',
            'anchor_date': 'linked UDS calendar date', 'actual_MRI_acquisition_date_available': False,
            'actual_MRI_cognition_acquisition_interval_verified': False,
            'age_source': 'NACCAGE at selected linked UDS visit'},
        'adjustment': {'performed_by_preparer': False, 'performed_by': 'shared runner once on full frozen K8 frame',
            'fit_population': 'CN across all eligible APOE groups within this cohort',
            'cognition_factors': ['Age', 'SexMale', 'Education'], 'MRI_factors': ['Age', 'SexMale', 'ICV'],
            'sex_coding': 'SexMale = 1 for Male, 0 for Female',
            'formula': 'y_adjusted = y - (X - mean_full_frozen_cohort(X)) @ beta_slopes; OLS intercept included when fitting',
            'K4_rule': 'subset the same corrected columns and same participants; never refit covariate models for K4',
            'manual_sign_flip': False, 'orientation': 'CN/AD-labelled measurement-model initialization in existing inference package'},
        'scope': 'authorized new common-panel validation; does not overwrite original manuscript cohorts or analyses',
        'model_fits_permutations_resampling_performed': False,
    }
    dump_json(out / 'manifest.json', manifest)
    verification = {'frozen_n_1886': len(final) == 1886, 'expected_nine_cells_match': True if 'cells' in expected else 'not checked',
        'unique_source_participants': bool(final.NACCID.is_unique), 'unique_pseudonyms': bool(remote.PTID.is_unique),
        'all_K8_finite': bool(finite_matrix(remote, EVENTS).all()), 'covariates_complete': True,
        'K4_exact_nested_subset': K4 == EVENTS[:4], 'CSV_round_trip_exact': True,
        'HPC_file_has_no_original_ID_or_visit_date': not bool(set(remote.columns) & {'NACCID','NACCVNUM','UDS_date'}),
        'all_pair_by_APOE_diagnosis_cells_supported': bool((counts > 0).all().all()),
        'minimum_pair_cell_n': int(counts.to_numpy().min()), 'no_model_fit': True}
    dump_json(out / 'verification.json', verification)
    for p in out.iterdir():
        if p.is_file():
            os.chmod(p, 0o600)
    print(json.dumps({'cohort':'nacc', 'n':len(final), 'counts':table, 'hpc_input':str(remote_path),
                      'input_sha256':sha256(remote_path), 'minimum_pair_cell_n':int(counts.to_numpy().min()), 'fits_run':0}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
