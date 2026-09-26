#!/usr/bin/env python3
"""Frozen shared MRI/cognition validation. Never selects panels from outcomes."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'vendor'))
from concord import compare, prepare_data
from concord.core import CompareSpec, permute_labels, scheme_definitions
from concord._pyebm import verify_pyebm

EVENTS = ['MMSE', 'MEM', 'Hippocampus', 'Entorhinal', 'Fusiform', 'MiddleTemporal', 'Precuneus', 'Ventricles']
PANELS = {'MRIcog4': EVENTS[:4], 'MRIcog8': EVENTS}
GROUPS = ('e2', 'e33', 'e4')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False, default=str) + '\n')
    tmp.replace(path)

def correct(raw):
    """One cohort-wide correction, independent of genotype and selected panel."""
    x = raw.copy()
    assert x.PTID.is_unique and x[EVENTS + ['Age', 'Education', 'ICV']].apply(np.isfinite).all().all()
    if set(x.Sex) != {'Male', 'Female'}:
        raise ValueError('Expected both Male and Female, without missing Sex')
    x['SexMale'] = x.Sex.eq('Male').astype(float)
    cn = x.Diagnosis.eq('CN').to_numpy()
    log = {}
    for name, events, covs in [('cognition', EVENTS[:2], ['Age', 'SexMale', 'Education']),
                               ('MRI', EVENTS[2:], ['Age', 'SexMale', 'ICV'])]:
        z = x[covs].to_numpy(float)
        design = np.column_stack([np.ones(cn.sum()), z[cn]])
        if np.linalg.matrix_rank(design) != len(covs) + 1:
            raise ValueError('CN adjustment design is rank deficient')
        center = z.mean(axis=0)
        for event in events:
            beta, _, _, _ = np.linalg.lstsq(design, x.loc[cn, event].to_numpy(float), rcond=None)
            x[event] = x[event].to_numpy(float) - (z - center) @ beta[1:]
            log[event] = dict(block=name, covariates=covs, control_n=int(cn.sum()),
                              full_frame_center=center.tolist(), intercept=float(beta[0]), slopes=beta[1:].tolist())
    return x[['PTID', 'APOE', 'Diagnosis'] + EVENTS], log

def read_campaign(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    for rel, expected in config['frozen_files'].items():
        if sha(path.parent / rel) != expected:
            raise ValueError(f'Frozen file changed: {rel}')
    return config, path.parent

def data_for(config, base, cohort, mode, rep):
    c = config['cohorts'][cohort]
    frame = pd.read_csv(base / c['adjusted_file'], float_precision='round_trip')
    assert len(frame) == c['n'] and frame.PTID.is_unique and not frame[EVENTS].isna().any().any()
    if mode == 'apoe':
        return frame, GROUPS, c['apoe_seed'], c['adjusted_sha256']
    plans = np.load(base / c['draw_file'])
    a, b = plans['a'][rep], plans['b'][rep]
    assert not set(a) & set(b)
    assert frame.iloc[np.r_[a,b]].APOE.eq('e33').all()
    pa, pb = frame.iloc[a].copy(), frame.iloc[b].copy()
    for part, counts in ((pa, c['n2_counts']['a']), (pb, c['n2_counts']['b'])):
        assert part.Diagnosis.value_counts().reindex(c['labels'], fill_value=0).tolist() == counts
    pa['APOE'], pb['APOE'] = 'pA', 'pB'
    data = pd.concat([pa, pb], ignore_index=True)
    draw_hash = hashlib.sha256(np.r_[a,b].astype('<i8').tobytes() + c['adjusted_sha256'].encode()).hexdigest()
    return data, ('pA', 'pB'), c['n2_fit_seed_base'] + rep, draw_hash

def validate(payload, mode, estimator, budget):
    assert payload['status'] == 'ok', payload['status']
    assert payload['fits']['error'] == payload['fits']['timeout'] == 0, payload['fits']
    expected = {'unrestricted', 'diagnosis'} if estimator == 'standard' else {'diagnosis'}
    if mode == 'apoe':
        expected = {'diagnosis_pair_e2-e33', 'diagnosis_pair_e2-e4', 'diagnosis_pair_e33-e4'}
    assert set(payload['tests']) == expected
    assert payload['decisions']['stop_when_decided'] is False
    assert payload['spec']['consensus'] == 'repaired' and payload['spec']['stability_resamples'] == 0
    assert payload['fits']['ok'] == 1 + len(expected) * budget
    for key, t in payload['tests'].items():
        assert t['budget'] == t['completed'] == budget and t['failed'] == 0
        assert len(t['pairs']) == 1
        if mode == 'apoe':
            assert set(t['pairs']) == {key.removeprefix('diagnosis_pair_')}
        for cell in t['pairs'].values():
            p = (cell['exceedances'] + 1) / (budget + 1)
            assert cell['p'] == p and cell['p_bounds'] == [p,p]
            assert np.isclose(cell['threshold'], .05 / (3 if mode == 'apoe' else 1))
            assert cell['reject'] == (p <= .05 / (3 if mode == 'apoe' else 1))

def compact(payload):
    """Only deterministic scientific fields, for serial/parallel checks."""
    d = {key: payload[key] for key in ['groups', 'biomarkers', 'orderings', 'positions', 'distances',
                                      'effective_sample_size', 'composition', 'fits']}
    d['tests'] = {k: {name: t[name] for name in ['budget','completed','failed','pairs','max']}
                  for k,t in payload['tests'].items()}
    return d

def operator_audit(frame, labels, seed):
    data, groups, events = prepare_data(frame, group_column='APOE', labels=labels, biomarkers=EVENTS, group_order=GROUPS)
    spec = CompareSpec(group_column='APOE', group_values=groups, labels=tuple(labels), biomarkers=events, seed=seed)
    defs = scheme_definitions(spec, ['diagnosis_pair'])
    for name, definition in defs.items():
        for index in (0,1,7):
            shuffled = permute_labels(data, spec, name, definition, index)
            allowed = data.APOE.isin(definition['groups'])
            assert shuffled.loc[~allowed].equals(data.loc[~allowed])
            assert shuffled.drop(columns='APOE').equals(data.drop(columns='APOE'))
            assert pd.crosstab(shuffled.APOE, shuffled.Diagnosis).equals(pd.crosstab(data.APOE,data.Diagnosis))
    return {'third_group_fixed': True, 'diagnosis_counts_preserved': True, 'biomarkers_and_rows_unchanged': True}

def run_one(args, config, base, cohort, mode, panel, estimator, rep, budget, workers, directory):
    data, groups, seed, draw_hash = data_for(config, base, cohort, mode, rep)
    out = Path(directory) / cohort / mode / panel / estimator / f'rep_{rep:04d}.json'
    identity = dict(cohort=cohort, mode=mode, panel=panel, estimator=estimator, replicate=rep,
                    B=budget, seed=seed, n=len(data), draw_sha256=draw_hash, campaign_sha256=sha(args.config))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if out.exists():
            old = json.loads(out.read_text())
            assert old['analysis'] == identity
            validate(old['result'], mode, estimator, budget)
            return old
        schemes = ('diagnosis_pair',) if mode == 'apoe' else (('unrestricted','diagnosis') if estimator == 'standard' else ('diagnosis',))
        start = time.monotonic()
        result = compare(data, group_column='APOE', labels=config['cohorts'][cohort]['labels'], biomarkers=PANELS[panel], group_order=groups,
                         estimator=estimator, consensus='repaired', schemes=schemes, B=budget, alpha=.05,
                         rule='le', stop_when_decided=False, stability_resamples=0, seed=seed, workers=workers,
                         fit_timeout_s=900, fast_likelihood=True, verbose=True, progress_every=30)
        payload = result.to_dict()
        value = dict(analysis=identity, result=payload, wall_seconds=time.monotonic()-start,
                     job_id=os.environ.get('SLURM_JOB_ID'), array_id=os.environ.get('SLURM_ARRAY_TASK_ID'))
        # Keep unsuccessful fits for diagnosis, never reclassify or change denominator.
        try:
            validate(payload, mode, estimator, budget)
        except Exception:
            save(out.with_name(out.stem + '_FAILED_' + datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '.json'), value)
            raise
        save(out, value)
        print(json.dumps({'completed':str(out), 'elapsed_s':value['wall_seconds'], 'fit_counts':payload['fits']}), flush=True)
        return value

def preflight(args, config, base):
    c = config['cohorts'][args.cohort]
    raw = pd.read_csv(base / c['raw_file'], float_precision='round_trip')
    rebuilt, _ = correct(raw)
    stored = pd.read_csv(base / c['adjusted_file'], float_precision='round_trip')
    # Cross-platform linalg need not produce bit-identical coefficients.
    assert np.allclose(rebuilt[EVENTS], stored[EVENTS], rtol=1e-9, atol=1e-7)
    audit = operator_audit(stored, c['labels'], c['apoe_seed'])
    rows = []
    tasks = [('n2',panel,estimator) for panel in PANELS for estimator in ('standard','invariant_min')]
    tasks.append(('apoe','MRIcog8','invariant_min'))
    for mode,panel,estimator in tasks:
        one = run_one(args,config,base,args.cohort,mode,panel,estimator,0,8,1,args.out / 'serial')
        many = run_one(args,config,base,args.cohort,mode,panel,estimator,0,8,args.workers,args.out / 'parallel')
        assert compact(one['result']) == compact(many['result']), (mode,panel,estimator)
        rows.append(dict(mode=mode,panel=panel,estimator=estimator,serial_seconds=one['wall_seconds'],
                         parallel_seconds=many['wall_seconds'],fits=many['result']['fits']['ok']))
    save(args.out / args.cohort / 'preflight_passed.json', dict(passed=True, operator_audit=audit,
         serial_parallel_B8_match=True, frozen_input_verified=True, timings=rows, pyebm=verify_pyebm(), campaign_sha256=sha(args.config)))

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=HERE/'campaign.json')
    ap.add_argument('--mode',choices=['preflight','n2','apoe','input-check'],required=True)
    ap.add_argument('--cohort',choices=['adni','nacc'],required=True)
    ap.add_argument('--panel',choices=PANELS,default='MRIcog8')
    ap.add_argument('--estimator',choices=['standard','invariant_min'],default='invariant_min')
    ap.add_argument('--start',type=int,default=0)
    ap.add_argument('--stop',type=int,default=1)
    ap.add_argument('--workers',type=int,default=8)
    ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--preflight-root',type=Path)
    args = ap.parse_args()
    config, base = read_campaign(args.config)
    if args.mode == 'input-check':
        frame,_,_,_=data_for(config,base,args.cohort,'apoe',0)
        prepare_data(frame,group_column='APOE',labels=config['cohorts'][args.cohort]['labels'],biomarkers=EVENTS,group_order=GROUPS)
        print(json.dumps({'cohort':args.cohort,'n':len(frame),'frozen_files_verified':True,'pyebm':verify_pyebm()}))
        return
    if args.mode == 'preflight':
        preflight(args,config,base)
        return
    if args.preflight_root is None:
        raise ValueError('Production requires --preflight-root')
    gate = json.loads((args.preflight_root / args.cohort / 'preflight_passed.json').read_text())
    assert gate['passed'] and gate['campaign_sha256'] == sha(args.config)
    if args.mode == 'apoe':
        assert args.panel == 'MRIcog8' and args.estimator == 'invariant_min'
        run_one(args,config,base,args.cohort,'apoe','MRIcog8','invariant_min',0,config['B'],args.workers,args.out)
    else:
        assert 0 <= args.start < args.stop <= config['R']
        for rep in range(args.start,args.stop):
            run_one(args,config,base,args.cohort,'n2',args.panel,args.estimator,rep,config['B'],args.workers,args.out)

if __name__ == '__main__':
    main()
