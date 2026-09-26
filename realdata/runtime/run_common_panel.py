#!/usr/bin/env python3
"""Append-only common-panel validation; formal launch requires a frozen gate.

N2 supplies ordering-response estimates only. APOE supplies full599 D-pair tests.
Preflight uses its own directory and fit/permutation seed namespace.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
import fcntl
import importlib.util
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import time
import numpy as np
from common import (GROUPS, LABELS, ENGINES, DESIGNS, PAIR_CODES, read_campaign, read_frame,
                    n2_frame, pair_shuffle, bootstrap_frame, file_sha, frame_sha, label_sha, canonical, save)
from engine_adapter import fit_task, ENGINE_DIR
from dependency_v2 import verify_pyebm
HERE = Path(__file__).resolve().parent
PREFLIGHT_OFFSET = 1000000000
# Engine files read from ENGINE_DIR (simulation/scripts); recorded under the names of the
# vendor/current copies that ran, so the hashes in EXECUTION_GATE.json still apply.
ENGINE_FILES = ('dependency_v2.py', 'engine_v2.py', 'fast_likelihood_v2.py', 'invariant_engine_v2.py', 'kde_engine_v2.py')


def sources():
    value = {str(p.relative_to(HERE)):file_sha(p) for p in sorted(HERE.rglob('*.py')) if '__pycache__' not in str(p)}
    value.update({f'vendor/current/{name}':file_sha(ENGINE_DIR/name) for name in ENGINE_FILES})
    return value


def environment():
    value = verify_pyebm()
    kde = importlib.util.find_spec('kde_ebm')
    if kde is None or kde.origin is None:
        raise RuntimeError('kde_ebm is required, including its GMM implementation')
    root = Path(kde.origin).parent
    value['kde_ebm_source_sha256'] = {str(p.relative_to(root)):file_sha(p) for p in sorted(root.rglob('*.py'))}
    return value


@contextmanager
def unit(folder, identity):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/'lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path = folder/'manifest.json'
        identity = identity | {'identity_sha256':canonical(identity)}
        if path.exists():
            if json.loads(path.read_text()) != identity:
                raise ValueError('Existing output identity differs; do not mix revisions')
        else:
            save(path, identity)
        yield


def task(data, cfg, engine, key, fit_seed, **details):
    biomarkers = details.pop('biomarkers', cfg['biomarkers'])
    return {'data':data[['PTID','APOE','Diagnosis']+list(biomarkers)].copy(), 'engine':engine,
            'key':key, 'biomarkers':list(biomarkers), 'reference':cfg['reference'],
            'fit_seed':int(fit_seed),'fit_timeout_s':cfg.get('fit_timeout_s',900), **details}


def run_records(folder, tasks, pool):
    expected = {t['key']:t for t in tasks}
    if len(expected)!=len(tasks):
        raise ValueError('Duplicate planned task key')
    records = {}
    path = folder/'fit_records.jsonl'
    if path.exists():
        with path.open('rb') as stream:
            for raw in stream:
                if not raw.endswith(b'\n'):
                    raise ValueError('Unterminated append-only record: preserve and review before recovery')
                row = json.loads(raw)
                key = row['key']
                if key not in expected or key in records:
                    raise ValueError('Duplicate or unplanned existing raw record')
                wanted = expected[key]
                if row['input_frame_sha256'] != frame_sha(wanted['data']) or row['permuted_labels_sha256'] != label_sha(wanted['data'].APOE):
                    raise ValueError('Existing raw record input/label mismatch')
                if any(row.get(k)!=v for k,v in wanted.items() if k!='data'):
                    raise ValueError('Existing raw record specification mismatch')
                records[key] = row
    # Failed records are retained, not silently retried or removed.
    pending = [t for t in tasks if t['key'] not in records]
    with path.open('a',buffering=1) as output:
        def accept(row):
            output.write(json.dumps(row,allow_nan=False)+'\n');output.flush();os.fsync(output.fileno())
            records[row['key']]=row
            if len(records)%30==0:
                save(folder/'progress.json',{'recorded':len(records),'planned':len(tasks),'failed':sum(r['status']!='ok' for r in records.values())})
        if pool is None:
            for one in pending:
                accept(fit_task(one))
        else:
            futures = [pool.submit(fit_task,one) for one in pending]
            for future in as_completed(futures):
                accept(future.result())
    if set(records)!=set(expected):
        raise ValueError('Raw-record set is incomplete')
    return records


def identity(args,cfg,env,mode,panel,rep,preflight):
    return {'cohort':args.cohort,'mode':mode,'panel':panel,'replicate':rep,'preflight':preflight,
            'campaign_sha256':file_sha(args.config),'source_sha256':sources(),'environment':env,
            'reference':cfg['reference'],'B':8 if preflight else cfg['B'],'R':cfg['R']}


def run_n2(args,cfg,base,frame,env,rep,pool,preflight=False):
    folder=args.out/args.cohort/'n2'/f'rep_{rep:04d}'
    tasks=[];draws={}
    fit_seed=cfg['cohorts'][args.cohort]['fit_seed_base']+rep+(PREFLIGHT_OFFSET if preflight else 0)
    for design in DESIGNS:
        data,draw=n2_frame(frame,cfg,args.cohort,rep,design,base)
        draws[design]=draw
        for engine in ENGINES:
            tasks.append(task(data,cfg,engine,f'{design}/{engine}',fit_seed,kind='observed',design=design,replicate=rep))
    if draws['matched']['selected_roster_sha256']!=draws['different']['selected_roster_sha256']:
        raise ValueError('Composition conditions are not paired on the same participants')
    with unit(folder,identity(args,cfg,env,'n2','primary',rep,preflight)|{'draws':draws}):
        records=run_records(folder,tasks,pool)
        result={'complete':all(r['status']=='ok' for r in records.values()),'expected_records':8,
                'failed_records':[k for k,r in records.items() if r['status']!='ok'],
                'draws':draws,'estimates':{k:{name:r.get(name) for name in ('status','orderings','pair_distances','diagnostics')} for k,r in records.items()},
                'raw_fit_records_sha256':file_sha(folder/'fit_records.jsonl')}
        save(folder/'result.json',result)
    return result


def encoded(frame,biomarkers):
    data=frame[['PTID','APOE','Diagnosis']+list(biomarkers)].copy()
    data['APOE']=data.APOE.map({g:i for i,g in enumerate(GROUPS)}).astype(int)
    return data


def run_apoe(args,cfg,frame,env,panel,pool,preflight=False):
    biomarkers=cfg['panels'][panel];data=encoded(frame,biomarkers)
    seed=cfg['cohorts'][args.cohort]['apoe_seed']+(PREFLIGHT_OFFSET if preflight else 0)
    budget=8 if preflight else 599
    tasks=[task(data,cfg,'concord','observed',seed,kind='observed',biomarkers=biomarkers)]
    streams={}
    for pair_index,(a,b) in enumerate(PAIR_CODES):
        import hashlib
        stream=hashlib.sha256()
        for pid in range(budget):
            shuffled=pair_shuffle(data,seed,pair_index,pid)
            stream.update(np.asarray(shuffled.APOE,dtype='<i8').tobytes())
            tasks.append(task(shuffled,cfg,'concord',f'pair_{pair_index}/{pid}',seed,kind='permutation',
                              pair_index=pair_index,perm_id=pid,biomarkers=biomarkers))
        streams[f'{a}-{b}']=stream.hexdigest()
    folder=args.out/args.cohort/'apoe'/panel
    with unit(folder,identity(args,cfg,env,'apoe',panel,0,preflight)|{'data_sha256':frame_sha(data),'permutation_streams_sha256':streams,'seed':seed}):
        records=run_records(folder,tasks,pool)
        observed=records['observed'];tests={}
        for pair_index,(a,b) in enumerate(PAIR_CODES):
            key=f'{a}-{b}';good=[records[f'pair_{pair_index}/{pid}'] for pid in range(budget) if records[f'pair_{pair_index}/{pid}']['status']=='ok']
            complete=len(good)==budget and observed['status']=='ok'
            exceedances=None
            if observed['status']=='ok':
                v=observed['pair_distances'][key]
                exceedances=sum(r['pair_distances'][key]>v+1e-9 or abs(r['pair_distances'][key]-v)<=1e-9 for r in good)
            lo=None if exceedances is None else (1+exceedances)/(budget+1)
            hi=None if exceedances is None else min(1.,(1+exceedances+budget-len(good))/(budget+1))
            p=lo if complete else None
            tests[key]={'groups':[GROUPS[a],GROUPS[b]],'budget':budget,'completed':len(good),'failed':budget-len(good),
                        'exceedances':exceedances,'p':p,'p_bounds':[lo,hi], 'adjusted_p':None if p is None else min(1.,3*p),
                        'reject':None if p is None else 3*(1+exceedances)*20<=budget+1,'threshold':.05/3}
        result={'complete':all(r['status']=='ok' for r in records.values()),'expected_records':1+3*budget,'observed':observed,
                'failed_records':[k for k,r in records.items() if r['status']!='ok'],'tests':tests,'reference':cfg['reference'],
                'permutation_streams_sha256':streams,'raw_fit_records_sha256':file_sha(folder/'fit_records.jsonl')}
        save(folder/'result.json',result)
    return result


def run_bootstrap(args,cfg,frame,env,rep,pool,preflight=False):
    data=encoded(frame,cfg['biomarkers'])
    seed=cfg['cohorts'][args.cohort]['bootstrap_seed']+rep+(PREFLIGHT_OFFSET if preflight else 0)
    data,draw=bootstrap_frame(data,seed,cfg['biomarkers'])
    folder=args.out/args.cohort/'bootstrap'/f'rep_{rep:04d}'
    one=task(data,cfg,'concord','bootstrap',seed,kind='observed',replicate=rep)
    with unit(folder,identity(args,cfg,env,'bootstrap','primary',rep,preflight)|draw):
        records=run_records(folder,[one],pool)
        result={'complete':records['bootstrap']['status']=='ok','draw':draw,'estimate':records['bootstrap'],
                'raw_fit_records_sha256':file_sha(folder/'fit_records.jsonl')}
        save(folder/'result.json',result)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--mode',choices=['preflight','n2','apoe','bootstrap','environment'],required=True)
    ap.add_argument('--cohort',choices=['adni','nacc'],required=True);ap.add_argument('--panel')
    ap.add_argument('--start',type=int,default=0);ap.add_argument('--stop',type=int,default=1)
    ap.add_argument('--workers',type=int,default=1);ap.add_argument('--gate',type=Path)
    args=ap.parse_args();cfg,base=read_campaign(args.config);env=environment()
    if args.mode=='environment':
        save(args.out,{'environment':env,'source_sha256':sources(),'campaign_sha256':file_sha(args.config)})
        return
    if args.workers<1 or args.workers>16:
        raise ValueError('Workers must be in 1..16')
    frame=read_frame(cfg,base,args.cohort)
    if args.mode!='preflight':
        if args.gate is None:
            raise ValueError('Formal execution requires --gate')
        gate=json.loads(args.gate.read_text())
        if gate.get('formal_launch_approved') is not True or gate.get('preflight_passed') is not True:
            raise ValueError('Formal launch and completed preflight approval required')
        if gate['campaign_sha256']!=file_sha(args.config) or gate['source_sha256']!=sources() or gate['environment_versions']!=env['versions'] or gate['kde_ebm_source_sha256']!=env['kde_ebm_source_sha256']:
            raise ValueError('Frozen source, campaign, or numeric environment changed')
        if not 0<=args.start<args.stop<=200:
            raise ValueError('Invalid fixed replicate range')
    if args.mode=='preflight':
        args.out=args.out/'preflight'
    start=time.monotonic()
    pool=ProcessPoolExecutor(max_workers=args.workers,mp_context=get_context('spawn')) if args.workers>1 else None
    try:
        if args.mode=='preflight':
            # No production output is read or written here.
            results=[run_n2(args,cfg,base,frame,env,0,pool,True)]
            results += [run_apoe(args,cfg,frame,env,p,pool,True) for p in cfg['panels']]
            results += [run_bootstrap(args,cfg,frame,env,0,pool,True)]
            if not all(r['complete'] for r in results):
                raise RuntimeError('Preflight fit failure retained in raw records')
            save(args.out/args.cohort/'preflight_passed.json',{'passed':True,'preflight_B':8,'formal_outputs_used':False,
                 'campaign_sha256':file_sha(args.config),'source_sha256':sources(),'environment':env,'seconds':time.monotonic()-start})
        elif args.mode=='apoe':
            if args.panel not in cfg['panels']:
                raise ValueError('Select a frozen APOE panel')
            if not run_apoe(args,cfg,frame,env,args.panel,pool)['complete']:
                raise RuntimeError('APOE failure retained; original B599 denominator preserved')
        else:
            for rep in range(args.start,args.stop):
                result=run_n2(args,cfg,base,frame,env,rep,pool) if args.mode=='n2' else run_bootstrap(args,cfg,frame,env,rep,pool)
                if not result['complete']:
                    raise RuntimeError(f'{args.mode} replicate {rep}: failed fit retained')
                print(json.dumps({'mode':args.mode,'cohort':args.cohort,'replicate':rep,'complete':True}),flush=True)
    finally:
        if pool is not None:
            pool.shutdown(wait=True,cancel_futures=True)


if __name__=='__main__':
    main()
