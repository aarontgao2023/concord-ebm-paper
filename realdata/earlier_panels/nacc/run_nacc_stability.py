#!/usr/bin/env python3
"""NACC amyloid PET population: 500 bootstrap refits (Supplementary Fig. 3, NACC panel).

Each refit resamples participants with replacement within the nine genotype-by-diagnosis cells
(concord.core.resample_within_strata), keeps the covariate adjustment of run_nacc.build_frame
fixed and refits CONCORD from scratch (pooled mixture, weights and continued search). Because the
cell sizes are fixed, the common proportions and weights are fixed too. Every fit record is kept.

  --reference min      CONCORD with the smallest-proportion common proportions (used in the paper)
  --reference pooled   CONCORD with pooled proportions (not reported)

Seed 20260917; R = 500. Before the refits the script checks that resampling keeps every cell
size and that serial and parallel refits agree. With --baseline (an earlier observed CONCORD fit
of this population, cell_concord_<reference>.json) it also checks that the observed ordering is
reproduced. Requires the concord-ebm package (version 0.1.1, commit 40fc046) and pyebm 2.0.3.

  python run_nacc_stability.py --table <derived>/nacc_p5_deid.csv --reference min \
      --out $CONCORD_WORK_DIR/earlier/nacc_stability/min

Outputs (aggregate): observed.json, records/rep_XXXX.json (one refit each; orderings and
diagnostics), bootstrap_orderings.csv, position_probabilities.csv, pair_precedence.csv,
summary.json, validation.json. export_stability_source.py collects them.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
import concord
from concord import core

HERE = Path(__file__).resolve().parent
GROUPS = ('e2', 'e33', 'e4')
R = 500


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, default=core._jsonable) + '\n')
    tmp.replace(path)


def invariant(record):
    return {k: record.get(k) for k in ('kind', 'scheme', 'index', 'status', 'orderings', 'distances')} | {
        'ess': record.get('diagnostics', {}).get('effective_sample_size')}


def run_tasks(frame, spec, tasks, workers):
    runner = core._Runner(frame, spec, {}, workers)
    try:
        yield from runner.map(tasks)
    finally:
        runner.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--reference', choices=('min', 'pooled'), required=True,
                        help='common-proportion rule of CONCORD (min is used in the paper)')
    parser.add_argument('--table', type=Path, required=True, help='nacc_p5_deid.csv from build_cohort_nacc.py')
    parser.add_argument('--nacc-script', type=Path, default=HERE / 'run_nacc.py',
                        help='run_nacc.py that defines build_frame (default: next to this script)')
    parser.add_argument('--baseline', type=Path, help='optional observed fit to reproduce (cell_concord_<reference>.json)')
    parser.add_argument('--plan', type=Path, help='optional protocol document whose hash is recorded')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    TABLE, SOURCE = args.table, args.nacc_script
    baseline_path = args.baseline
    baseline = json.loads(baseline_path.read_text()) if baseline_path else None
    sys.path.insert(0, str(SOURCE.parent))
    module_spec = importlib.util.spec_from_file_location('frozen_nacc_preprocessing', SOURCE)
    module = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(module)
    raw = pd.read_csv(TABLE)
    adjusted, correction = module.build_frame(raw)
    del raw
    frame, groups, events = core.prepare_data(adjusted, group_column='APOE',
        labels=module.LABELS, biomarkers=module.EVENTS, group_order=GROUPS)
    assert len(frame) == 1200 and len(events) == 5
    spec = core.CompareSpec(group_column='APOE', group_values=groups, labels=tuple(module.LABELS),
        biomarkers=events, estimator='invariant_' + args.reference, consensus='repaired',
        fast_likelihood=True, fit_timeout_s=900, seed=20260917)
    core.verify_pyebm()
    counts = core._composition(frame, spec)['counts']
    if baseline is not None:
        assert counts == baseline['composition']['counts']
    manifest = {'analysis': 'NACC diagnosis-by-genotype conditional resampling stability',
        'created_utc': datetime.now(timezone.utc).isoformat(), 'reference': args.reference,
        'R': R, 'sample_indices': [0, R-1], 'spec': asdict(spec),
        'fixed_covariate_correction': True, 'fixed_diagnostic_counts': counts,
        'reference_recomputed_but_counts_unchanged': True,
        'refit_pooled_measurement_model_each_resample': True,
        'input_sha256': sha(TABLE), 'preprocessing_sha256': sha(SOURCE),
        'baseline_sha256': sha(baseline_path) if baseline_path else None, 'runner_sha256': sha(__file__),
        'plan_sha256': sha(args.plan) if args.plan else None, 'workers': args.workers,
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'pandas': pd.__version__, 'concord': concord.__version__},
        'package_sources': {p.name: sha(p) for p in Path(concord.__file__).parent.glob('*.py')}}
    mp = args.out / 'manifest.json'
    if mp.exists():
        previous = json.loads(mp.read_text())
        for field in ('reference','R','spec','input_sha256','preprocessing_sha256','baseline_sha256',
                      'runner_sha256','plan_sha256','environment','package_sources'):
            assert previous[field] == json.loads(json.dumps(manifest[field], default=core._jsonable)), field
    else:
        write(mp, manifest)
    observed = next(run_tasks(frame, spec, [('observed','observed',0)], 1))
    assert observed['status'] == 'ok', observed
    observed_names = {g: [events[i] for i in observed['orderings'][gi]] for gi,g in enumerate(groups)}
    if baseline is not None:
        assert observed_names == baseline['orderings'], 'Observed point order must reproduce the baseline fit'
        for i,pair in enumerate(spec.pairs):
            assert np.isclose(observed['distances'][i], baseline['distances']['pairs'][spec.pair_name(pair)], atol=1e-12)
        for gi,g in enumerate(groups):
            assert np.isclose(observed['diagnostics']['effective_sample_size'][str(gi)], baseline['effective_sample_size'][g])
    write(args.out/'observed.json', observed)
    # Native resampling must leave every genotype-by-diagnosis cell size unchanged.
    for index in (0,1,499):
        sample = core.resample_within_strata(frame, spec, index)
        assert core._composition(sample, spec)['counts'] == counts
        assert len(sample) == len(frame)
    tasks = [('stability','stability',i) for i in range(2)]
    serial = sorted(run_tasks(frame, spec, tasks, 1), key=lambda r:r['index'])
    parallel = sorted(run_tasks(frame, spec, tasks, args.workers), key=lambda r:r['index'])
    assert all(r['status']=='ok' for r in serial + parallel)
    assert [invariant(r) for r in serial] == [invariant(r) for r in parallel]
    write(args.out/'preflight_validation.json', {'passed':True, 'observed_baseline_equal':True if baseline is not None else 'not checked',
        'serial_parallel_indices':[0,1], 'serial':[invariant(r) for r in serial],
        'parallel':[invariant(r) for r in parallel], 'cell_counts_preserved_indices':[0,1,499]})
    for record in parallel:
        path = args.out/'records'/f"rep_{record['index']:04d}.json"
        if path.exists():
            assert invariant(json.loads(path.read_text())) == invariant(record)
        else:
            write(path, record)
    todo=[]
    for i in range(R):
        path=args.out/'records'/f'rep_{i:04d}.json'
        if not path.exists(): todo.append(('stability','stability',i))
    completed=R-len(todo)
    for record in run_tasks(frame,spec,todo,args.workers):
        write(args.out/'records'/f"rep_{record['index']:04d}.json",record)
        completed+=1
        if completed%20==0: print(f'{args.reference}: {completed}/{R} saved',flush=True)
    records=[json.loads((args.out/'records'/f'rep_{i:04d}.json').read_text()) for i in range(R)]
    failed=[r['index'] for r in records if r['status']!='ok']
    write(args.out/'fit_completion.json',{'requested':R,'completed':R-len(failed),'failed_indices':failed})
    if failed: raise RuntimeError(f'Unresolved original resamples retained: {failed}')
    orders=np.asarray([r['orderings'] for r in records],dtype=int)
    assert orders.shape==(R,3,5)
    assert all(sorted(row.tolist())==list(range(5)) for row in orders.reshape(-1,5))
    positions=np.argsort(orders,axis=2)
    pos_rows=[];order_rows=[];precedence=[]
    for gi,g in enumerate(groups):
        for ei,event in enumerate(events):
            for pos in range(5):
                count=int(np.sum(positions[:,gi,ei]==pos))
                pos_rows.append(dict(cohort='NACC',reference=args.reference,group=g,event=event,
                                     position=pos+1,count=count,R=R,probability=count/R))
            for r in range(R):
                order_rows.append(dict(reference=args.reference,replicate=r,group=g,event=event,
                                       position=int(positions[r,gi,ei])+1))
        for a in range(5):
            for b in range(a+1,5):
                count=int(np.sum(positions[:,gi,a]<positions[:,gi,b]))
                precedence.append(dict(reference=args.reference,group=g,event_a=events[a],event_b=events[b],
                                       a_before_b=count,R=R,frequency=count/R))
    pd.DataFrame(pos_rows).to_csv(args.out/'position_probabilities.csv',index=False)
    pd.DataFrame(order_rows).to_csv(args.out/'bootstrap_orderings.csv',index=False)
    pd.DataFrame(precedence).to_csv(args.out/'pair_precedence.csv',index=False)
    summary=core._stability_summary(records, observed['orderings'], groups, events, R)
    centiloid=events.index('CENTILOID')
    summary['centiloid']={g:{'first':float(np.mean(positions[:,gi,centiloid]==0)),
                             'last':float(np.mean(positions[:,gi,centiloid]==4)),
                             'median_position':float(np.median(positions[:,gi,centiloid]+1))}
                         for gi,g in enumerate(groups)}
    summary['observed_orderings']=observed_names
    write(args.out/'summary.json',summary)
    write(args.out/'validation.json',{'passed':True,'R':R,'complete':R,'failed':0,
        'observed_baseline_equal':True if baseline is not None else 'not checked','serial_parallel_equal':True,
        'position_probability_sums_equal_one':bool(np.allclose(pd.DataFrame(pos_rows).groupby(['group','event']).probability.sum(),1)),
        'outputs_sha256':{p.name:sha(p) for p in args.out.glob('*.csv')}})
    print(json.dumps(summary['centiloid']),flush=True)


if __name__=='__main__':
    main()
