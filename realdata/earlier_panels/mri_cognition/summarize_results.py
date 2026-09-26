#!/usr/bin/env python3
"""Read-only production status and aggregate results for the frozen campaign.

Usage:
  python summarize_results.py --root /path/to/production --config campaign.json
  python summarize_results.py --root /path/to/production --json-out summary.json

JSON goes to stdout unless --json-out is given; a short report goes to stderr.
No inference package is imported and no participant-level table is parsed.
Frozen files are hashed; NPZ draw indices are used only to verify draw hashes.
Incomplete rates always retain denominator 200. B=8 preflight is never accepted.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

R, B, ALPHA = 200, 599, 0.05
COHORTS = ('adni', 'nacc')
PANELS = ('MRIcog4', 'MRIcog8')
ESTIMATORS = ('standard', 'invariant_min')
APOE_PAIRS = ('e2-e33', 'e2-e4', 'e33-e4')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def near(a, b):
    return isinstance(a, (int, float)) and not isinstance(a, bool) and math.isfinite(a) and math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-12)


def campaign_context(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text())
    errors = []
    for field, expected in [('R', R), ('B', B), ('alpha', ALPHA), ('consensus', 'repaired'),
                            ('stop_when_decided', False), ('stability_resamples', 0),
                            ('actual_APOE_panel', 'MRIcog8')]:
        if config.get(field) != expected:
            errors.append(f'campaign.{field} differs from production contract')
    if set(config.get('cohorts', {})) != set(COHORTS) or set(config.get('panels', {})) != set(PANELS):
        errors.append('campaign cohort/panel set differs from production contract')
    files = []
    for rel, expected in config.get('frozen_files', {}).items():
        p = path.parent / rel
        status = 'missing' if not p.is_file() else ('ok' if sha(p) == expected else 'hash_mismatch')
        files.append({'path': rel, 'status': status})
        if status != 'ok':
            errors.append(f'frozen file {status}: {rel}')
    if not files:
        errors.append('campaign has no frozen file manifest')
    draw_hashes = {}
    for cohort in COHORTS:
        c = config.get('cohorts', {}).get(cohort, {})
        try:
            # These are row indices, never exported or printed.
            with np.load(path.parent / c['draw_file'], allow_pickle=False) as plans:
                a, b = plans['a'], plans['b']
                if len(a) != R or len(b) != R:
                    raise ValueError('draw array does not have 200 rows')
                if a.shape[1] != sum(c['n2_counts']['a']) or b.shape[1] != sum(c['n2_counts']['b']):
                    raise ValueError('draw sizes differ from configured counts')
                for rep in range(R):
                    both = np.r_[a[rep], b[rep]]
                    if len(set(map(int, both))) != len(both) or both.min() < 0 or both.max() >= c['n']:
                        raise ValueError('draw indices are duplicated or out of range')
                    draw_hashes[(cohort, rep)] = hashlib.sha256(both.astype('<i8').tobytes() + c['adjusted_sha256'].encode()).hexdigest()
            adjusted = path.parent / c['adjusted_file']
            if not adjusted.is_file() or sha(adjusted) != c['adjusted_sha256']:
                errors.append(f'{cohort}: adjusted input SHA mismatch')
        except (OSError, ValueError, KeyError, IndexError) as exc:
            errors.append(f'{cohort}: draw manifest unavailable or invalid ({type(exc).__name__})')
    return {'config': config, 'config_path': path, 'campaign_sha256': sha(path),
            'errors': errors, 'frozen_files': files, 'draw_hashes': draw_hashes}


def quality(result):
    d = result.get('observed_diagnostics') or {}
    consensus = d.get('consensus') or []
    return {
        'quality_flags': d.get('quality_flags', []),
        'composition_flags': (result.get('composition') or {}).get('flags', []),
        'parameters_finite': d.get('parameters_finite'),
        'event_centers_finite': d.get('event_centers_finite'),
        'optimizer_calls': d.get('optimizer_calls'), 'optimizer_failures': d.get('optimizer_failures'),
        'warnings_count': d.get('warnings_count'),
        'consensus_local_optimal': [x.get('local_optimal') for x in consensus],
        'consensus_termination': [x.get('termination') for x in consensus],
    }


def validate(value, ctx, cohort, mode, panel, estimator, rep):
    errors = []
    def check(condition, message):
        if not condition:
            errors.append(message)
    try:
        config = ctx['config']; c = config['cohorts'][cohort]
        a, r = value['analysis'], value['result']
        seed = c['apoe_seed'] if mode == 'apoe' else c['n2_fit_seed_base'] + rep
        expected_n = c['n'] if mode == 'apoe' else sum(c['n2_counts']['a']) + sum(c['n2_counts']['b'])
        draw = c['adjusted_sha256'] if mode == 'apoe' else ctx['draw_hashes'].get((cohort, rep))
        for field, expected in dict(cohort=cohort, mode=mode, panel=panel, estimator=estimator,
                replicate=rep, B=B, seed=seed, n=expected_n, draw_sha256=draw,
                campaign_sha256=ctx['campaign_sha256']).items():
            check(a.get(field) == expected, f'analysis.{field} mismatch')
        check(draw is not None, 'expected draw hash unavailable')
        check(not ctx['errors'], 'campaign frozen-file validation failed')
        expected_tests = ({'diagnosis_pair_' + p for p in APOE_PAIRS} if mode == 'apoe' else
                          ({'unrestricted','diagnosis'} if estimator == 'standard' else {'diagnosis'}))
        groups = ['e2','e33','e4'] if mode == 'apoe' else ['pA','pB']
        check(r.get('status') == 'ok', 'result status is not ok')
        check(r.get('groups') == groups, 'result groups mismatch')
        check(r.get('biomarkers') == config['panels'][panel], 'event panel mismatch')
        check(set(r.get('tests', {})) == expected_tests, 'test/operator set mismatch')
        spec = r.get('spec', {}); dec = r.get('decisions', {})
        for field, expected in [('B',B),('estimator',estimator),('consensus','repaired'),('seed',seed),
                                ('stability_resamples',0),('labels',c['labels']),('group_values',groups),
                                ('biomarkers',config['panels'][panel])]:
            check(spec.get(field) == expected, f'spec.{field} mismatch')
        check(dec.get('stop_when_decided') is False, 'early stopping was enabled or unspecified')
        check(dec.get('rule') == 'le' and near(dec.get('alpha'), ALPHA), 'decision rule/alpha mismatch')
        fits = r.get('fits', {})
        check(fits.get('error') == 0 and fits.get('timeout') == 0, 'failed/timeout fit present')
        check(fits.get('ok') == 1 + len(expected_tests) * B, 'full expected fit count not reached')
        threshold = ALPHA / (3 if mode == 'apoe' else 1)
        for name in expected_tests:
            t = r.get('tests', {}).get(name, {})
            check(t.get('budget') == B and t.get('completed') == B and t.get('failed') == 0,
                  f'{name}: not a complete failure-free B599 test')
            pair = name.removeprefix('diagnosis_pair_') if mode == 'apoe' else 'pA-pB'
            check(set(t.get('pairs', {})) == {pair}, f'{name}: pair set mismatch')
            cell = t.get('pairs', {}).get(pair, {})
            ex = cell.get('exceedances')
            valid_ex = isinstance(ex,int) and not isinstance(ex,bool) and 0 <= ex <= B
            check(valid_ex, f'{name}: invalid inclusive exceedance count')
            if valid_ex:
                p = (ex + 1) / (B + 1)
                bounds = cell.get('p_bounds', [])
                check(near(cell.get('p'), p), f'{name}: plus-one p mismatch')
                check(len(bounds) == 2 and all(near(x,p) for x in bounds), f'{name}: p bounds not point-valued')
                check(near(cell.get('threshold'), threshold), f'{name}: rejection threshold mismatch')
                check(cell.get('reject') is (p <= threshold), f'{name}: rejection decision mismatch')
            d = r.get('distances', {}).get('pairs', {}).get(pair)
            check(isinstance(d,(int,float)) and math.isfinite(d) and 0 <= d <= 1,
                  f'{name}: missing/nonfinite/out-of-range order distance')
        ess = r.get('effective_sample_size', {})
        check(set(ess) == set(groups), 'ESS group set mismatch')
        check(all(isinstance(x,(int,float)) and math.isfinite(x) and x > 0 for x in ess.values()), 'invalid ESS')
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        errors.append(f'malformed result schema ({type(exc).__name__})')
    return errors


def read_unit(root, ctx, cohort, mode, panel, estimator, rep):
    directory = root / cohort / mode / panel / estimator
    path = directory / f'rep_{rep:04d}.json'
    failed = sorted(directory.glob(f'rep_{rep:04d}_FAILED_*.json'))
    unit = {'cohort':cohort,'mode':mode,'panel':panel,'estimator':estimator,'replicate':rep,
            'path':str(path.relative_to(root)), 'failed_attempt_files':len(failed), 'errors':[]}
    if not path.is_file():
        unit['status'] = 'failed' if failed else 'missing'
        return unit
    try:
        value = json.loads(path.read_text())
    except (OSError,ValueError) as exc:
        unit.update(status='invalid',errors=[f'canonical JSON unreadable ({type(exc).__name__})'])
        return unit
    unit['reported_draw_sha256'] = value.get('analysis', {}).get('draw_sha256')
    unit['errors'] = validate(value,ctx,cohort,mode,panel,estimator,rep)
    unit['status'] = 'invalid' if unit['errors'] else 'complete'
    if unit['status'] == 'complete':
        unit['result'] = value['result']
    return unit


def summarize(root, config_path, issue_limit=30):
    root = Path(root).resolve(); ctx = campaign_context(config_path)
    units = [read_unit(root,ctx,c,'n2',p,e,rep) for c in COHORTS for p in PANELS
             for e in ESTIMATORS for rep in range(R)]
    actual = [read_unit(root,ctx,c,'apoe','MRIcog8','invariant_min',0) for c in COHORTS]
    statuses = ('complete','missing','failed','invalid')
    arms = []
    for c in COHORTS:
        for p in PANELS:
            for e,scheme in [('standard','unrestricted'),('standard','diagnosis'),('invariant_min','diagnosis')]:
                u = [x for x in units if (x['cohort'],x['panel'],x['estimator']) == (c,p,e)]
                counts = {s:sum(x['status']==s for x in u) for s in statuses}
                successes = [x for x in u if x['status']=='complete']
                rejected = sum(x['result']['tests'][scheme]['pairs']['pA-pB']['reject'] for x in successes)
                unknown = R-counts['complete']
                flags = Counter(str(f) for x in successes for f in quality(x['result'])['quality_flags'])
                arms.append({'cohort':c,'panel':p,'method':'Standard' if e=='standard' else 'CONCORD',
                    'estimator':e,'operator':'U-all' if scheme=='unrestricted' else 'D-all','search':'repaired',
                    'planned_R':R, **counts,'unknown':unknown,'reject_count_validated':rejected,
                    'reject_rate_lower':rejected/R,'reject_rate_upper':(rejected+unknown)/R,
                    'final_rejection_rate':rejected/R if unknown==0 else None,
                    'rate_status':'final' if unknown==0 else 'incomplete_bounds_only',
                    'failed_attempt_files':sum(x['failed_attempt_files'] for x in u),
                    'completed_files_with_quality_flags':sum(bool(quality(x['result'])['quality_flags']) for x in successes),
                    'quality_flag_counts':dict(flags)})
    paired = {'expected_draw_blocks':len(COHORTS)*R,'all_four_files_complete':0,'partially_present':0,
              'no_canonical_files_present':0,'hash_conflict_blocks':0,'hash_conflict_details':[]}
    for c in COHORTS:
        for rep in range(R):
            block=[x for x in units if x['cohort']==c and x['replicate']==rep]
            hashes={x['reported_draw_sha256'] for x in block if x.get('reported_draw_sha256')}
            n=sum(x['status']=='complete' for x in block)
            if n==4: paired['all_four_files_complete']+=1
            elif hashes: paired['partially_present']+=1
            else: paired['no_canonical_files_present']+=1
            if len(hashes)>1:
                paired['hash_conflict_blocks']+=1
                if len(paired['hash_conflict_details']) < issue_limit:
                    paired['hash_conflict_details'].append({'cohort':c,'replicate':rep,'distinct_hashes':len(hashes)})
    apoe=[]
    for u in actual:
        row={k:u[k] for k in ['cohort','panel','estimator','status','path','failed_attempt_files','errors']}
        row['operator']='D-pair';row['search']='repaired';row['pairs']=[]
        if u['status']=='complete':
            r=u['result']; row.update(n=ctx['config']['cohorts'][u['cohort']]['n'],
                effective_sample_size=r['effective_sample_size'],quality=quality(r),
                fits=r['fits'],composition=r.get('composition'),orderings=r.get('orderings'))
            for pair in APOE_PAIRS:
                t=r['tests']['diagnosis_pair_'+pair]; x=t['pairs'][pair]
                row['pairs'].append({'pair':pair,'distance':r['distances']['pairs'][pair],
                    'B':B,'completed':t['completed'],'inclusive_exceedances':x['exceedances'],
                    'p_raw':x['p'],'p_bonferroni':min(1.0,3*x['p']),'reject':x['reject'],
                    'raw_threshold':x['threshold']})
        apoe.append(row)
    counts={s:sum(u['status']==s for u in units+actual) for s in statuses}
    expected={u['path'] for u in units+actual}
    extras=[]
    # Limit discovery to production cohort directories; do not recursively ingest
    # serial/parallel preflight directories from a parent results folder.
    for c in COHORTS:
        for p in (root/c).glob('**/rep_*.json'):
            if re.fullmatch(r'rep_\d+\.json',p.name) and str(p.relative_to(root)) not in expected:
                extras.append(str(p.relative_to(root)))
    issues=[{k:u[k] for k in ['cohort','mode','panel','estimator','replicate','path','status','errors']}
            for u in units+actual if u['status'] in ('invalid','failed')]
    complete = counts['complete']==1602 and not ctx['errors'] and paired['hash_conflict_blocks']==0
    return {'campaign':ctx['config'].get('name'), 'campaign_sha256':ctx['campaign_sha256'],
        'config_path':str(ctx['config_path']),'results_root':str(root),'production_complete':complete,
        'fixed_contract':{'R_per_N2_arm':R,'B_per_test':B,'N2_expected_files':1600,
            'N2_arms':12,'actual_APOE_expected_files':2,'actual_APOE_pairs_per_file':3,
            'unknown_denominator_policy':'missing, failed and invalid remain in denominator 200',
            'rate_bounds_meaning':'worst-case completion bounds, not statistical confidence intervals',
            'B8_preflight_accepted':False},
        'campaign_integrity':{'valid':not ctx['errors'],'errors':ctx['errors'],'frozen_files':ctx['frozen_files'],
            'expected_N2_draw_hashes_recomputed':len(ctx['draw_hashes'])},
        'file_status':{'expected':1602,**counts,'historical_FAILED_attempt_files':sum(u['failed_attempt_files'] for u in units+actual),
            'recovered_complete_files_with_FAILED_history':sum(u['status']=='complete' and u['failed_attempt_files']>0 for u in units+actual)},
        'N2':arms,'paired_draw_validation':paired,'actual_APOE':apoe,
        'issue_count':len(issues),'issues':issues[:issue_limit],'issues_truncated':len(issues)>issue_limit,
        'unexpected_canonical_file_count':len(extras),'unexpected_canonical_files':extras[:issue_limit],
        'participant_records_exported':False,'models_run':False}


def short_text(summary):
    f=summary['file_status']
    lines=[f"Campaign {summary['campaign']}: {'COMPLETE' if summary['production_complete'] else 'INCOMPLETE'}",
        f"Files {f['complete']}/{f['expected']} complete; {f['missing']} missing; {f['failed']} FAILED-only; {f['invalid']} invalid.",
        f"Historical FAILED files: {f['historical_FAILED_attempt_files']}; recovered canonical outputs: {f['recovered_complete_files_with_FAILED_history']}.",
        f"Campaign/input integrity: {'PASS' if summary['campaign_integrity']['valid'] else 'FAIL'}; paired-draw hash conflicts: {summary['paired_draw_validation']['hash_conflict_blocks']}."]
    for a in summary['N2']:
        rate=(f"final {a['final_rejection_rate']:.3f}" if a['final_rejection_rate'] is not None else
              f"bounds [{a['reject_rate_lower']:.3f}, {a['reject_rate_upper']:.3f}]")
        lines.append(f"N2 {a['cohort']} {a['panel']} {a['method']} {a['operator']}: {a['complete']}/200 complete; rejects {a['reject_count_validated']}/200; {rate}; unknown={a['unknown']}.")
    for a in summary['actual_APOE']:
        if a['status']!='complete':
            lines.append(f"APOE {a['cohort']} MRIcog8: {a['status']}.")
        else:
            vals='; '.join(f"{x['pair']} d={x['distance']:.4f}, p={x['p_raw']:.6f}, adjusted={x['p_bonferroni']:.6f}" for x in a['pairs'])
            lines.append(f"APOE {a['cohort']} MRIcog8: {vals}; flags={a['quality']['quality_flags']}.")
    if summary['issue_count']:
        lines.append(f"Validation/FAILED issues: {summary['issue_count']} (details in JSON).")
    if summary['unexpected_canonical_file_count']:
        lines.append(f"Unexpected canonical files excluded: {summary['unexpected_canonical_file_count']}.")
    return '\n'.join(lines)+'\n'


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True,help='Production output directory containing adni/ and nacc/.')
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('campaign.json'))
    parser.add_argument('--json-out',type=Path,help='Optional aggregate JSON path; otherwise stdout.')
    parser.add_argument('--text-out',type=Path,help='Optional short-text path; text is always also sent to stderr.')
    parser.add_argument('--issue-limit',type=int,default=30)
    parser.add_argument('--require-complete',action='store_true',help='Exit 1 while production is incomplete; default is 0 for a readable status report.')
    args=parser.parse_args(argv)
    if args.issue_limit < 0:
        parser.error('--issue-limit must be nonnegative')
    try:
        result=summarize(args.root,args.config,args.issue_limit)
    except (OSError,ValueError,KeyError) as exc:
        print(f'Unable to summarize campaign: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'; report=short_text(result)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True,exist_ok=True);args.json_out.write_text(encoded)
    else:
        sys.stdout.write(encoded)
    sys.stderr.write(report)
    if args.text_out:
        args.text_out.parent.mkdir(parents=True,exist_ok=True);args.text_out.write_text(report)
    return int(args.require_complete and not result['production_complete'])


if __name__=='__main__':
    raise SystemExit(main())
