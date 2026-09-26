#!/usr/bin/env python3
"""MRI-cognition populations: aggregate the stored results (Fig. 4e, four MRI-cognition rows).

No model, participant bootstrap, new draw or label permutation is run. The only resampling is a
paired bootstrap of the 200 stored draws for the paired differences. Outputs contain aggregate
per-draw results, never participants. n2_rejection_rates.csv holds the Fig. 4e rates: rows
Standard_U (separately fitted DEBM with continued search, unrestricted permutation) and
CONCORD_D (CONCORD, within-diagnosis permutation); Standard_D is not reported. Wilson 95%
intervals use z = NormalDist().inv_cdf(0.975).

All 1,602 result files must validate first (1,600 draws from mode n2 and the two APOE fits
from mode apoe of run_shared_panel.py; summarize_results.py holds the validator).

  python aggregate_completed_results.py --root $CONCORD_WORK_DIR/mri_cognition/results \
    --config campaign.json --out $CONCORD_WORK_DIR/earlier/mri_cognition_summary

Exit 1 if any required output is missing/failed/invalid; no summary is written.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
import sys

import numpy as np

# Operational validator only: this imports no CONCORD/model implementation.
import summarize_results as validator

R = 200
B = 599
BOOTSTRAP_R = 10000
SEED_ROOT = 20260920
COHORTS = ('adni', 'nacc')
PANELS = ('MRIcog4', 'MRIcog8')
ESTIMATORS = ('standard', 'invariant_min')
ARMS = (('Standard_U', 'standard', 'unrestricted'),
        ('Standard_D', 'standard', 'diagnosis'),
        ('CONCORD_D', 'invariant_min', 'diagnosis'))
METHOD = {'standard': 'Standard', 'invariant_min': 'CONCORD'}


class IncompleteCampaign(RuntimeError):
    def __init__(self, gate):
        self.gate = gate
        super().__init__('All 1602 required full-budget outputs must validate before aggregation.')


def json_text(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def endpoint_spec(family, **parts):
    endpoint = json_text({'family': family, **parts})
    digest = hashlib.sha256(endpoint.encode('utf-8')).digest()
    words = [int.from_bytes(digest[i:i + 4], 'little') for i in range(0, 32, 4)]
    return {'endpoint': endpoint, 'endpoint_sha256': digest.hex(),
            'seed_entropy': [SEED_ROOT] + words}


def paired_bootstrap(left, right, endpoint):
    """Percentile CI of mean(left-right), retaining paired outer draws."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if left.shape != (R,) or right.shape != (R,) or not np.isfinite(left).all() or not np.isfinite(right).all():
        raise ValueError('A paired endpoint requires exactly 200 finite observations per side.')
    delta = left - right
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(endpoint['seed_entropy'])))
    means = np.empty(BOOTSTRAP_R, dtype=float)
    for start in range(0, BOOTSTRAP_R, 500):
        stop = min(start + 500, BOOTSTRAP_R)
        index = rng.integers(0, R, size=(stop - start, R), dtype=np.int64)
        means[start:stop] = delta[index].mean(axis=1)
    low, high = np.quantile(means, [0.025, 0.975], method='linear')
    return {'R': R, 'mean_left': float(left.mean()), 'mean_right': float(right.mean()),
            'mean_difference': float(delta.mean()), 'sd_paired_difference': float(delta.std(ddof=1)),
            'ci95_low': float(low), 'ci95_high': float(high), 'bootstrap_resamples': BOOTSTRAP_R,
            'endpoint': endpoint['endpoint'], 'endpoint_sha256': endpoint['endpoint_sha256']}


def wilson95(rejects, n=R):
    if not isinstance(rejects, (int, np.integer)) or not 0 <= rejects <= n or n <= 0:
        raise ValueError('Invalid binomial count.')
    z = NormalDist().inv_cdf(0.975)
    p = rejects / n
    divisor = 1 + z*z/n
    center = (p + z*z/(2*n)) / divisor
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / divisor
    return max(0.0, center-half), min(1.0, center+half)


def expected_keys():
    for c in COHORTS:
        for p in PANELS:
            for e in ESTIMATORS:
                for rep in range(R):
                    yield (c, 'n2', p, e, rep)
        yield (c, 'apoe', 'MRIcog8', 'invariant_min', 0)


def require_complete(root, config_path):
    """Validate and hash exact bytes once; do not summarize partial endpoints."""
    root = Path(root).resolve()
    ctx = validator.campaign_context(config_path)
    records, sources, problems = {}, [], []
    status = Counter()
    failed_history = 0
    for key in expected_keys():
        c, mode, panel, est, rep = key
        relative = Path(c) / mode / panel / est / f'rep_{rep:04d}.json'
        path = root / relative
        history = list(path.parent.glob(f'rep_{rep:04d}_FAILED_*.json'))
        failed_history += len(history)
        errors = []
        state = 'complete'
        if not path.is_file():
            state = 'failed' if history else 'missing'
        else:
            try:
                blob = path.read_bytes()
                value = json.loads(blob)
                errors = validator.validate(value, ctx, c, mode, panel, est, rep)
                if errors:
                    state = 'invalid'
            except (OSError, ValueError, TypeError) as exc:
                state = 'invalid'; errors = [f'Unreadable result ({type(exc).__name__})']
        status[state] += 1
        if state != 'complete':
            if len(problems) < 30:
                problems.append({'path': str(relative), 'status': state, 'errors': errors})
            continue
        source = {'cohort': c, 'mode': mode, 'panel': panel, 'estimator': est, 'replicate': rep,
                  'path': str(relative), 'sha256': hashlib.sha256(blob).hexdigest(),
                  'bytes': len(blob), 'failed_attempt_files': len(history)}
        records[key] = {'analysis': value['analysis'], 'result': value['result'], 'source': source}
        sources.append(source)
    paired_ok = 0
    if len(records) == 1602:
        for c in COHORTS:
            for rep in range(R):
                hashes = {records[(c, 'n2', p, e, rep)]['analysis']['draw_sha256']
                          for p in PANELS for e in ESTIMATORS}
                if hashes != {ctx['draw_hashes'][(c, rep)]}:
                    problems.append({'cohort': c, 'replicate': rep, 'status': 'paired_hash_mismatch'})
                else:
                    paired_ok += 1
    gate = {'expected_files': 1602, **{s: status[s] for s in ['complete', 'missing', 'failed', 'invalid']},
            'campaign_integrity_valid': not ctx['errors'], 'campaign_errors': ctx['errors'],
            'paired_draw_blocks_verified': paired_ok, 'expected_paired_draw_blocks': 400,
            'historical_FAILED_attempt_files': failed_history, 'problems_first30': problems[:30]}
    if status['complete'] != 1602 or ctx['errors'] or paired_ok != 400 or problems:
        raise IncompleteCampaign(gate)
    # Detect a changing config during collection, without touching source files.
    if validator.sha(config_path) != ctx['campaign_sha256']:
        raise IncompleteCampaign(gate | {'config_changed_while_reading': True})
    return records, sources, ctx, gate


def scientific_tables(records, config):
    """All prespecified endpoints; no ranking, selection or quality exclusion."""
    if set(records) != set(expected_keys()):
        raise ValueError('Scientific aggregation requires the exact 1602-file task inventory.')
    estimation, inference, means, reductions, rates, rate_differences = [], [], [], [], [], []
    actual_pairs, actual_groups, actual_detail, quality_rows, endpoints = [], [], [], [], []
    distances, decisions = {}, {}
    for c in COHORTS:
        for panel in PANELS:
            for est in ESTIMATORS:
                distance = []
                ess_a, ess_b = [], []
                for rep in range(R):
                    record = records[(c, 'n2', panel, est, rep)]
                    a, result, source = record['analysis'], record['result'], record['source']
                    d = float(result['distances']['pairs']['pA-pB'])
                    q = validator.quality(result)
                    distance.append(d)
                    ess_a.append(float(result['effective_sample_size']['pA']))
                    ess_b.append(float(result['effective_sample_size']['pB']))
                    estimation.append({'cohort': c, 'panel': panel, 'K': len(config['panels'][panel]),
                        'replicate': rep, 'method': METHOD[est], 'estimator': est, 'search': 'repaired',
                        'n_drawn': a['n'], 'n_pA': sum(config['cohorts'][c]['n2_counts']['a']),
                        'n_pB': sum(config['cohorts'][c]['n2_counts']['b']), 'distance': d,
                        'ESS_pA': result['effective_sample_size']['pA'], 'ESS_pB': result['effective_sample_size']['pB'],
                        'ordering_pA_json': json_text(result['orderings']['pA']),
                        'ordering_pB_json': json_text(result['orderings']['pB']),
                        'draw_sha256': a['draw_sha256'], 'source_file': source['path'], 'source_sha256': source['sha256'],
                        'quality_flags_json': json_text(q['quality_flags']), 'quality_json': json_text(q)})
                    quality_rows.append({'cohort': c, 'mode': 'n2', 'panel': panel, 'estimator': est,
                        'replicate': rep, 'quality_flags_json': json_text(q['quality_flags']),
                        'quality_json': json_text(q), 'source_file': source['path']})
                    for arm, arm_est, scheme in ARMS:
                        if arm_est != est:
                            continue
                        cell = result['tests'][scheme]['pairs']['pA-pB']
                        decisions.setdefault((c, panel, arm), []).append(int(cell['reject']))
                        inference.append({'cohort': c, 'panel': panel, 'K': len(config['panels'][panel]),
                            'replicate': rep, 'arm': arm, 'method': METHOD[est], 'estimator': est,
                            'operator': 'U-all' if scheme == 'unrestricted' else 'D-all', 'search': 'repaired',
                            'B': B, 'completed_permutations': B, 'inclusive_exceedances': cell['exceedances'],
                            'p_raw': cell['p'], 'threshold': cell['threshold'], 'reject': int(cell['reject']),
                            'distance': d, 'draw_sha256': a['draw_sha256'],
                            'source_file': source['path'], 'source_sha256': source['sha256']})
                vector = np.asarray(distance, dtype=float)
                distances[(c, panel, est)] = vector
                means.append({'cohort': c, 'panel': panel, 'K': len(config['panels'][panel]),
                    'method': METHOD[est], 'estimator': est, 'search': 'repaired', 'R': R,
                    'mean_distance': float(vector.mean()), 'sd_distance': float(vector.std(ddof=1)),
                    'mean_ESS_pA': float(np.mean(ess_a)), 'sd_ESS_pA': float(np.std(ess_a,ddof=1)),
                    'mean_ESS_pB': float(np.mean(ess_b)), 'sd_ESS_pB': float(np.std(ess_b,ddof=1)),
                    'standard_observed_order_counting': 'one per draw shared by U-all and D-all' if est == 'standard' else 'one per draw'})
            ep = endpoint_spec('paired_distance_reduction', cohort=c, panel=panel, contrast='Standard_minus_CONCORD')
            endpoints.append(ep)
            reductions.append({'cohort': c, 'panel': panel, 'K': len(config['panels'][panel]),
                'contrast': 'Standard_minus_CONCORD', 'positive_means': 'smaller CONCORD distance',
                **paired_bootstrap(distances[(c, panel, 'standard')], distances[(c, panel, 'invariant_min')], ep)})
            for arm, est, scheme in ARMS:
                vector = np.asarray(decisions[(c, panel, arm)], dtype=float)
                if len(vector) != R:
                    raise ValueError('A rejection endpoint lost its 200-draw denominator.')
                count = int(vector.sum()); low, high = wilson95(count)
                rates.append({'cohort': c, 'panel': panel, 'K': len(config['panels'][panel]),
                    'arm': arm, 'method': METHOD[est], 'operator': 'U-all' if scheme == 'unrestricted' else 'D-all',
                    'search': 'repaired', 'R': R, 'rejects': count, 'nonrejects': R-count,
                    'rejection_rate': count/R, 'wilson95_low': low, 'wilson95_high': high, 'B': B})
            ep = endpoint_spec('paired_rate_difference', cohort=c, panel=panel, contrast='Standard_D_minus_CONCORD_D')
            endpoints.append(ep)
            left, right = decisions[(c, panel, 'Standard_D')], decisions[(c, panel, 'CONCORD_D')]
            rate_differences.append({'cohort': c, 'comparison': 'method_at_fixed_panel', 'panel': panel,
                'arm': 'D-all', 'contrast': 'Standard_D_minus_CONCORD_D',
                'positive_means': 'higher Standard D-all rejection rate',
                **discordance(left, right), **paired_bootstrap(left, right, ep)})
        for arm, _, _ in ARMS:
            ep = endpoint_spec('paired_rate_difference', cohort=c, arm=arm, contrast='MRIcog8_minus_MRIcog4')
            endpoints.append(ep)
            left, right = decisions[(c, 'MRIcog8', arm)], decisions[(c, 'MRIcog4', arm)]
            rate_differences.append({'cohort': c, 'comparison': 'panel_within_arm', 'panel': 'MRIcog8_minus_MRIcog4',
                'arm': arm, 'contrast': 'MRIcog8_minus_MRIcog4',
                'positive_means': 'higher K8 rejection rate',
                **discordance(left, right), **paired_bootstrap(left, right, ep)})
        record = records[(c, 'apoe', 'MRIcog8', 'invariant_min', 0)]
        result, source = record['result'], record['source']
        q = validator.quality(result)
        ref = result['composition']['min_reference']
        detail = {'cohort': c, 'panel': 'MRIcog8', 'n': config['cohorts'][c]['n'],
            'operator': 'D-pair', 'estimator': 'invariant_min', 'search': 'repaired',
            'reference': 'sparsest', 'reference_weights': ref, 'composition': result['composition'],
            'ESS': result['effective_sample_size'], 'orderings': result['orderings'],
            'positions': result.get('positions'), 'quality': q,
            'observed_diagnostics': result.get('observed_diagnostics'), 'fits': result['fits'],
            'source_file': source['path'], 'source_sha256': source['sha256'], 'pairs': []}
        for group in ('e2', 'e33', 'e4'):
            actual_groups.append({'cohort': c, 'panel': 'MRIcog8', 'group': group,
                'n': sum(config['cohorts'][c]['counts'][group].values()),
                'diagnosis_counts_json': json_text(config['cohorts'][c]['counts'][group]),
                'reference': 'sparsest', 'reference_weights_json': json_text(ref),
                'ESS': result['effective_sample_size'][group],
                'ordering_json': json_text(result['orderings'][group]),
                'quality_json': json_text(q), 'source_file': source['path'], 'source_sha256': source['sha256']})
        for pair in validator.APOE_PAIRS:
            cell = result['tests']['diagnosis_pair_' + pair]['pairs'][pair]
            row = {'cohort': c, 'panel': 'MRIcog8', 'K': 8, 'n': config['cohorts'][c]['n'],
                'pair': pair, 'estimator': 'invariant_min', 'operator': 'D-pair', 'search': 'repaired',
                'reference': 'sparsest', 'reference_weights_json': json_text(ref), 'B': B,
                'inclusive_exceedances': cell['exceedances'], 'p_raw': cell['p'],
                'p_bonferroni': min(1., 3*cell['p']), 'raw_threshold': cell['threshold'], 'reject': int(cell['reject']),
                'distance': result['distances']['pairs'][pair],
                'ESS_all_groups_json': json_text(result['effective_sample_size']),
                'orderings_all_groups_json': json_text(result['orderings']),
                'quality_flags_json': json_text(q['quality_flags']), 'quality_json': json_text(q),
                'source_file': source['path'], 'source_sha256': source['sha256']}
            actual_pairs.append(row); detail['pairs'].append(row)
        actual_detail.append(detail)
        quality_rows.append({'cohort': c, 'mode': 'apoe', 'panel': 'MRIcog8', 'estimator': 'invariant_min',
            'replicate': 0, 'quality_flags_json': json_text(q['quality_flags']),
            'quality_json': json_text(q), 'source_file': source['path']})
    tables = {'n2_estimation_per_draw': estimation, 'n2_inference_per_draw': inference,
        'n2_distance_summary': means, 'n2_paired_distance_reduction': reductions,
        'n2_rejection_rates': rates, 'n2_paired_rate_differences': rate_differences,
        'actual_APOE_pairs': actual_pairs, 'actual_APOE_groups': actual_groups, 'quality_by_result': quality_rows}
    expected_rows = [1600,2400,8,4,12,10,6,6,1602]
    if [len(v) for v in tables.values()] != expected_rows or len(endpoints) != 14:
        raise AssertionError('Prespecified output inventory is incomplete.')
    return tables, actual_detail, endpoints


def discordance(left, right):
    left, right = np.asarray(left, dtype=int), np.asarray(right, dtype=int)
    return {'left_rejects': int(left.sum()), 'right_rejects': int(right.sum()),
            'both_reject': int(((left==1)&(right==1)).sum()),
            'left_only_reject': int(((left==1)&(right==0)).sum()),
            'right_only_reject': int(((left==0)&(right==1)).sum()),
            'neither_reject': int(((left==0)&(right==0)).sum())}


def write_json(path, value):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    temp.replace(path)


def aggregate(root, config_path, out, protocol_path):
    records, sources, ctx, gate = require_complete(root, config_path)
    tables, actual_detail, endpoints = scientific_tables(records, ctx['config'])
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=True)
    outputs = []
    for name, rows in tables.items():
        p = out / (name+'.csv'); write_csv(p, rows); outputs.append(p)
    p = out/'source_files.csv'; write_csv(p,sources); outputs.append(p)
    for name, value in [('actual_APOE_detail',actual_detail),('completion_gate',gate),
                        ('summary', {'tables':tables,'bootstrap_endpoints':endpoints})]:
        p=out/(name+'.json');write_json(p,value);outputs.append(p)
    notes = ('# Complete shared-panel aggregation\n\n'
        'All 1,602 canonical results passed full-budget, frozen-input and paired-draw validation before summary. '
        'The output retains every prespecified endpoint, including unfavourable or null contrasts.\n\n'
        'The separately fitted DEBM (Standard) rows with unrestricted and within-diagnosis permutation share a single '
        'observed ordering per draw; estimation summaries therefore count it once. SD uses n−1. '
        'Each inference cell has denominator 200 and a two-sided Wilson 95% interval. '
        'Paired differences subtract the right condition from the left as explicitly named in each row.\n\n'
        'All 14 paired intervals use 10,000 resamples of the same 200 outer-draw index pairs, sampled with replacement. '
        'The interval is the linear-interpolated 2.5th and 97.5th percentile of bootstrap means. '
        'These are pointwise draw/Monte Carlo uncertainty intervals conditional on the prepared empirical pool, '
        'not simultaneous intervals or participant-population uncertainty intervals. No cross-cohort pairing is used.\n\n'
        'Endpoint streams use the exact canonical JSON string recorded in manifest.json: SHA256 UTF-8 digest, '
        'split into eight unsigned little-endian 32-bit integers, appended to seed root 20260920, then '
        'NumPy SeedSequence and PCG64. No Python hash(), outcome-dependent seed or stream selection is used.\n\n'
        'K4/K8 fix people and preparation but also change event identities and information; the contrast is not '
        'a pure causal event-count experiment. N2 rates concern the constructed conditional null. '
        'Actual APOE results retain all three pairs, their sparsest reference, ESS, orders and quality diagnostics.\n\n'
        'Per-draw exports contain aggregate fit/inference values and draw hashes, never participant records, '
        'membership indices or biomarker rows. No scientific model was fitted during this aggregation.\n')
    p=out/'README.md';p.write_text(notes);outputs.append(p)
    protocol_path=Path(protocol_path) if protocol_path else None
    manifest={'campaign':ctx['config']['name'],'campaign_sha256':ctx['campaign_sha256'],
        'config_path':str(Path(config_path).resolve()),'source_root':str(Path(root).resolve()),
        'script_sha256':validator.sha(__file__),'validator_sha256':validator.sha(validator.__file__),
        'protocol':{'path':protocol_path.name,'sha256':validator.sha(protocol_path)} if protocol_path and protocol_path.is_file() else None,
        'complete':True,'gate':gate,'source_file_count':len(sources),'source_file_hashes':'source_files.csv',
        'output_rows':{k:len(v) for k,v in tables.items()},'frozen_files':ctx['config']['frozen_files'],
        'panels':ctx['config']['panels'],'cohort_config':ctx['config']['cohorts'],
        'SD':'sample standard deviation, ddof=1','Wilson':'two-sided score interval, z=NormalDist().inv_cdf(0.975)',
        'paired_bootstrap':{'resamples':BOOTSTRAP_R,'draw_count':R,'seed_root':SEED_ROOT,
            'sample_unit':'paired outer-draw indices 0..199, with replacement, within cohort only',
            'statistic':'mean(left-right); orientation explicitly named in endpoint row',
            'CI':'percentile, NumPy quantile [0.025,0.975], method=linear',
            'generator':'numpy.random.Generator(PCG64(SeedSequence(seed_entropy)))',
            'endpoint_encoding':'canonical JSON sorted keys, compact separators, UTF-8; SHA256 split into 8 little-endian uint32 words',
            'seed_entropy':'[20260920] + eight endpoint digest words','batch_size':500,
            'pointwise_not_simultaneous':True,'conditional_on_prepared_empirical_pool':True,'endpoints':endpoints},
        'quality_policy':'all validated finite outputs retained; flags reported, not used to exclude endpoints',
        'manual_sign_selection':False,'new_subject_resamples':False,'new_label_permutations':False,'new_model_fits':False,
        'participant_records_exported':False,'software':{'python':sys.version.split()[0],'numpy':np.__version__},
        'outputs':{p.name:{'sha256':validator.sha(p),'bytes':p.stat().st_size} for p in outputs}}
    write_json(out/'manifest.json',manifest)
    return manifest


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True,help='Full production result directory containing adni/ and nacc/.')
    parser.add_argument('--config',type=Path,default=Path(__file__).with_name('campaign.json'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--protocol',type=Path,default=None,help='optional protocol document whose hash is recorded')
    args=parser.parse_args(argv)
    try:
        manifest=aggregate(args.root,args.config,args.out,args.protocol)
    except IncompleteCampaign as exc:
        print(json.dumps({'aggregation':'blocked_incomplete','gate':exc.gate},indent=2),file=sys.stderr)
        return 1
    except (OSError,ValueError,KeyError,AssertionError) as exc:
        print(f'Aggregation failed: {type(exc).__name__}: {exc}',file=sys.stderr)
        return 2
    print(json.dumps({'aggregation':'complete','source_files':1602,'output_rows':manifest['output_rows'],
                      'paired_bootstrap_endpoints':14,'out':str(args.out.resolve()),'models_run':0}))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
