"""Resumable simulation runner. Append every fit; derive decisions from fixed-B bounds.

The immutable JSON config specifies a named development/confirmation run. Each
chunk owns one directory and retains failed fits. A stopped or failed permutation
has an unknown outcome in the original planned budget, never a new denominator.
"""
from __future__ import annotations

import os
for _thread_var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_thread_var] = '1'

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import fcntl
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import signal
import sys
import time

import numpy as np

from design_v2 import resolve, simulate
from engine_v2 import EngineConfig, fit_orderings
from dependency_v2 import verify_pyebm
from oracle_engine_v2 import fit_oracle_orderings
from oracle_estimated_prior_v2 import fit_oracle_estimated_prior_orderings
from invariant_engine_v2 import fit_invariant_orderings, VARIANTS as INVARIANT_VARIANTS
from kde_engine_v2 import fit_kde_orderings, KDE_VARIANTS
from fast_likelihood_v2 import fast_likelihood_context, FAST_LIKELIHOOD_VERSION

PAIRS = ((0, 1), (0, 2), (1, 2))
GROUPS = ('e2', 'e33', 'e4')
STOP_REQUESTED = False


def utc():
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def kendall(a, b):
    ra, rb = np.argsort(a), np.argsort(b)
    n = len(a)
    return sum((ra[i]-ra[j])*(rb[i]-rb[j]) < 0
               for i in range(n) for j in range(i+1, n)) / (n*(n-1)/2)


def pair_taus(orders):
    return [float(kendall(orders[a], orders[b])) for a, b in PAIRS]


def decision_bounds(gt, eq, successes, budget, alpha=0.05):
    """Exact decisions or None over every completion of unknown planned fits."""
    if not 0 <= successes <= budget:
        raise ValueError('Invalid successful/planned fit count')
    unknown = budget-successes
    alpha = Fraction(str(alpha))
    result = {'alpha': float(alpha)}
    for rule in ('strict', 'le'):
        threshold = alpha / 3
        compare = (lambda n, d: n*threshold.denominator < d*threshold.numerator) \
            if rule == 'strict' else \
            (lambda n, d: n*threshold.denominator <= d*threshold.numerator)
        answers = []
        for g, e in zip(gt, eq):
            lower = 1+int(g)+int(e)
            upper = lower+unknown
            if compare(upper, budget+1):
                answers.append(True)
            elif not compare(lower, budget+1):
                answers.append(False)
            else:
                answers.append(None)
        result[rule] = {'pair_reject': answers}
    return result


def single_decision(exceedances, successes, budget, alpha, rule):
    threshold = Fraction(str(alpha))
    compare = (lambda n: n*threshold.denominator < (budget+1)*threshold.numerator) \
        if rule == 'strict' else \
        (lambda n: n*threshold.denominator <= (budget+1)*threshold.numerator)
    lower = 1+exceedances
    if compare(lower+budget-successes):
        return True
    if not compare(lower):
        return False
    return None


def compact_diagnostics(d):
    optimizer = d.get('optimizer', [])
    consensus = d.get('consensus', [])
    return {k: d[k] for k in ('quality_flags', 'error', 'elapsed_seconds',
                              'parameters_finite', 'event_centers_finite') if k in d} | {
        'optimizer_calls': len(optimizer),
        'optimizer_unsuccessful': sum(not x['success'] for x in optimizer),
        'optimizer_status_counts': dict(Counter(str(x['status']) for x in optimizer)),
        'consensus': consensus,
    }


class FitTimeout(TimeoutError):
    pass


def timeout_handler(signum, frame):
    raise FitTimeout('Single fit exceeded fixed timeout')


def request_stop(signum, frame):
    global STOP_REQUESTED
    STOP_REQUESTED = True


def permute(df, truth, seed, name, pid, definition):
    tag = int.from_bytes(hashlib.sha256(name.encode()).digest()[:4], 'little')
    rng = np.random.default_rng(np.random.SeedSequence([seed, tag, pid, 61000000]))
    labels = df.APOE.to_numpy().copy()
    stratify = definition.get('stratify', 'diagnosis')
    if stratify == 'none':
        strata = np.zeros(len(df), dtype=int)
    elif stratify == 'diagnosis':
        strata = df.Diagnosis.to_numpy()
    elif stratify == 'oracle_dx_stage':
        # DX is part of the estimator input; conditioning only on k is insufficient.
        strata = np.asarray([f'{dx}:{k}' for dx, k in
                             zip(df.Diagnosis, truth['latent_stage'])])
    elif stratify == 'diagnosis_observed_count':
        # Missingness is part of the fitted input too: condition on Diagnosis and on how many
        # biomarkers a subject has observed (the amount of information, not which ones).
        observed = df[list(truth['biomarker_names'])].notna().sum(axis=1).to_numpy()
        strata = np.asarray([f'{dx}:{k}' for dx, k in zip(df.Diagnosis, observed)])
    else:
        raise ValueError(f'Unknown stratification: {stratify}')
    allowed = definition.get('groups', [0, 1, 2])
    for stratum in np.unique(strata):
        indexes = np.flatnonzero((strata == stratum) & np.isin(labels, allowed))
        labels[indexes] = rng.permutation(labels[indexes])
    result = df.copy()
    result['APOE'] = labels
    return result


def fit_job(job):
    config, cell_spec, seed, engine, scheme, pid = job
    started = time.monotonic()
    base = {'kind': 'observed' if scheme is None else 'permutation', 'seed': seed,
            'engine': engine, 'scheme': scheme, 'perm_id': pid,
            'likelihood_implementation': (FAST_LIKELIHOOD_VERSION if config.get('fast_likelihood', False)
                                          else 'upstream-frozen-normal')
                                          if engine in ('original', 'repaired') or engine in INVARIANT_VARIANTS else 'oracle_not_applicable'}
    truth = {}
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(config.get('fit_timeout_s', 300))
    try:
        design = resolve(cell_spec['name'], **cell_spec.get('overrides', {}))
        df, truth = simulate(design, seed)
        orders = [truth['group_orderings'][g] for g in GROUPS]
        td = pair_taus(orders)
        truth.update(pair_null=[x == 0 for x in td], pair_distances=td,
                     h1_distance=max(td))
        if scheme is not None:
            if engine.startswith('oracle_'):
                raise ValueError('Oracle engines are mechanism diagnostics only')
            df = permute(df, truth, seed, scheme, pid, config['schemes'][scheme])
        mode = 'repaired' if (engine in INVARIANT_VARIANTS or engine in KDE_VARIANTS) else engine.rsplit('_', 1)[-1]
        ec = EngineConfig(mode=mode, audit_original_neighbors=scheme is None,
                          expected_events=len(truth['biomarker_names']),
                          biomarker_names=tuple(truth['biomarker_names']))
        with fast_likelihood_context(config.get('fast_likelihood', False) and
                                     (engine in ('original', 'repaired') or engine in INVARIANT_VARIANTS)):
            if engine in INVARIANT_VARIANTS:
                result = fit_invariant_orderings(df, ec, variant=engine)
            elif engine in KDE_VARIANTS:
                result = fit_kde_orderings(df, ec, variant=KDE_VARIANTS[engine])
            elif engine.startswith('oracle_estimated_'):
                result = fit_oracle_estimated_prior_orderings(df, truth, engine_config=ec)
            elif engine.startswith('oracle_'):
                prior = {'oracle_equal': 'equal_prior',
                         'oracle_known': 'known_group_prior'}[engine.rsplit('_', 1)[0]]
                result = fit_oracle_orderings(df, truth, prior_mode=prior, engine_config=ec)
            else:
                result = fit_orderings(df, ec)
        base.update(result.to_dict())
        base['taus'] = pair_taus(result.orderings) if result.ok else None
        if scheme is None:
            base['truth'] = truth
            base['distance_to_truth'] = [float(kendall(a, b)) for a, b in
                                         zip(result.orderings, orders)] if result.ok else None
        else:
            base['diagnostics'] = compact_diagnostics(base['diagnostics'])
    except Exception as exc:
        base.update(status='error', orderings=None, taus=None,
                    diagnostics={'error': f'{type(exc).__name__}: {exc}'})
        if scheme is None:
            base['truth'] = truth
    finally:
        signal.alarm(0)
    base['fit_seconds'] = time.monotonic()-started
    return base


def fit_paired_job(job):
    """One dataset/permutation, two consensus engines, no cross-label reuse."""
    from paired_engine_v2 import fit_paired_orderings
    config, cell_spec, seed, _, scheme, pid = job
    started = time.monotonic()
    truth = {}
    common = {'kind': 'observed' if scheme is None else 'permutation',
              'seed': seed, 'scheme': scheme, 'perm_id': pid,
              'likelihood_implementation': FAST_LIKELIHOOD_VERSION if config.get('fast_likelihood', False)
                                           else 'upstream-frozen-normal'}
    shared_id = hashlib.sha256(json.dumps([config['run_id'], cell_spec, seed, scheme, pid],
                                          sort_keys=True).encode()).hexdigest()
    records = []
    signal.signal(signal.SIGALRM, timeout_handler)
    total_budget = 2*config.get('fit_timeout_s', 300)
    signal.alarm(total_budget)
    try:
        design = resolve(cell_spec['name'], **cell_spec.get('overrides', {}))
        df, truth = simulate(design, seed)
        orders = [truth['group_orderings'][g] for g in GROUPS]
        td = pair_taus(orders)
        truth.update(pair_null=[x == 0 for x in td], pair_distances=td,
                     h1_distance=max(td))
        if scheme is not None:
            df = permute(df, truth, seed, scheme, pid, config['schemes'][scheme])
        ec = EngineConfig(mode='original', audit_original_neighbors=scheme is None,
                          expected_events=len(truth['biomarker_names']),
                          biomarker_names=tuple(truth['biomarker_names']))
        with fast_likelihood_context(config.get('fast_likelihood', False)):
            results = fit_paired_orderings(df, base_config=ec, abort_exceptions=(FitTimeout,))
        for engine in ('original', 'repaired'):
            result = results[engine]
            record = common | {'engine': engine} | result.to_dict()
            record['taus'] = pair_taus(result.orderings) if result.ok else None
            if scheme is None:
                record['truth'] = truth
                record['distance_to_truth'] = [float(kendall(a, b)) for a, b in
                                                zip(result.orderings, orders)] if result.ok else None
            else:
                record['diagnostics'] = compact_diagnostics(record['diagnostics'])
            records.append(record)
    except Exception as exc:
        completed = getattr(exc, 'completed_results', {})
        records = []
        for engine in ('original', 'repaired'):
            result = completed.get(engine)
            if result is not None:
                record = common | {'engine': engine} | result.to_dict()
                record['taus'] = pair_taus(result.orderings) if result.ok else None
                if scheme is None:
                    record['distance_to_truth'] = [float(kendall(a, b)) for a, b in
                                                  zip(result.orderings, orders)] if result.ok else None
                else:
                    record['diagnostics'] = compact_diagnostics(record['diagnostics'])
            else:
                record = common | {'engine': engine, 'status': 'error', 'orderings': None,
                                   'taus': None, 'diagnostics': {
                                       'error': f'{type(exc).__name__}: {exc}',
                                       'paired_total_timeout_s': total_budget}}
            if scheme is None:
                record['truth'] = truth
            records.append(record)
    finally:
        signal.alarm(0)
    elapsed = time.monotonic()-started
    for record in records:
        record['shared_fit_id'] = shared_id
        record['fit_seconds'] = elapsed
        record['diagnostics'].update(shared_fit_id=shared_id, paired_standard=True,
                                      shared_compute_seconds=elapsed,
                                      paired_total_timeout_s=total_budget)
    return records


def execute_job(job):
    return fit_paired_job(job) if job[3] == '__paired__' else [fit_job(job)]


def collapse_paired_jobs(jobs, config):
    if not config.get('paired_standard', False):
        return jobs
    if not {'original', 'repaired'} <= set(config['engines']):
        raise ValueError('paired_standard requires both original and repaired engines')
    collapsed, seen = [], set()
    for cfg, cell, seed, engine, scheme, pid in jobs:
        if engine not in ('original', 'repaired'):
            collapsed.append((cfg, cell, seed, engine, scheme, pid))
            continue
        task_key = (json.dumps(cfg, sort_keys=True), json.dumps(cell, sort_keys=True),
                    seed, scheme, pid)
        if task_key not in seen:
            collapsed.append((cfg, cell, seed, '__paired__', scheme, pid))
            seen.add(task_key)
    return collapsed


def replay_consistent(old, new):
    # Timeout/other failure records are retained, never replaced by a successful retry.
    if old['status'] == 'ok' and new['status'] == 'ok':
        if any(old.get(field) != new.get(field) for field in ('orderings', 'taus')):
            raise ValueError(f'Deterministic paired replay disagrees with saved result: {key(old)}')
    return {'key': list(key(old)), 'old_status': old['status'], 'new_status': new['status'],
            'policy': 'keep_original_record', 'utc': utc()}


def key(record):
    return (record['seed'], record['engine'], record['scheme'], record['perm_id'])


def source_hashes():
    root = Path(__file__).resolve().parent
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.glob('*.py')) if p.name in {
                'run_v2.py', 'design_v2.py', 'engine_v2.py',
                'dependency_v2.py', 'oracle_engine_v2.py', 'oracle_estimated_prior_v2.py', 'paired_engine_v2.py', 'fast_likelihood_v2.py', 'invariant_engine_v2.py', 'kde_engine_v2.py'}}


def load_events(path):
    """Recover only a torn final append; preserve its bytes for audit/replay."""
    path = Path(path)
    records = {}
    if not path.exists():
        return records
    raw = path.read_bytes()
    lines = raw.splitlines(keepends=True)
    offset = 0
    for i, line in enumerate(lines):
        if not line.strip():
            offset += len(line)
            continue
        try:
            r = json.loads(line)
        except (json.JSONDecodeError, UnicodeDecodeError):
            if i != len(lines)-1 or line.endswith(b'\n'):
                raise ValueError(f'Corrupted fit log before a recoverable torn tail: {path}')
            quarantine = path.with_name(path.name + '.torn-' + hashlib.sha256(line).hexdigest()[:12])
            quarantine.write_bytes(line)
            with path.open('r+b') as stream:
                stream.truncate(offset)
                stream.flush()
                os.fsync(stream.fileno())
            break
        k = key(r)
        if k in records:
            raise ValueError(f'Duplicate fit record: {k}')
        records[k] = r
        offset += len(line)
        if i == len(lines)-1 and not line.endswith(b'\n'):
            with path.open('ab') as stream:
                stream.write(b'\n')
                stream.flush()
                os.fsync(stream.fileno())
    return records


def state_for(records, seed, engine, scheme, budget, alpha):
    observed = records.get((seed, engine, None, None))
    if observed is None or observed['status'] != 'ok':
        return None
    rs = [r for k, r in records.items() if k[:3] == (seed, engine, scheme)]
    good = [r for r in rs if r['status'] == 'ok']
    o = np.asarray(observed['taus'])
    gt = np.zeros(3, dtype=int)
    eq = np.zeros(3, dtype=int)
    mgt = meq = 0
    for r in good:
        t = np.asarray(r['taus'])
        gt += t > o+1e-9
        eq += np.abs(t-o) <= 1e-9
        mgt += int(t.max() > o.max()+1e-9)
        meq += int(abs(t.max()-o.max()) <= 1e-9)
    decisions = decision_bounds(gt, eq, len(good), budget, alpha)
    for rule in ('strict', 'le'):
        decisions[rule]['max_reject'] = single_decision(
            mgt+meq, len(good), budget, alpha, rule)
    return {'requested_nperm': budget, 'nperm': len(good), 'attempted': len(rs),
            'failed': len(rs)-len(good), 'complete': len(good) == budget,
            'gt': gt.tolist(), 'eq': eq.tolist(), 'maxgt': mgt, 'maxeq': meq,
            'exact_decision': decisions, 'perm_ids': {r['perm_id'] for r in rs}}


def selected_full_p(config):
    selected = set(config.get('full_p_seeds', []))
    selected.update(range(config['base_seed'], config['base_seed']+config.get('full_p_first_n', 0)))
    return selected


def resolved_for(state, definition, rule):
    # Pair-only schemes expose their matched endpoint; global max is optional.
    indexes = [definition['tested_pair_index']] if 'tested_pair_index' in definition else range(3)
    decision = state['exact_decision'][rule]
    resolved = all(decision['pair_reject'][i] is not None for i in indexes)
    if definition.get('require_max', False):
        resolved = resolved and decision['max_reject'] is not None
    return resolved


def make_rows(config, cell_spec, seeds, records):
    rows = []
    for seed in seeds:
        for engine in config['engines']:
            observed = records.get((seed, engine, None, None))
            if observed is None:
                continue
            schemes = {}
            for scheme in config.get('schemes', {}):
                state = state_for(records, seed, engine, scheme,
                                  config.get('bperm', 0), config.get('alpha', .05))
                if state is not None:
                    state.pop('perm_ids')
                    state['definition'] = config['schemes'][scheme]
                    schemes[scheme] = state
            rows.append({'run_id': config['run_id'], 'phase': config['phase'],
                         'cell': cell_spec.get('id', cell_spec['name']), 'engine': engine,
                         'seed': seed, 'status': observed['status'],
                         'truth': observed.get('truth', {}), 'schemes': schemes,
                         'observed': {k: observed.get(k) for k in
                                      ('orderings', 'taus', 'distance_to_truth', 'fit_seconds')},
                         'diagnostics': observed['diagnostics'] | {
                             'sequential_stopping': bool(config.get('schemes')),
                             'full_p_selected': seed in selected_full_p(config)}})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--cell-index', type=int, default=0)
    ap.add_argument('--chunk', type=int, default=0)
    ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--time-budget-s', type=float, default=12000)
    args = ap.parse_args()
    config = json.loads(args.config.read_text())
    if config['phase'] not in ('development', 'confirmation', 'diagnostic'):
        raise ValueError('Must declare development/confirmation/diagnostic phase')
    if not 0 <= args.chunk < args.nchunks or args.workers < 1:
        raise ValueError('Invalid chunk/workers')
    cell = config['cells'][args.cell_index]
    seeds = list(range(config['base_seed'], config['base_seed']+config['datasets']))[args.chunk::args.nchunks]
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'runner.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        identity = {'config': config, 'cell_index': args.cell_index, 'chunk': args.chunk,
                    'nchunks': args.nchunks, 'source_sha256': source_hashes()}
        signature = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
        manifest_path = args.output/'manifest.json'
        current_environment = verify_pyebm()
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            if manifest['signature'] != signature:
                raise ValueError('Source/config changed; choose a new run output directory')
            old_environment = manifest['environment']
            for field in ('versions', 'source_sha256', 'wheel_sha256'):
                if old_environment.get(field) != current_environment.get(field):
                    raise ValueError(f'Numerical environment changed ({field}); use a new run')
        else:
            atomic_json(manifest_path, identity | {'signature': signature, 'created_utc': utc(),
                                                   'environment': current_environment, 'seeds': seeds})
        event_path = args.output/'fit_records.jsonl'
        records = load_events(event_path)
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGUSR1):
            signal.signal(sig, request_stop)
        started = time.monotonic()
        full_p = selected_full_p(config)
        rule = config.get('rule', 'strict')
        if rule not in ('strict', 'le'):
            raise ValueError('Unknown rejection rule')

        def stopping():
            stop_file = Path(config.get('stop_file', str(args.output.parent/'STOP')))
            return STOP_REQUESTED or stop_file.exists() or time.monotonic()-started > args.time_budget_s

        def checkpoint(done=False):
            rows = make_rows(config, cell, seeds, records)
            temporary = args.output/'rows.jsonl.tmp'
            with temporary.open('w') as stream:
                for row in rows:
                    stream.write(json.dumps(row, allow_nan=False)+'\n')
            os.replace(temporary, args.output/'rows.jsonl')
            atomic_json(args.output/'progress.json', {'updated_utc': utc(), 'done': done,
                        'stopped': stopping(), 'fit_records': len(records), 'rows': len(rows),
                        'expected_rows': len(seeds)*len(config['engines']),
                        'computational_jobs': len({r.get('shared_fit_id', str(key(r)))
                                                   for r in records.values()}),
                        'statuses': dict(Counter(r['status'] for r in records.values())),
                        'elapsed_seconds': time.monotonic()-started})

        done = False
        with event_path.open('a', buffering=1) as stream, get_context('spawn').Pool(args.workers) as pool:
            while not stopping():
                jobs = []
                for seed in seeds:
                    for engine in config['engines']:
                        if (seed, engine, None, None) not in records:
                            jobs.append((config, cell, seed, engine, None, None))
                if not jobs:
                    for seed in seeds:
                        for engine in config['engines']:
                            for scheme in config.get('schemes', {}):
                                state = state_for(records, seed, engine, scheme,
                                                  config['bperm'], config.get('alpha', .05))
                                if state is None or state['attempted'] >= config['bperm']:
                                    continue
                                resolved = resolved_for(state, config['schemes'][scheme], rule)
                                if resolved and seed not in full_p:
                                    continue
                                missing = [i for i in range(config['bperm']) if i not in state['perm_ids']]
                                jobs.extend((config, cell, seed, engine, scheme, i)
                                            for i in missing[:config.get('permutation_batch', 8)])
                if not jobs:
                    done = True
                    break
                # Bound work before the next stop check; completion order never selects seeds.
                jobs = collapse_paired_jobs(jobs, config)[:max(args.workers*2, 8)]
                for record_batch in pool.imap_unordered(execute_job, jobs, chunksize=1):
                    # Validate every replay before appending any half of this pair.
                    replay_audits = {}
                    batch_keys = [key(record) for record in record_batch]
                    if len(batch_keys) != len(set(batch_keys)):
                        raise ValueError('Worker returned duplicate keys within one batch')
                    for candidate in record_batch:
                        candidate_key = key(candidate)
                        if candidate_key in records:
                            if not config.get('paired_standard', False):
                                raise ValueError(f'Duplicate newly completed fit {candidate_key}')
                            replay_audits[candidate_key] = replay_consistent(records[candidate_key], candidate)
                    for record in record_batch:
                        k = key(record)
                        if k in records:
                            if not config.get('paired_standard', False):
                                raise ValueError(f'Duplicate newly completed fit {k}')
                            audit = replay_audits[k]
                            with (args.output/'replay_audit.jsonl').open('a') as replay_stream:
                                replay_stream.write(json.dumps(audit)+'\n')
                            continue
                        stream.write(json.dumps(record, allow_nan=False)+'\n')
                        stream.flush()
                        os.fsync(stream.fileno())
                        records[k] = record
                        if len(records) % max(1, args.workers) == 0:
                            checkpoint()
                checkpoint()
                print(json.dumps({'utc': utc(), 'cell': cell.get('id', cell['name']),
                                  'chunk': args.chunk, 'fits_saved': len(records)}), flush=True)
        checkpoint(done)
        print(json.dumps({'done': done, 'output': str(args.output), 'fits': len(records)}), flush=True)


if __name__ == '__main__':
    main()
