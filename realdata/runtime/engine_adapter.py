"""New-namespace adapter: frozen engines plus an explicit fixed diagnosis target.

DEBM mixtures use the all-group original .01 stopping rule and continued consensus.
Each task clears the mixture cache and fits afresh, so every successful shared-score
record carries its own convergence evidence. No cross-process cache claim is needed.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
from contextlib import contextmanager
from pathlib import Path
import signal
import sys
import time
import traceback
import numpy as np
HERE = Path(__file__).resolve().parent
# The engines are the files in simulation/scripts of this repository, byte-identical to the
# vendor/current copies that ran (their hashes are listed in EXECUTION_GATE.json).
# CONCORD_ENGINE_DIR overrides the location when the runtime is copied elsewhere.
ENGINE_DIR = Path(os.environ.get('CONCORD_ENGINE_DIR', HERE.parents[1] / 'simulation' / 'scripts')).resolve()
sys.path.insert(0, str(ENGINE_DIR))
from engine_v2 import EngineConfig, fit_orderings
import invariant_engine_v2 as invariant
from fast_likelihood_v2 import fast_likelihood_context
from kde_engine_v2 import fit_kde_orderings, _seed_from
from allgroup_engine import allgroup_stopping
from common import reference_check, frame_sha, label_sha, kendall


@contextmanager
def explicit_reference(reference):
    target = reference_check(reference)
    original = invariant.composition_weights
    def weights(diag, gv, mask, ref):
        groups, out, ess = np.unique(gv), {}, {}
        for g in groups:
            d = diag[(gv == g) & mask]
            comp = np.array([(d == c).mean() for c in invariant.CODES])
            if (comp <= 0).any():
                raise ValueError('Fixed-reference diagnosis lacks support in a group')
            w = np.zeros(len(d))
            for j, code in enumerate(invariant.CODES):
                w[d == code] = target[j] / comp[j]
            out[g] = w
            ess[str(g)] = float(w.sum()**2 / np.square(w).sum())
        return groups, out, ess
    invariant.composition_weights = weights
    try:
        yield
    finally:
        invariant.composition_weights = original


def timeout_handler(signum, frame):
    raise TimeoutError('Fit exceeded the fixed task timeout')


def json_evidence(value):
    """Keep nonfinite diagnostic evidence explicit without writing invalid JSON."""
    if isinstance(value, dict):
        return {str(k):json_evidence(v) for k,v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_evidence(v) for v in value]
    if isinstance(value, np.generic):
        value=value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return {'nonfinite_float':repr(value)}
    return value


def fit_task(task):
    data, engine, biomarkers, reference = task['data'], task['engine'], task['biomarkers'], task['reference']
    groups = tuple(sorted(data.APOE.unique().tolist()))
    config = EngineConfig(mode='repaired', expected_events=len(biomarkers), group_column='APOE',
                          group_values=groups, labels=('CN','MCI','AD'), biomarker_names=tuple(biomarkers),
                          audit_original_neighbors=True)
    record = {k:v for k,v in task.items() if k not in ('data',)}
    record.update(input_frame_sha256=frame_sha(data), permuted_labels_sha256=label_sha(data.APOE),
                  group_values=list(groups), group_diagnosis_counts={str(g):data.loc[data.APOE.eq(g),'Diagnosis'].value_counts().reindex(['CN','MCI','AD'],fill_value=0).astype(int).tolist() for g in groups})
    start = time.monotonic()
    np.random.seed(int(task['fit_seed']) % (2**32))
    invariant._CACHE.clear()
    measurement = None
    old_signal = signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(task['fit_timeout_s'])
    try:
        if engine == 'likelihood_ebm':
            record['effective_engine_seed'] = _seed_from(data[list(biomarkers)].to_numpy(float), data.APOE.to_numpy())
            record['rng_policy'] = 'Frozen likelihood engine derives its actual RNG seed from biomarker matrix and group labels; requested fit_seed is not used by this engine'
            result = fit_kde_orderings(data, config, variant='per_group', greedy_n_init=5, greedy_n_iter=500)
            record['measurement_rule'] = 'kde_ebm GMM mixture implementation; separate group fits'
        else:
            record['effective_engine_seed'] = int(task['fit_seed']) % (2**32)
            record['rng_policy'] = 'NumPy global RNG initialized from the task fit_seed before the frozen DEBM engine'
            with allgroup_stopping(detailed=task['kind']=='observed') as measurement, fast_likelihood_context(True):
                if engine == 'separate':
                    result = fit_orderings(data, config)
                elif engine == 'shared':
                    result = invariant.fit_invariant_orderings(data, config, variant='shared', use_cache=False)
                elif engine == 'concord':
                    with explicit_reference(reference):
                        result = invariant.fit_invariant_orderings(data, config, variant='invariant_min', use_cache=False)
                    result.diagnostics.update(variant='explicit_reference', reference_composition=list(reference), reference_fixed_under_permutation=True)
                else:
                    raise ValueError(f'Unknown engine: {engine}')
            record['measurement'] = measurement.summary()
            if result.ok and measurement.termination != 'all_groups_at_upstream_tolerance':
                raise ValueError('Successful DEBM fit lacks a fresh all-group convergence certificate')
            if result.ok:
                certificate = record['measurement']
                final = certificate.get('final_mixing_mean_abs_change_by_group')
                if final is None:
                    final = certificate['trajectory'][-1]['mixing_mean_abs_change_by_group']
                if not np.isfinite(final).all() or (np.asarray(final)<0).any() or np.max(final)>=.01:
                    raise ValueError('Nonfinite or nonconverged final all-group mixing residual')
                if not 1 <= certificate['outer_iterations'] <= 100 or result.diagnostics.get('parameters_finite') is not True:
                    raise ValueError('Invalid all-group iteration/parameter certificate')
        record.update(result.to_dict())
        if result.ok:
            k = len(biomarkers)
            if len(record['orderings']) != len(groups) or any(sorted(o)!=list(range(k)) for o in record['orderings']):
                raise ValueError('Invalid event/group ordering dimensions')
            record['pair_distances'] = {f'{groups[i]}-{groups[j]}':kendall(record['orderings'][i],record['orderings'][j]) for i in range(len(groups)) for j in range(i+1,len(groups))}
    except Exception as exc:
        record.update(status='error', error=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc())
        if measurement is not None:
            record['failed_measurement'] = measurement.summary()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_signal)
        invariant._CACHE.clear()
    record['fit_seconds'] = time.monotonic()-start
    return json_evidence(record)
