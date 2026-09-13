"""Two unchanged DEBM pipelines sharing one identical-input mixture evaluation.

The original engine runs first through engine_v2.fit_orderings. Its complete
mixture outputs are copied before consensus. The repaired engine then executes
the same full pipeline, replaying that mixture output only after exact equality
of every do_mixturemodel input. Consensus, event-center normalization, staging,
model compilation and result validation execute normally for each engine.

There is no persistent cache and no reuse between datasets or permutations.
The cache exists only within one call on one DataFrame, guarded by the same
process-local patch lock as engine_v2. All monkeypatches restore on exceptions,
including caller-designated timeouts. No scientific stopping rule is changed.
"""
from __future__ import annotations

from collections import namedtuple
from copy import copy, deepcopy
from dataclasses import replace
from hashlib import sha256
from time import monotonic
import warnings

import numpy as np
import pandas as pd

try:
    from .dependency_v2 import verify_pyebm
    from .engine_v2 import (EngineConfig, FitResult, _PATCH_LOCK,
                            _OptimizationRecorder, fit_orderings)
except ImportError:
    from dependency_v2 import verify_pyebm
    from engine_v2 import (EngineConfig, FitResult, _PATCH_LOCK,
                           _OptimizationRecorder, fit_orderings)


PAIRED_ENGINE_VERSION = "same-input-mixture-replay-v2.1"


class _AbortPair(BaseException):
    """Pass a caller timeout through the stable wrapper's Exception handler."""
    def __init__(self, cause):
        self.cause = cause


def _snapshot(value):
    """Copy dynamic upstream parameter classes by values, not class identity.

    pyebm mutates attributes of freshly created namedtuple *classes* instead of
    instances. deepcopy(class) returns the original class and is not a snapshot.
    """
    if isinstance(value, type) and hasattr(value, "_fields"):
        fields = {name: _snapshot(member) for name, member in vars(value).items()
                  if not name.startswith("_") and not callable(member)}
        return {"parameter_class_name": value.__name__,
                "parameter_field_names": tuple(value._fields), "parameter_fields": fields}
    if isinstance(value, np.ndarray):
        return value.copy()
    if isinstance(value, pd.Series):
        return value.copy(deep=True)
    if isinstance(value, dict):
        return {key: _snapshot(member) for key, member in value.items()}
    if isinstance(value, (tuple, list)):
        return type(value)(_snapshot(member) for member in value)
    if value is None or isinstance(value, (str, int, float, bool, np.generic)):
        return deepcopy(value)
    raise TypeError(f"Unsupported mixture-input snapshot type: {type(value).__name__}")


def _identical(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, np.ndarray):
        if left.dtype != right.dtype or left.shape != right.shape:
            return False
        return bool(np.array_equal(left, right, equal_nan=True)) if left.dtype.kind in "fc" else bool(
            np.array_equal(left, right))
    if isinstance(left, pd.Series):
        return bool(left.equals(right))
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(_identical(left[key], right[key]) for key in left)
    if isinstance(left, (list, tuple)):
        return len(left) == len(right) and all(_identical(a, b) for a, b in zip(left, right))
    if isinstance(left, (float, np.floating)) and np.isnan(left):
        return bool(np.isnan(right))
    return bool(left == right)


def _clone_mixture_output(output):
    parameters, *arrays = output
    if not isinstance(parameters, type) or tuple(parameters._fields) != ("Control", "Disease", "Mixing"):
        raise TypeError("Unexpected upstream BiomarkerParams representation")
    clone = namedtuple(parameters.__name__, parameters._fields)
    for field in parameters._fields:
        setattr(clone, field, deepcopy(getattr(parameters, field)))
    return (clone, *(deepcopy(array) for array in arrays))


def _numeric_record(value):
    """JSON-safe details; explicit nonfinite values are represented by null."""
    if isinstance(value, np.ndarray):
        return _numeric_record(value.tolist())
    if isinstance(value, (list, tuple)):
        return [_numeric_record(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    raise TypeError(f"Unsupported numeric model output type: {type(value).__name__}")


def _array_digest(value):
    array = np.ascontiguousarray(value)
    digest = sha256()
    digest.update(str((array.dtype.str, array.shape)).encode())
    digest.update(array.tobytes())
    return digest.hexdigest()


def _emit_warnings(records):
    for item in records:
        warnings.warn_explicit(item.message, item.category, item.filename, item.lineno,
                               source=item.source)


def _annotate_results(results, state, started):
    total_seconds = monotonic() - started
    for mode, result in results.items():
        result.diagnostics["paired_engine"] = {
            "version": PAIRED_ENGINE_VERSION, "mode": mode, "mixture_source_mode": "original",
            "same_input_equality_required": True, "reuse_scope": "this function call only",
            "mixture_compute_calls": state["compute_calls"], "mixture_replay_calls": state["replay_calls"],
            "optimizer_records_replayed": mode == "repaired" and bool(state["replay_calls"]),
            "mixture_seconds": state["mixture_seconds"], "total_wall_seconds": total_seconds,
            "cost_accounting": "two statistical results share this one job; count total_wall_seconds once",
        }
        if mode in state["artifacts"]:
            result.diagnostics["model_artifacts"] = state["artifacts"][mode]


def fit_paired_orderings(
        df, base_config: EngineConfig | None = None, *,
        abort_exceptions: tuple[type[BaseException], ...] = (TimeoutError,)) -> dict[str, FitResult]:
    """Return original/repaired FitResults using one same-DF mixture fit.

    All config fields except mode apply equally to both engines. The first mode
    is always original. Standard diagnostics match independent fit_orderings
    calls except elapsed seconds; extra paired/model-artifact fields identify
    shared computation. optimizer records are replayed *diagnostic evidence*,
    not evidence that the repaired pass ran the optimizer a second time.

    The caller owns the joint job's timeout. Include its timeout exception class
    in abort_exceptions if it does not inherit TimeoutError. Such exceptions are
    re-raised after all patches restore, rather than swallowed by fit_orderings;
    its completed_results attribute retains any fully completed original result.
    Ordinary mixture failures are replayed; an engine-specific downstream error
    does not prevent the other engine from running from a successful mixture.
    """
    verify_pyebm()
    from pyebm import core_utilities as cu, debm
    from pyebm.mixture_model import gaussian_mixture_model as gmm

    if not isinstance(abort_exceptions, tuple) or any(
            not isinstance(item, type) or not issubclass(item, BaseException) for item in abort_exceptions):
        raise TypeError("abort_exceptions must be a tuple of exception classes")
    config = base_config or EngineConfig()
    state = {"mode": "original", "input": None, "output": None, "exception": None,
             "warnings": [], "optimizer": [], "compute_calls": 0, "replay_calls": 0,
             "mixture_seconds": 0., "artifacts": {}, "posterior_sha256": {}}
    results = {}
    started = monotonic()
    with _PATCH_LOCK:
        original_mixture, original_fit = cu.do_mixturemodel, debm.fit

        def mixture_call(*args, **kwargs):
            inputs = _snapshot((args, kwargs))
            if not isinstance(gmm.opt, _OptimizationRecorder):
                raise RuntimeError("Paired mixture requires the stable optimizer recorder")
            if state["mode"] == "original":
                if state["compute_calls"]:
                    raise RuntimeError("Paired engine only supports the one-mixture default pipeline")
                state["compute_calls"] += 1
                state["input"] = inputs
                optimizer_start = len(gmm.opt.records)
                mixture_started = monotonic()
                recorded = []
                try:
                    with warnings.catch_warnings(record=True) as recorded:
                        try:
                            output = original_mixture(*args, **kwargs)
                            state["output"] = _clone_mixture_output(output)
                            state["posterior_sha256"]["original"] = _array_digest(output[1])
                        except BaseException as exc:
                            state["exception"] = exc
                            raise
                finally:
                    state["mixture_seconds"] = monotonic() - mixture_started
                    state["optimizer"] = deepcopy(gmm.opt.records[optimizer_start:])
                    state["warnings"] = list(recorded)
                    _emit_warnings(recorded)
                return output

            if state["replay_calls"] or state["input"] is None:
                raise RuntimeError("No single matching mixture fit is available in this paired call")
            if not _identical(state["input"], inputs):
                raise RuntimeError("Mixture inputs differ: replay across changed data or labels is prohibited")
            state["replay_calls"] += 1
            gmm.opt.records.extend(deepcopy(state["optimizer"]))
            _emit_warnings(state["warnings"])
            if state["exception"] is not None:
                raise copy(state["exception"]).with_traceback(None)
            output = _clone_mixture_output(state["output"])
            state["posterior_sha256"]["repaired"] = _array_digest(output[1])
            return output

        def pipeline_call(*args, **kwargs):
            try:
                output = original_fit(*args, **kwargs)
            except abort_exceptions as exc:
                raise _AbortPair(exc) from exc
            model = output[0]
            # Capture upstream model objects while their parameter classes still
            # belong to this pass. Preserve their model-space order explicitly.
            state["artifacts"][state["mode"]] = {
                "model_biomarkers": list(model.BiomarkerList),
                "event_centers": _numeric_record(model.EventCenters),
                "biomarker_parameters": [
                    {field: _numeric_record(getattr(parameters, field))
                     for field in ("Control", "Disease", "Mixing")}
                    for parameters in model.BiomarkerParameters],
                "posterior_sha256": state["posterior_sha256"].get(state["mode"]),
                "nonfinite_encoding": "null; see standard finite/quality diagnostics",
            }
            return output

        try:
            cu.do_mixturemodel, debm.fit = mixture_call, pipeline_call
            for mode in ("original", "repaired"):
                state["mode"] = mode
                results[mode] = fit_orderings(df, replace(config, mode=mode))
        except _AbortPair as abort:
            _annotate_results(results, state, started)
            abort.cause.completed_results = results
            raise abort.cause.with_traceback(None)
        finally:
            cu.do_mixturemodel, debm.fit = original_mixture, original_fit

    _annotate_results(results, state, started)
    return results
