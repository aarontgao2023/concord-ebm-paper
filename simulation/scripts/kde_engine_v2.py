"""Cross-method arm: the original (likelihood-based) event-based model of Fonteijn et al. (2012) /
Young et al. (2014) as implemented in the UCL `kde_ebm` package (GMM mixture models, greedy ascent
for the maximum-likelihood sequence).

Two variants, mirroring the DEBM arms of the main study:
  kde_gmm         per-group: mixture models fitted on the group's own CN and AD subjects
                  (kde_ebm default, fit_all_gmm_models with y<2), maximum-likelihood ordering
                  from greedy ascent on the group's subjects.
  kde_gmm_pooled  mixture models fitted once on all subjects' CN and AD, per-group ordering
                  from greedy ascent on the group's subjects (the pooled-measurement analogue).

Complete-data cells only: kde_ebm's missing-value imputation raises on 2-D input in this
version, so the arm is run where nothing is missing (the `_CD` cells). The KDE mixture path of
kde_ebm is never called (its awkde dependency is stubbed on the cluster). Greedy ascent is
kde_ebm's own (n_init restarts x n_iter proposals); the MCMC stage is not needed for a point
ordering and is not run. Randomness: numpy's global RNG is seeded from a hash of the data
matrix and labels, so a given (dataset, relabelling) always yields the same fit.
"""
from __future__ import annotations

import hashlib
import os
import warnings
from time import monotonic

import numpy as np

os.environ.setdefault("TQDM_DISABLE", "1")

try:
    from .engine_v2 import EngineConfig, FitResult, validate_orderings
except ImportError:
    from engine_v2 import EngineConfig, FitResult, validate_orderings

KDE_ENGINE_VERSION = "kde-ebm-arm-v1.0"
KDE_VARIANTS = {"kde_gmm": "per_group", "kde_gmm_pooled": "pooled"}
DIAG_CODE = {"CN": 0, "AD": 1, "MCI": 2, 1: 0, 3: 1, 2: 2, "1": 0, "3": 1, "2": 2}


def _seed_from(X, labels):
    h = hashlib.sha256(np.ascontiguousarray(np.nan_to_num(X, nan=-9.0)).tobytes())
    h.update(np.asarray(labels).astype(str).tobytes())
    return int.from_bytes(h.digest()[:4], "little")


def _best_greedy(prob_mat, n_iter, n_init):
    from kde_ebm.mcmc.mcmc import greedy_ascent_creation
    runs = greedy_ascent_creation(prob_mat, n_iter, n_init)
    finals = [runs[k][-1] for k in runs]
    best = max(finals, key=lambda o: o.score)
    return np.asarray(best.ordering, dtype=int), float(best.score), [float(o.score) for o in finals]


def fit_kde_orderings(df, engine_config: EngineConfig | None = None, variant: str = "per_group",
                      greedy_n_init: int = 5, greedy_n_iter: int = 500) -> FitResult:
    config = engine_config or EngineConfig()
    from kde_ebm import mixture_model
    started = monotonic()
    metadata = {"PTID", "Diagnosis", "EXAMDATE", config.group_column}
    names = list(config.biomarker_names) if config.biomarker_names is not None else [
        c for c in df.columns if c not in metadata]
    X_all = df[names].to_numpy(dtype=float)
    if np.isnan(X_all).any():
        return FitResult(None, "error", {"engine_version": KDE_ENGINE_VERSION, "variant": variant,
                                         "error": "ValueError: kde_ebm arm requires complete data"})
    y_all = np.asarray([DIAG_CODE[d] for d in df["Diagnosis"].tolist()], dtype=int)
    g_all = df[config.group_column].to_numpy()
    groups = list(config.group_values)
    if tuple(sorted(np.unique(g_all))) != tuple(groups):
        raise ValueError(f"Expected groups {groups}; found {sorted(np.unique(g_all))}")
    diagnostics = {"engine_version": KDE_ENGINE_VERSION, "variant": variant,
                   "greedy_n_init": greedy_n_init, "greedy_n_iter": greedy_n_iter,
                   "biomarkers": names, "greedy_scores": {}, "quality_flags": []}
    orderings, status = None, "error"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            np.random.seed(_seed_from(X_all, g_all))
            pooled_models = mixture_model.fit_all_gmm_models(X_all, y_all) if variant == "pooled" else None
            orderings = []
            for g in groups:
                rows = g_all == g
                X, y = X_all[rows], y_all[rows]
                models = pooled_models if variant == "pooled" else mixture_model.fit_all_gmm_models(X, y)
                prob_mat = mixture_model.get_prob_mat(X, models)
                if not np.isfinite(prob_mat).all():
                    raise ValueError("non-finite mixture probabilities")
                order, score, scores = _best_greedy(prob_mat, greedy_n_iter, greedy_n_init)
                orderings.append(order)
                diagnostics["greedy_scores"][str(g)] = {"best": score, "restarts": scores}
            orderings = validate_orderings(orderings, len(names), len(groups))
            status = "ok"
        except Exception as exc:  # noqa: BLE001  — recorded, never silently dropped
            diagnostics["error"] = f"{type(exc).__name__}: {exc}"
            orderings = None
    diagnostics["seconds"] = monotonic() - started
    return FitResult(orderings, status, diagnostics)
