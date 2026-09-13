"""Known-component, data-estimated-prior DEBM mechanism diagnostic.

For each group/biomarker, estimate the abnormal mixing probability q using only
that group's nonmissing X and its known Gaussian component densities. No true
ordering, individual stage, diagnosis, diagnosis counts, or stage law enters q
or initialization. This removes the true-order initialization advantage of the
known_group_prior oracle, but still uses unavailable component-density truth.
It is a mechanism arm, not an inference candidate.

The reduced log likelihood is sum log((1-q) + q*exp(log(f1/f0))). Its first and
second derivatives are sum (f1-f0)/mix and -sum ((f1-f0)/mix)**2. Consequently a
score root is a global maximum; a boundary score determines a boundary maximum.
All fits use fixed q bounds [1e-6, 1-1e-6]. Flat or empty likelihoods receive a
declared q=.5 fallback and an identification flag, never an invented estimate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from time import monotonic
import warnings

import numpy as np
from scipy.optimize import brentq

try:
    from .dependency_v2 import verify_pyebm
    from .engine_v2 import (ENGINE_VERSION, EngineConfig, FitResult,
                            _temporary_engine, map_orderings)
    from .oracle_engine_v2 import oracle_posteriors
except ImportError:
    from dependency_v2 import verify_pyebm
    from engine_v2 import (ENGINE_VERSION, EngineConfig, FitResult,
                           _temporary_engine, map_orderings)
    from oracle_engine_v2 import oracle_posteriors


ORACLE_ESTIMATED_PRIOR_VERSION = "oracle-components-estimated-prior-v2.1"
Q_LOWER, Q_UPPER = 1e-6, 1 - 1e-6


@dataclass(frozen=True)
class PriorFit:
    """JSON-serializable scalar prior estimate and identification diagnostics."""
    abnormal_prior: float
    status: str
    identified: bool
    n_observed: int
    n_missing: int
    score_lower: float | None
    score_upper: float | None
    score_at_estimate: float | None
    curvature_at_estimate: float | None
    relative_loglik_at_estimate: float | None
    root_iterations: int = 0

    def to_dict(self) -> dict:
        return asdict(self)


def _component_log_bayes_factors(values, component_parameters):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("Prior estimation requires a one-dimensional measurement vector")
    if np.isinf(values).any():
        raise ValueError("Measurements contain infinities; missing values must be NaN")
    mean0 = float(component_parameters["normal_mean"])
    mean1 = float(component_parameters["abnormal_mean"])
    sd = float(component_parameters["sd"])
    if not np.isfinite([mean0, mean1, sd]).all() or sd <= 0:
        raise ValueError("Component means and positive standard deviation must be finite")
    observed = values[~np.isnan(values)]
    # Equal component SDs, as in design_v2. Algebra avoids subtracting two
    # quadratics/log densities; extreme but finite log ratios are supported.
    log_bayes_factors = (mean1 - mean0) * (observed - (mean1 + mean0) / 2) / sd**2
    if not np.isfinite(log_bayes_factors).all():
        raise FloatingPointError("Nonfinite component log Bayes factors")
    return log_bayes_factors, int(len(values) - len(observed)), mean0 == mean1


def _scaled_components(log_bayes_factors):
    """Divide both densities by their pointwise maximum, avoiding overflow."""
    positive = log_bayes_factors >= 0
    f0 = np.exp(-np.maximum(log_bayes_factors, 0))
    f1 = np.exp(np.minimum(log_bayes_factors, 0))
    difference = np.empty_like(log_bayes_factors)
    # expm1 retains information when the two densities are almost equal.
    difference[positive] = -np.expm1(-log_bayes_factors[positive])
    difference[~positive] = np.expm1(log_bayes_factors[~positive])
    return f0, f1, difference


def mixture_loglik_derivatives(q: float, log_bayes_factors) -> tuple[float, float, float]:
    """Reduced log likelihood, score and curvature at interior q.

    The likelihood differs from the ordinary mixture log likelihood by the
    q-independent sum of normal-component log densities. Missing entries are
    not accepted here; estimate_abnormal_prior removes them explicitly.
    """
    q = float(q)
    ratios = np.asarray(log_bayes_factors, dtype=float)
    if not np.isfinite(q) or not 0 < q < 1:
        raise ValueError("q must be finite and strictly between zero and one")
    if ratios.ndim != 1 or not np.isfinite(ratios).all():
        raise ValueError("Log Bayes factors must be a finite one-dimensional vector")
    f0, f1, difference = _scaled_components(ratios)
    scaled_mixture = (1 - q) * f0 + q * f1
    individual_scores = difference / scaled_mixture
    loglik = np.logaddexp(np.log1p(-q), np.log(q) + ratios).sum()
    return float(loglik), float(individual_scores.sum()), -float(individual_scores @ individual_scores)


def estimate_abnormal_prior(values, component_parameters: dict) -> PriorFit:
    """Globally maximize the fixed-component observed-X mixture likelihood.

    q describes the biomarker's observed measurements within this group. Under
    informative missingness it need not be the whole group's abnormal fraction.
    q=.5 for empty/flat likelihoods is a deterministic fallback, flagged as
    unidentified. Boundary estimates are identified constrained maxima.
    """
    ratios, missing, identical_components = _component_log_bayes_factors(values, component_parameters)
    observed = len(ratios)
    if not observed:
        return PriorFit(.5, "no_observations", False, 0, missing,
                        None, None, None, None, None)
    if np.all(ratios == 0):
        status = "identical_components" if identical_components else "flat_observed_likelihood"
        return PriorFit(.5, status, False, observed, missing, 0., 0., 0., 0., 0.)
    f0, f1, difference = _scaled_components(ratios)

    def score(q):
        return float(np.sum(difference / ((1 - q) * f0 + q * f1)))

    lower_score, upper_score = score(Q_LOWER), score(Q_UPPER)
    iterations = 0
    if lower_score <= 0:
        q, status = Q_LOWER, "lower_boundary"
    elif upper_score >= 0:
        q, status = Q_UPPER, "upper_boundary"
    else:
        q, solver = brentq(score, Q_LOWER, Q_UPPER, xtol=5e-15,
                           rtol=8 * np.finfo(float).eps, full_output=True)
        if not solver.converged:
            raise RuntimeError("Bounded concave prior score did not converge")
        iterations, status = int(solver.iterations), "interior"
    loglik, fitted_score, curvature = mixture_loglik_derivatives(q, ratios)
    if not np.isfinite([q, loglik, fitted_score, curvature]).all() or curvature > 0:
        raise FloatingPointError("Invalid fixed-component prior likelihood diagnostics")
    return PriorFit(float(q), status, True, observed, missing, lower_score,
                    upper_score, fitted_score, curvature, loglik, iterations)


def fit_oracle_estimated_prior_orderings(
        df, truth: dict, engine_config: EngineConfig | None = None) -> FitResult:
    """Run original/repaired Mallows using known densities and X-estimated q.

    The only truth keys accessed are biomarker_names, group_order (label names),
    and component_parameters. In particular no Diagnosis column is required.
    Output indices follow engine_config.biomarker_names when supplied. Callers
    must align the known component parameters with the actual simulated labels;
    do not use group-specific oracle parameters as a permutation inference arm.
    """
    config = engine_config or EngineConfig(mode="repaired")
    provenance = verify_pyebm()
    from pyebm.central_ordering.generalized_mallows import weighted_mallows

    names = list(truth["biomarker_names"])
    group_names = list(truth["group_order"])
    expected_names = list(config.biomarker_names) if config.biomarker_names is not None else names
    if config.expected_events is not None and len(names) != config.expected_events:
        raise ValueError("Truth biomarker count differs from expected_events")
    if (len(names) != len(set(names)) or len(expected_names) != len(set(expected_names))
            or set(names) != set(expected_names)):
        raise ValueError("Truth and requested biomarker names must match uniquely")
    if len(group_names) != len(config.group_values) or len(set(group_names)) != len(group_names):
        raise ValueError("Truth group labels do not match engine group count uniquely")
    if tuple(sorted(df[config.group_column].unique())) != config.group_values:
        raise ValueError("Observed groups do not match engine configuration")
    if not set(names) <= set(df.columns):
        raise ValueError("A truth biomarker is absent from the input data")

    diagnostics = {
        "engine_version": ENGINE_VERSION, "oracle_version": ORACLE_ESTIMATED_PRIOR_VERSION,
        "mode": config.mode, "oracle_prior_mode": "estimated_group_prior",
        "arm": "oracle-components-estimated-group-prior", "inference_eligible": False,
        "provenance": provenance, "biomarkers": expected_names, "posterior_biomarkers": names,
        "truth_order_used_in_prior_and_initialization": False,
        "individual_latent_stage_used": False, "diagnosis_used_in_prior_estimation": False,
        "truth_dx_counts_used": False,
        "prior_conditioning": "group and biomarker observed X only; fixed known component densities",
        "prior_objective": "sum log((1-q)*f_normal(X)+q*f_abnormal(X)) over nonmissing X",
        "prior_bounds": [Q_LOWER, Q_UPPER],
        "unidentified_prior_fallback": .5,
        "posterior_covariance": "univariate component marginal densities, including in CORR designs",
        "missingness_conditioning": "estimate among this biomarker's observed values; no individual DX or missingness-pattern model",
        "abnormal_priors": {}, "normal_mixing_params": {}, "initial_orderings": {},
        "prior_fits": {}, "prior_fit_calls": 0, "prior_root_solver_calls": 0,
        "prior_boundary_fits": 0, "prior_unidentified_fits": 0,
        "quality_flags": [], "optimizer_calls": 0, "optimizer_failures": 0,
    }
    started = monotonic()
    status, orders = "error", None
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        try:
            raw_orders, event_centers = [], []
            with _temporary_engine(config, diagnostics):
                for code, group in zip(config.group_values, group_names):
                    values = df.loc[df[config.group_column] == code, names].to_numpy(dtype=float)
                    components = truth["component_parameters"][group]
                    fits = [estimate_abnormal_prior(values[:, i], components[name])
                            for i, name in enumerate(names)]
                    priors = np.asarray([fit.abnormal_prior for fit in fits])
                    diagnostics["prior_fits"][group] = {name: fit.to_dict() for name, fit in zip(names, fits)}
                    diagnostics["abnormal_priors"][group] = priors.tolist()
                    diagnostics["normal_mixing_params"][group] = (1 - priors).tolist()
                    diagnostics["initial_orderings"][group] = np.argsort(1 - priors).tolist()
                    diagnostics["prior_fit_calls"] += len(fits)
                    for name, fit in zip(names, fits):
                        diagnostics["prior_root_solver_calls"] += int(fit.status == "interior")
                        diagnostics["prior_boundary_fits"] += int(fit.status.endswith("_boundary"))
                        diagnostics["prior_unidentified_fits"] += int(not fit.identified)
                        if not fit.identified:
                            diagnostics["quality_flags"].append(f"prior_{fit.status}:{group}:{name}")
                    posterior = oracle_posteriors(values, components, names, priors)
                    order, centers, _ = weighted_mallows.fitMallows(posterior, 1 - priors)
                    raw_orders.append(order)
                    event_centers.append(centers)
            try:
                orders = map_orderings(raw_orders, names, expected_names, len(config.group_values))
                status = "ok"
            except (ValueError, TypeError, OverflowError) as exc:
                status = "invalid"
                diagnostics["error"] = f"{type(exc).__name__}: {exc}"
            diagnostics["event_centers_finite"] = bool(np.isfinite(event_centers).all())
            if not diagnostics["event_centers_finite"]:
                diagnostics["quality_flags"].append("nonfinite_event_centers")
        except Exception as exc:
            diagnostics["error"] = f"{type(exc).__name__}: {exc}"
    if any(item.get("termination") == "evaluation_limit" for item in diagnostics.get("consensus", [])):
        diagnostics["quality_flags"].append("consensus_evaluation_limit")
    diagnostics["warnings_count"] = len(captured)
    diagnostics["warnings"] = sorted({f"{item.category.__name__}: {item.message}" for item in captured})[:30]
    diagnostics["seconds"] = monotonic() - started
    return FitResult(orders, status, diagnostics)
