"""Oracle component-density arms for simulation mechanism diagnostics only.

Both arms use known univariate normal/abnormal Gaussian components, with missing
measurements retained as NaN. They differ in abnormal-state prior q:

* equal_prior: q=0.5 for every biomarker and group;
* known_group_prior: q is the population abnormal probability conditional on the
  group's realized diagnosis counts and the true group ordering/stage laws.

Known-group priors use true-order information, including through the upstream
initialization argsort(1-q). These are oracle performance diagnostics, not
deployable tests. Neither arm reads individual latent stages or assigns abnormal
states from those stages. No mixture-model optimization is performed.
"""
from __future__ import annotations

from time import monotonic
from typing import Literal
import warnings

import numpy as np
from scipy.special import expit, logit
from scipy.stats import beta

try:
    from .dependency_v2 import verify_pyebm
    from .engine_v2 import (ENGINE_VERSION, EngineConfig, FitResult,
                            _temporary_engine, map_orderings, validate_orderings)
except ImportError:
    from dependency_v2 import verify_pyebm
    from engine_v2 import (ENGINE_VERSION, EngineConfig, FitResult,
                           _temporary_engine, map_orderings, validate_orderings)


ORACLE_VERSION = "oracle-components-v2.1"
PRIOR_MODES = ("equal_prior", "known_group_prior")


def analytic_group_priors(truth: dict) -> dict[str, np.ndarray]:
    """Pr(event abnormal), integrating Beta stages over realized group DX counts.

    For zero-based event rank r and k=round(I*u), Pr(k>r) is
    Beta.sf((r+0.5)/I, a, b). Half-integer rounding ties have probability zero.
    IID diagnosis designs therefore use a prior conditional on the realized
    diagnosis composition, not the unconditional pooled diagnosis probability.
    Priors do not condition on an individual's DX, missingness pattern, or k.
    """
    events = len(truth["biomarker_names"])
    probabilities = {}
    for group in truth["group_order"]:
        order = validate_orderings([truth["group_orderings"][group]], events, 1)[0]
        ranks = np.argsort(order)
        thresholds = (ranks + 0.5) / events
        counts = truth["group_dx_counts"][group]
        if any(not np.isfinite(value) or value < 0 for value in counts.values()):
            raise ValueError("Diagnosis counts must be finite and nonnegative")
        total = sum(counts.values())
        if total <= 0:
            raise ValueError(f"Group {group} has no subjects")
        probability = np.zeros(events)
        for diagnosis, count in counts.items():
            a, b = truth["stage_beta_parameters"][group][diagnosis]
            if not np.isfinite([a, b]).all() or a <= 0 or b <= 0:
                raise ValueError("Beta shapes must be finite and positive")
            probability += count / total * beta.sf(thresholds, a, b)
        if not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1)):
            raise FloatingPointError("Invalid analytic abnormal probabilities")
        probabilities[group] = probability
    return probabilities


def oracle_posteriors(values: np.ndarray, component_parameters: dict,
                      biomarker_names: list[str], abnormal_priors: np.ndarray) -> np.ndarray:
    """Known-component posterior probabilities; q is the abnormal prior.

    pyebm's mixing_params instead represent the normal prior 1-q. For equal
    component standard deviations, log(f_abnormal/f_normal) is linear in x.
    Computing expit(logit(q)+log Bayes factor) avoids density underflow.
    """
    values = np.asarray(values, dtype=float)
    priors = np.asarray(abnormal_priors, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(biomarker_names):
        raise ValueError("Measurements do not match the biomarker names")
    if np.isinf(values).any():
        raise ValueError("Measurements contain infinities; missing values must be NaN")
    if priors.shape != (len(biomarker_names),) or not np.isfinite(priors).all() or np.any(
            (priors < 0) | (priors > 1)):
        raise ValueError("Abnormal priors must be finite probabilities")
    means0 = np.asarray([component_parameters[name]["normal_mean"] for name in biomarker_names])
    means1 = np.asarray([component_parameters[name]["abnormal_mean"] for name in biomarker_names])
    scales = np.asarray([component_parameters[name]["sd"] for name in biomarker_names])
    if not np.isfinite([means0, means1, scales]).all() or np.any(scales <= 0):
        raise ValueError("Component means and positive standard deviations must be finite")
    log_bayes_factor = (means1 - means0) * (values - (means1 + means0) / 2) / scales**2
    return expit(log_bayes_factor + logit(priors))


def fit_oracle_orderings(
        df, truth: dict,
        prior_mode: Literal["equal_prior", "known_group_prior"] = "equal_prior",
        engine_config: EngineConfig | None = None) -> FitResult:
    """Use original/repaired Mallows on oracle posteriors; return FitResult.

    The caller supplies the actual simulated labels and aligned truth. Do not use
    this function as a permutation-test arm: known-group priors and group-specific
    component parameters are not estimable from an arbitrary relabeling's data.
    """
    if prior_mode not in PRIOR_MODES:
        raise ValueError(f"Unknown oracle prior mode: {prior_mode!r}")
    config = engine_config or EngineConfig(mode="repaired")
    provenance = verify_pyebm()
    from pyebm.central_ordering.generalized_mallows import weighted_mallows

    names = list(truth["biomarker_names"])
    expected_names = list(config.biomarker_names) if config.biomarker_names is not None else names
    if config.expected_events is not None and len(names) != config.expected_events:
        raise ValueError("Truth biomarker count differs from expected_events")
    if len(truth["group_order"]) != len(config.group_values):
        raise ValueError("Truth group count differs from engine configuration")
    if tuple(sorted(df[config.group_column].unique())) != config.group_values:
        raise ValueError("Observed groups do not match engine configuration")
    if not set(names) <= set(df.columns):
        raise ValueError("A truth biomarker is absent from the input data")
    priors = (analytic_group_priors(truth) if prior_mode == "known_group_prior"
              else {group: np.full(len(names), 0.5) for group in truth["group_order"]})
    diagnostics = {
        "engine_version": ENGINE_VERSION, "oracle_version": ORACLE_VERSION,
        "mode": config.mode, "oracle_prior_mode": prior_mode,
        "arm": f"oracle-components-{'equal-prior' if prior_mode == 'equal_prior' else 'known-group-prior'}",
        "inference_eligible": False, "provenance": provenance,
        "biomarkers": expected_names, "posterior_biomarkers": names,
        "truth_order_used_in_prior_and_initialization": prior_mode == "known_group_prior",
        "individual_latent_stage_used": False,
        "prior_conditioning": "constant q=0.5" if prior_mode == "equal_prior" else
        "group and realized diagnosis counts; Beta stage laws integrated analytically",
        "posterior_covariance": "univariate component marginal densities, including in CORR designs",
        "missingness_conditioning": "prior does not condition on individual missingness or diagnosis",
        "abnormal_priors": {group: value.tolist() for group, value in priors.items()},
        "normal_mixing_params": {group: (1 - value).tolist() for group, value in priors.items()},
        "initial_orderings": {group: np.argsort(1 - value).tolist() for group, value in priors.items()},
        "quality_flags": [], "optimizer_calls": 0, "optimizer_failures": 0,
    }
    started = monotonic()
    status, orders = "error", None
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always")
        try:
            raw_orders, event_centers = [], []
            with _temporary_engine(config, diagnostics):
                for code, group in zip(config.group_values, truth["group_order"]):
                    rows = df[config.group_column] == code
                    if not rows.any():
                        raise ValueError(f"No subjects in group {group}")
                    for diagnosis, count in truth["group_dx_counts"][group].items():
                        if int((rows & (df["Diagnosis"] == diagnosis)).sum()) != count:
                            raise ValueError("Input group-by-diagnosis counts differ from truth")
                    posterior = oracle_posteriors(df.loc[rows, names].to_numpy(),
                                                   truth["component_parameters"][group],
                                                   names, priors[group])
                    order, centers, _ = weighted_mallows.fitMallows(posterior, 1 - priors[group])
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
    if any(item.get("termination") == "evaluation_limit"
           for item in diagnostics.get("consensus", [])):
        diagnostics["quality_flags"].append("consensus_evaluation_limit")
    diagnostics["warnings_count"] = len(captured)
    diagnostics["warnings"] = sorted({f"{item.category.__name__}: {item.message}"
                                      for item in captured})[:30]
    diagnostics["seconds"] = monotonic() - started
    return FitResult(orders, status, diagnostics)
