"""Development benchmark: exact objective, full trajectory and model comparison.

Examples (the pinned pyebm environment must already be available):
  python scripts/v2/benchmark_fast_likelihood_v2.py --micro-calls 300
  python scripts/v2/benchmark_fast_likelihood_v2.py --full-fit --events 4 14 \
      --modes original repaired --output /tmp/fast_likelihood_bridge.json

Output is JSON. It is a development equivalence/throughput check, not a
confirmation simulation. No code or existing run is changed by this command.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
from functools import wraps
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
import warnings

import numpy as np
import scipy.stats

from design_v2 import DesignConfig, resolve, simulate
from engine_v2 import EngineConfig, _PATCH_LOCK, fit_orderings
from fast_likelihood_v2 import (FAST_LIKELIHOOD_VERSION, calculate_likelihood_gmm,
                                fast_likelihood_context)


def assert_exact(left, right, path="root"):
    """Exact numeric equality (including NaN), not an allclose tolerance."""
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        a, b = np.asarray(left), np.asarray(right)
        if a.dtype != b.dtype or a.shape != b.shape or not np.array_equal(a, b, equal_nan=True):
            raise AssertionError(f"Array mismatch at {path}")
    elif isinstance(left, dict) and isinstance(right, dict):
        if left.keys() != right.keys():
            raise AssertionError(f"Dictionary keys differ at {path}")
        for key in left:
            assert_exact(left[key], right[key], f"{path}.{key}")
    elif isinstance(left, (tuple, list)) and isinstance(right, (tuple, list)):
        if len(left) != len(right):
            raise AssertionError(f"Sequence lengths differ at {path}")
        for index, (a, b) in enumerate(zip(left, right)):
            assert_exact(a, b, f"{path}[{index}]")
    elif isinstance(left, (float, np.floating)) and isinstance(right, (float, np.floating)):
        if left != right and not (np.isnan(left) and np.isnan(right)):
            raise AssertionError(f"Numeric mismatch at {path}: {left!r} != {right!r}")
    elif left != right:
        raise AssertionError(f"Value mismatch at {path}: {left!r} != {right!r}")


class _TrajectoryRecorder:
    def __init__(self, original):
        self.original = original
        self.calls = []

    def __getattr__(self, name):
        return getattr(self.original, name)

    def minimize(self, fun, x0, *args, **kwargs):
        call = {"objective_name": fun.__name__, "x0": np.asarray(x0).copy(),
                "method": kwargs.get("method"), "options": deepcopy(kwargs.get("options")),
                "bounds": deepcopy(kwargs.get("bounds")), "evaluations": []}
        self.calls.append(call)

        @wraps(fun)
        def objective(x, *fargs, **fkwargs):
            value = fun(x, *fargs, **fkwargs)
            call["evaluations"].append((np.asarray(x).copy(), deepcopy(value)))
            return value

        result = self.original.minimize(objective, x0, *args, **kwargs)
        call["result"] = {name: deepcopy(getattr(result, name)) for name in
                          ("x", "fun", "jac", "nit", "nfev", "njev", "success", "status", "message")
                          if hasattr(result, name)}
        return result


@contextmanager
def capture_pipeline():
    """Capture all SLSQP evaluations and raw outputs without changing numerics."""
    from pyebm import core_utilities as cu
    from pyebm.mixture_model import gaussian_mixture_model as gmm
    with _PATCH_LOCK:
        previous_opt, previous_mixture, previous_compile = gmm.opt, cu.do_mixturemodel, cu.compile_model_output
        recorder = _TrajectoryRecorder(previous_opt)
        captured = {"optimizer": recorder.calls, "mixture_outputs": [], "models": []}

        def mixture(*args, **kwargs):
            output = previous_mixture(*args, **kwargs)
            captured["mixture_outputs"].append({
                "parameters": {name: deepcopy(getattr(output[0], name)) for name in ("Control", "Disease", "Mixing")},
                "probabilities_and_likelihoods": deepcopy(output[1:]),
            })
            return output

        def compile_model(*args, **kwargs):
            model = previous_compile(*args, **kwargs)
            captured["models"].append({"biomarkers": list(model.BiomarkerList),
                                       "orderings": deepcopy(model.MeanCentralOrdering),
                                       "event_centers": deepcopy(model.EventCenters)})
            return model

        try:
            gmm.opt, cu.do_mixturemodel, cu.compile_model_output = recorder, mixture, compile_model
            yield captured
        finally:
            gmm.opt, cu.do_mixturemodel, cu.compile_model_output = previous_opt, previous_mixture, previous_compile


def full_fit_equivalence(frame, config):
    results, captures, times = {}, {}, {}
    for enabled, label in ((False, "frozen"), (True, "direct")):
        began = perf_counter()
        with fast_likelihood_context(enabled), capture_pipeline() as captured:
            results[label] = fit_orderings(frame, config)
        captures[label] = captured
        times[label] = perf_counter() - began
    for label, result in results.items():
        if not result.ok:
            raise AssertionError(f"{label} full fit failed: {result.diagnostics}")
    for label in results:
        captures[label]["fit_result"] = results[label].to_dict()
        captures[label]["fit_result"]["diagnostics"] = {
            key: value for key, value in results[label].diagnostics.items() if key != "seconds"}
    assert_exact(captures["frozen"], captures["direct"])
    calls = captures["direct"]["optimizer"]
    trajectory_hash = sha256()
    for call in calls:
        for x, value in call["evaluations"]:
            trajectory_hash.update(x.tobytes())
            trajectory_hash.update(np.asarray(value).tobytes())
    return {"events": len(results["direct"].orderings[0]), "subjects": len(frame),
            "mode": config.mode, "full_result_exact_equal": True,
            "all_objective_evaluations_exact_equal": True, "all_optimizer_results_exact_equal": True,
            "all_mixture_outputs_exact_equal": True, "event_centers_exact_equal": True,
            "optimizer_calls": len(calls), "objective_evaluations": sum(len(call["evaluations"]) for call in calls),
            "trajectory_sha256": trajectory_hash.hexdigest(),
            "frozen_seconds": times["frozen"], "direct_seconds": times["direct"],
            "saving_fraction": 1 - times["direct"] / times["frozen"],
            "speedup": times["frozen"] / times["direct"],
            "versions": results["direct"].diagnostics["provenance"]["versions"]}


def microbenchmark(calls=300, seed=831771):
    from pyebm.mixture_model import gaussian_mixture_model as gmm
    rng = np.random.default_rng(seed)
    data = rng.normal(size=971)
    data[rng.random(len(data)) < .2] = np.nan
    groups = rng.integers(0, 3, size=len(data))
    parameters = [np.array([rng.normal(), np.exp(rng.uniform(-3, 1)), rng.normal()+2,
                            np.exp(rng.uniform(-3, 1)), rng.uniform()]) for _ in range(calls)]
    rows = []
    for label, mixing in (("single", []), ("grouped", [.3, .5, .7])):
        args = (data, ["APOE"], [groups], mixing)
        # Check the scalar objective before timing either implementation.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            old = np.asarray([gmm.calculate_likelihood_gmm(p, *args) for p in parameters])
            new = np.asarray([calculate_likelihood_gmm(p, *args) for p in parameters])
        assert_exact(old, new)
        elapsed = []
        for function in (gmm.calculate_likelihood_gmm, calculate_likelihood_gmm):
            began = perf_counter()
            for param in parameters:
                function(param, *args)
            elapsed.append(perf_counter()-began)
        rows.append({"mode": label, "calls": calls, "exact_equal": True,
                     "frozen_seconds": elapsed[0], "direct_seconds": elapsed[1],
                     "saving_fraction": 1-elapsed[1]/elapsed[0], "speedup": elapsed[0]/elapsed[1]})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--micro-calls", type=int, default=300)
    parser.add_argument("--full-fit", action="store_true")
    parser.add_argument("--events", type=int, nargs="+", default=[4, 14])
    parser.add_argument("--modes", choices=["original", "repaired"], nargs="+", default=["original"])
    parser.add_argument("--seed", type=int, default=31009402)
    parser.add_argument("--permuted", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.micro_calls < 1:
        parser.error("--micro-calls must be positive")
    report = {"phase": "development", "implementation": FAST_LIKELIHOOD_VERSION,
              "microbenchmark": microbenchmark(args.micro_calls), "full_fits": []}
    if args.full_fit:
        for events in args.events:
            if not 2 <= events <= len(DesignConfig().biomarker_names):
                parser.error("--events must be between 2 and 14")
            design = resolve("REF_H0", biomarker_names=DesignConfig().biomarker_names[:events])
            frame, truth = simulate(design, args.seed)
            if args.permuted:
                frame = frame.copy()
                frame["APOE"] = np.random.default_rng(args.seed+90000000).permutation(frame.APOE.to_numpy())
            for mode in args.modes:
                config = EngineConfig(mode=mode, expected_events=events, audit_original_neighbors=True,
                                      biomarker_names=tuple(truth["biomarker_names"]))
                report["full_fits"].append(full_fit_equivalence(frame, config) | {
                    "seed": args.seed, "permuted": args.permuted})
    report["all_checks_passed"] = True
    payload = json.dumps(report, indent=2, allow_nan=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload+"\n")
    print(payload)


if __name__ == "__main__":
    main()
