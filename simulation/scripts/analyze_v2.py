#!/usr/bin/env python3
"""Analyze EBM simulation JSONL or explicitly selected legacy CSV runs.

Examples (run from the project root)::

    python scripts/v2/analyze_v2.py --legacy-dir results/by_cell \
      --run-id legacy_B99 legacy_B599 --engine pyebm_v1 --all-cells \
      --rules strict le --out-dir results/v2_legacy_review
    python scripts/v2/analyze_v2.py --jsonl results/new_run/*.jsonl \
      --run-id simulation_v2 --engine pyebm --cell REF STRESS_STAGE \
      --rules le --out-dir results/new_run/analysis

Identity is (run_id, engine, cell, seed); duplicate identities are errors, never
silently dropped. JSONL fit-record streams must NOT be supplied as row inputs.

Expected row fields: run_id, engine, cell, seed, status, truth, observed, schemes.
truth.pair_null contains three booleans in PAIRS order. observed may include
taus/orderings/distance_to_truth. A scheme records requested_nperm, nperm, gt[3],
eq[3], complete, and optional maxgt/maxeq. An incomplete scheme may record
exact_decision={alpha:0.05, strict:{pair_reject:[bool|null]*3,max_reject:...},
le:{...}}. Decisions can also be certified from worst/best case bounds on all
unfinished or failed permutations. Only complete, full-budget rows contribute
permutation p-values; early-stop decisions remain eligible for rejection rates.
For sequential runs, diagnostics.sequential_stopping=true and
diagnostics.full_p_selected identify the preselected full-budget subset. Only
that subset enters p-distribution summaries. Mechanism-only rows with schemes={}
are supported; paired order-recovery comparisons across selected engines use
matching seeds within a run and cell and appear in mechanism_comparisons.csv.
scheme_comparisons retains legacy U/D, new unrestricted/diagnosis, and
diagnosis/oracle_dx_stage contrasts separately. --primary-comparison chooses
which contrast is also exposed through the backward-compatible paired field.
paired.p_shift remains median-of-dataset-medians. mean_pair_p_shift instead uses
one mean across the three pairs per preselected full-budget dataset, then a mean
and seed-level bootstrap/MCSE across datasets. Its per-pair and max-p components
are reported separately and also require full_p_selected=true.

Intervals quantify outer-dataset Monte Carlo uncertainty. Paired effects are
bootstrapped at dataset level. No empirical power-size calibration or equivalence
claim is made by this program. Partial-null error rates are diagnostics, not a
proof of strong FWER control for a procedure justified only under the global null.
Uses the Python standard library only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import random
import statistics
import sys

PAIRS = ("e2_vs_e33", "e2_vs_e4", "e33_vs_e4")
RULES = ("strict", "le")
IDENTITY = ("run_id", "engine", "cell", "seed")
SCHEME_COMPARISONS = {
    "stratified_minus_unstratified": ("unstratified", "stratified"),
    "diagnosis_minus_unrestricted": ("unrestricted", "diagnosis"),
    "oracle_dx_stage_minus_diagnosis": ("diagnosis", "oracle_dx_stage"),
}


def integer(value, label):
    """Reject nonintegral counts instead of silently truncating them."""
    if isinstance(value, bool):
        raise ValueError(f"{label} must be an integer count, not bool")
    try:
        numeric = Fraction(str(value))
        number = int(numeric)
        if numeric != number or number < 0:
            raise ValueError
    except (ValueError, TypeError, OverflowError):
        raise ValueError(f"{label} must be a nonnegative integer: {value!r}") from None
    return number


def alpha_fraction(alpha):
    value = Fraction(str(alpha))
    if not 0 < value < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    return value


def rejection_rank(budget, alpha=0.05, comparisons=3, rule="le"):
    """Largest rejecting integer numerator in (1 + exceedances)/(B + 1)."""
    budget = integer(budget, "budget")
    if rule not in RULES or comparisons < 1:
        raise ValueError("unsupported decision rule or comparisons")
    a = alpha_fraction(alpha)
    numerator = a.numerator * (budget + 1)
    denominator = a.denominator * comparisons
    rank, remainder = divmod(numerator, denominator)
    if rule == "strict" and remainder == 0:
        rank -= 1
    return max(0, rank)


def rejects(exceedances, budget, alpha=0.05, comparisons=3, rule="le"):
    count = integer(exceedances, "exceedances")
    if count > budget:
        raise ValueError("exceedances cannot exceed the permutation budget")
    return 1 + count <= rejection_rank(budget, alpha, comparisons, rule)


def bounded_decision(exceedances, nperm, budget, alpha, comparisons, rule):
    """Certify a decision only when every possible completion gives it."""
    if budget < nperm:
        raise ValueError("nperm exceeds requested_nperm")
    best = rejects(exceedances, budget, alpha, comparisons, rule)
    worst = rejects(exceedances + budget - nperm, budget, alpha, comparisons, rule)
    return best if best == worst else None


def wilson(k, n, z=1.959963984540054):
    if not n:
        return None
    p = k / n
    denominator = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


def binomial_summary(values, attempted=None):
    observed = [int(v) for v in values if v is not None]
    if any(v not in (0, 1) for v in observed):
        raise ValueError("binary endpoint contains a nonbinary value")
    n, k = len(observed), sum(observed)
    total = len(values) if attempted is None else attempted
    return {"n_evaluable": n, "n_unevaluable": total - n, "reject_count": k,
            "rate": k / n if n else None,
            "mcse": math.sqrt((k / n) * (1 - k / n) / n) if n else None,
            "wilson_ci95": wilson(k, n),
            "all_attempted_rate_bounds": [k / total, (k + total - n) / total]
            if total else None}


def percentile(sorted_values, probability):
    position = (len(sorted_values) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    return (sorted_values[lower] * (upper - position)
            + sorted_values[upper] * (position - lower)) if lower != upper else sorted_values[lower]


def bootstrap(values, statistic, repeats, seed, context):
    if not values or repeats <= 0:
        return None
    digest = hashlib.sha256(f"{seed}:{context}".encode()).digest()
    rng = random.Random(int.from_bytes(digest[:8], "big"))
    estimates = sorted(statistic(rng.choices(values, k=len(values))) for _ in range(repeats))
    return [percentile(estimates, 0.025), percentile(estimates, 0.975)]


def continuous_summary(values, repeats, seed, context, use_median=False):
    clean = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    statistic = statistics.median if use_median else statistics.mean
    return {"n_datasets": len(clean), "estimate": statistic(clean) if clean else None,
            "estimand": "median_of_dataset_medians" if use_median else "mean_dataset_difference",
            "bootstrap_ci95": bootstrap(clean, statistic, repeats, seed, context),
            "bootstrap_degenerate_sample": bool(clean) and min(clean) == max(clean),
            "interval_caution": "A constant observed sample yields a degenerate percentile interval; it does not prove equivalence."
            if clean and min(clean) == max(clean) else None,
            "mcse_mean": statistics.stdev(clean) / math.sqrt(len(clean))
            if len(clean) > 1 and not use_median else None}


def _legacy_float(row, key):
    value = row.get(key)
    return float(value) if value not in (None, "") else None


def load_legacy(directory):
    paths = sorted(Path(directory).glob("*.csv"))
    if not paths:
        raise ValueError(f"no CSV files in {directory}")
    rows = []
    for path in paths:
        with path.open(newline="") as handle:
            source_rows = list(csv.DictReader(handle))
        if not source_rows:
            continue
        for line, old in enumerate(source_rows, 2):
            cell = old["cell"]
            # The single B599 failure has 598 successful fits; filename/run
            # determines its intended budget, not a pooled median denominator.
            intended = 599 if "B599" in path.stem else max(
                integer(old[f"{s}_nperm"], "legacy nperm") for s in ("unstratified", "stratified"))
            run_id = f"legacy_B{intended}"
            distance = float(old.get("h1_d", 0))
            target = float(cell[4:]) if cell.startswith("PWR_") else 0.0
            truth = {"pair_null": [True, distance == 0, distance == 0],
                     "h1_distance": distance, "target_distance": target,
                     "target_distance_realized": round(target * 91) / 91,
                     "distance_tolerance": 0.000051,
                     "source": "legacy_cell_design"}
            row = {"run_id": run_id, "engine": "pyebm_v1", "cell": cell,
                   "seed": integer(old["seed"], "seed"), "status": "ok",
                   "truth": truth, "observed": {"taus": [_legacy_float(old, f"obs_{p}") for p in PAIRS]},
                   "schemes": {}, "diagnostics": {"legacy": True, "n_mci_e2": integer(old["n_mci_e2"], "n_mci_e2")},
                   "source": f"{path.resolve()}:{line}"}
            for scheme in ("unstratified", "stratified"):
                n = integer(old[f"{scheme}_nperm"], "nperm")
                entry = {"requested_nperm": intended, "nperm": n,
                         "complete": n == intended, "attempted": intended,
                         "gt": [integer(old[f"{scheme}_ngt_{p}"], "gt") for p in PAIRS],
                         "eq": [integer(old[f"{scheme}_neq_{p}"], "eq") for p in PAIRS],
                         "legacy_recorded_pge": [_legacy_float(old, f"{scheme}_pge_{p}") for p in PAIRS],
                         "legacy_recorded_reject_ge": integer(old[f"{scheme}_reject_ge"], "legacy rejection")}
                for key in ("maxgt", "maxeq"):
                    if old.get(f"{scheme}_{key}") not in (None, ""):
                        entry[key] = integer(old[f"{scheme}_{key}"], key)
                row["schemes"][scheme] = entry
            rows.append(row)
    return rows


def load_jsonl(paths):
    rows = []
    for name in paths:
        path = Path(name)
        with path.open() as handle:
            for line, text in enumerate(handle, 1):
                if not text.strip():
                    continue
                try:
                    row = json.loads(text)
                except json.JSONDecodeError as exc:
                    raise ValueError(f"{path}:{line}: {exc}") from exc
                if not isinstance(row, dict):
                    raise ValueError(f"{path}:{line}: expected a dataset object")
                row = dict(row)
                row["source"] = f"{path.resolve()}:{line}"
                rows.append(row)
    return rows


def validate_identities(rows):
    seen = {}
    for row in rows:
        missing = [key for key in IDENTITY if key not in row]
        if missing:
            raise ValueError(f"missing identity fields {missing}: {row.get('source')}")
        if any(not isinstance(row[key], str) or not row[key] for key in IDENTITY[:3]):
            raise ValueError(f"run_id, engine, and cell must be nonempty strings: {row.get('source')}")
        if type(row["seed"]) is not int or row["seed"] < 0:
            raise ValueError(f"seed must be a nonnegative integer: {row.get('source')}")
        if not isinstance(row.get("schemes", {}), dict):
            raise ValueError(f"schemes must be an object: {row.get('source')}")
        identity = tuple(row[key] for key in IDENTITY)
        if identity in seen:
            raise ValueError(f"duplicate dataset identity {identity}: {seen[identity]} and {row.get('source')}")
        seen[identity] = row.get("source")


def null_pairs(row):
    values = row.get("truth", {}).get("pair_null")
    if values is None:
        return None
    if not isinstance(values, list) or len(values) != len(PAIRS) or any(type(x) is not bool for x in values):
        raise ValueError(f"truth.pair_null must contain three booleans: {row.get('source')}")
    return values


def any_known(values):
    if not values:
        return None
    if any(v is True for v in values):
        return True
    return False if all(v is False for v in values) else None


def evaluate_scheme(entry, alpha, rule):
    """Return certified decisions and full-budget p-values for one dataset."""
    empty = {"pair_reject": [None] * 3, "max_reject": None,
             "pge": None, "pmax": None, "nperm": None, "budget": None,
             "complete": False, "decision_source": "unavailable", "tested_pair_index": None}
    if not entry:
        return empty
    result = dict(empty)
    tested = entry.get("definition", {}).get("tested_pair_index")
    if tested is not None and (type(tested) is not int or tested not in range(3)):
        raise ValueError("definition.tested_pair_index must be 0, 1, 2, or null")
    result["tested_pair_index"] = tested
    n = integer(entry.get("nperm", 0), "nperm")
    budget = integer(entry.get("requested_nperm", n), "requested_nperm")
    if budget < n:
        raise ValueError("requested_nperm is smaller than nperm")
    complete = entry.get("complete", n == budget)
    if type(complete) is not bool or complete != (n == budget):
        raise ValueError("complete must agree with nperm == requested_nperm")
    result.update(nperm=n, budget=budget, complete=complete)
    gt, eq = entry.get("gt"), entry.get("eq")
    if gt is not None or eq is not None:
        if not isinstance(gt, list) or not isinstance(eq, list) or len(gt) != 3 or len(eq) != 3:
            raise ValueError("gt and eq must both contain three counts")
        counts = [integer(g, "gt") + integer(e, "eq") for g, e in zip(gt, eq)]
        if any(k > n for k in counts):
            raise ValueError("gt + eq exceeds nperm")
        result["pair_reject"] = [bounded_decision(k, n, budget, alpha, 3, rule) for k in counts]
        result["decision_source"] = "full_counts" if complete else "completion_bounds"
        if complete and budget > 0:
            result["pge"] = [(1 + k) / (1 + n) for k in counts]
    if "maxgt" in entry or "maxeq" in entry:
        count = integer(entry.get("maxgt"), "maxgt") + integer(entry.get("maxeq"), "maxeq")
        if count > n:
            raise ValueError("maxgt + maxeq exceeds nperm")
        result["max_reject"] = bounded_decision(count, n, budget, alpha, 1, rule)
        if complete and budget > 0:
            result["pmax"] = (1 + count) / (1 + n)
    recorded = entry.get("exact_decision", {})
    if rule in recorded and "alpha" not in recorded:
        raise ValueError("recorded exact decisions must specify their alpha")
    if rule in recorded and alpha_fraction(recorded["alpha"]) == alpha_fraction(alpha):
        declaration = recorded[rule]
        values = declaration.get("pair_reject", [None] * 3)
        if len(values) != 3 or any(v is not None and type(v) is not bool for v in values):
            raise ValueError("exact_decision.pair_reject must contain bool/null values")
        for index, value in enumerate(values):
            derived = result["pair_reject"][index]
            if value is not None and derived is not None and value != derived:
                raise ValueError("recorded exact decision conflicts with count bounds")
            if value is not None:
                result["pair_reject"][index] = value
        max_value = declaration.get("max_reject")
        if max_value is not None:
            if type(max_value) is not bool:
                raise ValueError("max_reject must be bool/null")
            if result["max_reject"] is not None and result["max_reject"] != max_value:
                raise ValueError("recorded max decision conflicts with count bounds")
            result["max_reject"] = max_value
        result["decision_source"] += "+recorded_exact"
    if tested is not None:
        # A label permutation of one pair does not test the other comparisons
        # or the global max, even if their raw statistics were also recorded.
        result["pair_reject"] = [value if i == tested else None for i, value in enumerate(result["pair_reject"])]
        if result["pge"] is not None:
            result["pge"] = [value if i == tested else None for i, value in enumerate(result["pge"])]
        result["max_reject"] = result["pmax"] = None
    return result


def endpoints(row, evaluated):
    pairs = evaluated["pair_reject"]
    null = null_pairs(row)
    output = {"global_reject": any_known(pairs), "max_global_reject": evaluated["max_reject"]}
    output.update({f"pair_{pair}": value for pair, value in zip(PAIRS, pairs)})
    output["any_true_alternative_reject"] = any_known([value for value, is_null in zip(pairs, null) if not is_null]) if null else None
    output["any_true_null_reject"] = any_known([value for value, is_null in zip(pairs, null) if is_null]) if null else None
    if evaluated.get("tested_pair_index") is not None:
        for endpoint in ("global_reject", "max_global_reject", "any_true_alternative_reject", "any_true_null_reject"):
            output[endpoint] = None
    return output


def endpoint_applicable(row, endpoint, evaluated=None):
    tested = evaluated.get("tested_pair_index") if evaluated else None
    if tested is not None:
        return endpoint == f"pair_{PAIRS[tested]}"
    null = null_pairs(row)
    if null is None:
        return True  # Unknown truth is an unevaluable result, not an absent hypothesis.
    if endpoint == "any_true_alternative_reject":
        return not all(null)
    if endpoint == "any_true_null_reject":
        return any(null)
    return True


def _numeric_summary(values):
    numbers = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    return {"n": len(numbers), "mean": statistics.mean(numbers) if numbers else None,
            "median": statistics.median(numbers) if numbers else None,
            "minimum": min(numbers) if numbers else None, "maximum": max(numbers) if numbers else None}


def p_summary_eligible(row):
    """Do not select a p-distribution sample by whether a sequential run finished."""
    diagnostics = row.get("diagnostics", {})
    return not diagnostics.get("sequential_stopping", False) or diagnostics.get("full_p_selected") is True


def compare_schemes(rows, evaluations, outcomes, scheme_a, scheme_b,
                    context, context_base, bootstrap_repeats, bootstrap_seed):
    """Pair by dataset, with explicit applicability and full-p selection counts.

    The historical p_shift remains median(median_pair(p_b-p_a)). New mean-p
    endpoints require explicit preselection even for a nonsequential input.
    """
    a_values, b_values = evaluations[scheme_a], evaluations[scheme_b]
    direction = f"{scheme_b} minus {scheme_a}"
    metrics = {}
    for endpoint in outcomes[scheme_a][0]:
        applicable = [(a[endpoint], b[endpoint])
                      for row, a, b, av, bv in zip(rows, outcomes[scheme_a], outcomes[scheme_b], a_values, b_values)
                      if endpoint_applicable(row, endpoint, av) and endpoint_applicable(row, endpoint, bv)]
        known = [(a, b) for a, b in applicable if a is not None and b is not None]
        effect = continuous_summary([int(b)-int(a) for a, b in known], bootstrap_repeats,
                                    bootstrap_seed, context+":"+endpoint)
        effect.update(direction=direction, n_applicable=len(applicable),
                      n_inapplicable=len(rows)-len(applicable), n_unevaluable=len(applicable)-len(known),
                      paired_table=dict(Counter(f"{int(a)}{int(b)}" for a, b in known)),
                      paired_table_scheme_order=[scheme_a, scheme_b])
        metrics[endpoint] = effect

    def p_effect(endpoint, getter, applicability, *, legacy=False, use_median=False):
        values = []
        n_applicable = n_selected = n_not_selected = 0
        for row, a, b in zip(rows, a_values, b_values):
            if not applicability(a, b):
                continue
            n_applicable += 1
            selected = p_summary_eligible(row) if legacy else row.get("diagnostics", {}).get("full_p_selected") is True
            if not selected:
                n_not_selected += 1
                continue
            n_selected += 1
            # evaluate_scheme exposes p values only at each complete successful
            # planned budget. Failed observed fits supply no p values.
            value = getter(a, b) if a["complete"] and b["complete"] else None
            if value is not None:
                values.append(value)
        seed_context = context_base+":p_shift" if legacy else context_base+":"+endpoint
        effect = continuous_summary(values, bootstrap_repeats, bootstrap_seed, seed_context, use_median=use_median)
        effect.update(direction=direction, n_applicable=n_applicable,
                      n_inapplicable=len(rows)-n_applicable, n_selected=n_selected,
                      n_not_selected=n_not_selected, n_unevaluable=n_selected-len(values),
                      selection="historical p-summary eligibility; complete budgets only" if legacy
                                else "full_p_selected=true and both schemes complete at their planned budgets",
                      bootstrap_unit="simulated dataset (seed)")
        metrics[endpoint] = effect
        return effect

    joint = lambda a, b: a["tested_pair_index"] is None and b["tested_pair_index"] is None
    def pair_deltas(a, b):
        if a["pge"] is None or b["pge"] is None or any(v is None for v in a["pge"]+b["pge"]):
            return None
        return [y-x for x, y in zip(a["pge"], b["pge"])]
    def reduce_pairs(a, b, statistic):
        differences = pair_deltas(a, b)
        return statistic(differences) if differences is not None else None
    median = p_effect("p_shift", lambda a, b: reduce_pairs(a, b, statistics.median), joint,
                      legacy=True, use_median=True)
    median["direction"] += "; median across three pairs within each dataset first"
    mean = p_effect("mean_pair_p_shift", lambda a, b: reduce_pairs(a, b, statistics.mean), joint)
    mean.update(estimand="mean_of_dataset_mean_pair_p_differences",
                direction=direction+"; mean across three pairs within each dataset, then mean across datasets")
    for index, pair in enumerate(PAIRS):
        def applicable(a, b, i=index):
            return all(value["tested_pair_index"] in (None, i) for value in (a, b))
        def difference(a, b, i=index):
            if a["pge"] is None or b["pge"] is None or a["pge"][i] is None or b["pge"][i] is None:
                return None
            return b["pge"][i]-a["pge"][i]
        p_effect(f"pair_p_shift_{pair}", difference, applicable)
    p_effect("max_p_shift", lambda a, b: b["pmax"]-a["pmax"]
             if a["pmax"] is not None and b["pmax"] is not None else None, joint)
    return {"scheme_a": scheme_a, "scheme_b": scheme_b, "direction": direction,
            "n_rows": len(rows), "metrics": metrics}


def oracle_engine(name):
    return "oracle" in name.lower()


def truth_distances(row):
    if row.get("status") != "ok":
        return None
    values = row.get("observed", {}).get("distance_to_truth")
    if isinstance(values, dict):
        values = [values.get(group) for group in ("e2", "e33", "e4")]
    if not isinstance(values, list) or len(values) != 3:
        return None
    if any(value is None or not math.isfinite(float(value)) for value in values):
        return None
    return [float(value) for value in values]


def mechanism_comparisons(rows, engine_order, bootstrap_repeats=2000, bootstrap_seed=20260908):
    """Paired recovery errors by seed; oracle arms remain diagnostic references."""
    grouped = defaultdict(lambda: defaultdict(dict))
    for row in rows:
        grouped[(row["run_id"], row["cell"])][row["engine"]][row["seed"]] = row
    output = []
    for (run_id, cell), engines in sorted(grouped.items()):
        ordered = [engine for engine in dict.fromkeys(engine_order) if engine in engines]
        for engine_a, engine_b in combinations(ordered, 2):
            a, b = engines[engine_a], engines[engine_b]
            shared = sorted(a.keys() & b.keys())
            observations = []
            for seed in shared:
                first, second = truth_distances(a[seed]), truth_distances(b[seed])
                if first is not None and second is not None:
                    observations.append((first + [statistics.mean(first)], second + [statistics.mean(second)]))
            if not observations:
                continue  # Legacy CSVs do not contain order recovery errors.
            diagnostic = oracle_engine(engine_a) or oracle_engine(engine_b)
            for index, group in enumerate(("e2", "e33", "e4", "group_mean")):
                first = [value[0][index] for value in observations]
                second = [value[1][index] for value in observations]
                differences = [y - x for x, y in zip(first, second)]
                effect = continuous_summary(differences, bootstrap_repeats, bootstrap_seed,
                    f"mechanism:{run_id}:{cell}:{engine_a}:{engine_b}:{group}")
                output.append({"run_id": run_id, "cell": cell, "engine_a": engine_a,
                    "engine_b": engine_b, "group": group,
                    "comparison_role": "oracle_mechanism_diagnostic" if diagnostic else "algorithm_mechanism_comparison",
                    "oracle_is_inference_candidate": False,
                    "n_seed_union": len(a.keys() | b.keys()), "n_shared_seeds": len(shared),
                    "n_missing_or_failed_pairs": len(a.keys() | b.keys()) - len(observations),
                    "mean_distance_a": statistics.mean(first), "mean_distance_b": statistics.mean(second),
                    "paired_difference": effect,
                    "direction": "engine_b minus engine_a; negative values mean improved order recovery",
                    "interval_scope": "pointwise dataset-bootstrap interval; mechanism evidence, not a validity test"})
    return output


def summarize_pairwise_families(rows, evaluations, alpha, rule):
    families = sorted({entry.get("definition", {}).get("family") for row in rows
                       for entry in row.get("schemes", {}).values()
                       if entry.get("definition", {}).get("family")})
    result = {}
    for family in families:
        assembled, member_details = [], []
        for row_index, row in enumerate(rows):
            members = {}
            for name, entry in row.get("schemes", {}).items():
                definition = entry.get("definition", {})
                if definition.get("family") != family:
                    continue
                index = definition.get("tested_pair_index")
                if type(index) is not int or index not in range(3):
                    raise ValueError(f"family {family!r} member {name!r} lacks a valid tested_pair_index")
                if index in members:
                    raise ValueError(f"family {family!r} has duplicate members for pair {index}: {row.get('source')}")
                members[index] = (name, evaluations[name][row_index])
            values = [members.get(i, (None, {}))[1] for i in range(3)]
            decisions = [value.get("pair_reject", [None] * 3)[i] for i, value in enumerate(values)]
            complete = len(members) == 3 and all(value.get("pge") is not None for value in values)
            pge = [value["pge"][i] for i, value in enumerate(values)] if complete else None
            assembled.append({"pair_reject": decisions, "max_reject": None, "pge": pge, "complete": complete})
            member_details.append({PAIRS[i]: {"scheme": name, "nperm": value["nperm"], "budget": value["budget"]}
                                   for i, (name, value) in members.items()})
        outcome_rows = [endpoints(row, value) for row, value in zip(rows, assembled)]
        metrics = {}
        for endpoint in outcome_rows[0]:
            applicable = [outcome[endpoint] for row, outcome in zip(rows, outcome_rows)
                          if endpoint != "max_global_reject" and endpoint_applicable(row, endpoint)]
            metrics[endpoint] = binomial_summary(applicable)
            metrics[endpoint]["n_inapplicable"] = len(rows) - len(applicable)
        family_result = {"scope": "bonferroni_family_of_three_separate_pairwise_permutations", "alpha": float(alpha_fraction(alpha)),
            "rule": rule, "metrics": metrics, "n_rows": len(rows),
            "n_full_p_datasets": sum(value["complete"] for value in assembled),
            "n_p_summary_datasets": sum(value["complete"] and p_summary_eligible(row) for row, value in zip(rows, assembled)),
            "member_scheme_names": {pair: sorted({detail[pair]["scheme"] for detail in member_details if pair in detail}) for pair in PAIRS},
            "per_pair_budget_counts": {pair: dict(Counter(str(detail[pair]["budget"]) for detail in member_details if pair in detail)) for pair in PAIRS},
            "full_p_summary": {pair: _numeric_summary([value["pge"][i] for row, value in zip(rows, assembled)
                                 if value["pge"] is not None and p_summary_eligible(row)]) for i, pair in enumerate(PAIRS)},
            "notes": ["Each pair uses its own successful counts and planned budget; Bonferroni threshold is alpha/3.",
                      "No global-max test is assembled from separate pairwise permutation distributions.",
                      "Partial-null error rates remain diagnostics until the pairwise procedure's assumptions are established."]}
        result[family] = family_result
    return result


def summarize_group(rows, alpha=0.05, rule="le", bootstrap_repeats=2000, bootstrap_seed=20260908,
                    primary_comparison="auto"):
    first = rows[0]
    context_base = ":".join(str(first[k]) for k in IDENTITY[:3])
    context = context_base + ":" + rule
    truth_counts = Counter("unknown" if null_pairs(row) is None else
                           "global_null" if all(null_pairs(row)) else
                           "partial_null" if any(null_pairs(row)) else "all_alternative" for row in rows)
    summary = {key: first[key] for key in IDENTITY[:3]}
    summary.update(rule=rule, alpha=float(alpha_fraction(alpha)), n_rows=len(rows),
                   status_counts=dict(Counter(str(row.get("status", "unknown")) for row in rows)),
                   truth_counts=dict(truth_counts), schemes={}, paired={}, scheme_comparisons={}, primary_comparison=None)
    summary["engine_role"] = "oracle_mechanism_diagnostic_only" if oracle_engine(first["engine"]) else "non_oracle_engine"
    distances = [row.get("truth", {}).get("h1_distance") for row in rows]
    summary["realized_distance"] = _numeric_summary(distances)
    summary["realized_distance"]["counts"] = dict(Counter(str(d) for d in distances if d is not None))
    assessed = []
    for row in rows:
        truth = row.get("truth", {})
        target, actual = truth.get("target_distance_realized"), truth.get("h1_distance")
        if target is not None and actual is not None:
            assessed.append(abs(float(target) - float(actual)) <= float(truth.get("distance_tolerance", 1e-9)))
    summary["target_distance_check"] = {"n_assessable": len(assessed), "n_achieved": sum(assessed)}
    def observed_component(row, field, index, key):
        value = row.get("observed", {}).get(field)
        if isinstance(value, dict):
            return value.get(key)
        if isinstance(value, list) and index < len(value):
            return value[index]
        return None
    summary["observed_statistics"] = {
        "taus": {pair: _numeric_summary([observed_component(row, "taus", i, pair) for row in rows])
                 for i, pair in enumerate(PAIRS)},
        "distance_to_truth": {group: _numeric_summary([observed_component(row, "distance_to_truth", i, group) for row in rows])
                              for i, group in enumerate(("e2", "e33", "e4"))}}
    summary["interpretation"] = {
        "global_endpoint": "FWER" if truth_counts.get("global_null") == len(rows) else "omnibus_rejection_probability",
        "partial_null": "diagnostic_only; global-null permutation validity does not establish strong FWER control",
        "confidence_intervals": "outer-dataset Monte Carlo uncertainty; conditional on the specified DGP and procedure",
        "p_values": "only complete successful permutation budgets; sequential runs additionally require preselected full_p_selected=true",
        "failure_denominator": "rates use evaluable rows; all-attempted bounds and unevaluable counts are reported"}
    scheme_names = sorted({name for row in rows for name in row.get("schemes", {})})
    evaluations, outcomes = {}, {}
    for name in scheme_names:
        values = []
        for row in rows:
            try:
                entry = row.get("schemes", {}).get(name)
                # Failed observed fits do not become negative rejection outcomes.
                failed = str(row.get("status", "ok")).lower().startswith(("fail", "error"))
                values.append(evaluate_scheme(None if failed else entry, alpha, rule))
            except ValueError as exc:
                raise ValueError(f"{row.get('source')} [{name}]: {exc}") from exc
        declared_targets = {row.get("schemes", {}).get(name, {}).get("definition", {}).get("tested_pair_index")
                            for row in rows}
        declared_targets.discard(None)
        if len(declared_targets) > 1:
            raise ValueError(f"scheme {name!r} changes tested_pair_index within one run/engine/cell")
        if declared_targets:
            target = next(iter(declared_targets))
            for value in values:
                value["tested_pair_index"] = target
        evaluations[name] = values
        outcomes[name] = [endpoints(row, value) for row, value in zip(rows, values)]
        metric_names = list(outcomes[name][0])
        metrics = {}
        for metric in metric_names:
            applicable = [value[metric] for row, value, evaluated in zip(rows, outcomes[name], values)
                          if endpoint_applicable(row, metric, evaluated)]
            metrics[metric] = binomial_summary(applicable)
            metrics[metric]["n_inapplicable"] = len(rows) - len(applicable)
        result = {"metrics": metrics,
                  "n_full_p_datasets": sum(value["pge"] is not None for value in values),
                  "n_p_summary_datasets": sum(value["pge"] is not None and p_summary_eligible(row) for row, value in zip(rows, values)),
                  "n_max_available": sum(value["max_reject"] is not None for value in values),
                  "nperm_counts": dict(Counter(str(value["nperm"]) for value in values)),
                  "requested_nperm_counts": dict(Counter(str(value["budget"]) for value in values)),
                  "decision_source_counts": dict(Counter(value["decision_source"] for value in values)),
                  "attainable_rank_bounds": []}
        tested_indices = {value["tested_pair_index"] for value in values if value["tested_pair_index"] is not None}
        if len(tested_indices) > 1:
            raise ValueError(f"scheme {name!r} changes tested_pair_index within one run/engine/cell")
        result["scope"] = "single_pair_only" if tested_indices else "joint_label_permutation"
        result["tested_pair_index"] = next(iter(tested_indices), None)
        for budget in sorted({value["budget"] for value in values if value["budget"] is not None}):
            rank = rejection_rank(budget, alpha, 3, rule)
            result["attainable_rank_bounds"].append({"budget": budget, "pair_rank": rank,
                "per_pair_upper_bound": rank / (budget + 1),
                "bonferroni_fwer_upper_bound": min(1, 3 * rank / (budget + 1)),
                "note": "upper bound under valid exchangeability; ties and dependence can make size smaller"})
        result["full_p_summary"] = {pair: _numeric_summary([value["pge"][i] for row, value in zip(rows, values)
                                                            if value["pge"] is not None and p_summary_eligible(row)]) for i, pair in enumerate(PAIRS)}
        result["full_max_p_summary"] = _numeric_summary([value["pmax"] for row, value in zip(rows, values)
                                                          if value["pmax"] is not None and p_summary_eligible(row)])
        summary["schemes"][name] = result
    summary["pairwise_families"] = summarize_pairwise_families(rows, evaluations, alpha, rule)
    for comparison, (scheme_a, scheme_b) in SCHEME_COMPARISONS.items():
        if {scheme_a, scheme_b}.issubset(scheme_names):
            # Keep the old bootstrap RNG contexts byte-for-byte for legacy summaries.
            suffix = "" if comparison == "stratified_minus_unstratified" else ":"+comparison
            summary["scheme_comparisons"][comparison] = compare_schemes(
                rows, evaluations, outcomes, scheme_a, scheme_b, context+suffix, context_base+suffix,
                bootstrap_repeats, bootstrap_seed)
    if primary_comparison == "auto":
        selected = next((name for name in SCHEME_COMPARISONS if name in summary["scheme_comparisons"]), None)
    elif primary_comparison not in summary["scheme_comparisons"]:
        raise ValueError(f"Requested primary comparison {primary_comparison!r} unavailable in {context_base}")
    else:
        selected = primary_comparison
    if selected is not None:
        summary["primary_comparison"] = selected
        summary["paired"] = summary["scheme_comparisons"][selected]["metrics"]
    # Legacy audit preserves the historical per-successful-B policy explicitly;
    # it must not be confused with complete planned-budget p-value inference.
    if all(row.get("diagnostics", {}).get("legacy") for row in rows):
        legacy = {"policy": "historical denominator is successful nperm; audit only, including incomplete planned budgets", "schemes": {}}
        shifts = []
        for name in scheme_names:
            effective = [[rejects(g + e, entry["nperm"], alpha, 3, rule) for g, e in zip(entry["gt"], entry["eq"])]
                         for row in rows for entry in [row["schemes"][name]]]
            legacy["schemes"][name] = {"global_reject": binomial_summary([any(values) for values in effective]),
                "pair_reject_count": sum(sum(values) for values in effective), "n_pairs": len(rows) * 3,
                "recorded_strict_flag_mismatches": sum(any(values) != bool(row["schemes"][name]["legacy_recorded_reject_ge"])
                                                        for values, row in zip(effective, rows)) if rule == "strict" else None}
        for row in rows:
            a, b = row["schemes"].get("unstratified"), row["schemes"].get("stratified")
            if a and b:
                shifts.append(statistics.median(y - x for x, y in zip(a["legacy_recorded_pge"], b["legacy_recorded_pge"])))
        legacy["recorded_p_shift_dataset_median"] = continuous_summary(shifts, bootstrap_repeats, bootstrap_seed,
                                                                        context_base + ":legacy_p_shift", use_median=True)
        summary["legacy_audit"] = legacy
    return summary


def write_outputs(output, directory):
    destination = Path(directory)
    destination.mkdir(parents=True, exist_ok=True)
    with (destination / "analysis.json").open("w") as handle:
        json.dump(output, handle, indent=2, allow_nan=False)
        handle.write("\n")
    metric_rows, paired_rows, family_rows, scheme_comparison_rows = [], [], [], []
    for group in output["groups"]:
        identity = {key: group[key] for key in ("run_id", "engine", "cell", "rule", "alpha")}
        for scheme, detail in group["schemes"].items():
            for endpoint, metric in detail["metrics"].items():
                ci = metric["wilson_ci95"] or [None, None]
                bounds = metric["all_attempted_rate_bounds"] or [None, None]
                metric_rows.append(identity | {"scheme": scheme, "endpoint": endpoint,
                    "n_rows": group["n_rows"], "n_evaluable": metric["n_evaluable"],
                    "n_unevaluable": metric["n_unevaluable"], "n_inapplicable": metric["n_inapplicable"],
                    "reject_count": metric["reject_count"],
                    "rate": metric["rate"], "mcse": metric["mcse"], "ci95_low": ci[0], "ci95_high": ci[1],
                    "all_attempted_low": bounds[0], "all_attempted_high": bounds[1],
                    "n_full_p_datasets": detail["n_full_p_datasets"],
                    "n_p_summary_datasets": detail["n_p_summary_datasets"],
                    "realized_distance_mean": group["realized_distance"]["mean"]})
        for endpoint, metric in group["paired"].items():
            ci = metric["bootstrap_ci95"] or [None, None]
            paired_rows.append(identity | {"endpoint": endpoint, "n_datasets": metric["n_datasets"],
                "estimate": metric["estimate"], "estimand": metric["estimand"],
                "ci95_low": ci[0], "ci95_high": ci[1], "mcse_mean": metric["mcse_mean"],
                "direction": metric["direction"]})
        for comparison, detail in group.get("scheme_comparisons", {}).items():
            for endpoint, metric in detail["metrics"].items():
                ci = metric["bootstrap_ci95"] or [None, None]
                scheme_comparison_rows.append(identity | {"comparison": comparison,
                    "scheme_a": detail["scheme_a"], "scheme_b": detail["scheme_b"],
                    "is_primary": group.get("primary_comparison") == comparison,
                    "endpoint": endpoint, "n_datasets": metric["n_datasets"],
                    "n_applicable": metric["n_applicable"], "n_inapplicable": metric["n_inapplicable"],
                    "n_unevaluable": metric["n_unevaluable"], "n_selected": metric.get("n_selected"),
                    "n_not_selected": metric.get("n_not_selected"),
                    "estimate": metric["estimate"], "estimand": metric["estimand"],
                    "ci95_low": ci[0], "ci95_high": ci[1], "mcse_mean": metric["mcse_mean"],
                    "direction": metric["direction"]})
        for family, detail in group.get("pairwise_families", {}).items():
            for endpoint, metric in detail["metrics"].items():
                ci = metric["wilson_ci95"] or [None, None]
                bounds = metric["all_attempted_rate_bounds"] or [None, None]
                family_rows.append(identity | {"family": family, "endpoint": endpoint,
                    "n_rows": group["n_rows"], "n_evaluable": metric["n_evaluable"],
                    "n_unevaluable": metric["n_unevaluable"], "n_inapplicable": metric["n_inapplicable"],
                    "reject_count": metric["reject_count"], "rate": metric["rate"], "mcse": metric["mcse"],
                    "ci95_low": ci[0], "ci95_high": ci[1], "all_attempted_low": bounds[0], "all_attempted_high": bounds[1],
                    "n_full_p_datasets": detail["n_full_p_datasets"], "n_p_summary_datasets": detail["n_p_summary_datasets"]})
    for filename, rows in (("summary_metrics.csv", metric_rows), ("paired_effects.csv", paired_rows),
                           ("pairwise_family_metrics.csv", family_rows),
                           ("scheme_comparisons.csv", scheme_comparison_rows)):
        with (destination / filename).open("w", newline="") as handle:
            fallback = ["run_id", "engine", "cell", "rule", "alpha", "endpoint"]
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else fallback)
            writer.writeheader()
            writer.writerows(rows)
    comparison_rows = []
    for comparison in output.get("mechanism_comparisons", []):
        effect = comparison["paired_difference"]
        ci = effect["bootstrap_ci95"] or [None, None]
        comparison_rows.append({key: value for key, value in comparison.items() if key != "paired_difference"} |
                               {"n_evaluable_pairs": effect["n_datasets"], "paired_mean_difference": effect["estimate"],
                                "ci95_low": ci[0], "ci95_high": ci[1], "mcse_mean": effect["mcse_mean"]})
    with (destination / "mechanism_comparisons.csv").open("w", newline="") as handle:
        fallback = ["run_id", "cell", "engine_a", "engine_b", "group", "comparison_role", "n_evaluable_pairs", "paired_mean_difference"]
        writer = csv.DictWriter(handle, fieldnames=list(comparison_rows[0]) if comparison_rows else fallback)
        writer.writeheader()
        writer.writerows(comparison_rows)
    return destination / "analysis.json"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--jsonl", nargs="+", default=[], help="dataset-level JSONL inputs (not fit records)")
    parser.add_argument("--legacy-dir", help="directory containing the selected legacy per-cell CSV family")
    parser.add_argument("--run-id", nargs="+", required=True, help="explicit included run IDs")
    parser.add_argument("--engine", nargs="+", required=True, help="explicit included engine IDs")
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--cell", nargs="+", help="explicit included cell names")
    selection.add_argument("--all-cells", action="store_true", help="intentionally include every cell in selected runs/engines")
    parser.add_argument("--rules", nargs="+", choices=RULES, default=["le"])
    parser.add_argument("--alpha", default="0.05", help="family-wise alpha, parsed as an exact decimal fraction")
    parser.add_argument("--bootstrap", type=int, default=2000, help="dataset bootstrap replicates; zero disables intervals")
    parser.add_argument("--bootstrap-seed", type=int, default=20260908)
    parser.add_argument("--primary-comparison", choices=["auto", *SCHEME_COMPARISONS], default="auto",
                        help="comparison exposed through paired; all available contrasts remain in scheme_comparisons")
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args(argv)
    try:
        alpha_fraction(args.alpha)
        if not args.jsonl and not args.legacy_dir:
            raise ValueError("provide --jsonl and/or --legacy-dir")
        if args.bootstrap < 0:
            raise ValueError("bootstrap replicates cannot be negative")
        rows = load_jsonl(args.jsonl)
        if args.legacy_dir:
            rows.extend(load_legacy(args.legacy_dir))
        validate_identities(rows)
        for field, requested in (("run_id", args.run_id), ("engine", args.engine), ("cell", args.cell or [])):
            absent = set(requested) - {row[field] for row in rows}
            if absent:
                raise ValueError(f"requested {field} values absent from inputs: {sorted(absent)}")
        chosen = [row for row in rows if row["run_id"] in args.run_id and row["engine"] in args.engine
                  and (args.all_cells or row["cell"] in args.cell)]
        if not chosen:
            raise ValueError("no datasets match all explicit filters")
        grouped = defaultdict(list)
        for row in chosen:
            grouped[tuple(row[key] for key in IDENTITY[:3])].append(row)
        summaries = [summarize_group(sorted(group, key=lambda row: str(row["seed"])), args.alpha, rule,
                                    args.bootstrap, args.bootstrap_seed, args.primary_comparison)
                     for _, group in sorted(grouped.items()) for rule in dict.fromkeys(args.rules)]
        output = {"schema_version": "analysis_v2.2", "pair_order": list(PAIRS),
                  "inputs": {"jsonl": [str(Path(path).resolve()) for path in args.jsonl],
                             "legacy_dir": str(Path(args.legacy_dir).resolve()) if args.legacy_dir else None},
                  "selection": {"run_id": args.run_id, "engine": args.engine, "cell": args.cell,
                                "all_cells": args.all_cells, "n_loaded": len(rows), "n_selected": len(chosen)},
                  "primary_comparison_requested": args.primary_comparison,
                  "bootstrap": {"replicates": args.bootstrap, "seed": args.bootstrap_seed,
                                "unit": "simulated dataset", "method": "percentile"},
                  "notes": ["No silent deduplication, pooling across runs, or empirical size matching.",
                            "Comparing ratios to an attainable bound does not establish nominal validity.",
                            "P-shift is descriptive supporting evidence, not a substitute for error-rate calibration.",
                            "Shared seeds across cells may couple simulations; do not pool them as independent repetitions."],
                  "groups": summaries,
                  "mechanism_comparisons": mechanism_comparisons(chosen, args.engine, args.bootstrap, args.bootstrap_seed)}
        path = write_outputs(output, args.out_dir)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"analysis error: {exc}\n")
    print(f"Analyzed {len(chosen)} dataset rows in {len(grouped)} explicit run/engine/cell groups: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
