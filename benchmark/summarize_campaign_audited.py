"""Summarize the paired benchmark (Fig. 3b-d, Supplementary Table 2); no model is fitted.

Reads the per-dataset result directories written by runtime/run_dataset.py and the
reused results in runtime/inputs/historical_reused_rows.jsonl, and writes per-dataset
endpoints, rates with Wilson intervals (rates.csv) and paired bootstrap differences
between methods (paired_differences.csv) into a NEW output directory, only after every
table has been written successfully.

The published tables were produced with the four record-audit files of the original
run (--record-audit, --transport-audit, --historical-audit, --native-binding), which
check every result file against hashes taken when the results were returned from the
cluster. Without them the same tables are computed from the result directories alone;
rates.csv and paired_differences.csv do not depend on the audit files.
"""
from pathlib import Path
import argparse
import collections
import csv
import hashlib
import json
import math
import os
import platform
import tempfile

import numpy as np

# SHA-256 of runtime/FROZEN_PACKAGE.json as published in this repository.
CLOSURE = "f7589961c1d458b36a1f4a2a165ffc46f14d25c19468bfd3a5d6635ded2bd8b8"
# SHA-256 of FROZEN_PACKAGE.json of the executed runtime, which the historical audit file names.
EXECUTED_CLOSURE = "99052c6b8ca24729f08aa87da7293b3cd409217382eec3782da3f578827a8904"
ARMS = ("repaired", "shared", "invariant_min")
NAMES = {"repaired": "Separately fitted DEBM", "shared": "Pooled-score DEBM", "invariant_min": "CONCORD"}
GROUPS = ("e2", "e33", "e4")
PAIRS = ("e2_vs_e33", "e2_vs_e4", "e33_vs_e4")
ENDPOINTS = ("any_pair", "pair_1", "pair_2", "pair_3", "truly_affected_pair", "unchanged_pair")
CONTINUOUS = tuple(f"truth_error_g{i}" for i in (1, 2, 3)) + tuple(f"pair_distance_{i}" for i in (1, 2, 3))
CONTRASTS = (("shared", "repaired"), ("invariant_min", "shared"), ("invariant_min", "repaired"))
BOOTSTRAPS = 10000


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def strict_count(value, name):
    require(type(value) is int and value >= 0, f"invalid {name}")
    return value


def decision_summary(exceedances, successes, attempted, failed):
    """Inclusive Bonferroni-3 decisions for the fixed 599-slot stream."""
    for name, value in (("successes", successes), ("attempted", attempted), ("failed", failed)):
        strict_count(value, name)
    require(successes + failed == attempted <= 599, "permutation denominator mismatch")
    require(len(exceedances) == 3, "expected three exceedance counts")
    e = [strict_count(v, "exceedance") for v in exceedances]
    require(all(v <= successes for v in e), "exceedances exceed successes")
    lower_num = [1 + v for v in e]
    upper_num = [1 + v + 599 - successes for v in e]
    decisions = [True if hi <= 10 else False if lo > 10 else None
                 for lo, hi in zip(lower_num, upper_num)]
    require(all(v is not None for v in decisions), "unresolved fixed-budget decision")
    return {"exceedances": e, "decisions": decisions,
            "successful_permutations": successes, "attempted_permutations": attempted,
            "failed_permutations": failed, "full_tail": successes == 599,
            "raw_p_lower": [v / 600 for v in lower_num],
            "raw_p_upper": [v / 600 for v in upper_num],
            "adjusted_p_lower": [min(1., v / 200) for v in lower_num],
            "adjusted_p_upper": [min(1., v / 200) for v in upper_num]}


def endpoints(decisions, affected):
    require(len(decisions) == len(affected) == 3, "three pairs required")
    require(all(type(v) is bool for v in decisions + affected), "decisions/truth must be Boolean")
    out = {"any_pair": any(decisions), **{f"pair_{i+1}": v for i, v in enumerate(decisions)}}
    # Empty truth-defined sets are NA, never false. Under H0 all pairs are
    # unchanged, so unchanged_pair equals any_pair explicitly.
    out["truly_affected_pair"] = any(v for v, t in zip(decisions, affected) if t) if any(affected) else None
    out["unchanged_pair"] = any(v for v, t in zip(decisions, affected) if not t) if not all(affected) else None
    return out


def wilson(k, n):
    require(n > 0 and 0 <= k <= n, "invalid rate denominator")
    z = 1.959963984540054
    p, d = k / n, 1 + z * z / n
    mid = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [mid - half, mid + half]


def paired_ci(delta, seed):
    """10,000 paired percentile bootstrap resamples of datasets."""
    delta = np.asarray(delta, dtype=float)
    require(len(delta) > 0 and np.isfinite(delta).all(), "invalid paired differences")
    rng, means = np.random.default_rng(seed), []
    for _ in range(100):
        means.extend(delta[rng.integers(len(delta), size=(100, len(delta)))].mean(axis=1).tolist())
    return np.quantile(means, [.025, .975]).tolist()


def write_csv(path, rows, fields=None):
    require(bool(rows), f"empty planned table {path}")
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def view_arm(record, historical=False):
    if historical:
        require(record["status"] == "ok", "historical observed fit failed")
        s = record["schemes"]["diagnosis"]
        require(s["requested_nperm"] == 599 and s["definition"]["stratify"] == "diagnosis", "wrong historical test")
        e = [a + b for a, b in zip(s["gt"], s["eq"])]
        state = decision_summary(e, s["nperm"], s["attempted"], s["failed"])
        require(state["decisions"] == s["exact_decision"]["le"]["pair_reject"], "historical decision mismatch")
    else:
        state = decision_summary(record["exceedances"], record["computed_permutations"],
                                 record["attempted_permutations"], record["failed_permutations"])
        require(state["decisions"] == record["pair_reject"] and record["all_decisions_resolved"] is True,
                "new decision mismatch")
        require(record["full_permutation_budget_complete"] == state["full_tail"], "wrong full-tail flag")
        for field in ("adjusted_p_lower", "adjusted_p_upper"):
            require(np.allclose(state[field], record[field], rtol=0, atol=1e-12), f"wrong {field}")
        if state["full_tail"]:
            require(np.allclose(record["adjusted_p"], state["adjusted_p_lower"], rtol=0, atol=1e-12), "wrong exact P")
        else:
            require(record["adjusted_p"] is None, "truncated tail misrepresented as point P")
        require(record["observed"]["status"] == "ok", "new observed fit failed")
    obs = record["observed"]
    vals = list(obs["distance_to_truth"]) + list(obs["taus"])
    require(len(vals) == 6 and all(np.isfinite(v) and 0 <= v <= 1 for v in vals), "invalid continuous endpoints")
    return state, {k: float(v) for k, v in zip(CONTINUOUS, vals)}


def verify_runtime(runtime):
    require(sha(runtime / "FROZEN_PACKAGE.json") == CLOSURE, "wrong frozen closure")
    frozen = read(runtime / "FROZEN_PACKAGE.json")
    for name, expected in frozen["files"].items():
        require(sha(runtime / name) == expected, f"changed frozen file: {name}")
    require(sha(runtime / "EXECUTION_GATE.json") == frozen["execution_gate_sha256"], "changed execution gate")


def local_witness(folder):
    """Hashes of one result directory, used when no record audit is supplied."""
    log = folder / "fit_records.jsonl"
    return {"result.json": sha(folder / "result.json"), "manifest.json": sha(folder / "manifest.json"),
            "truth.json": sha(folder / "truth.json"),
            "manifest_identity_sha256": read(folder / "manifest.json")["identity_sha256"],
            "raw_fit_log_sha256": sha(log) if log.exists() else None}


def verify_gates(runtime, record_path, transport_path, history_path, binding_path):
    verify_runtime(runtime)
    record, transport, history = map(read, (record_path, transport_path, history_path))
    require(record["status"] == "RECORD_LAYER_PASS" and record["replay_data_and_all599_label_arrays"] is True
            and record["limited_smoke_check"] is False and not record["issues"] and not record["unexpected_result_paths"],
            "complete record-layer audit has not passed")
    require(all(record["counts"][k] == n for k, n in (("datasets", 7000), ("new_arms", 18000), ("reused_arms", 3000))),
            "record-audit scope incomplete")
    require(transport["status"] == "TRANSPORT_PASS" and transport["result_files"] == 7000 and transport["raw_record_files"] == 7000,
            "transport gate incomplete")
    for key in ("issues", "unexpected_local_files", "missing_local_files", "remote_changed_during_hash", "remote_symlinks", "local_symlinks", "incomplete_summary_paths"):
        require(not transport[key], f"transport problem: {key}")
    require(history["status"] == "HISTORICAL_ARCHIVE_JOIN_PASS" and not history["issues"] and not history["missing_compact_keys"],
            "historical archive gate incomplete")
    require(history["compact_index_sha256"] == sha(runtime / "inputs/historical_reused_rows.jsonl") and history["frozen_package_sha256"] == EXECUTED_CLOSURE,
            "historical gate names different frozen input")
    binding = read(binding_path)
    require(binding["status"] == "NATIVE_AUDIT_RETURN_BINDING_PASS" and binding["matched_dataset_files"] == 28000,
            "native audit has not been bound to the full local return")
    require(binding["counts"] == record["counts"] and binding["adverse_fit_records"] == record["adverse_fit_records"],
            "native binding refers to a different record audit")
    for target in (record_path, history_path):
        require(any(Path(path).name == target.name and value == sha(target)
                    for path, value in binding["evidence_sha256"].items()), "binding source identity mismatch")
    require(any(Path(path).name == "remote_transport_manifest.json" and value == transport["remote_manifest_sha256"]
                for path, value in binding["evidence_sha256"].items()), "binding refers to a different transfer")
    evidence = {}
    for row in record["evidence"]:
        key = (row["campaign"], row["cell"], row["seed"])
        require(key not in evidence, "duplicate audit evidence")
        evidence[key] = row
    require(len(evidence) == 7000, "record evidence incomplete")
    return record, evidence


def summarize_rows(cfg, rows, continuous):
    """Table derivation, separately callable by synthetic fixtures."""
    grouped, cg = collections.defaultdict(dict), collections.defaultdict(dict)
    for target, source in ((grouped, rows), (cg, continuous)):
        for row in source:
            cell, key = (row["campaign"], row["cell"]), (row["seed"], row["arm"])
            require(key not in target[cell], "duplicate estimator-dataset row")
            target[cell][key] = row
    rates, contrasts, summaries, continuous_contrasts = [], [], [], []
    for campaign in ("power", "robustness"):
        for ci, spec in enumerate(cfg[campaign]):
            cell = (campaign, spec["id"])
            seeds = list(range(spec["base_seed"], spec["base_seed"] + spec["datasets"]))
            expected = {(seed, arm) for seed in seeds for arm in ARMS}
            require(set(grouped[cell]) == set(cg[cell]) == expected, "missing/unplanned rows or changed denominator")
            eps = [ep for ep in ENDPOINTS if grouped[cell][(seeds[0], ARMS[0])][ep] is not None]
            for ep in ENDPOINTS:
                require(all((r[ep] is not None) == (ep in eps) for r in grouped[cell].values()), "endpoint applicability changed within cell")
            for arm in ARMS:
                for ep in eps:
                    vals = [grouped[cell][(seed, arm)][ep] for seed in seeds]
                    require(all(type(v) is bool for v in vals), "non-Boolean endpoint")
                    k, n = sum(vals), len(vals)
                    lo, hi = wilson(k, n)
                    rates.append({"campaign": campaign, "cell": spec["id"], "arm": arm, "endpoint": ep,
                                  "rejections": k, "planned_R": n, "rate": k / n, "wilson_low": lo, "wilson_high": hi})
                for ep in CONTINUOUS:
                    vals = np.array([cg[cell][(seed, arm)][ep] for seed in seeds])
                    require(np.isfinite(vals).all(), "invalid continuous values")
                    q1, med, q3 = np.quantile(vals, [.25, .5, .75])
                    summaries.append({"campaign": campaign, "cell": spec["id"], "arm": arm, "endpoint": ep,
                                      "planned_R": len(seeds), "mean": float(vals.mean()), "median": float(med),
                                      "q25": float(q1), "q75": float(q3)})
            for mi, (left, right) in enumerate(CONTRASTS):
                common = {"campaign": campaign, "cell": spec["id"], "left_arm": left, "right_arm": right,
                          "paired_R": len(seeds), "bootstrap_replicates": BOOTSTRAPS}
                for ep in eps:
                    ei = ENDPOINTS.index(ep)
                    d = np.array([int(grouped[cell][(seed, left)][ep]) - int(grouped[cell][(seed, right)][ep]) for seed in seeds])
                    bs_seed = 92600000 + ci * 1000 + mi * 100 + ei
                    lo, hi = paired_ci(d, bs_seed)
                    contrasts.append(common | {"endpoint": ep, "difference": float(d.mean()), "paired_bootstrap_low": lo,
                                               "paired_bootstrap_high": hi, "bootstrap_seed": bs_seed})
                for ei, ep in enumerate(CONTINUOUS):
                    d = np.array([cg[cell][(seed, left)][ep] - cg[cell][(seed, right)][ep] for seed in seeds])
                    bs_seed = 92700000 + (100000 if campaign == "robustness" else 0) + ci * 1000 + mi * 100 + ei
                    lo, hi = paired_ci(d, bs_seed)
                    continuous_contrasts.append(common | {"endpoint": ep, "difference": float(d.mean()), "paired_bootstrap_low": lo,
                                                          "paired_bootstrap_high": hi, "bootstrap_seed": bs_seed})
    require(set(grouped) == set(cg) == {(c, s["id"]) for c in ("power", "robustness") for s in cfg[c]}, "extra cells")
    return rates, contrasts, summaries, continuous_contrasts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runtime", type=Path, default=Path(__file__).resolve().parent / "runtime",
                    help="benchmark runtime directory (default: runtime/ next to this script)")
    ap.add_argument("--results", type=Path, required=True,
                    help="directory given as --output to run_batch.py / run_dataset.py")
    ap.add_argument("--output", type=Path, required=True, help="new directory for the summary tables")
    for name in ("record-audit", "transport-audit", "historical-audit", "native-binding"):
        ap.add_argument("--" + name, type=Path, help="record-audit file of the original run (all four or none)")
    a = ap.parse_args()
    require(not a.output.exists(), "use a new output directory; existing results are never overwritten")
    audit_files = (a.record_audit, a.transport_audit, a.historical_audit, a.native_binding)
    audited = all(p is not None for p in audit_files)
    require(audited or all(p is None for p in audit_files), "give all four audit files or none")
    if audited:
        audit, evidence = verify_gates(a.runtime, *audit_files)
    else:
        verify_runtime(a.runtime)
        audit, evidence = {"adverse_fit_records": []}, None
    cfg = read(a.runtime / "campaign.json")
    require(cfg["bperm"] == 599 and cfg["alpha"] == .05 and cfg["rule"] == "inclusive_bonferroni", "wrong scientific settings")
    hist = {}
    for line in (a.runtime / "inputs/historical_reused_rows.jsonl").open():
        row = json.loads(line)
        key = (row["run_id"], row["cell"], row["seed"], row["engine"])
        require(key not in hist, "duplicate historical row")
        hist[key] = row
    require(len(hist) == 5000, "historical index scope mismatch")
    adverse = collections.Counter((r["campaign"], r["cell"], r["seed"], r["arm"]) for r in audit["adverse_fit_records"])
    caps = collections.Counter((r["campaign"], r["cell"], r["seed"], r["arm"]) for r in audit["adverse_fit_records"] if r["cap"])
    rows, pair_rows, continuous, sensitivity, provenance, used_hist = [], [], [], [], [], set()
    planned = {(c, s["id"], s["base_seed"] + i) for c in ("power", "robustness") for s in cfg[c] for i in range(s["datasets"])}
    require(not audited or planned == set(evidence), "audit evidence differs from planned datasets")
    result_paths = {str(p.relative_to(a.results)) for p in a.results.glob("*/*/seed_*/result.json")}
    require(result_paths == {f"{c}/{cell}/seed_{seed}/result.json" for c, cell, seed in planned}, "missing/extra result paths")
    for campaign in ("power", "robustness"):
        for spec in cfg[campaign]:
            for seed in range(spec["base_seed"], spec["base_seed"] + spec["datasets"]):
                location = (campaign, spec["id"], seed)
                folder = a.results / campaign / spec["id"] / f"seed_{seed}"
                witness = evidence[location] if audited else local_witness(folder)
                for name in ("result.json", "manifest.json", "truth.json"):
                    require(sha(folder / name) == witness[name], f"post-audit changed bytes: {location}/{name}")
                result, truth = read(folder / "result.json"), read(folder / "truth.json")
                require(result["complete"] is True and result["unresolved_arms"] == [] and result["benchmark"] is False and result["budget"] == 599,
                        "incomplete/nonformal result")
                require(set(result["arms"]) == set(spec["new_arms"]), "new-arm roster mismatch")
                affected = [truth["pair_truth"][name]["inversions"] > 0 for name in PAIRS]
                for arm in ARMS:
                    new = arm in spec["new_arms"]
                    if new:
                        r = result["arms"][arm]
                    else:
                        key = (spec["historical_run"], spec["historical_cell"], seed, arm)
                        r = hist[key]
                        used_hist.add(key)
                    state, metrics = view_arm(r, historical=not new)
                    loc = {"campaign": campaign, "cell": spec["id"], "seed": seed, "arm": arm}
                    source = "new_allgroup01" if new else "certified_historical_reuse"
                    ep = endpoints(state["decisions"], affected)
                    key = location + (arm,)
                    if new and audited:
                        require(adverse[key] == state["failed_permutations"], "failure inventory disagrees with arm summary")
                    common = loc | {"source": source, "planned_permutations": 599,
                                    **{k: state[k] for k in ("successful_permutations", "attempted_permutations", "failed_permutations", "full_tail")},
                                    "outer_cap_records": caps[key] if new and audited else None,
                                    "outer_cap_inventory_scope": ("new raw records" if audited else "not audited") if new else "not present in compact historical index"}
                    rows.append(common | ep | {f"pair_{i+1}_true_affected": v for i, v in enumerate(affected)})
                    continuous.append(loc | {"source": source} | metrics)
                    for j, name in enumerate(PAIRS):
                        pair_rows.append(common | {"pair": name, "true_affected": affected[j], "reject": state["decisions"][j],
                                                   "exceedances": state["exceedances"][j],
                                                   **{k: state[k][j] for k in ("raw_p_lower", "raw_p_upper", "adjusted_p_lower", "adjusted_p_upper")},
                                                   "raw_p_exact": state["raw_p_lower"][j] if state["full_tail"] else None,
                                                   "adjusted_p_exact": state["adjusted_p_lower"][j] if state["full_tail"] else None})
                if campaign == "power" and spec["id"] != "GAUSSIAN_H0":
                    key = (spec["historical_run"], spec["historical_cell"], seed, "repaired")
                    r = hist[key]
                    used_hist.add(key)
                    state, metrics = view_arm(r, historical=True)
                    sensitivity.append({"campaign": campaign, "cell": spec["id"], "seed": seed, "arm": "historical_repaired",
                                        "historical_run": spec["historical_run"], "historical_cell": spec["historical_cell"],
                                        "role": "first-group stopping; secondary implementation sensitivity only",
                                        **{k: state[k] for k in ("successful_permutations", "attempted_permutations", "failed_permutations", "full_tail")},
                                        **endpoints(state["decisions"], affected), **metrics,
                                        **{f"pair_{j+1}_{bound}": state[bound][j] for j in range(3) for bound in ("adjusted_p_lower", "adjusted_p_upper")}})
                provenance.append({"campaign": campaign, "cell": spec["id"], "seed": seed,
                                   **{k: witness[k] for k in ("result.json", "manifest.json", "truth.json", "manifest_identity_sha256", "raw_fit_log_sha256")}})
    require(len(rows) == 21000 and len(pair_rows) == 63000 and len(continuous) == 21000 and len(sensitivity) == 2000,
            "wrong final table dimensions")
    require(used_hist == set(hist), "historical roster not exhausted exactly as primary reuse plus separate sensitivity")
    tables = summarize_rows(cfg, rows, continuous)
    require([len(t) for t in tables] == [189, 189, 198, 198], "summary endpoint counts differ from frozen grid")
    # historical_first_group_stopping.csv: separately fitted DEBM with the pyebm default stopping rule
    # on the same datasets (earlier simulation runs); excluded from rates.csv and paired_differences.csv.
    files = {"per_dataset_endpoints.csv": rows, "pair_decisions_and_bounds.csv": pair_rows,
             "truth_error_and_distance.csv": continuous, "historical_first_group_stopping.csv": sensitivity,
             "rates.csv": tables[0], "paired_differences.csv": tables[1], "continuous_summary.csv": tables[2],
             "paired_continuous_differences.csv": tables[3]}
    if audited:
        files["audited_fit_inventory.csv"] = [{"campaign": key.split("/")[0], "cell": key.split("/")[1], **counts}
                                              for key, counts in audit["by_cell"].items()]
    files["dataset_provenance.csv"] = provenance
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".summary-staging-", dir=a.output.parent) as temp:
        stage = Path(temp) / "published"
        stage.mkdir()
        for name, values in files.items():
            write_csv(stage / name, values)
        manifest = {"status": "SUMMARY_READY_FOR_INDEPENDENT_RECONCILIATION" if audited else "SUMMARY_WITHOUT_RECORD_AUDIT",
                    "scope": ("audited" if audited else "unaudited") + " fixed-budget decisions; no fitting",
                    "reporting_environment": {"python": platform.python_version(), "numpy": np.__version__},
                    "counts": {name: len(values) for name, values in files.items()},
                    "files": {name: sha(stage / name) for name in files},
                    "source_sha256": sha(__file__), "frozen_closure_sha256": CLOSURE,
                    "campaign_sha256": sha(a.runtime / "campaign.json"),
                    "historical_index_sha256": sha(a.runtime / "inputs/historical_reused_rows.jsonl"),
                    "audit_sources": {p.name: sha(p) for p in audit_files} if audited else {},
                    "arm_labels": NAMES, "groups": list(GROUPS), "pairs": list(PAIRS),
                    "endpoint_definitions": {"any_pair": "at least one of the three pairs rejected",
                                             "truly_affected_pair": "at least one truly changed pair rejected; NA under H0",
                                             "unchanged_pair": "at least one unchanged pair rejected; equals any_pair under H0; NA if no unchanged pair exists",
                                             "truth_error_g1_g2_g3": "normalized Kendall inversion fraction versus each group's true order",
                                             "pair_distance_1_2_3": "normalized Kendall inversion fraction between estimated group orders"},
                    "uncertainty": {"rates": "two-sided 95% Wilson intervals, no cross-endpoint multiplicity adjustment",
                                    "contrasts": "left minus right; 10,000 paired bootstrap dataset resamples; 2.5/97.5 percentile interval; NumPy default linear quantiles",
                                    "unit": "independent simulated dataset within a cell; original R retained; no pooling across cells",
                                    "decision_bootstrap_seed": "92600000 + cell_index*1000 + contrast_index*100 + endpoint_index; cell_index counts within each campaign",
                                    "continuous_bootstrap_seed": "92700000 + (100000 for robustness) + cell_index*1000 + contrast_index*100 + endpoint_index"},
                    "cautions": ["Truncated tails have bounds and exact decisions, never point P values.",
                                 "Successful, attempted and failed permutation counts are distinct; failures retain unknown slots inside B=599.",
                                 "Historical first-group-stopping rows are excluded from all primary rates and contrasts.",
                                 "Missing compact warning/cap fields remain unknown; consult audited_fit_inventory.csv scope.",
                                 "This manifest is not the independent reconciliation gate."]}
        (stage / "summary_manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
        require(not a.output.exists(), "output appeared while reporting; refusing replacement")
        os.rename(stage, a.output)
    print(json.dumps({"status": manifest["status"], "output": str(a.output), "counts": manifest["counts"]}))


if __name__ == "__main__":
    main()
