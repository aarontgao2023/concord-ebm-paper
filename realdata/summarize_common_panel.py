#!/usr/bin/env python3
"""Six-event panel, step 5: aggregate tables from the stored fit records (no fits, no network).

Reads the 804 result units (10,792 fit records) that runtime/run_common_panel.py wrote under
<formal-root>/{adni,nacc}/{n2,apoe,bootstrap}/ and writes, to a new folder, the aggregate tables
behind Figs 1c, 4b-d and 5, Tables 2-3, Supplementary Tables 4-5 and Supplementary Fig. 2:
composition_*.csv (controlled composition experiment, with 10,000 paired resamples of the 200
draws, seed SeedSequence([20260923, 88000000, cohort]), cohort 0 = ADNI, 1 = NACC), APOE_*.csv
(observed orderings and the pairwise within-diagnosis permutation tests), bootstrap_orderings.csv,
event_pair_order_summary.csv, rank_*.csv and reference_counts_weights_ESS.csv. No output holds
participant-level values.

  python summarize_common_panel.py --formal-root $CONCORD_WORK_DIR/formal \\
      --new-output $CONCORD_WORK_DIR/summary

Every unit, record count and test is checked before any table is written. The original run was
additionally gated on an independent audit report and a file-transfer manifest; pass
--formal-audit, --transport-report and --auditor to apply those checks as well.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import tempfile
import numpy as np

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
CAMPAIGN_SHA = "dc5511321ad49f4ad8a8a97abc9a738d3f955ed8e81cc174840e6d724cee5a95"
COHORTS = ["adni", "nacc"]
EVENTS = ["ABETA", "TAU", "PTAU", "MEM", "EXF", "LAN"]
GROUPS = ["e2", "e33", "e4"]
DX = ["CN", "MCI", "AD"]
ENGINES = ["likelihood_ebm", "separate", "shared", "concord"]
DESIGNS = ["matched", "different"]
PANELS = {"CSFcog6": EVENTS, "omit_ptau": [e for e in EVENTS if e != "PTAU"]}
PAIR_CODES = list(itertools.combinations(range(3), 2))
EVENT_PAIRS = list(itertools.combinations(range(6), 2))
CATEGORIES = ["concordant_change", "opposing_change", "ADNI_specific_change", "NACC_specific_change", "no_change_in_either"]
CI_DRAWS = 10000
CI_SEED_COMPONENTS = {cohort: [20260923, 88000000, i] for i, cohort in enumerate(COHORTS)}


def require(condition, message):
    if not bool(condition):
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def js(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False)


def planned_units():
    return {f"{c}/{kind}/{tail}": {"cohort": c, "mode": kind, "tail": tail,
                                  "records": {"n2": 8, "apoe": 1798, "bootstrap": 1}[kind]}
            for c in COHORTS for kind in ("n2", "apoe", "bootstrap")
            for tail in (list(PANELS) if kind == "apoe" else [f"rep_{i:04d}" for i in range(200)])}


def inventory_sources(formal_root, runtime=RUNTIME):
    """Default check: exactly the planned formal units, no preflight output, frozen campaign."""
    root = Path(formal_root).resolve()
    require(not any("preflight" in x.lower() or x == "candidate_only" for x in root.parts), "Preflight cannot supply formal summaries")
    expected = planned_units()
    actual = {str(p.parent.relative_to(root)) for name in ("manifest.json", "result.json", "fit_records.jsonl") for p in root.rglob(name)}
    require(actual == set(expected), "Returned unit inventory differs from the planned set")
    for unit in expected:
        for name in ("manifest.json", "result.json", "fit_records.jsonl"):
            require((root / unit / name).is_file(), "A planned unit lacks one of its three files")
    require(sha(Path(runtime)/"campaign.json") == CAMPAIGN_SHA, "Frozen campaign changed")
    return read(Path(runtime)/"campaign.json"), {
        "formal_root": str(root), "formal_units": 804, "formal_fit_records": 10792,
        "external_audit_and_transport_gates": "not applied (no --formal-audit given)"}


def validated_sources(formal_root, audit_path, transport_path, runtime=RUNTIME, auditor=None):
    """Fail before reading scientific records if either full-scope gate is absent."""
    root = Path(formal_root).resolve()
    require(not any("preflight" in x.lower() or x == "candidate_only" for x in root.parts), "Preflight cannot supply formal summaries")
    audit, transport = read(audit_path), read(transport_path)
    require(audit.get("audit_status") == "PASS_FORMAL_COMPLETE" and audit.get("formal_complete") is True and audit.get("mode") == "formal", "Full formal audit PASS required; partial/fixture evidence rejected")
    require(audit.get("issues") == [] and audit.get("missing_units") == [], "Audit has unresolved issues/missing units")
    require(audit.get("campaign_sha256") == CAMPAIGN_SHA and audit.get("preflight_records_in_formal_denominator") == 0, "Audit campaign/preflight scope mismatch")
    require(Path(audit["outputs_root"]).resolve() == root, "Audit belongs to different returned outputs")
    require(audit["audit_source_sha256"] == sha(auditor), "Audit source changed after recorded PASS; rerun audit")
    require(audit["counts"].get("result_units_checked") == 804 and audit["counts"].get("raw_records_checked") == 10792, "Audit did not cover 804 units/10792 records")
    require(audit["unit_inventory"].get("absent_units") == 0 and audit["unit_inventory"].get("expected_units") == 804, "Audit inventory incomplete")
    require(transport.get("status") == "PASS" and transport.get("mismatches") == [] and transport.get("unmapped") == [], "Transport PASS without mismatches/unmapped entries required")
    require(transport.get("formal_files", 0) >= 2412 and transport.get("remote_gate_exists") is True, "Transport evidence does not cover formal results and production gate")
    remote_manifest_path = Path(transport["remote_manifest"])
    if not remote_manifest_path.is_absolute():
        remote_manifest_path = Path(transport_path).resolve().parent / remote_manifest_path
    remote = read(remote_manifest_path)
    require(isinstance(remote.get("files"), dict), "Transport remote manifest lacks file fingerprints")
    if "remote_manifest_sha256" in transport:
        require(transport["remote_manifest_sha256"] == sha(remote_manifest_path), "Transport manifest hash mismatch")
    expected = planned_units()
    proof = {entry["unit"]: entry for entry in audit["evidence"]}
    require(len(audit["evidence"]) == 804 and set(proof) == set(expected), "Audited unit identities not exactly planned")
    actual = {str(p.parent.relative_to(root)) for name in ("manifest.json", "result.json", "fit_records.jsonl") for p in root.rglob(name)}
    require(actual == set(expected), "Returned unit inventory differs from audited expected set")
    for unit in expected:
        require(proof[unit]["raw_records"] == expected[unit]["records"], "Wrong per-unit raw record count")
        for name, field in [("manifest.json", "manifest_sha256"), ("fit_records.jsonl", "records_sha256"), ("result.json", "result_sha256")]:
            file = root / unit / name; relative = f"formal/{unit}/{name}"
            require(relative in remote["files"], "Formal core file absent from remote transport manifest")
            h = sha(file); witness = remote["files"][relative]
            require(h == proof[unit][field] == witness["sha256"] and file.stat().st_size == witness["bytes"], "File changed since formal audit or lossless return")
    require(sha(Path(runtime)/"campaign.json") == CAMPAIGN_SHA, "Frozen campaign changed")
    return read(Path(runtime)/"campaign.json"), {
        "audit_path": str(Path(audit_path).resolve()), "audit_sha256": sha(audit_path),
        "transport_path": str(Path(transport_path).resolve()), "transport_sha256": sha(transport_path),
        "remote_manifest_path": str(remote_manifest_path), "remote_manifest_sha256": sha(remote_manifest_path),
        "formal_root": str(root), "verified_core_files": 2412, "formal_units": 804, "formal_fit_records": 10792}


def positions(order):
    require(sorted(order) == list(range(len(order))), "Invalid event permutation")
    result = np.empty(len(order), dtype=int)
    for rank, event in enumerate(order, 1):
        result[event] = rank
    return result


def distance(a, b):
    x, y = positions(a), positions(b)
    return sum((x[i] < x[j]) != (y[i] < y[j]) for i, j in itertools.combinations(range(len(a)), 2)) / (len(a)*(len(a)-1)/2)


def describe(values):
    x = np.asarray(values, dtype=float)
    require(x.ndim == 1 and len(x) > 0 and np.isfinite(x).all(), "Missing values cannot be silently removed")
    lo, median, hi = np.quantile(x, [.025, .5, .975], method="linear")
    return {"n": len(x), "mean": float(x.mean()), "sd": float(x.std(ddof=1)) if len(x)>1 else 0.,
            "median": float(median), "distribution_p025": float(lo), "distribution_p975": float(hi)}


def paired_mean_intervals(matrix, cohort, draws=CI_DRAWS):
    """Same row-resample indices for every planned contrast within each cohort."""
    x = np.asarray(matrix, dtype=float)
    require(x.ndim == 2 and np.isfinite(x).all(), "Paired bootstrap requires complete common rows")
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(CI_SEED_COMPONENTS[cohort])))
    means = np.empty((draws, x.shape[1]))
    for start in range(0, draws, 256):
        stop = min(start + 256, draws)
        indices = rng.integers(0, len(x), size=(stop-start, len(x)))
        means[start:stop] = x[indices].mean(axis=1)
    interval = np.quantile(means, [.025, .975], axis=0, method="linear")
    return interval[0], interval[1]


def composition_tables(orders, r_expected=200, ci_draws=CI_DRAWS):
    tables = {name: [] for name in ["composition_fits", "composition_event_positions", "composition_distance_summary", "composition_paired_contrasts", "composition_paired_contrast_summary"]}
    for cohort in COHORTS:
        array = np.asarray(orders[cohort], dtype=int)
        require(array.shape == (r_expected, 2, 4, 2, 6), "Composition order tensor incomplete")
        ds = np.zeros((r_expected, 2, 4))
        for rep, di, ei in itertools.product(range(r_expected), range(2), range(4)):
            a, b = array[rep, di, ei]
            ds[rep, di, ei] = distance(a, b)
            base = {"cohort": cohort, "replicate": rep, "design": DESIGNS[di], "engine": ENGINES[ei]}
            tables["composition_fits"].append({**base, "planned_replicates": r_expected, "group_A_n": 160, "group_B_n": 160,
                "status": "ok", "normalized_kendall": ds[rep, di, ei], "group_A_event_order": js([EVENTS[i] for i in a]), "group_B_event_order": js([EVENTS[i] for i in b])})
            for gi, order in enumerate((a, b)):
                for event, rank in zip(EVENTS, positions(order)):
                    tables["composition_event_positions"].append({**base, "group": "A" if gi == 0 else "B", "event": event, "position_1_based": int(rank)})
        series, descriptions = [], []
        for di, ei in itertools.product(range(2), range(4)):
            series.append(ds[:, di, ei]); descriptions.append({"kind": "absolute_distance", "design": DESIGNS[di], "engine_left": ENGINES[ei], "engine_right": ""})
        for ei in range(4):
            series.append(ds[:, 1, ei] - ds[:, 0, ei]); descriptions.append({"kind": "different_minus_matched", "design": "paired_designs", "engine_left": ENGINES[ei], "engine_right": ENGINES[ei]})
        for di in range(2):
            for a, b in itertools.combinations(range(4), 2):
                series.append(ds[:, di, a] - ds[:, di, b]); descriptions.append({"kind": "engine_left_minus_right", "design": DESIGNS[di], "engine_left": ENGINES[a], "engine_right": ENGINES[b]})
        matrix = np.column_stack(series); low, high = paired_mean_intervals(matrix, cohort, ci_draws)
        for i, (values, meta) in enumerate(zip(series, descriptions)):
            row = {"cohort": cohort, **meta, "planned_replicates": r_expected, "complete_replicates": r_expected,
                   **describe(values), "mean_ci_low": float(low[i]), "mean_ci_high": float(high[i]),
                   "ci_method": "paired-replicate bootstrap, pointwise 95% percentile", "bootstrap_draws": ci_draws,
                   "seed_components": js(CI_SEED_COMPONENTS[cohort]), "scope": "conditional_on_fixed_empirical_source_pool"}
            if meta["kind"] == "absolute_distance":
                tables["composition_distance_summary"].append(row)
            else:
                tables["composition_paired_contrast_summary"].append(row)
                for rep, value in enumerate(values):
                    tables["composition_paired_contrasts"].append({"cohort": cohort, "replicate": rep, **meta, "paired_difference": float(value)})
    return tables


def sign_distribution(values, observed):
    x = np.asarray(values)
    return {"negative_count": int((x<0).sum()), "zero_count": int((x==0).sum()), "positive_count": int((x>0).sum()),
            "negative_fraction": float((x<0).mean()), "zero_fraction": float((x==0).mean()), "positive_fraction": float((x>0).mean()),
            "observed_direction_support": float((np.sign(x)==np.sign(observed)).mean())}


def category(adni_change, nacc_change):
    a, b = np.sign(adni_change), np.sign(nacc_change)
    if a and b:
        return "concordant_change" if a == b else "opposing_change"
    if a:
        return "ADNI_specific_change"
    if b:
        return "NACC_specific_change"
    return "no_change_in_either"


def stability_tables(observed, boot_orders, b_expected=200):
    names = ["bootstrap_orderings", "bootstrap_event_positions", "rank_distribution", "rank_summary",
             "bootstrap_genotype_position_contrasts", "genotype_position_contrast_summary",
             "bootstrap_event_pair_orders", "event_pair_order_summary", "bootstrap_genotype_event_pair_changes",
             "genotype_event_pair_change_summary", "cross_cohort_position_contrasts", "cross_cohort_event_pair_changes", "cross_cohort_category_counts"]
    tables = {name: [] for name in names}
    position_stats, event_stats = {}, {}
    for cohort in COHORTS:
        obs_orders = np.asarray(observed[cohort], int); boots = np.asarray(boot_orders[cohort], int)
        require(obs_orders.shape == (3,6) and boots.shape == (b_expected,3,6), "Primary-six stability tensor incomplete")
        obs = np.stack([positions(order) for order in obs_orders])
        ranks = np.array([[positions(order) for order in replicate] for replicate in boots])
        for rep, gi in itertools.product(range(b_expected), range(3)):
            base = {"cohort": cohort, "panel": "CSFcog6", "replicate": rep, "genotype": GROUPS[gi]}
            tables["bootstrap_orderings"].append({**base, "event_order": js([EVENTS[i] for i in boots[rep,gi]]), "status": "ok"})
            for ei, event in enumerate(EVENTS):
                tables["bootstrap_event_positions"].append({**base, "event": event, "position_1_based": int(ranks[rep,gi,ei])})
        for gi, ei in itertools.product(range(3), range(6)):
            values = ranks[:,gi,ei]; base = {"cohort": cohort, "panel": "CSFcog6", "genotype": GROUPS[gi], "event": EVENTS[ei]}
            tables["rank_summary"].append({**base, "observed_position": int(obs[gi,ei]), "planned_bootstraps": b_expected, **describe(values),
                                           "scope": "conditional_rank_stability_fixed_preprocessing"})
            for rank in range(1,7):
                count = int((values==rank).sum())
                tables["rank_distribution"].append({**base, "rank": rank, "count": count, "denominator": b_expected, "fraction": count/b_expected})
        for a,b in PAIR_CODES:
            for ei,event in enumerate(EVENTS):
                delta = ranks[:,b,ei]-ranks[:,a,ei]; obsdelta = int(obs[b,ei]-obs[a,ei])
                base = {"cohort": cohort, "panel": "CSFcog6", "group_A": GROUPS[a], "group_B": GROUPS[b], "event": event}
                for rep,value in enumerate(delta):
                    tables["bootstrap_genotype_position_contrasts"].append({**base, "replicate": rep, "position_B_minus_A": int(value)})
                row = {**base, "observed_position_B_minus_A": obsdelta, "planned_bootstraps": b_expected, **describe(delta), **sign_distribution(delta,obsdelta)}
                tables["genotype_position_contrast_summary"].append(row); position_stats[cohort,a,b,ei] = row
        for i,j in EVENT_PAIRS:
            before = ranks[:,:,i] < ranks[:,:,j]; obsbefore = obs[:,i] < obs[:,j]
            for gi in range(3):
                base = {"cohort": cohort, "panel": "CSFcog6", "genotype": GROUPS[gi], "event_A": EVENTS[i], "event_B": EVENTS[j]}
                tables["event_pair_order_summary"].append({**base, "observed_A_before_B": int(obsbefore[gi]), "A_before_B_count": int(before[:,gi].sum()),
                    "B_before_A_count": int((~before[:,gi]).sum()), "denominator": b_expected, "A_before_B_fraction": float(before[:,gi].mean())})
                for rep,value in enumerate(before[:,gi]):
                    tables["bootstrap_event_pair_orders"].append({**base, "replicate":rep, "A_before_B":int(value)})
            for a,b in PAIR_CODES:
                delta = before[:,b].astype(int)-before[:,a].astype(int); obsdelta=int(obsbefore[b])-int(obsbefore[a])
                base = {"cohort":cohort,"panel":"CSFcog6","group_A":GROUPS[a],"group_B":GROUPS[b],"event_A":EVENTS[i],"event_B":EVENTS[j]}
                for rep,value in enumerate(delta):
                    tables["bootstrap_genotype_event_pair_changes"].append({**base,"replicate":rep,"before_indicator_B_minus_A":int(value)})
                row = {**base,"observed_before_in_group_A":int(obsbefore[a]),"observed_before_in_group_B":int(obsbefore[b]),
                       "observed_before_indicator_B_minus_A":obsdelta,"planned_bootstraps":b_expected,**sign_distribution(delta,obsdelta),
                       "support_kind":"descriptive_bootstrap_frequency_not_event_pair_P"}
                tables["genotype_event_pair_change_summary"].append(row); event_stats[cohort,a,b,i,j]=row
    for a,b in PAIR_CODES:
        for ei,event in enumerate(EVENTS):
            left,right=position_stats["adni",a,b,ei],position_stats["nacc",a,b,ei]
            row={"group_A":GROUPS[a],"group_B":GROUPS[b],"event":event,
                 "ADNI_observed_position_B_minus_A":left["observed_position_B_minus_A"],"NACC_observed_position_B_minus_A":right["observed_position_B_minus_A"],
                 "category":category(left["observed_position_B_minus_A"],right["observed_position_B_minus_A"])}
            for prefix,source in [("ADNI",left),("NACC",right)]:
                row.update({prefix+"_"+key:source[key] for key in ["planned_bootstraps","negative_fraction","zero_fraction","positive_fraction","observed_direction_support","distribution_p025","distribution_p975"]})
            tables["cross_cohort_position_contrasts"].append(row)
        for i,j in EVENT_PAIRS:
            left,right=event_stats["adni",a,b,i,j],event_stats["nacc",a,b,i,j]
            row={"group_A":GROUPS[a],"group_B":GROUPS[b],"event_A":EVENTS[i],"event_B":EVENTS[j],
                 "ADNI_observed_change":left["observed_before_indicator_B_minus_A"],"NACC_observed_change":right["observed_before_indicator_B_minus_A"],
                 "category":category(left["observed_before_indicator_B_minus_A"],right["observed_before_indicator_B_minus_A"])}
            for prefix,source in [("ADNI",left),("NACC",right)]:
                row.update({prefix+"_"+key:source[key] for key in ["planned_bootstraps","negative_count","zero_count","positive_count","negative_fraction","zero_fraction","positive_fraction","observed_direction_support"]})
            row["cross_cohort_scope"]="same_named_events_and_genotype_contrast; cohort_bootstrap_indices_not_paired"
            tables["cross_cohort_event_pair_changes"].append(row)
        for cat in CATEGORIES:
            count=sum(row["category"]==cat and row["group_A"]==GROUPS[a] and row["group_B"]==GROUPS[b] for row in tables["cross_cohort_event_pair_changes"])
            tables["cross_cohort_category_counts"].append({"group_A":GROUPS[a],"group_B":GROUPS[b],"category":cat,"event_pair_count":count,"planned_event_pairs":15})
    return tables


def inference_row(cohort,panel,pair_index,observed,permutations,saved_test):
    a,b=PAIR_CODES[pair_index]; key=f"{a}-{b}"
    values=[row["pair_distances"][key] for row in permutations]
    require(len(values)==599 and all(row["status"]=="ok" for row in permutations), "Incomplete permutation family; denominator cannot change")
    observed_d=observed["pair_distances"][key]
    e=sum(value>=observed_d for value in values); p=(1+e)/600; adj=min(1.,3*p)
    require(saved_test["budget"]==599 and saved_test["completed"]==599 and saved_test["failed"]==0, "Saved full-tail budget mismatch")
    require(saved_test["exceedances"]==e and saved_test["p"]==p and saved_test["adjusted_p"]==adj and saved_test["p_bounds"]==[p,p], "Independent inference arithmetic differs")
    require(saved_test["reject"]==(e<=9), "Inclusive Bonferroni rejection boundary differs")
    return {"cohort":cohort,"panel":panel,"analysis_role":"primary" if panel=="CSFcog6" else "same_person_sensitivity",
            "group_A":GROUPS[a],"group_B":GROUPS[b],"K":len(PANELS[panel]),"normalized_kendall":observed_d,
            "permutations_planned":599,"permutations_completed":599,"failed_permutations":0,"inclusive_exceedances":e,
            "raw_P":p,"Bonferroni_factor":3,"adjusted_P":adj,"reject_alpha_0_05":e<=9,
            "family":"three_APOE_pairs_within_this_cohort_and_panel","threshold":"0.05/3"}


def load_formal(root,cfg):
    """Called only after validated_sources passes. Never reads preflight results."""
    tables={name:[] for name in ["fit_record_ledger","campaign_coverage","APOE_pair_tests","APOE_observed_orderings","APOE_observed_positions","APOE_permutation_records"]}
    comp={c:np.full((200,2,4,2,6),-1,int) for c in COHORTS}
    boots={c:np.full((200,3,6),-1,int) for c in COHORTS}; obsprimary={}
    totals=Counter()
    for unit,info in planned_units().items():
        cohort,mode,tail=info["cohort"],info["mode"],info["tail"]
        folder=Path(root)/unit; manifest=read(folder/"manifest.json"); result=read(folder/"result.json")
        require(manifest["preflight"] is False and manifest["B"]==599 and result["complete"] is True, "Incomplete/preflight unit cannot be summarised")
        records={}
        with (folder/"fit_records.jsonl").open() as handle:
            for line in handle:
                require(line.endswith("\n"),"Truncated source record")
                row=json.loads(line)
                require(row["key"] not in records,"Duplicate fit key")
                records[row["key"]]=row
        require(len(records)==info["records"] and all(row["status"]=="ok" for row in records.values()),"Failed/missing fit retained; summary blocked rather than dropped")
        for key,row in records.items():
            d=row["diagnostics"]; groups=row["group_values"]
            require(len(row["orderings"])==len(groups),"Group order count")
            for a,b in itertools.combinations(range(len(groups)),2):
                require(distance(row["orderings"][a],row["orderings"][b])==row["pair_distances"][f"{groups[a]}-{groups[b]}"],"Distance-to-order reconciliation")
            m=row.get("measurement",{}); opt=d.get("optimizer",[])
            tables["fit_record_ledger"].append({"unit":unit,"cohort":cohort,"mode":mode,"record_key":key,"engine":row["engine"],
                "kind":row["kind"],"status":row["status"],"fit_seed":row["fit_seed"],"input_frame_sha256":row["input_frame_sha256"],
                "labels_sha256":row["permuted_labels_sha256"],"measurement_termination":m.get("termination","not_recorded_for_likelihood_engine"),
                "outer_iterations":m.get("outer_iterations",""),"optimizer_calls":len(opt) if row["engine"]!="likelihood_ebm" else "not_recorded",
                "optimizer_failures":d.get("optimizer_failures","not_recorded"),"quality_flags":js(d.get("quality_flags",[])),
                "warnings_count":d.get("warnings_count","not_recorded"),"fit_seconds":row["fit_seconds"]})
        totals[cohort,mode,"units"]+=1; totals[cohort,mode,"records"]+=len(records)
        if mode=="n2":
            rep=int(tail[4:]); expected={f"{de}/{en}" for de in DESIGNS for en in ENGINES}
            require(set(records)==expected,"Composition record identities incomplete")
            for di,de in enumerate(DESIGNS):
                for ei,en in enumerate(ENGINES):
                    comp[cohort][rep,di,ei]=records[f"{de}/{en}"]["orderings"]
        elif mode=="bootstrap":
            require(set(records)=={"bootstrap"},"Bootstrap record identity")
            boots[cohort][int(tail[4:])]=records["bootstrap"]["orderings"]
        else:
            panel=tail; names=PANELS[panel]; obs=records["observed"]
            require(set(records)=={"observed"}|{f"pair_{pi}/{p}" for pi in range(3) for p in range(599)},"APOE expected pair/permutation indices")
            require(obs["biomarkers"]==names,"APOE event mapping")
            if panel=="CSFcog6":
                obsprimary[cohort]=obs["orderings"]
            for gi,g in enumerate(GROUPS):
                expected_counts=cfg["cohorts"][cohort]["counts"][g]
                require(obs["group_diagnosis_counts"][str(gi)]==expected_counts,"APOE panel changed population cells")
                base={"cohort":cohort,"panel":panel,"analysis_role":"primary" if panel=="CSFcog6" else "same_person_sensitivity","genotype":g}
                tables["APOE_observed_orderings"].append({**base,"n":sum(expected_counts),"event_order":js([names[i] for i in obs["orderings"][gi]])})
                for event,rank in zip(names,positions(obs["orderings"][gi])):
                    tables["APOE_observed_positions"].append({**base,"event":event,"position_1_based":int(rank)})
            for pi,(a,b) in enumerate(PAIR_CODES):
                permutations=[records[f"pair_{pi}/{p}"] for p in range(599)]
                tables["APOE_pair_tests"].append(inference_row(cohort,panel,pi,obs,permutations,result["tests"][f"{a}-{b}"]))
                for pid,row in enumerate(permutations):
                    tables["APOE_permutation_records"].append({"cohort":cohort,"panel":panel,"analysis_role":"primary" if panel=="CSFcog6" else "same_person_sensitivity",
                        "shuffled_group_A":GROUPS[a],"shuffled_group_B":GROUPS[b],"permutation_index":pid,"status":row["status"],
                        "distance_e2_e33":row["pair_distances"]["0-1"],"distance_e2_e4":row["pair_distances"]["0-2"],"distance_e33_e4":row["pair_distances"]["1-2"],
                        "order_e2":js([names[i] for i in row["orderings"][0]]),"order_e33":js([names[i] for i in row["orderings"][1]]),"order_e4":js([names[i] for i in row["orderings"][2]])})
    for cohort in COHORTS:
        for mode,expected_units,expected_records in [("n2",200,1600),("apoe",2,3596),("bootstrap",200,200)]:
            require(totals[cohort,mode,"units"]==expected_units and totals[cohort,mode,"records"]==expected_records,"Campaign coverage mismatch")
            tables["campaign_coverage"].append({"cohort":cohort,"mode":mode,"planned_units":expected_units,"completed_units":totals[cohort,mode,"units"],
                "planned_fit_records":expected_records,"recorded_fit_records":totals[cohort,mode,"records"],"failed_records":0,"unresolved_caps_or_diagnostic_flags":0,"preflight_records_included":0})
    return tables,comp,obsprimary,boots


def reference_tables(cfg):
    rows,draws=[],[]
    q=np.array(cfg["reference"])
    for cohort in COHORTS:
        entries=[("actual_APOE",g,cfg["cohorts"][cohort]["counts"][g]) for g in GROUPS]
        entries += [(design,"A" if gi==0 else "B",counts) for design,pair in [("matched",[[80,48,32],[80,48,32]]),("different",[[112,32,16],[48,64,48]])] for gi,counts in enumerate(pair)]
        for design,g,counts in entries:
            cc=np.array(counts); n=int(cc.sum()); prop=cc/n; w=q/prop; ess=float(np.dot(cc,w)**2/np.dot(cc,w*w))
            for j,d in enumerate(DX):
                rows.append({"cohort":cohort,"application":design,"group":g,"diagnosis":d,"cell_n":int(cc[j]),"group_n":n,"original_proportion":float(prop[j]),
                             "reference_proportion":float(q[j]),"standardization_weight":float(w[j]),"group_ESS":ess,"ESS_over_n":ess/n})
        for j,d in enumerate(DX):
            available=cfg["cohorts"][cohort]["counts"]["e33"][j]; selected=[160,96,64][j]
            draws.append({"cohort":cohort,"source_genotype":"e33","diagnosis":d,"source_pool_n":available,"selected_per_replicate":selected,
                          "selected_source_fraction":selected/available,"same_selected_people_both_designs":True,"replicates_planned":200})
    return {"reference_counts_weights_ESS":rows,"composition_sampling_support":draws}


def expected_row_counts(r=200,b=200):
    return {"composition_fits":16*r,"composition_event_positions":192*r,"composition_distance_summary":16,"composition_paired_contrasts":32*r,
        "composition_paired_contrast_summary":32,"bootstrap_orderings":6*b,"bootstrap_event_positions":36*b,"rank_distribution":216,"rank_summary":36,
        "bootstrap_genotype_position_contrasts":36*b,"genotype_position_contrast_summary":36,"bootstrap_event_pair_orders":90*b,"event_pair_order_summary":90,
        "bootstrap_genotype_event_pair_changes":90*b,"genotype_event_pair_change_summary":90,"cross_cohort_position_contrasts":18,
        "cross_cohort_event_pair_changes":45,"cross_cohort_category_counts":15,"fit_record_ledger":10792,"campaign_coverage":6,"APOE_pair_tests":12,
        "APOE_observed_orderings":12,"APOE_observed_positions":66,"APOE_permutation_records":7188,"reference_counts_weights_ESS":42,"composition_sampling_support":6}


def reconcile(tables,r=200,b=200,partial=False):
    expected=expected_row_counts(r,b)
    require(partial or set(tables)==set(expected),"Missing/unplanned summary table")
    for name,rows in tables.items():
        require(len(rows)==expected[name],f"Unexpected row denominator: {name}")
    if "rank_distribution" in tables:
        sums=Counter()
        for row in tables["rank_distribution"]:
            require(row["denominator"]==b and 1<=row["rank"]<=6,"Rank denominator/scale")
            sums[row["cohort"],row["genotype"],row["event"]]+=row["count"]
        require(len(sums)==36 and set(sums.values())=={b},"Rank distributions do not sum to B200")
    if "genotype_event_pair_change_summary" in tables:
        for row in tables["genotype_event_pair_change_summary"]:
            require(row["negative_count"]+row["zero_count"]+row["positive_count"]==b,"Event-pair directional support denominator")
    if "cross_cohort_category_counts" in tables:
        sums=Counter()
        for row in tables["cross_cohort_category_counts"]:
            sums[row["group_A"],row["group_B"]]+=row["event_pair_count"]
        require(len(sums)==3 and set(sums.values())=={15},"Cross-cohort categories omit event pairs")
    if "composition_paired_contrasts" in tables:
        original={(x["cohort"],x["replicate"],x["design"],x["engine"]):x["normalized_kendall"] for x in tables["composition_fits"]}
        for row in tables["composition_paired_contrasts"]:
            c,rep,left,right=row["cohort"],row["replicate"],row["engine_left"],row["engine_right"]
            direct=original[c,rep,"different",left]-original[c,rep,"matched",left] if row["kind"]=="different_minus_matched" else original[c,rep,row["design"],left]-original[c,rep,row["design"],right]
            require(direct==row["paired_difference"],"Independent paired contrast reconciliation")
    if "APOE_pair_tests" in tables:
        for row in tables["APOE_pair_tests"]:
            require(row["permutations_completed"]==row["permutations_planned"]==599 and row["adjusted_P"]==min(1.,3*((row["inclusive_exceedances"]+1)/600)),"Full-tail table denominator arithmetic")
    return {"status":"INTERNAL_AGGREGATE_CHECKS_PASS","table_row_counts":{k:len(v) for k,v in tables.items()},"total_output_rows":sum(map(len,tables.values()))}


README = """# Six-event panel: aggregate tables

Produced only after every planned unit was found (804 units; 10,792 fits). When
an external audit report and transport manifest were supplied, their hashes are
bound in summary_manifest.json. No preflight output enters these tables. Failed,
missing or diagnostically unresolved records block publication; none is dropped
or replaced. This release has the original R=200 and B=599 denominators.
The summary builder's row-count/arithmetic checks are internal consistency
checks.

Composition event positions are 1-based; both group orderings and all four
models are retained. different_minus_matched is the difference in their
between-group normalized Kendall distances, not a difference in P values.
engine_left_minus_right follows the column names exactly. Mean intervals use
10,000 paired-replicate bootstrap resamples, the same replicate indices for all
planned comparisons within each cohort, and pointwise 2.5/97.5 percentiles.
Seeds are SeedSequence([20260923,88000000,cohort_index]), index 0=ADNI, 1=NACC.
These intervals describe the drawing experiment conditional on its empirical
source pool; they are not cohort-population sampling intervals or simultaneous
confidence bands. Distribution percentiles separately describe replicate spread.

APOE tests retain all three pairwise within-diagnosis permutation tests (599
relabelings each), inclusive exceedances, plus-one P and Bonferroni factor three
separately within each cohort/panel.
CSFcog6 is primary; omit_ptau is a same-person sensitivity and is never selected
in place of the primary by significance. Their event positions use their own K.

Rank distributions use the 200 bootstrap refits of the six-event panel (resampled
within the nine genotype-by-diagnosis cells), with preprocessing and common proportions fixed. Their 2.5/97.5 percentiles describe
the empirical bootstrap distribution; they are not another 10,000 model fits.
Genotype position contrast is rank(group B)-rank(group A): positive means later
in group B. An event-pair indicator is 1 if event A precedes event B; its genotype
contrast is indicator(group B)-indicator(group A), taking -1,0,+1. Every one of
the 15 named event pairs is retained for every genotype contrast and both cohorts.
Observed-direction support is the fraction of bootstrap contrasts with the same
sign as the observed contrast; for an observed zero it is the zero fraction.
These supports are descriptive frequencies, never event-level P values.

Cross-cohort categories compare observed contrast signs for the same named events
and genotype groups. Concordant/opposing require nonzero changes in both cohorts;
cohort-specific means one observed contrast is zero, not evidence of statistical
effect heterogeneity. Both-zero is reported explicitly. Support stays separate
by cohort. Bootstrap replicate numbers across cohorts are never paired and no
joint replication P value is computed. Concordant estimates or two small global
P values alone do not establish biological replication; nonrejection is not
equivalence. Acquisition, assay and population-selection limitations remain;
these descriptive tables do not remove them.

All output consists of aggregate source data and fit-level analysis identities,
not participant-level profiles or original participant identifiers.
"""


def publish_new(output,tables,evidence,cfg):
    output=Path(output).resolve()
    require(not output.exists(),"Refusing to overwrite an existing summary release")
    output.parent.mkdir(parents=True,exist_ok=True)
    stage=Path(tempfile.mkdtemp(prefix="."+output.name+".staging-",dir=output.parent))
    try:
        reconciliation=reconcile(tables)
        hashes={}
        for name,rows in tables.items():
            path=stage/(name+".csv")
            fields=list(rows[0])
            require(all(set(row)==set(fields) for row in rows),"Inconsistent CSV fields")
            with path.open("w",newline="") as handle:
                writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
            with path.open(newline="") as handle:
                require(sum(1 for _ in csv.DictReader(handle))==len(rows),"CSV round-trip row count")
            hashes[path.name]={"sha256":sha(path),"rows":len(rows)}
        (stage/"README.md").write_text(README)
        manifest={"status":"SUMMARY_READY_FOR_INDEPENDENT_RECONCILIATION","created_utc":datetime.now(timezone.utc).isoformat(),"source_evidence":evidence,
            "summary_script_sha256":sha(__file__),"campaign_sha256":CAMPAIGN_SHA,"reference":cfg["reference"],
            "paired_mean_bootstrap_draws":CI_DRAWS,"paired_mean_bootstrap_seed_components":CI_SEED_COMPONENTS,
            "reconciliation":reconciliation,"output_files":hashes,"new_model_fits":0,"new_permutation_tests":0,
            "cohorts_pooled":False,"preflight_records_used":0,"failed_records_dropped":0,
            "rank_support_is_event_pair_P":False,"primary_selected_by_P":False}
        (stage/"summary_manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
        require(not output.exists(),"Output appeared concurrently; refusing replacement")
        stage.rename(output)
    except Exception as exc:
        (stage/"FAILED_BUILD.json").write_text(json.dumps({"status":"UNPUBLISHED_FAILED_BUILD","error":str(exc)})+"\n")
        raise


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--formal-root",type=Path,required=True,help="folder written by run_common_panel.py (formal mode)")
    ap.add_argument("--formal-audit",type=Path,help="optional audit report; requires --transport-report and --auditor")
    ap.add_argument("--transport-report",type=Path,help="optional file-transfer manifest")
    ap.add_argument("--auditor",type=Path,help="optional audit script whose hash the audit report records")
    ap.add_argument("--runtime",type=Path,default=RUNTIME)
    ap.add_argument("--new-output",type=Path,required=True,help="new folder for the tables (must not exist)")
    args=ap.parse_args(argv)
    if args.formal_audit is None:
        cfg,evidence=inventory_sources(args.formal_root,args.runtime)
    else:
        if args.transport_report is None or args.auditor is None:
            ap.error("--formal-audit requires --transport-report and --auditor")
        cfg,evidence=validated_sources(args.formal_root,args.formal_audit,args.transport_report,args.runtime,args.auditor)
    tables,comp,observed,boots=load_formal(args.formal_root,cfg)
    tables.update(composition_tables(comp));tables.update(stability_tables(observed,boots));tables.update(reference_tables(cfg))
    publish_new(args.new_output,tables,evidence,cfg)
    print(json.dumps({"status":"SUMMARY_READY_FOR_INDEPENDENT_RECONCILIATION","output":str(args.new_output),"tables":len(tables),"new_fits":0}))


if __name__=="__main__":
    main()
