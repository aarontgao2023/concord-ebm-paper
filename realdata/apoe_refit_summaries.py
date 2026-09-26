#!/usr/bin/env python3
"""APOE bootstrap refits: summaries for Fig. 5e-g, the paired-refit statistic of the Results and
Supplementary Fig. 3 (standard library only).

No model is fitted and no permutation is run. Every quantity is computed from stored fit-level
orderings:

  six-event panel (from summarize_common_panel.py, folder --summary-dir)
    APOE_observed_orderings.csv, bootstrap_orderings.csv (200 refits per cohort, resampled within
    genotype-by-diagnosis cells; the three genotype orderings of a refit come from the same
    resample and the same pooled fit, so genotypes are paired within refit)
  earlier panels (Supplementary Fig. 3; optional)
    --adni-bootstrap  fig4_ADNI_bootstrap_orderings.csv, 500 refits of the ADNI 12-event panel
                      (earlier_panels/aggregate_earlier_panels.py)
    --nacc-bootstrap  r2_nacc_bootstrap_orderings.csv, 500 refits of the NACC amyloid PET panel
                      (earlier_panels/nacc/export_stability_source.py)
    --point-orders    optional table of observed orderings (columns cohort, reference, group,
                      method, event_order); fills the observed-rank column only

Outputs (folder --out):
  apoe_cp_block_by_group.csv          per genotype: S = number of the 9 (CSF event, cognitive
                                      score) pairs with the CSF event first, and related
                                      indicators, over the 200 refits
  apoe_cp_block_contrasts.csv         paired genotype differences within refit; the Results
                                      statistic is S(e4) - S(e3/e3) > 0 in 70% (ADNI) and 87%
                                      (NACC) of refits (statistic S_csf_before_cog_pairs_0to9,
                                      column frac_positive)
  apoe_cp_e4_prefix_frequencies.csv   Abeta42 < p-tau < total tau frequencies (Fig. 5f check)
  apoe_hist_amyloid_position.csv      amyloid event placed first, earlier panels (Supplementary
                                      Fig. 3 uses rows ADNI_K12/concord_min and NACC5_PET/min)

In the earlier-panel files, concord_min and min are CONCORD at the smallest-proportion rule,
standard_original is the separately fitted DEBM with the unmodified pyebm search and pooled is
CONCORD with pooled proportions (not reported). All statistics here were defined after the
observed orderings were seen; they are descriptive and are not prespecified endpoints.

  python apoe_refit_summaries.py --summary-dir $CONCORD_WORK_DIR/summary \\
      --adni-bootstrap $CONCORD_WORK_DIR/earlier/fig4_ADNI_bootstrap_orderings.csv \\
      --nacc-bootstrap $CONCORD_WORK_DIR/earlier/r2_nacc_bootstrap_orderings.csv \\
      --out $CONCORD_WORK_DIR/resummary
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path

CSF6 = ["ABETA", "TAU", "PTAU"]
COG = ["MEM", "EXF", "LAN"]
GENOS = ["e2", "e33", "e4"]
COHORTS = ["adni", "nacc"]
OUT = None   # set in main()


# ---------------------------------------------------------------- helpers
def pos(order):
    return {e: i + 1 for i, e in enumerate(order)}


def block_S(order, csf):
    """Number of (CSF, cognitive) cross pairs in which the CSF event precedes."""
    p = pos(order)
    return sum(1 for c in csf for k in COG if p[c] < p[k])


def block_S_nolan(order, csf):
    p = pos(order)
    return sum(1 for c in csf for k in ("MEM", "EXF") if p[c] < p[k])


def block_S_nomem(order, csf):
    p = pos(order)
    return sum(1 for c in csf for k in ("EXF", "LAN") if p[c] < p[k])


def csf_before_lan(order, csf):
    p = pos(order)
    return sum(1 for c in csf if p[c] < p["LAN"])


def mean_rank_gap(order, csf):
    p = pos(order)
    return sum(p[k] for k in COG) / 3.0 - sum(p[c] for c in csf) / len(csf)


def pct(vals, q):
    """Linear-interpolated percentile (same convention as numpy default)."""
    v = sorted(vals)
    if not v:
        return float("nan")
    h = (len(v) - 1) * q
    lo = math.floor(h)
    hi = math.ceil(h)
    return v[lo] + (v[hi] - v[lo]) * (h - lo)


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def write_csv(name, header, rows):
    path = os.path.join(OUT, name)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([round(x, 6) if isinstance(x, float) else x for x in r])
    return path


def load_long(path, arm_col, rep_col, grp_col, ev_col, pos_col):
    d = {}
    for r in csv.DictReader(open(path)):
        d.setdefault((r[arm_col], int(r[rep_col]), r[grp_col]), {})[r[ev_col]] = int(r[pos_col])
    return d


def main(argv=None):
    global OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-dir", type=Path, required=True,
                    help="folder with APOE_observed_orderings.csv and bootstrap_orderings.csv")
    ap.add_argument("--adni-bootstrap", type=Path, help="fig4_ADNI_bootstrap_orderings.csv (optional)")
    ap.add_argument("--nacc-bootstrap", type=Path, help="r2_nacc_bootstrap_orderings.csv (optional)")
    ap.add_argument("--point-orders", type=Path, help="observed orderings of the earlier panels (optional)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    OUT = str(args.out)
    CP = str(args.summary_dir)

    # ------------------------------------------------------------ load six-event panel
    obs = {}
    for r in csv.DictReader(open(os.path.join(CP, "APOE_observed_orderings.csv"))):
        obs[(r["cohort"], r["panel"], r["genotype"])] = json.loads(r["event_order"])

    boot = {}  # (cohort, replicate) -> {geno: order}
    for r in csv.DictReader(open(os.path.join(CP, "bootstrap_orderings.csv"))):
        assert r["status"] == "ok" and r["panel"] == "CSFcog6"
        boot.setdefault((r["cohort"], int(r["replicate"])), {})[r["genotype"]] = json.loads(r["event_order"])
    reps = {c: sorted(k[1] for k in boot if k[0] == c) for c in COHORTS}
    assert all(len(reps[c]) == 200 for c in COHORTS)
    assert all(len(boot[k]) == 3 for k in boot)

    # ------------------------------------------------------------ 1. per-group block summaries
    stats = {
        "S_csf_before_cog_pairs_0to9": lambda o: block_S(o, CSF6),
        "block_all_csf_first": lambda o: 1 if block_S(o, CSF6) == 9 else 0,
        "S_csf_before_MEM_EXF_0to6": lambda o: block_S_nolan(o, CSF6),
        "S_csf_before_EXF_LAN_0to6": lambda o: block_S_nomem(o, CSF6),
        "csf_events_before_LAN_0to3": lambda o: csf_before_lan(o, CSF6),
        "LAN_after_all_csf": lambda o: 1 if csf_before_lan(o, CSF6) == 3 else 0,
        "TAU_before_LAN": lambda o: 1 if pos(o)["TAU"] < pos(o)["LAN"] else 0,
        "mean_rank_gap_cog_minus_csf": lambda o: mean_rank_gap(o, CSF6),
        "ABETA_rank1": lambda o: 1 if o[0] == "ABETA" else 0,
        "MEM_rank6": lambda o: 1 if o[-1] == "MEM" else 0,
    }
    rows1 = []
    for c in COHORTS:
        for g in GENOS:
            o = obs[(c, "CSFcog6", g)]
            for sname, f in stats.items():
                vals = [f(boot[(c, b)][g]) for b in reps[c]]
                rows1.append([c, "CSFcog6", g, sname, f(o), len(vals), mean(vals),
                              pct(vals, 0.025), pct(vals, 0.975), min(vals), max(vals)])
    write_csv("apoe_cp_block_by_group.csv",
              ["cohort", "panel", "genotype", "statistic", "observed", "bootstrap_n",
               "bootstrap_mean", "bootstrap_p025", "bootstrap_p975", "bootstrap_min", "bootstrap_max"], rows1)

    # ------------------------------------------------------------ 2. paired genotype contrasts (B - A)
    contrast_stats = {k: stats[k] for k in [
        "S_csf_before_cog_pairs_0to9", "S_csf_before_MEM_EXF_0to6", "S_csf_before_EXF_LAN_0to6",
        "csf_events_before_LAN_0to3", "LAN_after_all_csf",
        "mean_rank_gap_cog_minus_csf", "block_all_csf_first", "TAU_before_LAN", "ABETA_rank1"]}
    rows2 = []
    for c in COHORTS:
        for gA, gB in [("e33", "e4"), ("e2", "e33"), ("e2", "e4")]:
            for sname, f in contrast_stats.items():
                od = f(obs[(c, "CSFcog6", gB)]) - f(obs[(c, "CSFcog6", gA)])
                d = [f(boot[(c, b)][gB]) - f(boot[(c, b)][gA]) for b in reps[c]]
                n = len(d)
                rows2.append([c, gA, gB, sname, od, n, mean(d), pct(d, 0.025), pct(d, 0.975),
                              sum(1 for x in d if x > 0) / n, sum(1 for x in d if x == 0) / n,
                              sum(1 for x in d if x < 0) / n])
    write_csv("apoe_cp_block_contrasts.csv",
              ["cohort", "group_A", "group_B", "statistic", "observed_B_minus_A", "bootstrap_n",
               "bootstrap_mean_B_minus_A", "bootstrap_p025", "bootstrap_p975",
               "frac_positive", "frac_zero", "frac_negative"], rows2)

    # ------------------------------------------------------------ 3. CSF-order frequencies (Fig. 5f)
    rows12 = []
    for c in COHORTS:
        for g in GENOS:
            o = obs[(c, "CSFcog6", g)]
            f_prefix = lambda x: 1 if x[:3] == ["ABETA", "PTAU", "TAU"] else 0
            f_csf_internal = lambda x: 1 if pos(x)["ABETA"] < pos(x)["PTAU"] < pos(x)["TAU"] else 0
            vals_p = [f_prefix(boot[(c, b)][g]) for b in reps[c]]
            vals_i = [f_csf_internal(boot[(c, b)][g]) for b in reps[c]]
            rows12.append([c, g, f_prefix(o), mean(vals_p), f_csf_internal(o), mean(vals_i), len(vals_p)])
    write_csv("apoe_cp_e4_prefix_frequencies.csv",
              ["cohort", "genotype", "observed_first3_is_ABETA_PTAU_TAU", "boot_frac_first3_is_ABETA_PTAU_TAU",
               "observed_csf_internal_ABETA_lt_PTAU_lt_TAU", "boot_frac_csf_internal_ABETA_lt_PTAU_lt_TAU", "bootstrap_n"], rows12)
    print("six-event panel summaries written")

    # ------------------------------------------------------------ 4. earlier panels: amyloid position
    if args.adni_bootstrap is None or args.nacc_bootstrap is None:
        print("earlier-panel inputs not given; apoe_hist_amyloid_position.csv not written")
        return 0
    adni_b = load_long(args.adni_bootstrap, "arm", "replicate", "group", "event", "position")
    nacc_b = load_long(args.nacc_bootstrap, "reference", "replicate", "group", "event", "position")
    # position base of each file (both are written 1-based); positions below are made 0-based
    nacc_min_pos = min(min(v.values()) for v in nacc_b.values())
    adni_min_pos = min(min(v.values()) for v in adni_b.values())

    point = {}
    if args.point_orders is not None:
        for r in csv.DictReader(open(args.point_orders)):
            if r["method"] == "exact":
                point[(r["cohort"], r["reference"], r["group"])] = json.loads(r["event_order"])

    rows8 = []
    for cohort, data, amy, K, base, fits in [
            ("ADNI_K12", adni_b, "ABETA", 12, adni_min_pos, ["concord_min", "standard_original"]),
            ("NACC5_PET", nacc_b, "CENTILOID", 5, nacc_min_pos, ["min", "pooled"])]:
        for fit in fits:
            rlist = sorted({k[1] for k in data if k[0] == fit})
            for g in GENOS:
                vals = [(data[(fit, r, g)][amy] - base) for r in rlist]  # 0-based
                nr = [v / (K - 1) for v in vals]
                ref_key = {"concord_min": "min", "min": "min", "pooled": "pooled"}.get(fit)
                obs_nr = ""
                if ref_key is not None and point:
                    ck = "adni" if cohort.startswith("ADNI") else "nacc"
                    oo = point[(ck, ref_key, g)]
                    obs_nr = oo.index(amy) / (K - 1)
                rows8.append([cohort, fit, g, amy, K, len(vals), obs_nr, mean(nr), pct(nr, 0.025), pct(nr, 0.975),
                              sum(1 for v in vals if v == 0) / len(vals),
                              sum(1 for v in vals if v == K - 1) / len(vals),
                              sum(1 for v in vals if v >= (K - 1) / 2) / len(vals)])
    write_csv("apoe_hist_amyloid_position.csv",
              ["cohort_panel", "arm_or_reference", "genotype", "amyloid_event", "K", "resamples",
               "observed_normalized_rank_0first_1last", "boot_mean_normalized_rank", "boot_p025", "boot_p975",
               "frac_amyloid_first", "frac_amyloid_last", "frac_amyloid_in_later_half"], rows8)
    print("earlier-panel amyloid positions written; position base ADNI", adni_min_pos, "NACC", nacc_min_pos)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
