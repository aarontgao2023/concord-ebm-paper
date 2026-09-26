#!/usr/bin/env python3
"""Supplementary Tables 1-2 (simulation) of the CONCORD manuscript.

Reads only simulation-only aggregate files in figures/inputs/simulation/ (no fitting):
  * S1_event_parameters.csv      the 14 simulated biomarkers (= design_v2.BIOMARKERS)
  * group_stage_parameters.csv   group sizes and stage distributions (= design_v2.REFERENCE_COUNTS
                                 and BASE_STAGE_BETA)
  * rates.csv                    paired benchmark rates

Writes into figures/output/tables/ (CONCORD_OUTPUT_DIR/tables):
  SuppTable1_generator.{tex,csv}
  SuppTable2_benchmark.{tex,csv}
  _check_make_si_sim_tables.tex   (compile check; \\input's both tables)

Captions and table notes are edited in the manuscript source (tables/*.tex there). The caption
and note text below is a copy of the submitted version, so that the output can be compared with
the manuscript file by file.

Run:  python figures/tables/make_si_sim_tables.py
"""
from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np

FIGDIR = Path(__file__).resolve().parents[1]  # figures/


def _env_dir(name: str, default: Path) -> Path:
    v = os.environ.get(name)
    return Path(v).expanduser() if v else default


SIM_INPUTS = _env_dir("CONCORD_FIG_INPUTS", FIGDIR / "inputs") / "simulation"
OUT = _env_dir("CONCORD_OUTPUT_DIR", FIGDIR / "output") / "tables"
EVENTS_CSV = SIM_INPUTS / "S1_event_parameters.csv"
STAGE_CSV = SIM_INPUTS / "group_stage_parameters.csv"
RATES_CSV = SIM_INPUTS / "rates.csv"

Z95 = 1.959963984540054  # two-sided 95% normal quantile for Wilson intervals


# --------------------------------------------------------------------------- helpers
def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def note(text: str) -> list[str]:
    """Table note below the tabular, as in the manuscript."""
    return [r"\par\smallskip", r"{\footnotesize\raggedright " + text + r"\par}"]


def write_csv(stem: str, header: list[str], rows: list[list[str]]) -> Path:
    path = OUT / f"{stem}.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    return path


def write_tex(stem: str, text: str) -> Path:
    path = OUT / f"{stem}.tex"
    assert all(ord(c) < 128 for c in text), f"non-ASCII character in {stem}.tex"
    path.write_text(text, encoding="utf-8")
    return path


def thousands(n: int) -> str:
    return f"{n:,}"


def wilson(k: int, n: int, z: float = Z95) -> tuple[float, float]:
    p = k / n
    den = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return float(centre - half), float(centre + half)


# =========================================================================== Table 1
# Generator code names in input order (design_v2.BIOMARKERS). The manuscript numbers the simulated
# biomarkers 1-14 in this order and does not use these names: 1-5 are modeled on CSF markers, 6-7 on
# cognitive scores and 8-14 on imaging measures.
SIM_CODES = ["ABETA", "PTAU", "TAU", "NG", "NFL", "ADAS13", "MMSE", "Hippocampus", "Entorhinal", "MidTemp",
             "Fusiform", "WholeBrain", "Ventricles", "Precuneus"]
MODALITY = {"CSF": "CSF", "COG": "Cognitive", "IMG": "Imaging"}   # column "Modeled on"

# Caption titles and notes: copies of the manuscript text (edited there, not here).
TITLE_S1 = "Simulated biomarkers"
NOTE_S1 = (
    "Each biomarker defines one event. Component AUC, area under the curve separating measurements "
    "before and after the event. Observation probability, probability that a measurement is observed "
    "in simulations with missing measurements, based on the biomarker availability of the ADNI "
    "comparison whose group composition the simulations follow. Values are the inputs of the "
    "simulation code. Group sizes, stage distributions and designs are given in Supplementary "
    "Methods 1. AUC, area under the curve; CSF, cerebrospinal fluid."
)

EXPECTED_EVENTS = {  # code: (component AUC, observation probability)
    "ABETA": (0.88, 0.735), "ADAS13": (0.92, 0.990), "Precuneus": (0.70, 0.990),
    "NG": (0.68, 0.275), "NFL": (0.72, 0.288),
}
EXPECTED_COUNTS = {"g1": (57, 6, 12), "g2": (244, 66, 101), "g3": (110, 156, 219)}
EXPECTED_BETA = {"CN": (0.6, 9.0), "MCI": (2.0, 2.0), "AD": (9.0, 0.6)}


def beta_str(a: float, b: float) -> str:
    if float(a).is_integer() and float(b).is_integer():
        return f"Beta({int(a)}, {int(b)})"
    return f"Beta({a:.1f}, {b:.1f})"


def make_table1() -> dict:
    events = read_csv(EVENTS_CSV)
    stage = read_csv(STAGE_CSV)

    # ---- checks on events
    codes = [r["event"] for r in events]
    assert len(codes) == 14 and len(set(codes)) == 14, codes
    assert codes == SIM_CODES, codes
    ev = {r["event"]: r for r in events}
    for code, (auc, obs) in EXPECTED_EVENTS.items():
        assert abs(float(ev[code]["component_AUC"]) - auc) < 1e-12, (code, ev[code])
        assert abs(float(ev[code]["observation_probability"]) - obs) < 1e-12, (code, ev[code])
    assert {r["modality"] for r in events} <= set(MODALITY)

    # ---- checks on groups and stage laws
    counts: dict[str, dict[str, int]] = {}
    beta: dict[str, tuple[float, float]] = {}
    for r in stage:
        counts.setdefault(r["group"], {})[r["diagnosis"]] = int(r["n"])
        ab = (float(r["stage_beta_a"]), float(r["stage_beta_b"]))
        assert beta.setdefault(r["diagnosis"], ab) == ab, "stage law differs across groups"
    for g, exp in EXPECTED_COUNTS.items():
        assert tuple(counts[g][d] for d in ("CN", "MCI", "AD")) == exp, (g, counts[g])
    for d, exp in EXPECTED_BETA.items():
        assert beta[d] == exp, (d, beta[d])
    group_n = {g: sum(counts[g].values()) for g in counts}
    assert group_n == {"g1": 75, "g2": 411, "g3": 485}, group_n
    n_total = sum(group_n.values())
    assert n_total == 971, n_total
    K = len(events)

    # ---- rows
    tex_rows, csv_rows = [], []
    for i, r in enumerate(events, start=1):
        name = f"Biomarker {i}"
        mod = MODALITY[r["modality"]]
        auc = f"{float(r['component_AUC']):.2f}"
        obs = f"{float(r['observation_probability']):.3f}"
        tex_rows.append(f"{name} & {mod} & {auc} & {obs} \\\\")
        csv_rows.append([name, mod, auc, obs])

    tex = "\n".join([
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{{TITLE_S1}}}",
        r"\label{supptab:generator}",
        r"\smallskip",
        r"\small",
        r"\begin{tabular}{@{}llrr@{}}",
        r"\toprule",
        r"Biomarker & Modeled on & Component AUC & Observation probability \\",
        r"\midrule",
        *tex_rows,
        r"\bottomrule",
        r"\end{tabular}",
        *note(NOTE_S1),
        r"\end{table}",
        "",
    ])
    header = ["Biomarker", "Modeled on", "Component AUC", "Observation probability"]
    p_tex = write_tex("SuppTable1_generator", tex)
    p_csv = write_csv("SuppTable1_generator", header, csv_rows)
    return {"paths": [p_tex, p_csv], "rows": csv_rows, "K": K, "n": n_total,
            "group_n": group_n, "beta": beta}


# =========================================================================== Table 2
METHODS = [  # model label in rates.csv -> column label
    ("repaired", "Separately fitted DEBM"),
    ("shared", "Pooled-score DEBM"),
    ("invariant_min", "CONCORD"),
]
MEASUREMENTS = {
    "gauss": ("Independent Gaussian noise", "Independent Gaussian noise"),
    "corr": ("Gaussian noise correlated within modality", "Gaussian noise correlated within modality"),
    "t5": ("Heavy-tailed noise (t distribution, 5 d.f.)", r"Heavy-tailed noise ($t$ distribution, 5 d.f.)"),
}
CHANGE = {
    "H0": "None",
    "EXACT27_G1": "27 of 91 pairs reversed, CN-heavy group",
    "EXACT27_G3": "27 of 91 pairs reversed, AD-heavy group",
    "SINGLE13_G3": "One event moved 13, AD-heavy group",
    "DISJOINT6_G1": "Six adjacent swaps, CN-heavy group",
    "SINGLE6_G1": "One event moved 6, CN-heavy group",
}
CONDITIONS = [  # (cell, measurement key, change key)
    ("GAUSSIAN_H0", "gauss", "H0"),
    ("EXACT27_G1", "gauss", "EXACT27_G1"),
    ("EXACT27_G3", "gauss", "EXACT27_G3"),
    ("SINGLE13_G3", "gauss", "SINGLE13_G3"),
    ("DISJOINT6_G1", "gauss", "DISJOINT6_G1"),
    ("CORR_H0", "corr", "H0"),
    ("CORR_SINGLE6_G1", "corr", "SINGLE6_G1"),
    ("CORR_SINGLE13_G3", "corr", "SINGLE13_G3"),
    ("T5_H0", "t5", "H0"),
    ("T5_SINGLE6_G1", "t5", "SINGLE6_G1"),
    ("T5_SINGLE13_G3", "t5", "SINGLE13_G3"),
]
EXPECTED_RATES = {  # cell: (separately fitted, pooled-score, CONCORD) in %
    "GAUSSIAN_H0": (3.4, 2.8, 4.2),
    "CORR_H0": (3.7, 3.3, 2.6),
    "T5_H0": (4.7, 1.9, 2.5),
    "EXACT27_G3": (13.6, 35.0, 54.2),
    "SINGLE13_G3": (8.2, 44.4, 54.2),
    "CORR_SINGLE13_G3": (8.8, 37.6, 51.0),
    "EXACT27_G1": (8.2, 18.4, 17.6),
}


def make_table2() -> dict:
    with RATES_CSV.open(newline="", encoding="utf-8") as fh:
        rates = list(csv.DictReader(fh))
    idx = {(r["cell"], r["arm"], r["endpoint"]): r for r in rates}
    assert {r["cell"] for r in rates} == {c for c, _, _ in CONDITIONS}
    assert {r["arm"] for r in rates} == {a for a, _ in METHODS}

    tex_rows, csv_rows, printed = [], [], {}
    current_meas = None
    for cell, meas, change in CONDITIONS:
        is_null = change == "H0"
        endpoint = "any_pair" if is_null else "truly_affected_pair"
        rate_label = "False positive" if is_null else "Detection"
        cells_txt, cells_tex, pct, Rs = [], [], [], set()
        for arm, _ in METHODS:
            r = idx[(cell, arm, endpoint)]
            k, R = int(r["rejections"]), int(r["planned_R"])
            Rs.add(R)
            assert abs(float(r["rate"]) - k / R) < 1e-12, (cell, arm)
            lo, hi = wilson(k, R)
            assert abs(lo - float(r["wilson_low"])) < 1e-9, (cell, arm, lo, r["wilson_low"])
            assert abs(hi - float(r["wilson_high"])) < 1e-9, (cell, arm, hi, r["wilson_high"])
            p = 100 * k / R
            pct.append(round(p, 1))
            cells_txt.append(f"{p:.1f} [{100 * lo:.1f}–{100 * hi:.1f}]")
            cells_tex.append(f"{p:.1f} [{100 * lo:.1f}--{100 * hi:.1f}]")
        assert len(Rs) == 1, (cell, Rs)  # same datasets for all methods
        R = Rs.pop()
        assert R == (1000 if is_null else 500), (cell, R)
        if cell in EXPECTED_RATES:
            assert tuple(pct) == EXPECTED_RATES[cell], (cell, pct, EXPECTED_RATES[cell])
        printed[cell] = (rate_label, R, cells_txt)

        meas_txt, meas_tex = MEASUREMENTS[meas]
        if meas != current_meas:
            if current_meas is not None:
                tex_rows.append(r"\addlinespace")
            tex_rows.append(f"\\multicolumn{{6}}{{@{{}}l}}{{{meas_tex}}} \\\\")
            current_meas = meas
        tex_rows.append(f"{CHANGE[change]} & {rate_label} & {thousands(R)} & "
                        + " & ".join(cells_tex) + r" \\")
        csv_rows.append([meas_txt, CHANGE[change], rate_label, thousands(R), *cells_txt])

    title = "False-positive and detection rates in the paired simulation benchmark"
    legend = (
        "Percentage of datasets in which at least one of the three group comparisons was rejected "
        "(conditions without a changed ordering) or at least one comparison involving the changed "
        "group was rejected (conditions with a changed ordering), with Wilson 95\\% intervals. All "
        "methods analyzed the same datasets with within-diagnosis permutation (599 relabelings; "
        "0.05/3 per comparison). Measurements were missing as in Supplementary Table 1. Changes are "
        "defined in Supplementary Methods 1. DEBM, discriminative event-based model; d.f., degrees "
        "of freedom."
    )

    def stack(*lines: str) -> str:
        return r"\begin{tabular}[b]{@{}r@{}}" + r"\\".join(lines) + r"\end{tabular}"

    tex = "\n".join([
        r"\begin{table}[htbp]",
        r"\centering",
        f"\\caption{{{title}}}",
        r"\label{supptab:benchmark}",
        r"\smallskip",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash\setlength{\leftskip}{1em}}Xlrrrr@{}}",
        r"\toprule",
        r" & & & \multicolumn{3}{c}{\% of datasets [Wilson 95\% interval]} \\",
        r"\cmidrule(l){4-6}",
        r"\multicolumn{1}{@{}l}{Changed sequence} & Rate & Datasets & "
        + stack("Separately", "fitted DEBM") + " & "
        + stack("Pooled-score", "DEBM") + r" & CONCORD \\",
        r"\midrule",
        *tex_rows,
        r"\bottomrule",
        r"\end{tabularx}",
        *note(legend),
        r"\end{table}",
        "",
    ])
    header = ["Measurement noise", "Changed sequence", "Rate", "Datasets",
              *[f"{label} (% [Wilson 95% interval])" for _, label in METHODS]]
    p_tex = write_tex("SuppTable2_benchmark", tex)
    p_csv = write_csv("SuppTable2_benchmark", header, csv_rows)
    return {"paths": [p_tex, p_csv], "rows": printed}


# =========================================================================== check doc
CHECK_TEX = r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{booktabs,tabularx,amsmath,array}
\usepackage[a4paper,margin=2.5cm]{geometry}
\renewcommand{\tablename}{Supplementary Table}
\begin{document}
\input{SuppTable1_generator.tex}
\input{SuppTable2_benchmark.tex}
\end{document}
"""


def main() -> None:
    missing = [p for p in (EVENTS_CSV, STAGE_CSV, RATES_CSV) if not p.is_file()]
    if missing:
        sys.exit("Input file(s) not found:\n" + "".join(f"  - {p}\n" for p in missing))
    OUT.mkdir(parents=True, exist_ok=True)
    t1 = make_table1()
    t2 = make_table2()
    (OUT / "_check_make_si_sim_tables.tex").write_text(CHECK_TEX, encoding="utf-8")

    print("Supplementary Table 1: K =", t1["K"], "n =", t1["n"], t1["group_n"],
          "beta =", {d: beta_str(*ab) for d, ab in t1["beta"].items()})
    for row in t1["rows"]:
        print("  ", " | ".join(row))
    print("Supplementary Table 2 (separately fitted / pooled-score / CONCORD):")
    for cell, (lab, R, cells) in t2["rows"].items():
        print(f"   {cell:17s} {lab:14s} R={R:5d}  " + " ; ".join(cells))
    for p in t1["paths"] + t2["paths"]:
        print("wrote", p)
    print("wrote", OUT / "_check_make_si_sim_tables.tex")


if __name__ == "__main__":
    main()
