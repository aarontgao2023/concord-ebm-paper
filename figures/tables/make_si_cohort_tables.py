"""Supplementary Tables 3-5 (cohorts and APOE) of the CONCORD manuscript.

Reads aggregate ADNI/NACC results (no participant-level data) that are not distributed with this
repository; place them in figures/inputs/realdata/ (see figures/README.md):
  - S5_selection_flow.csv              participants retained at each selection step
                                       (written by realdata/cohort/selection_flow.py)
  - reference_counts_weights_ESS.csv   CN/MCI/AD counts, weights and effective sample sizes
  - APOE_pair_tests.csv                pairwise within-diagnosis permutation tests

Writes into figures/output/tables/ (CONCORD_OUTPUT_DIR/tables):
  SuppTable3_selection.{tex,csv}
  SuppTable4_weights.{tex,csv}
  SuppTable5_without_ptau.{tex,csv}

Captions and table notes are edited in the manuscript source (tables/*.tex there). The caption
and note text below is a copy of the submitted version, so that the output can be compared with
the manuscript file by file.

Run:  python figures/tables/make_si_cohort_tables.py
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


REAL_INPUTS = _env_dir("CONCORD_FIG_INPUTS", FIGDIR / "inputs") / "realdata"
OUT = _env_dir("CONCORD_OUTPUT_DIR", FIGDIR / "output") / "tables"
SELECTION_CSV = REAL_INPUTS / "S5_selection_flow.csv"
WEIGHTS_CSV = REAL_INPUTS / "reference_counts_weights_ESS.csv"
PAIRS_CSV = REAL_INPUTS / "APOE_pair_tests.csv"

COHORTS = [("adni", "ADNI"), ("nacc", "NACC")]
GENOTYPES = ["e2", "e33", "e4"]
DIAGNOSES = ["CN", "MCI", "AD"]
GENO_TEX = {"e2": r"$\varepsilon$2", "e33": r"$\varepsilon$3/$\varepsilon$3", "e4": r"$\varepsilon$4"}
GENO_TXT = {"e2": "ε2", "e33": "ε3/ε3", "e4": "ε4"}
# Row order follows main Table 3.
PAIR_ORDER = [("e33", "e4"), ("e2", "e33"), ("e2", "e4")]

APOE_TEX = r"\textit{APOE}"


# ----------------------------------------------------------------- helpers
def read_csv_path(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def note(text: str) -> list[str]:
    """Table note below the tabular, as in the manuscript."""
    return [r"\par\smallskip", r"{\footnotesize\raggedright " + text + r"\par}"]


def close(a: float, b: float, tol: float) -> bool:
    return abs(float(a) - float(b)) <= tol


def fmt_int(n: int) -> str:
    return f"{int(n):,}"


def hdr(*lines: str, align: str = "r") -> str:
    """Multi-line column header without extra packages."""
    if len(lines) == 1:
        return lines[0]
    body = r"\\".join(lines)
    return rf"\begin{{tabular}}[b]{{@{{}}{align}@{{}}}}{body}\end{{tabular}}"


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def write_tex(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def caption(title: str) -> str:
    return rf"\caption{{{title}}}"


# ------------------------------------------------ Supplementary Table 3
def table_selection() -> dict:
    rows = read_csv_path(SELECTION_CSV)
    stages = [
        ("Supplied PHC CSF participants",
         r"ADSP-PHC-harmonized CSF data supplied",
         "ADSP-PHC-harmonized CSF data supplied"),
        ("Eligible CSF anchor diagnosis APOE",
         rf"Eligible CSF anchor visit with diagnosis and {APOE_TEX} genotype",
         "Eligible CSF anchor visit with diagnosis and APOE genotype"),
        ("Valid age sex education",
         r"Valid age, recorded sex and years of education",
         "Valid age, recorded sex and years of education"),
        ("Complete six events and covariates",
         r"Complete six-event CSF--cognition profile",
         "Complete six-event CSF–cognition profile"),
    ]
    counts = {}
    for r in rows:
        counts[(r["cohort"], r["stage"])] = int(r["participants"])
    expected = {
        "ADNI": [1249, 1231, 1227, 1203],
        "NACC": [2094, 1833, 1828, 1584],
    }
    table = {}
    for coh, exp in expected.items():
        got = [counts[(coh, s[0])] for s in stages]
        assert got == exp, (coh, got, exp)
        assert all(a >= b for a, b in zip(got, got[1:])), (coh, got)
        table[coh] = got
    assert len(rows) == 8, len(rows)

    def excl(coh, i):
        return "" if i == 0 else fmt_int(table[coh][i - 1] - table[coh][i])

    csv_header = ["Selection step",
                  "ADNI retained (n)", "ADNI excluded (n)",
                  "NACC retained (n)", "NACC excluded (n)"]
    csv_rows, tex_rows = [], []
    for i, (_, tex_lab, txt_lab) in enumerate(stages):
        cells = [fmt_int(table["ADNI"][i]), excl("ADNI", i),
                 fmt_int(table["NACC"][i]), excl("NACC", i)]
        csv_rows.append([txt_lab] + cells)
        tex_rows.append(" & ".join([tex_lab] + cells) + r" \\")

    # caption title and note: copies of the manuscript text (edited there, not here)
    title = "Selection of the ADNI and NACC analysis populations"
    legend = (
        r"Participants retained after each step and excluded at that step (Supplementary "
        r"Methods 2). ADSP-PHC, Alzheimer's Disease Sequencing Project Phenotype Harmonization "
        r"Consortium; CSF, cerebrospinal fluid."
    )
    tex = [
        r"\begin{table}[htbp]",
        r"\centering",
        caption(title),
        r"\label{supptab:selection}",
        r"\small",
        r"\begin{tabularx}{\linewidth}{@{}>{\raggedright\arraybackslash}Xrrrr@{}}",
        r"\toprule",
        r" & \multicolumn{2}{c}{ADNI} & \multicolumn{2}{c}{NACC} \\",
        r"\cmidrule(lr){2-3}\cmidrule(l){4-5}",
        " & ".join(["Selection step"] + [hdr("Retained", "($n$)"), hdr("Excluded", "($n$)")] * 2)
        + r" \\",
        r"\midrule",
        *tex_rows,
        r"\bottomrule",
        r"\end{tabularx}",
        *note(legend),
        r"\end{table}",
    ]
    write_tex(OUT / "SuppTable3_selection.tex", tex)
    write_csv(OUT / "SuppTable3_selection.csv", csv_header, csv_rows)
    return {
        "ADNI": table["ADNI"],
        "NACC": table["NACC"],
    }


# ------------------------------------------------ Supplementary Table 4
def table_weights() -> dict:
    rows = [r for r in read_csv_path(WEIGHTS_CSV) if r["application"] == "actual_APOE"]
    assert len(rows) == 18, len(rows)
    cell = {(r["cohort"], r["group"], r["diagnosis"]): r for r in rows}

    # Independent recomputation of the common proportions, weights and ESS.
    props = np.array([[[int(cell[(c, g, d)]["cell_n"]) for d in DIAGNOSES]
                       for g in GENOTYPES] for c, _ in COHORTS], dtype=float)
    ns = props.copy()
    props = props / props.sum(axis=2, keepdims=True)
    common = props.reshape(-1, 3).min(axis=0)
    common = common / common.sum()
    for (c, _), ci in zip(COHORTS, range(2)):
        for gi, g in enumerate(GENOTYPES):
            w = common / props[ci, gi]
            ess = (ns[ci, gi] * w).sum() ** 2 / (ns[ci, gi] * w ** 2).sum()
            for di, d in enumerate(DIAGNOSES):
                r = cell[(c, g, d)]
                assert int(r["group_n"]) == int(ns[ci, gi].sum())
                assert close(r["original_proportion"], props[ci, gi, di], 1e-12)
                assert close(r["reference_proportion"], common[di], 1e-12)
                assert close(r["standardization_weight"], w[di], 1e-10)
                assert close(r["group_ESS"], ess, 1e-8)

    # Expected numbers printed in the manuscript.
    assert [round(100 * x, 1) for x in common] == [48.9, 28.6, 22.5], common
    a = lambda c, g, d, k: float(cell[(c, g, d)][k])
    assert int(a("adni", "e2", "CN", "cell_n")) == 50
    assert round(100 * a("adni", "e2", "CN", "original_proportion"), 1) == 51.5
    assert round(a("adni", "e2", "CN", "standardization_weight"), 2) == 0.95
    assert int(a("adni", "e2", "MCI", "cell_n")) == 39
    assert round(a("adni", "e2", "MCI", "standardization_weight"), 2) == 0.71
    assert int(a("adni", "e2", "AD", "cell_n")) == 8
    assert round(a("adni", "e2", "AD", "standardization_weight"), 2) == 2.73
    wmax = max(float(r["standardization_weight"]) for r in rows)
    assert round(wmax, 2) == 2.73, wmax
    ess_exp = {"adni": [75.7, 464.0, 327.3], "nacc": [99.9, 635.5, 530.2]}
    for c, _ in COHORTS:
        got = [round(a(c, g, "CN", "group_ESS"), 1) for g in GENOTYPES]
        assert got == ess_exp[c], (c, got)
    assert sum(int(cell[("adni", g, "CN")]["group_n"]) for g in GENOTYPES) == 1203
    assert sum(int(cell[("nacc", g, "CN")]["group_n"]) for g in GENOTYPES) == 1584

    csv_header = ["Cohort", "Genotype", "Diagnosis", "n", "Proportion in group (%)",
                  "Common proportion (%)", "Weight", "Group n",
                  "Group effective sample size (n)"]
    csv_rows, tex_rows = [], []
    for ci, (c, clab) in enumerate(COHORTS):
        if ci > 0:
            tex_rows.append(r"\addlinespace")
        tex_rows.append(rf"\multicolumn{{8}}{{@{{}}l}}{{{clab}}} \\")
        for gi, g in enumerate(GENOTYPES):
            for di, d in enumerate(DIAGNOSES):
                r = cell[(c, g, d)]
                vals = [
                    fmt_int(r["cell_n"]),
                    f"{100 * float(r['original_proportion']):.1f}",
                    f"{100 * float(r['reference_proportion']):.1f}",
                    f"{float(r['standardization_weight']):.2f}",
                ]
                gn = fmt_int(r["group_n"])
                ess = f"{float(r['group_ESS']):.1f}"
                csv_rows.append([clab, GENO_TXT[g], d] + vals + [gn, ess])
                first = di == 0
                tex_rows.append(" & ".join(
                    [r"\quad " + GENO_TEX[g] if first else "", d] + vals
                    + ([gn, ess] if first else ["", ""])) + r" \\")

    # caption title and note: copies of the manuscript text (edited there, not here)
    title = rf"Weights and effective sample sizes in the {APOE_TEX} analyses"
    legend = (
        r"Weight = common proportion $\div$ proportion in the group. After weighting, every group "
        r"has 48.9\% CN, 28.6\% MCI and 22.5\% AD. Effective sample size = "
        r"$(\sum_j w_j)^2/\sum_j w_j^2$ over the participants $j$ of the group. CN, cognitively "
        r"normal; MCI, mild cognitive impairment; AD, Alzheimer's disease dementia."
    )
    tex = [
        r"\begin{table}[htbp]",
        r"\centering",
        caption(title),
        r"\label{supptab:weights}",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{@{}llrrrrrr@{}}",
        r"\toprule",
        " & ".join([
            hdr("Genotype", align="l"), hdr("Diagnosis", align="l"), hdr("$n$"),
            hdr("Proportion", r"in group (\%)"), hdr("Common", r"proportion (\%)"),
            hdr("Weight"), hdr("Group", "$n$"), hdr("Effective", "sample size"),
        ]) + r" \\",
        r"\midrule",
        *tex_rows,
        r"\bottomrule",
        r"\end{tabular}",
        *note(legend),
        r"\end{table}",
    ]
    write_tex(OUT / "SuppTable4_weights.tex", tex)
    write_csv(OUT / "SuppTable4_weights.csv", csv_header, csv_rows)
    group_n = {(c, g): int(cell[(c, g, "CN")]["group_n"]) for c, _ in COHORTS for g in GENOTYPES}
    return {
        "common_pct": [round(100 * float(x), 1) for x in common],
        "max_weight": round(wmax, 2),
        "ESS": {c: [round(a(c, g, "CN", "group_ESS"), 1) for g in GENOTYPES] for c, _ in COHORTS},
        "group_n": group_n,
    }


# ------------------------------------------------ Supplementary Table 5
def table_without_ptau(group_n: dict) -> dict:
    allrows = read_csv_path(PAIRS_CSV)
    panels = sorted({r["panel"] for r in allrows})
    assert panels == ["CSFcog6", "omit_ptau"], panels
    rows = [r for r in allrows if r["panel"] == "omit_ptau"]
    assert len(rows) == 6, len(rows)
    res = {}
    for r in rows:
        assert r["analysis_role"] == "same_person_sensitivity", r["analysis_role"]
        k = int(r["K"])
        assert k == 5
        n_pairs = k * (k - 1) // 2
        planned, done = int(r["permutations_planned"]), int(r["permutations_completed"])
        assert planned == done == 599 and int(r["failed_permutations"]) == 0
        e = int(r["inclusive_exceedances"])
        raw = (1 + e) / (done + 1)
        assert close(r["raw_P"], raw, 1e-12)
        assert int(r["Bonferroni_factor"]) == 3
        adj = min(1.0, 3 * raw)
        assert close(r["adjusted_P"], adj, 1e-9)
        rev = float(r["normalized_kendall"]) * n_pairs
        assert close(rev, round(rev), 1e-9), rev
        res[(r["cohort"], r["group_A"], r["group_B"])] = (int(round(rev)), e, adj, n_pairs)

    expected = {
        ("adni", "e2", "e33"): (1, 381, 1.000), ("adni", "e2", "e4"): (0, 599, 1.000),
        ("adni", "e33", "e4"): (1, 0, 0.005),
        ("nacc", "e2", "e33"): (0, 599, 1.000), ("nacc", "e2", "e4"): (1, 457, 1.000),
        ("nacc", "e33", "e4"): (1, 372, 1.000),
    }
    for key, (rev, e, p) in expected.items():
        got = res[key]
        assert got[0] == rev and got[1] == e and round(got[2], 3) == p, (key, got)

    csv_header = ["Comparison"]
    for _, clab in COHORTS:
        csv_header += [f"{clab} event pairs reversed (of 10)",
                       f"{clab} relabelings at least as extreme (of 599)",
                       f"{clab} adjusted P"]
    csv_rows, tex_rows = [], []
    for ga, gb in PAIR_ORDER:
        cells = []
        for c, _ in COHORTS:
            rev, e, adj, _ = res[(c, ga, gb)]
            cells += [str(rev), fmt_int(e), f"{adj:.3f}"]
        csv_rows.append([f"{GENO_TXT[ga]} vs {GENO_TXT[gb]}"] + cells)
        tex_rows.append(" & ".join([f"{GENO_TEX[ga]} vs {GENO_TEX[gb]}"] + cells) + r" \\")

    # group sizes (Table 2): the sensitivity analysis uses the same participants
    assert [group_n[("adni", g)] for g in GENOTYPES] == [97, 560, 546]
    assert [group_n[("nacc", g)] for g in GENOTYPES] == [143, 759, 682]
    # caption title and note: copies of the manuscript text (edited there, not here)
    title = rf"{APOE_TEX} comparisons without p-tau"
    legend = (
        r"Prespecified sensitivity analysis with five events (10 event pairs) on the participants of "
        r"Table 2 of the main text; otherwise as Table 3 of the main text."
    )
    colh = [hdr("Pairs reversed", "(of 10)"),
            hdr(r"Relabelings $\geq$", "observed (of 599)"),
            hdr("Adjusted", "$P$")]
    tex = [
        r"\begin{table}[htbp]",
        r"\centering",
        caption(title),
        r"\label{supptab:without-ptau}",
        r"\small",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}lrrrrrr@{}}",
        r"\toprule",
        r" & \multicolumn{3}{c}{ADNI} & \multicolumn{3}{c}{NACC} \\",
        r"\cmidrule(lr){2-4}\cmidrule(l){5-7}",
        " & ".join(["Comparison"] + colh + colh) + r" \\",
        r"\midrule",
        *tex_rows,
        r"\bottomrule",
        r"\end{tabular*}",
        *note(legend),
        r"\end{table}",
    ]
    write_tex(OUT / "SuppTable5_without_ptau.tex", tex)
    write_csv(OUT / "SuppTable5_without_ptau.csv", csv_header, csv_rows)
    return {f"{c}:{ga}-{gb}": res[(c, ga, gb)][:3] for c, _ in COHORTS for ga, gb in PAIR_ORDER}


def main() -> None:
    missing = [p for p in (SELECTION_CSV, WEIGHTS_CSV, PAIRS_CSV) if not p.is_file()]
    if missing:
        sys.exit("Supplementary Tables 3-5 need aggregate ADNI/NACC inputs that are not in this "
                 "repository.\n" + "".join(f"  - {p}\n" for p in missing)
                 + f"Place them in {REAL_INPUTS} (or set CONCORD_FIG_INPUTS); see figures/README.md.")
    OUT.mkdir(parents=True, exist_ok=True)
    s3 = table_selection()
    s4 = table_weights()
    s5 = table_without_ptau(s4["group_n"])
    print("[SuppTable3_selection]", s3)
    print("[SuppTable4_weights] common % =", s4["common_pct"], "max weight =", s4["max_weight"],
          "ESS =", s4["ESS"])
    print("[SuppTable5_without_ptau]", s5)
    print("[OK] all assertions passed; 3 tables written to", OUT)


if __name__ == "__main__":
    main()
