"""Main Tables 1-3 of the CONCORD manuscript.

Table 1  Event-based models compared in this study (qualitative; no input file).
Table 2  Participants in ADNI and NACC.
Table 3  Comparisons of APOE groups in ADNI and NACC.

Tables 2 and 3 read aggregate ADNI/NACC results (no participant-level data) that are not
distributed with this repository; place them in figures/inputs/realdata/ (see figures/README.md):
  Table2_demographics_aggregate.csv   group-level age, sex and education
  reference_counts_weights_ESS.csv    CN/MCI/AD counts, weights and effective sample sizes
  APOE_pair_tests.csv                 pairwise within-diagnosis permutation tests
Without them only Table 1 is written.

Writes <Stem>.tex (booktabs table environment) and <Stem>.csv (same cells, UTF-8) into
figures/output/tables/ (CONCORD_OUTPUT_DIR/tables). Key values are asserted against the numbers
printed in the manuscript.

Captions and table notes are edited in the manuscript source (tables/*.tex there). The caption
and note text below is a copy of the submitted version, so that the output can be compared with
the manuscript file by file.

Run:  python figures/tables/make_main_tables.py
"""
from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

FIGDIR = Path(__file__).resolve().parents[1]  # figures/


def _env_dir(name: str, default: Path) -> Path:
    v = os.environ.get(name)
    return Path(v).expanduser() if v else default


INPUTS = _env_dir("CONCORD_FIG_INPUTS", FIGDIR / "inputs")
REAL_INPUTS = INPUTS / "realdata"
OUT = _env_dir("CONCORD_OUTPUT_DIR", FIGDIR / "output") / "tables"
DEMOG = REAL_INPUTS / "Table2_demographics_aggregate.csv"
ESS = REAL_INPUTS / "reference_counts_weights_ESS.csv"
PAIRS = REAL_INPUTS / "APOE_pair_tests.csv"

EM = "---"  # empty cell (em dash) in LaTeX
EM_U = "—"  # empty cell in CSV

# Genotype labels: LaTeX and plain Unicode.
GENO_TEX = {
    "e2": r"$\varepsilon$2",
    "e33": r"$\varepsilon$3/$\varepsilon$3",
    "e4": r"$\varepsilon$4",
    "all": "All",
}
GENO_TXT = {"e2": "ε2", "e33": "ε3/ε3", "e4": "ε4", "all": "All"}
COHORTS = ("ADNI", "NACC")
GENOS = ("e2", "e33", "e4")

REF_PROP = {"CN": 0.489, "MCI": 0.286, "AD": 0.225}  # reported reference, one decimal in %


# --------------------------------------------------------------------------- helpers
def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def thousands(n: int) -> str:
    return f"{n:,}"


def n_pct(k: int, n: int) -> str:
    return f"{thousands(k)} ({100.0 * k / n:.1f})"


def mean_sd(m: float, s: float) -> str:
    return f"{m:.1f} ({s:.1f})"


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def hdr(lines: list[str], align: str = "r") -> str:
    """Stacked column header (nested tabular, bottom-aligned); plain LaTeX, no extra packages."""
    if len(lines) == 1:
        return lines[0]
    return rf"\begin{{tabular}}[b]{{@{{}}{align}@{{}}}}" + r" \\ ".join(lines) + r"\end{tabular}"


def write_csv(path: Path, header_rows: list[list[str]], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        for h in header_rows:
            w.writerow(h)
        for r in rows:
            w.writerow(r)


def write_tex(path: Path, body: str) -> None:
    path.write_text(body.strip() + "\n", encoding="utf-8")


def note(text: str) -> str:
    """Table note below the tabular, as in the manuscript."""
    return "\\par\\smallskip\n{\\footnotesize\\raggedright " + text + "\\par}"


# Caption titles and notes: copies of the manuscript text (edited there, not here).
NOTE_T1 = (r"Abnormality model, fitted in each group or once to the pooled sample of the compared groups. "
           r"CN/MCI/AD proportions, the proportions at which each group's ordering is estimated. "
           r"Permutation test, the test in current use (separately fitted DEBM) or used here (pooled-score "
           r"DEBM and CONCORD); Fig.~\ref{fig:testing} applies within-diagnosis permutation to all three "
           r"DEBM variants. EBM, event-based model; DEBM, discriminative EBM; CN, cognitively normal; MCI, "
           r"mild cognitive impairment; AD, Alzheimer's disease dementia.")
NOTE_T2 = (r"$\varepsilon$2, $\varepsilon$2 carriers without $\varepsilon$4; $\varepsilon$4, "
           r"$\varepsilon$4 carriers without $\varepsilon$2. Percentages are of the row total. Diagnostic "
           r"composition differed between genotypes in both cohorts ($\chi^2$ test, 4 d.f., $P<0.001$). "
           r"ADNI, Alzheimer's Disease Neuroimaging Initiative; NACC, National Alzheimer's Coordinating "
           r"Center; CN, cognitively normal; MCI, mild cognitive impairment; AD, Alzheimer's disease "
           r"dementia; SD, standard deviation; d.f., degrees of freedom.")
NOTE_T3 = (r"All groups were weighted to 48.9\% CN, 28.6\% MCI and 22.5\% AD. Pairs reversed, event pairs "
           r"ordered differently by the two groups. Relabelings $\geq$ observed ($E$), pairwise "
           r"within-diagnosis permutations that reversed at least as many pairs as observed. Adjusted $P$, "
           r"$(1+E)/600$ with Bonferroni correction for three comparisons. CN, cognitively normal; MCI, "
           r"mild cognitive impairment; AD, Alzheimer's disease dementia.")


# --------------------------------------------------------------------------- Table 1
def table1() -> dict:
    stem = "Table1_models"
    # (model, ordering, abnormality model, CN/MCI/AD proportions, permutation test); short single-line cells
    rows_tex = [
        ("Likelihood EBM", "Sequence likelihood", "Per group", "Observed", EM),
        ("Stage-aware EBM", "Bayesian, stage-aware", "Per group", "Observed", EM),
        ("Separately fitted DEBM", "Ranking consensus", "Per group", "Observed", "Unrestricted"),
        ("Pooled-score DEBM", "Ranking consensus", "Pooled sample", "Observed", "Within diagnosis"),
        ("CONCORD", "Ranking consensus", "Pooled sample", "Same for all groups", "Within diagnosis"),
    ]
    header = ["Model", "Ordering", "Abnormality model", "CN/MCI/AD proportions", "Permutation test"]

    def plain(s: str) -> str:
        return s.replace(EM, EM_U)

    rows_txt = [[plain(c) for c in r] for r in rows_tex]
    assert len(rows_tex) == 5
    assert [r[0] for r in rows_txt] == ["Likelihood EBM", "Stage-aware EBM", "Separately fitted DEBM",
                                        "Pooled-score DEBM", "CONCORD"]

    body_lines = [" & ".join(r) + r" \\" for r in rows_tex]
    tex = rf"""
\begin{{table}}[htbp]
\centering
\caption{{Event-based models compared in this study}}
\label{{tab:models}}
\footnotesize
\setlength{{\tabcolsep}}{{3pt}}
\renewcommand{{\arraystretch}}{{1.2}}
\begin{{tabular*}}{{\linewidth}}{{@{{\extracolsep{{\fill}}}}lllll@{{}}}}
\toprule
{" & ".join(header)} \\
\midrule
{chr(10).join(body_lines)}
\bottomrule
\end{{tabular*}}
{note(NOTE_T1)}
\end{{table}}
"""
    write_tex(OUT / f"{stem}.tex", tex)
    write_csv(OUT / f"{stem}.csv", [header], rows_txt)
    return {"stem": stem, "rows": rows_txt}


# --------------------------------------------------------------------------- Table 2
def table2() -> dict:
    stem = "Table2_participants"
    demog = read_csv(DEMOG)
    ess_rows = [r for r in read_csv(ESS) if r["application"] == "actual_APOE"]

    d = {(r["cohort"], r["genotype"]): r for r in demog}
    assert set(d) == {(c, g) for c in COHORTS for g in (*GENOS, "all")}, sorted(d)

    # --- expected counts (as printed in Table 2) ---
    expected = {
        ("ADNI", "e2"): (97, 50, 39, 8), ("ADNI", "e33"): (560, 222, 273, 65),
        ("ADNI", "e4"): (546, 98, 299, 149),
        ("NACC", "e2"): (143, 106, 15, 22), ("NACC", "e33"): (759, 480, 104, 175),
        ("NACC", "e4"): (682, 281, 94, 307),
    }
    for key, (n, cn, mci, ad) in expected.items():
        r = d[key]
        got = tuple(int(r[k]) for k in ("n", "CN", "MCI", "AD"))
        assert got == (n, cn, mci, ad), (key, got)
        assert cn + mci + ad == n, key
    for c in COHORTS:
        allr = d[(c, "all")]
        for k in ("n", "CN", "MCI", "AD", "female_n"):
            assert int(allr[k]) == sum(int(d[(c, g)][k]) for g in GENOS), (c, k)
    assert int(d[("ADNI", "all")]["n"]) == 1203 and int(d[("NACC", "all")]["n"]) == 1584
    # female percentage consistent with counts
    for key, r in d.items():
        assert close(100.0 * int(r["female_n"]) / int(r["n"]), float(r["female_pct"]), 0.051), key

    # --- effective sample size: recompute from cell counts and weights, compare with stored ---
    ess: dict[tuple[str, str], float] = {}
    for c in COHORTS:
        for g in GENOS:
            cells = [r for r in ess_rows if r["cohort"] == c.lower() and r["group"] == g]
            assert len(cells) == 3, (c, g)
            n_c = np.array([float(r["cell_n"]) for r in cells])
            w_c = np.array([float(r["standardization_weight"]) for r in cells])
            ess_calc = (n_c * w_c).sum() ** 2 / (n_c * w_c ** 2).sum()
            stored = {float(r["group_ESS"]) for r in cells}
            assert len(stored) == 1
            stored_v = stored.pop()
            assert close(ess_calc, stored_v, 1e-6), (c, g, ess_calc, stored_v)
            # cell counts agree with the demographics aggregate
            by_dx = {r["diagnosis"]: int(r["cell_n"]) for r in cells}
            assert (by_dx["CN"], by_dx["MCI"], by_dx["AD"]) == expected[(c, g)][1:], (c, g)
            # reference proportions are the reported 48.9 / 28.6 / 22.5 %
            for r in cells:
                assert close(float(r["reference_proportion"]), REF_PROP[r["diagnosis"]], 5e-4), r
                w_expect = float(r["reference_proportion"]) / float(r["original_proportion"])
                assert close(float(r["standardization_weight"]), w_expect, 1e-9), r
            ess[(c, g)] = stored_v
    ess_expected = {("ADNI", "e2"): 75.7, ("ADNI", "e33"): 464.0, ("ADNI", "e4"): 327.3,
                    ("NACC", "e2"): 99.9, ("NACC", "e33"): 635.5, ("NACC", "e4"): 530.2}
    for key, v in ess_expected.items():
        assert round(ess[key], 1) == v, (key, ess[key])

    # --- genotype x diagnosis chi-square (df 4) ---
    chi2: dict[str, float] = {}
    pval: dict[str, float] = {}
    for c in COHORTS:
        obs = np.array([[int(d[(c, g)][k]) for k in ("CN", "MCI", "AD")] for g in GENOS], float)
        exp = obs.sum(1, keepdims=True) * obs.sum(0, keepdims=True) / obs.sum()
        x = float(((obs - exp) ** 2 / exp).sum())
        df = (obs.shape[0] - 1) * (obs.shape[1] - 1)
        assert df == 4
        chi2[c] = x
        pval[c] = math.exp(-x / 2.0) * (1.0 + x / 2.0)  # chi-square survival function, 4 d.f.
    assert close(chi2["ADNI"], 104.2, 0.05), chi2
    assert close(chi2["NACC"], 111.5, 0.05), chi2
    assert pval["ADNI"] < 1e-3 and pval["NACC"] < 1e-3

    header = ["Cohort", "APOE group", "n", "CN, n (%)", "MCI, n (%)", "AD, n (%)",
              "Age, years, mean (SD)", "Female, n (%)", "Education, years, mean (SD)"]
    header_tex = [
        "Cohort", hdr([r"\textit{APOE}", "group"], "l"), r"$n$",
        r"CN, $n$ (\%)", r"MCI, $n$ (\%)", r"AD, $n$ (\%)",
        hdr(["Age, years,", "mean (SD)"]), hdr(["Female,", r"$n$ (\%)"]),
        hdr(["Education, years,", "mean (SD)"]),
    ]

    rows_txt: list[list[str]] = []
    tex_blocks: list[str] = []
    for ci, c in enumerate(COHORTS):
        block = []
        for gi, g in enumerate((*GENOS, "all")):
            r = d[(c, g)]
            n = int(r["n"])
            cells = [
                thousands(n),
                n_pct(int(r["CN"]), n), n_pct(int(r["MCI"]), n), n_pct(int(r["AD"]), n),
                mean_sd(float(r["age_mean"]), float(r["age_sd"])),
                f"{thousands(int(r['female_n']))} ({float(r['female_pct']):.1f})",
                mean_sd(float(r["education_mean"]), float(r["education_sd"])),
            ]
            rows_txt.append([c if gi == 0 else "", GENO_TXT[g]] + [x if x is not None else EM_U for x in cells])
            tex_cells = [c if gi == 0 else "", GENO_TEX[g]] + [x if x is not None else EM for x in cells]
            block.append(" & ".join(tex_cells) + r" \\")
        tex_blocks.append("\n".join(block))

    # the chi-square values (104.2 and 111.5) are asserted above and printed by main()
    tex = rf"""
\begin{{table}}[htbp]
\centering
\caption{{Participants in ADNI and NACC}}
\label{{tab:participants}}
\small
\setlength{{\tabcolsep}}{{3pt}}
\renewcommand{{\arraystretch}}{{1.15}}
\begin{{tabular*}}{{\linewidth}}{{@{{\extracolsep{{\fill}}}}l l *{{7}}{{r}}@{{}}}}
\toprule
{" & ".join(header_tex)} \\
\midrule
{tex_blocks[0]}
\midrule
{tex_blocks[1]}
\bottomrule
\end{{tabular*}}
{note(NOTE_T2)}
\end{{table}}
"""
    write_tex(OUT / f"{stem}.tex", tex)
    write_csv(OUT / f"{stem}.csv", [header], rows_txt)
    return {"stem": stem, "rows": rows_txt, "chi2": chi2, "p": pval, "ess": ess}


# --------------------------------------------------------------------------- Table 3
def table3() -> dict:
    stem = "Table3_apoe_tests"
    tests = [r for r in read_csv(PAIRS) if r["panel"] == "CSFcog6" and r["analysis_role"] == "primary"]
    assert len(tests) == 6, len(tests)
    comps = [("e33", "e4"), ("e2", "e33"), ("e2", "e4")]
    comp_tex = {("e33", "e4"): rf"{GENO_TEX['e33']} vs {GENO_TEX['e4']}",
                ("e2", "e33"): rf"{GENO_TEX['e2']} vs {GENO_TEX['e33']}",
                ("e2", "e4"): rf"{GENO_TEX['e2']} vs {GENO_TEX['e4']}"}
    comp_txt = {k: f"{GENO_TXT[k[0]]} vs {GENO_TXT[k[1]]}" for k in comps}

    expected = {  # reversed pairs, exceedances, adjusted P
        ("ADNI", ("e33", "e4")): (2, 0, 0.005), ("ADNI", ("e2", "e33")): (1, 506, 1.000),
        ("ADNI", ("e2", "e4")): (1, 184, 0.925),
        ("NACC", ("e33", "e4")): (3, 0, 0.005), ("NACC", ("e2", "e33")): (1, 525, 1.000),
        ("NACC", ("e2", "e4")): (4, 18, 0.095),
    }
    res: dict[tuple[str, tuple[str, str]], tuple[int, int, float]] = {}
    for r in tests:
        c = r["cohort"].upper()
        key = (c, (r["group_A"], r["group_B"]))
        assert key[1] in comps, key
        K = int(r["K"])
        assert K == 6
        n_pairs = K * (K - 1) // 2
        assert n_pairs == 15
        rev_f = float(r["normalized_kendall"]) * n_pairs
        rev = int(round(rev_f))
        assert close(rev_f, rev, 1e-9), (key, rev_f)
        B = int(r["permutations_planned"])
        assert B == 599 and int(r["permutations_completed"]) == 599 and int(r["failed_permutations"]) == 0
        E = int(r["inclusive_exceedances"])
        raw = (1 + E) / (B + 1)
        assert close(raw, float(r["raw_P"]), 1e-12), key
        assert int(r["Bonferroni_factor"]) == 3
        adj = min(1.0, 3 * raw)
        assert close(adj, float(r["adjusted_P"]), 1e-12), key
        res[key] = (rev, E, adj)
    for key, (rev, E, adj) in expected.items():
        assert res[key][0] == rev and res[key][1] == E and close(res[key][2], adj, 5e-4), (key, res[key])
    min_adj = min(1.0, 3 * 1 / 600)
    assert close(min_adj, 0.005, 1e-12)

    sub = ["Pairs reversed (of 15)", "Relabelings >= observed (of 599)", "Adjusted P"]
    sub_tex = [hdr(["Pairs reversed", "(of 15)"]),
               hdr([r"Relabelings $\geq$", "observed (of 599)"]),
               hdr(["Adjusted", "$P$"])]
    header_txt = ["Comparison"] + [f"{c}: {s}" for c in COHORTS for s in sub]

    rows_txt, rows_tex = [], []
    for k in comps:
        cells = []
        for c in COHORTS:
            rev, E, adj = res[(c, k)]
            cells += [str(rev), str(E), f"{adj:.3f}"]
        rows_txt.append([comp_txt[k]] + cells)
        rows_tex.append(" & ".join([comp_tex[k]] + cells) + r" \\")

    tex = rf"""
\begin{{table}}[htbp]
\centering
\caption{{Comparisons of \textit{{APOE}} groups in ADNI and NACC}}
\label{{tab:apoe-tests}}
\small
\setlength{{\tabcolsep}}{{4pt}}
\renewcommand{{\arraystretch}}{{1.15}}
\begin{{tabular*}}{{\linewidth}}{{@{{\extracolsep{{\fill}}}}l *{{6}}{{r}}@{{}}}}
\toprule
 & \multicolumn{{3}}{{c}}{{ADNI}} & \multicolumn{{3}}{{c}}{{NACC}} \\
\cmidrule(lr){{2-4}}\cmidrule(l){{5-7}}
Comparison & {" & ".join(sub_tex + sub_tex)} \\
\midrule
{chr(10).join(rows_tex)}
\bottomrule
\end{{tabular*}}
{note(NOTE_T3)}
\end{{table}}
"""
    write_tex(OUT / f"{stem}.tex", tex)
    write_csv(OUT / f"{stem}.csv", [header_txt], rows_txt)
    return {"stem": stem, "rows": rows_txt}


# --------------------------------------------------------------------------- main
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t1 = table1()
    print(f"[ok] {t1['stem']}: {len(t1['rows'])} models")
    missing = [p for p in (DEMOG, ESS, PAIRS) if not p.is_file()]
    if missing:
        sys.exit("Tables 2 and 3 need aggregate ADNI/NACC inputs that are not in this repository.\n"
                 + "".join(f"  - {p}\n" for p in missing)
                 + f"Place them in {REAL_INPUTS} (or set CONCORD_FIG_INPUTS); see figures/README.md.")
    t2 = table2()
    t3 = table3()
    print(f"[ok] {t2['stem']}: chi2 ADNI {t2['chi2']['ADNI']:.2f} (P={t2['p']['ADNI']:.1e}), "
          f"NACC {t2['chi2']['NACC']:.2f} (P={t2['p']['NACC']:.1e})")
    for r in t2["rows"]:
        print("    ", " | ".join(r))
    print(f"[ok] {t3['stem']}")
    for r in t3["rows"]:
        print("    ", " | ".join(r))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
