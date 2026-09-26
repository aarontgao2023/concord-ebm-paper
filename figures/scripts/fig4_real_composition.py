"""Figure 4 | Changing only the diagnostic composition of groups in ADNI and NACC.

180 x 150 mm, Nature Communications print size. Every plotted number is read from stored
aggregate files (no participant-level data, no model fits, no permutations). Panel a is a
schematic of the fixed design constants.

Inputs: aggregate ADNI/NACC results in figures/inputs/realdata/ (not distributed with this
repository; see figures/README.md):
  real_comp_gap.csv                     b  mean-position gap, matched and different composition
  composition_distance_summary.csv      c  discordant event pairs, matched composition
  real_comp_event_shifts.csv            d  per-event position shifts
  fig3_constructed_null_summary.csv     e  ADNI CSF-MRI-cognition (12 events), NACC amyloid PET-cognition
  precision_adni_rates.csv              e  ADNI CSF-MRI-cognition (8 and 5 events)
  n2_rejection_rates.csv                e  ADNI and NACC MRI-cognition (8 and 4 events)
Model labels in the input files: likelihood_ebm = likelihood EBM, separate = separately fitted
DEBM, shared = pooled-score DEBM, concord = CONCORD. In e, the keys standard_U / Standard_U are
separately fitted DEBM with unrestricted permutation and concord_D / CONCORD_D are CONCORD with
within-diagnosis permutation.
"""
from __future__ import annotations

import csv
import sys
from decimal import ROUND_HALF_UP, Decimal
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

S.apply_style()

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

STEM = "Fig4_real_composition"

F_GAP = S.REAL_INPUTS / "real_comp_gap.csv"
F_DIST = S.REAL_INPUTS / "composition_distance_summary.csv"
F_SHIFT = S.REAL_INPUTS / "real_comp_event_shifts.csv"
F_SD21 = S.REAL_INPUTS / "fig3_constructed_null_summary.csv"
F_PREC = S.REAL_INPUTS / "precision_adni_rates.csv"
F_SHARED = S.REAL_INPUTS / "n2_rejection_rates.csv"
S.require_inputs([F_GAP, F_DIST, F_SHIFT, F_SD21, F_PREC, F_SHARED], "Fig. 4")

ENGINE = {"likelihood_ebm": "likelihood", "separate": "separate", "shared": "pooled", "concord": "concord"}
MODELS = ["likelihood", "separate", "pooled", "concord"]  # CONCORD last (drawn on top)
COHORTS = [("adni", "ADNI"), ("nacc", "NACC")]
M = S.METHOD
INK = S.PAL["ink"]
GREY = S.PAL["grey_mid"]
ROWBAND = "#F3F3F3"
SEP = "#D9D9D9"

TOL_RATE = 0.002
TOL_GAP = 0.01


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def check(name, got, want, tol):
    assert abs(got - want) <= tol, f"{name}: got {got:.4f}, expected {want:.4f} (tol {tol})"


def write_csv(name, header, rows):
    S.SRC.mkdir(parents=True, exist_ok=True)
    with open(S.SRC / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def minus(s):
    return s.replace("-", "−")


def fmt_inc(x):
    """Signed increment (two decimals; three below 0.1), rounded half-up from the stored
    4-decimal value (0.565 -> +0.57, 0.845 -> +0.85), typographic minus."""
    q = Decimal("0.001") if abs(x) < 0.1 else Decimal("0.01")
    return minus(f"{Decimal(str(x)).quantize(q, rounding=ROUND_HALF_UP):+f}")


# ----------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------
# b: mean-position gap, matched vs different composition, and the paired increment
gap = {}
for r in read_csv(F_GAP):
    key = (r["cohort"], ENGINE[r["engine"]], r["design"])
    gap[key] = dict(
        gap=float(r["gap_mean_abs_mean_shift"]),
        lo=float(r["boot_ci95_low"]),
        hi=float(r["boot_ci95_high"]),
        dmm=float(r["diff_minus_matched_gap"]) if r["diff_minus_matched_gap"] else None,
        dmm_lo=float(r["dmm_ci95_low"]) if r["dmm_ci95_low"] else None,
        dmm_hi=float(r["dmm_ci95_high"]) if r["dmm_ci95_high"] else None,
    )
EXP_GAP = {  # (different, matched)
    ("adni", "likelihood"): (0.618, 0.053), ("adni", "separate"): (0.887, 0.137),
    ("adni", "pooled"): (0.435, 0.022), ("adni", "concord"): (0.025, 0.015),
    ("nacc", "likelihood"): (0.793, 0.125), ("nacc", "separate"): (1.478, 0.068),
    ("nacc", "pooled"): (0.907, 0.062), ("nacc", "concord"): (0.038, 0.057),
}
for (c, m), (d_exp, m_exp) in EXP_GAP.items():
    check(f"gap {c} {m} different", gap[(c, m, "different")]["gap"], d_exp, TOL_GAP)
    check(f"gap {c} {m} matched", gap[(c, m, "matched")]["gap"], m_exp, TOL_GAP)
# paired increment (different - matched), as stored; the printed labels are checked
EXP_INC = {
    ("adni", "likelihood"): "+0.57", ("adni", "separate"): "+0.75", ("adni", "pooled"): "+0.41",
    ("adni", "concord"): "+0.010", ("nacc", "likelihood"): "+0.67", ("nacc", "separate"): "+1.41",
    ("nacc", "pooled"): "+0.85", ("nacc", "concord"): "−0.018",
}
inc = {}
for (c, m), want in EXP_INC.items():
    g0, g1 = gap[(c, m, "matched")], gap[(c, m, "different")]
    check(f"increment {c} {m} = different - matched", g1["dmm"], g1["gap"] - g0["gap"], 0.001)
    inc[(c, m)] = g1["dmm"]
    assert fmt_inc(g1["dmm"]) == want, (c, m, fmt_inc(g1["dmm"]), want)
for m in ("likelihood", "separate", "pooled"):  # comparator increments: interval excludes 0
    assert all(gap[(c, m, "different")]["dmm_lo"] > 0 for c, _ in COHORTS)
assert all(gap[(c, "concord", "different")]["dmm_lo"] < 0 < gap[(c, "concord", "different")]["dmm_hi"]
           for c, _ in COHORTS)

# c: discordant event pairs between two groups of identical composition (matched design)
dist = {}
for r in read_csv(F_DIST):
    if r["kind"] != "absolute_distance" or r["design"] != "matched":
        continue
    assert int(r["complete_replicates"]) == 200
    dist[(r["cohort"], ENGINE[r["engine_left"]])] = dict(
        mean=float(r["mean"]), lo=float(r["mean_ci_low"]), hi=float(r["mean_ci_high"]), n=int(r["n"])
    )
EXP_DIST = {
    ("adni", "likelihood"): 0.346, ("adni", "separate"): 0.326, ("adni", "pooled"): 0.066,
    ("adni", "concord"): 0.068, ("nacc", "likelihood"): 0.512, ("nacc", "separate"): 0.396,
    ("nacc", "pooled"): 0.114, ("nacc", "concord"): 0.115,
}
for k, v in EXP_DIST.items():
    check(f"Kendall {k}", dist[k]["mean"], v, TOL_RATE)

# d: per-event position shift, CN-heavy - AD-heavy = -(B - A); B = AD-heavy, A = CN-heavy
shift = {}
for r in read_csv(F_SHIFT):
    if r["design"] != "different":
        continue
    b_minus_a = float(r["mean_shift_B_minus_A"])
    shift[(r["cohort"], ENGINE[r["engine"]], r["event"])] = dict(
        est=-b_minus_a, lo=-float(r["ci95_high"]), hi=-float(r["ci95_low"]), n=int(r["n_draws"])
    )
EXP_SHIFT = {
    ("MEM", "likelihood"): (1.15, 1.05), ("MEM", "separate"): (1.41, 2.33),
    ("MEM", "pooled"): (0.69, 1.79), ("MEM", "concord"): (-0.03, 0.07),
    ("PTAU", "likelihood"): (-0.81, -1.21), ("PTAU", "separate"): (-1.26, -1.97),
    ("PTAU", "pooled"): (-0.42, -0.65), ("PTAU", "concord"): (0.04, 0.01),
}
for (ev, m), (a_exp, n_exp) in EXP_SHIFT.items():
    check(f"shift adni {m} {ev}", shift[("adni", m, ev)]["est"], a_exp, TOL_GAP)
    check(f"shift nacc {m} {ev}", shift[("nacc", m, ev)]["est"], n_exp, TOL_GAP)
concord_max = {c: max(abs(shift[(c, "concord", e)]["est"]) for e in S.EVENT_ORDER) for c, _ in COHORTS}
assert max(concord_max.values()) <= 0.09 + 1e-9
assert all(shift[(c, "concord", e)]["lo"] <= 0 <= shift[(c, "concord", e)]["hi"]
           for c, _ in COHORTS for e in S.EVENT_ORDER)

# e: false positives when both groups are drawn from ε3/ε3 carriers (rows keyed N2 / n2 in the inputs)
sd21 = {(r["cohort"], r["null"], r["arm"]): r for r in read_csv(F_SD21)}
prec = {(r["panel"], r["stage"], r["arm"]): r for r in read_csv(F_PREC)}
shp = {(r["cohort"], r["panel"], r["arm"]): r for r in read_csv(F_SHARED)}


def cnt_sd21(cohort, arm):
    r = sd21[(cohort, "N2", arm)]
    return int(r["reject"]), int(r["R"]), float(r["rate"])


def cnt_prec(panel, arm):
    r = prec[(panel, "n2", arm)]
    return int(r["reject"]), int(r["planned"]), float(r["rate"])


def cnt_shared(cohort, panel, arm):
    r = shp[(cohort, panel, arm)]
    return int(r["rejects"]), int(r["R"]), float(r["rejection_rate"])


SETTINGS = [  # (cohort, modality, events, source, unrestricted counts, CONCORD counts, expected U, expected D)
    ("ADNI", "CSF–MRI–cognition", 12, F_SD21.name, cnt_sd21("ADNI", "standard_U"), cnt_sd21("ADNI", "concord_D"),
     0.730, 0.030),
    ("ADNI", "CSF–MRI–cognition", 8, F_PREC.name, cnt_prec("P8", "standard_U"), cnt_prec("P8", "concord_D"),
     0.515, 0.040),
    ("ADNI", "CSF–MRI–cognition", 5, F_PREC.name, cnt_prec("P5", "standard_U"), cnt_prec("P5", "concord_D"),
     0.145, 0.015),
    ("ADNI", "MRI–cognition", 8, F_SHARED.name, cnt_shared("adni", "MRIcog8", "Standard_U"),
     cnt_shared("adni", "MRIcog8", "CONCORD_D"), 0.615, 0.025),
    ("ADNI", "MRI–cognition", 4, F_SHARED.name, cnt_shared("adni", "MRIcog4", "Standard_U"),
     cnt_shared("adni", "MRIcog4", "CONCORD_D"), 0.140, 0.010),
    ("NACC", "MRI–cognition", 8, F_SHARED.name, cnt_shared("nacc", "MRIcog8", "Standard_U"),
     cnt_shared("nacc", "MRIcog8", "CONCORD_D"), 0.045, 0.025),
    ("NACC", "MRI–cognition", 4, F_SHARED.name, cnt_shared("nacc", "MRIcog4", "Standard_U"),
     cnt_shared("nacc", "MRIcog4", "CONCORD_D"), 0.010, 0.020),
    ("NACC", "Amyloid PET–cognition", 5, F_SD21.name, cnt_sd21("NACC", "standard_U"),
     cnt_sd21("NACC", "concord_D"), 0.005, 0.020),
]
fp = []
for cohort, modality, n_ev, src, (ku, nu, ru), (kd, nd, rd), eu, ed in SETTINGS:
    label = f"{modality}, {n_ev} events"
    check(f"FP {cohort} {label} unrestricted", ku / nu, eu, TOL_RATE)
    check(f"FP {cohort} {label} CONCORD", kd / nd, ed, TOL_RATE)
    check(f"FP {cohort} {label} stored rate U", ru, ku / nu, 1e-9)
    check(f"FP {cohort} {label} stored rate D", rd, kd / nd, 1e-9)
    fp.append(dict(cohort=cohort, label=label, tick=f"{modality}\n{n_ev} events", src=src, ku=ku, nu=nu,
                   kd=kd, nd=nd, ru=ku / nu, rd=kd / nd, wu=S.wilson(ku, nu), wd=S.wilson(kd, nd)))
fp.sort(key=lambda d: -d["ru"])
# cohorts must be contiguous after sorting for the bracket labels
order_cohorts = [d["cohort"] for d in fp]
assert order_cohorts == sorted(order_cohorts, key=lambda c: order_cohorts.index(c))
assert max(d["rd"] for d in fp) <= 0.040 + 1e-9

# ----------------------------------------------------------------------------------------
# Source data
# ----------------------------------------------------------------------------------------
DESIGN_A = [  # split, group, CN, MCI, AD
    ("Pool", "ε3/ε3 participants per draw", 160, 96, 64),
    ("Matched composition", "Group 1", 80, 48, 32),
    ("Matched composition", "Group 2", 80, 48, 32),
    ("Different composition", "CN-heavy group", 112, 32, 16),
    ("Different composition", "AD-heavy group", 48, 64, 48),
]
write_csv("Fig4a.csv", ["split", "group", "n_CN", "n_MCI", "n_AD", "n_total", "pct_CN", "draws_per_cohort",
                        "cohorts"],
          [(s, g, a, b, c, a + b + c, round(100 * a / (a + b + c), 1), 200, "ADNI; NACC (analyzed separately)")
           for s, g, a, b, c in DESIGN_A])
rows = []
for c, cl in COHORTS:
    for m in MODELS:
        for dsg in ("matched", "different"):
            g = gap[(c, m, dsg)]
            rows.append((cl, M[m]["label"], dsg, g["gap"], g["lo"], g["hi"],
                         g["dmm"] if g["dmm"] is not None else "", g["dmm_lo"] if g["dmm_lo"] is not None else "",
                         g["dmm_hi"] if g["dmm_hi"] is not None else ""))
write_csv("Fig4b.csv", ["cohort", "model", "design", "mean_position_gap", "boot_ci95_low", "boot_ci95_high",
                        "paired_increment_different_minus_matched", "paired_increment_ci95_low",
                        "paired_increment_ci95_high"], rows)
write_csv("Fig4c.csv", ["cohort", "model", "design", "n_draws", "discordant_event_pairs_pct", "ci95_low_pct",
                        "ci95_high_pct"],
          [(cl, M[m]["label"], "matched", dist[(c, m)]["n"], 100 * dist[(c, m)]["mean"], 100 * dist[(c, m)]["lo"],
            100 * dist[(c, m)]["hi"]) for c, cl in COHORTS for m in MODELS])
write_csv("Fig4d.csv", ["cohort", "model", "event", "n_draws", "position_shift_CNheavy_minus_ADheavy", "ci95_low",
                        "ci95_high"],
          [(cl, M[m]["label"], S.EVENT_LABEL[e], shift[(c, m, e)]["n"], round(shift[(c, m, e)]["est"], 6),
            round(shift[(c, m, e)]["lo"], 6), round(shift[(c, m, e)]["hi"], 6))
           for c, cl in COHORTS for m in MODELS for e in S.EVENT_ORDER])
rows = []
for d in fp:
    rows.append((d["cohort"], d["label"], "Separately fitted DEBM, unrestricted permutation", d["ku"], d["nu"],
                 100 * d["ru"], 100 * d["wu"][0], 100 * d["wu"][1], d["src"]))
    rows.append((d["cohort"], d["label"], "CONCORD, within-diagnosis permutation", d["kd"], d["nd"],
                 100 * d["rd"], 100 * d["wd"][0], 100 * d["wd"][1], d["src"]))
write_csv("Fig4e.csv", ["cohort", "setting", "analysis", "rejections", "comparisons", "false_positive_pct",
                        "wilson95_low_pct", "wilson95_high_pct", "source_file"], rows)

# ----------------------------------------------------------------------------------------
# Figure scaffold (all positions in mm from the top-left corner)
# ----------------------------------------------------------------------------------------
FW, FH = 180.0, 150.0
fig = plt.figure(figsize=(FW * S.MM, FH * S.MM))


def fx(x):
    return x / FW


def fy(y):
    return 1 - y / FH


L1 = 10.5         # letter/title baseline, top row (as Figs 2 and 3)
L2 = 69.0         # letter/title baseline, bottom row
TITLE_DX = 4.5    # title starts ~3 mm right of the letter (as in Figs 1-3 and 5)

S.method_legend(fig, MODELS, y_top=1.0)


def letter_title(ch, x, y, title):
    S.panel_label_mm(fig, ch, x=x, y_top=y)
    S.panel_title_mm(fig, title, x=x + TITLE_DX, y_top=y)


# ----------------------------------------------------------------------------------------
# a  Design schematic
# ----------------------------------------------------------------------------------------
A_X, A_Y, A_W, A_H = 1.0, 13.0, 39.5, 50.0
axa = S.add_axes_mm(fig, A_X, A_Y, A_W, A_H)
axa.set_xlim(0, A_W)
axa.set_ylim(A_H, 0)
axa.axis("off")
letter_title("a", 1.0, L1, "Design")

BAR_X0, BAR_W, BAR_H, SPLIT_GAP = 3.6, 35.4, 3.0, 1.2
HALF = (BAR_W - SPLIT_GAP) / 2
TXT_ON = S.DIAG_TEXT


def stack(ax, x0, y0, w, counts, seg_labels=False, pct_cn=False):
    tot = sum(counts)
    x = x0
    for d, n in zip(("CN", "MCI", "AD"), counts):
        seg = w * n / tot
        ax.add_patch(Rectangle((x, y0), seg, BAR_H, facecolor=S.DIAG[d], edgecolor="white", lw=0.4))
        if seg_labels:
            ax.text(x + seg / 2, y0 + BAR_H / 2 + 0.05, d, ha="center", va="center", fontsize=5.5, color=TXT_ON[d])
        x += seg
    if pct_cn:  # % CN inside the CN segment, as in Fig. 1c
        ax.text(x0 + 0.5, y0 + BAR_H / 2 + 0.05, f"{100 * counts[0] / tot:.0f}%", ha="left", va="center",
                fontsize=5.5, color=INK)
    ax.add_patch(Rectangle((x0, y0), w, BAR_H, facecolor="none", edgecolor=INK, lw=0.4))


axa.text(0.7, 0.2, "ε3/ε3 participants", fontsize=6, fontweight="bold", va="top", color=INK)
axa.text(0.7, 2.9, r"$\mathit{n}$ = 320 per draw, 200 draws", fontsize=6, va="top", color=GREY)
POOL_Y = 6.4
stack(axa, BAR_X0, POOL_Y, BAR_W, (160, 96, 64), seg_labels=True)


def split_block(ax, y_head, title, groups):
    ax.text(BAR_X0, y_head, title, fontsize=6, fontweight="bold", va="top", color=INK)
    yb = y_head + 3.0
    for i, (name, counts) in enumerate(groups):
        x0 = BAR_X0 + i * (HALF + SPLIT_GAP)
        stack(ax, x0, yb, HALF, counts, pct_cn=True)
        ax.text(x0 + HALF / 2, yb + BAR_H + 0.6, name, fontsize=6, ha="center", va="top", color=INK)
        ax.text(x0 + HALF / 2, yb + BAR_H + 3.1, " / ".join(str(c) for c in counts), fontsize=5.5,
                ha="center", va="top", color=GREY)
    return yb + BAR_H / 2


y_m = split_block(axa, 13.6, "Matched composition", [("Group 1", (80, 48, 32)), ("Group 2", (80, 48, 32))])
y_d = split_block(axa, 29.2, "Different composition",
                  [("CN-heavy group", (112, 32, 16)), ("AD-heavy group", (48, 64, 48))])
# tree connector: the same 320 participants split two ways
TREE_X = 1.4
axa.plot([TREE_X, TREE_X], [POOL_Y + BAR_H / 2, y_d], color=S.REF_LINE["color"], lw=0.6, solid_capstyle="butt")
axa.plot([TREE_X, BAR_X0 - 0.9], [POOL_Y + BAR_H / 2] * 2, color=S.REF_LINE["color"], lw=0.6)
for yy in (y_m, y_d):
    axa.annotate("", xy=(BAR_X0 - 0.3, yy), xytext=(TREE_X, yy),
                 arrowprops=dict(arrowstyle="-|>", lw=0.6, color=S.REF_LINE["color"], mutation_scale=4.5,
                                 shrinkA=0, shrinkB=0))
axa.text(0.7, A_H - 0.1, "Same participants in both splits;\nADNI and NACC analyzed separately",
         fontsize=5.5, va="bottom", color=GREY, linespacing=1.25)

# ----------------------------------------------------------------------------------------
# b  Mean-position gap (dumbbell with composition-driven increase) and c  discordance (bars),
#    shared rows
# ----------------------------------------------------------------------------------------
BC_TOP, BC_H = 14.0, 42.0
B_X, B_W = 69.0, 56.5
C_X, C_W = 131.5, 40.0
axb = S.add_axes_mm(fig, B_X, BC_TOP, B_W, BC_H)
axc = S.add_axes_mm(fig, C_X, BC_TOP, C_W, BC_H)
letter_title("b", 43.2, L1, "Matched vs different composition")
letter_title("c", 128.0, L1, "Matched composition")

BLOCK_OFF = {"adni": 0.0, "nacc": 5.0}
ROWS = {(c, m): BLOCK_OFF[c] + i for c, _ in COHORTS for i, m in enumerate(MODELS)}
HEAD_Y = {"adni": -0.9, "nacc": 4.1}
Y_LIM = (8.5, -1.3)
for ax in (axb, axc):
    ax.set_ylim(*Y_LIM)
    ax.axhline(HEAD_Y["nacc"], color=SEP, lw=0.5, zorder=0)
    ax.tick_params(axis="y", length=0, pad=4.5)
    ax.spines["left"].set_visible(False)
# CONCORD call-out rows (b: band runs under the row label)
B_BAND_X0 = -(B_X - 43.2) / B_W
for c, _ in COHORTS:
    S.concord_band(axb, ROWS[(c, "concord")], 0.84, xmin=B_BAND_X0, clip_on=False)
    S.concord_band(axc, ROWS[(c, "concord")], 0.84)

BAND_LW = 3.6  # composition-driven increase (pt); label sits just above the band
tr_above = offset_copy(axb.transData, fig=fig, y=BAND_LW / 2 + 0.9, units="points")
for c, _ in COHORTS:
    for m in MODELS:
        y = ROWS[(c, m)]
        col, mk = M[m]["color"], M[m]["marker"]
        is_c = m == "concord"
        msz = S.ms(m, 4.3 if is_c else 3.8)
        g0, g1 = gap[(c, m, "matched")], gap[(c, m, "different")]
        z = 10 if is_c else 3
        axb.plot([g0["gap"], g1["gap"]], [y, y], color=col, alpha=0.30, lw=BAND_LW, solid_capstyle="butt",
                 zorder=z)
        axb.plot([g1["lo"], g1["hi"]], [y, y], color=col, lw=0.9 if is_c else 0.6, zorder=z + 1,
                 solid_capstyle="butt")
        axb.plot(g1["gap"], y, marker=mk, ms=msz, mfc=col, mec="white", mew=0.4, ls="none", zorder=z + 3,
                 clip_on=False)
        # hollow = matched composition; for CONCORD the two points nearly coincide, so the ring is
        # transparent and drawn on top so that both remain visible
        axb.plot(g0["gap"], y, marker=mk, ms=msz, mfc="none" if is_c else "white", mec=col,
                 mew=0.9 if is_c else 0.8, ls="none", zorder=z + 4 if is_c else z + 2, clip_on=False)
        if is_c:
            axb.text(max(g1["hi"], g0["gap"]) + 0.05, y, fmt_inc(inc[(c, m)]), color=col, fontsize=6,
                     fontweight="bold", va="center", ha="left", zorder=20)
        else:
            axb.text((g0["gap"] + g1["gap"]) / 2, y, fmt_inc(inc[(c, m)]), transform=tr_above,
                     color=S.TEXT_COLOR[m], fontsize=6, va="baseline", ha="center", zorder=20)
B_XMAX = 1.9
axb.set_xlim(0, B_XMAX)
axb.set_xticks([0, 0.5, 1.0, 1.5])
axb.set_xticklabels(["0", "0.5", "1.0", "1.5"])
axb.set_xlabel(S.LABEL["gap"])
axb.set_yticks([ROWS[k] for k in ROWS])
axb.set_yticklabels([M[m]["label"] for c, _ in COHORTS for m in MODELS])
for tl, (c, m) in zip(axb.get_yticklabels(), ROWS):
    if m == "concord":
        tl.set_color(M["concord"]["color"])
        tl.set_fontweight("bold")
for c, cl in COHORTS:
    axb.text(-0.02, HEAD_Y[c], cl, transform=axb.get_yaxis_transform(), ha="right", va="center",
             fontsize=6.5, fontweight="bold", color=INK)
# encoding key inside b (empty upper right of the ADNI block)
KX, KT = 1.07, 0.085  # key handle x (data) and text offset
k_ys = [-0.95, -0.28, 0.39]
axb.plot(KX, k_ys[0], marker="o", ms=3.8, mfc="white", mec=GREY, mew=0.8, ls="none", clip_on=False)
axb.plot(KX, k_ys[1], marker="o", ms=3.8, mfc=GREY, mec="white", mew=0.4, ls="none", clip_on=False)
axb.plot([KX - 0.045, KX + 0.045], [k_ys[2], k_ys[2]], color=GREY, alpha=0.30, lw=BAND_LW, solid_capstyle="butt",
         clip_on=False)
for yk, txt in zip(k_ys, ["Matched composition", "Different composition", "Increase due to\ncomposition"]):
    axb.text(KX + KT, yk, txt, fontsize=6, color=INK, va="top" if "\n" in txt else "center", ha="left",
             linespacing=1.1, transform=offset_copy(axb.transData, fig=fig, y=2.1 if "\n" in txt else 0,
                                                     units="points"))

# c
for c, _ in COHORTS:
    for m in MODELS:
        y = ROWS[(c, m)]
        d = dist[(c, m)]
        is_c = m == "concord"
        axc.barh(y, 100 * d["mean"], height=0.62, color=S.FILL[m], edgecolor=M[m]["color"],
                 linewidth=0.6, zorder=3 if is_c else 2)
        axc.plot([100 * d["lo"], 100 * d["hi"]], [y, y], color=INK, lw=S.ERRBAR["elinewidth"], zorder=4,
                 solid_capstyle="butt")
        axc.text(100 * d["hi"] + 1.4, y, S.fmt_pct(100 * d["mean"]),
                 color=S.TEXT_COLOR[m], fontsize=6 if is_c else 5.5,
                 fontweight="bold" if is_c else "normal", va="center", ha="left")
axc.set_xlim(0, 60)
axc.set_xticks([0, 20, 40, 60])
axc.set_yticks([])
axc.set_xlabel(S.LABEL["distance"] + " (%)")

# ----------------------------------------------------------------------------------------
# d  Per-event position shifts, CN-heavy - AD-heavy group
# ----------------------------------------------------------------------------------------
D_TOP, D_BOT = L2 + 8.2, FH - 8.6
D_H = D_BOT - D_TOP
D_X1, D_W, D_GAP = 19.5, 39.5, 3.0
D_X2 = D_X1 + D_W + D_GAP
axd = [S.add_axes_mm(fig, D_X1, D_TOP, D_W, D_H), S.add_axes_mm(fig, D_X2, D_TOP, D_W, D_H)]
letter_title("d", 1.0, L2, "Each event, different composition")

GROUP_GAP = 0.35
EV_Y = {e: i + (GROUP_GAP if i >= 3 else 0.0) for i, e in enumerate(S.EVENT_ORDER)}
DODGE = 0.225
D_OFF = {m: (k - 1.5) * DODGE for k, m in enumerate(MODELS)}
D_YLIM = (5 + GROUP_GAP + 0.5, -0.5)
D_XLIM = (-2.45, 2.75)
for j, (c, cl) in enumerate(COHORTS):
    ax = axd[j]
    ax.set_ylim(*D_YLIM)
    ax.set_xlim(*D_XLIM)
    for e in S.EVENT_ORDER:
        ax.axhspan(EV_Y[e] - 0.45, EV_Y[e] + 0.45, color=ROWBAND, lw=0, zorder=0)
    ax.axvline(0, **S.REF_LINE)
    for m in MODELS:
        is_c = m == "concord"
        col, mk = M[m]["color"], M[m]["marker"]
        for e in S.EVENT_ORDER:
            s = shift[(c, m, e)]
            y = EV_Y[e] + D_OFF[m]
            z = 6 if is_c else 3
            ax.plot([s["lo"], s["hi"]], [y, y], color=col, lw=0.9 if is_c else 0.6, zorder=z,
                    solid_capstyle="butt")
            ax.plot(s["est"], y, marker=mk, ms=S.ms(m, 4.0 if is_c else 3.5), mfc=col, mec="white",
                    mew=0.35, ls="none", zorder=z + 1)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xticklabels(["−2", "−1", "0", "1", "2"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_yticks([EV_Y[e] for e in S.EVENT_ORDER])
    if j == 0:
        ax.set_yticklabels([S.EVENT_LABEL[e] for e in S.EVENT_ORDER])
    else:
        ax.set_yticklabels([])
    # cohort header line: cohort name (left) and the CONCORD call-out (right)
    x_left = D_X1 if j == 0 else D_X2
    fig.text(fx(x_left), fy(D_TOP - 1.3), cl, fontsize=6.5, fontweight="bold", color=INK, ha="left",
             va="baseline")
    fig.text(fx(x_left + D_W), fy(D_TOP - 1.3), f"CONCORD all |shift| ≤ {S.fmt_pos(concord_max[c])}",
             fontsize=6, fontweight="bold", color=M["concord"]["color"], ha="right", va="baseline")
# direction cue (positive = later in the CN-heavy group), once for both cohorts
fig.text(fx(D_X1 + D_W + D_GAP / 2), fy(L2 + 3.7), "← Earlier in CN-heavy group  ·  Later in CN-heavy group →",
         fontsize=6, color=GREY, ha="center", va="baseline")
fig.text(fx(D_X1 + D_W + D_GAP / 2), fy(D_BOT + 4.5), S.LABEL["shift"], ha="center", va="top", fontsize=7,
         color=INK)
# CSF / Cognition brackets left of the event labels
axl = axd[0]
BR_X = (7.4 - D_X1) / D_W  # axes fraction
for name, evs, col in [("CSF", S.EVENT_ORDER[:3], S.EVENT_COLOR["ABETA"]),
                       ("Cognition", S.EVENT_ORDER[3:], S.EVENT_COLOR["MEM"])]:
    y0, y1 = EV_Y[evs[0]] - 0.32, EV_Y[evs[-1]] + 0.32
    tr = axl.get_yaxis_transform()
    axl.plot([BR_X + 0.03, BR_X, BR_X, BR_X + 0.03], [y0, y0, y1, y1], transform=tr, color=col, lw=0.7,
             clip_on=False, solid_capstyle="butt")
    axl.text(BR_X - 0.03, (y0 + y1) / 2, name, transform=tr, rotation=90, ha="right", va="center",
             fontsize=6.5, color=col, fontweight="bold")

# ----------------------------------------------------------------------------------------
# e  False positives in constructed null comparisons (log10 axis, as Fig. 3a)
# ----------------------------------------------------------------------------------------
E_LET = 104.5
E_X = 133.5
E_W = FW - 3.8 - E_X
axe = S.add_axes_mm(fig, E_X, D_TOP, E_W, D_H)
letter_title("e", E_LET, L2, "Groups drawn from ε3/ε3 carriers")
E_GAP = 0.8  # white gap between the ADNI and NACC blocks; holds the CONCORD call-out
e_y = []
prev = None
off = 0.0
for i, d in enumerate(fp):
    if prev is not None and d["cohort"] != prev:
        off += E_GAP
    e_y.append(i + off)
    prev = d["cohort"]
col_u, col_d = M["separate"]["color"], M["concord"]["color"]
mk_u = M["separate"]["marker"]
ms_u, ms_d = S.ms("separate", 3.8), S.ms("concord", 4.2)
E_XLIM = (0.25, 100)
axe.set_xscale("log")
axe.axvline(5, **S.REF_LINE)
E_DODGE = 0.17


def e_interval(lo, hi, yy, col, lw, z):
    """Wilson interval in %; a lower bound below the axis minimum is drawn from the axis edge and
    marked there by a small left-pointing triangle (truncation cue)."""
    axe.plot([max(lo, E_XLIM[0]), hi], [yy, yy], color=col, lw=lw, zorder=z, solid_capstyle="butt")
    if lo < E_XLIM[0]:
        # arrowhead with its tip on the axis edge, pointing left (body inside the axes)
        axe.plot(E_XLIM[0], yy, marker=[(0, 0), (1, 0.55), (1, -0.55), (0, 0)], ms=4.6, color=col, mew=0,
                 ls="none", clip_on=False, zorder=z)
        return 1
    return 0


n_trunc = 0
for y, d in zip(e_y, fp):
    axe.axhspan(y - 0.42, y + 0.42, color=ROWBAND, lw=0, zorder=0)
    yu, yd = y - E_DODGE, y + E_DODGE
    n_trunc += e_interval(100 * d["wu"][0], 100 * d["wu"][1], yu, col_u, S.ERRBAR["elinewidth"], 3)
    axe.plot(100 * d["ru"], yu, marker=mk_u, ms=ms_u, mfc="white", mec=col_u, mew=0.8, ls="none", zorder=4)
    n_trunc += e_interval(100 * d["wd"][0], 100 * d["wd"][1], yd, col_d, 0.9, 5)
    axe.plot(100 * d["rd"], yd, marker="o", ms=ms_d, mfc=col_d, mec="white", mew=0.4, ls="none", zorder=6)
assert n_trunc == 1  # only NACC amyloid PET, separately fitted DEBM (1/200; Wilson lower bound 0.09%)
E_YLIM = (e_y[-1] + 0.55, -0.8)
axe.set_ylim(*E_YLIM)
axe.set_xlim(*E_XLIM)
axe.xaxis.set_major_locator(FixedLocator([0.5, 1, 2, 5, 10, 20, 50, 100]))
axe.xaxis.set_minor_locator(NullLocator())
axe.set_xticklabels(["0.5", "1", "2", "5", "10", "20", "50", "100"])
fig.text(fx(E_X + E_W / 2), fy(D_BOT + 4.5), "False-positive rate (%, log scale)", ha="center", va="top", fontsize=7,
         color=INK)
axe.set_yticks(e_y)
axe.set_yticklabels([d["tick"] for d in fp], fontsize=6, linespacing=1.05)
axe.tick_params(axis="y", length=0)
axe.spines["left"].set_visible(False)
axe.text(5 * 1.07, E_YLIM[1] + 0.02, "Nominal 5%", fontsize=6, color=GREY, ha="left", va="top")
top_u = fp[0]
axe.text(100 * top_u["wu"][0] / 1.12, e_y[0] - E_DODGE, S.fmt_pct(100 * top_u["ru"]), color=S.TEXT_COLOR["separate"],
         fontsize=6, fontweight="bold", ha="right", va="center")
# cohort brackets
EBR_X = (E_LET + 3.6 - E_X) / E_W
for coh in ("ADNI", "NACC"):
    ys = [y for y, d in zip(e_y, fp) if d["cohort"] == coh]
    y0, y1 = min(ys) - 0.34, max(ys) + 0.34
    tr = axe.get_yaxis_transform()
    axe.plot([EBR_X + 0.03, EBR_X, EBR_X, EBR_X + 0.03], [y0, y0, y1, y1], transform=tr, color=INK, lw=0.6,
             clip_on=False, solid_capstyle="butt")
    axe.text(EBR_X - 0.035, (y0 + y1) / 2, coh, transform=tr, rotation=90, ha="right", va="center",
             fontsize=6.5, fontweight="bold", color=INK)
# key (marker-only handles)
ekey = [Line2D([0], [0], ls="none", marker=mk_u, ms=ms_u, mfc="white", mec=col_u, mew=0.8),
        Line2D([0], [0], ls="none", marker="o", ms=ms_d, mfc=col_d, mec="white", mew=0.4)]
eleg = fig.legend(ekey, ["Separately fitted DEBM, unrestricted permutation", "CONCORD, within-diagnosis permutation"],
                  loc="upper left", bbox_to_anchor=(fx(E_LET + TITLE_DX - 0.4), fy(L2 + 1.5)), fontsize=6,
                  handlelength=1.0, handletextpad=0.5, labelspacing=0.3, frameon=False, borderaxespad=0, borderpad=0)
eleg.get_texts()[1].set_color(col_d)
eleg.get_texts()[1].set_fontweight("bold")
# CONCORD call-out in the white gap between the ADNI and NACC blocks, right-aligned to the axes edge
i_gap = [d["cohort"] for d in fp].index("NACC")
y_call = (e_y[i_gap - 1] + e_y[i_gap]) / 2
axe.text(1.0, y_call, f"CONCORD ≤ {S.fmt_pct(100 * max(d['rd'] for d in fp))}\nin all eight settings", color=col_d,
         fontsize=6, fontweight="bold", ha="right", va="center", linespacing=1.15,
         transform=axe.get_yaxis_transform())

paths = S.save(fig, STEM)
print("\n".join(str(p) for p in paths))
print("concord max |shift|:", concord_max)
print("increments:", {k: fmt_inc(v) for k, v in inc.items()})
