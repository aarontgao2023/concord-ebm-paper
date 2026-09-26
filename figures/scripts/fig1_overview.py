"""Fig. 1 | Overview of CONCORD (Nature Communications, 180 x 140 mm).

Top row (the method):
a  Separately fitted DEBM (grey lane) versus CONCORD (blue lane), step by step (conceptual).
b  Within-diagnosis permutation schematic (conceptual).
Bottom row (why it is needed):
c  Diagnostic composition (CN/MCI/AD) of the APOE groups in ADNI and NACC (the counts of Table 2)
   and of the three simulated groups (simulation design constants).
d  Three-event example: expected ranking loss of all six orderings under the CN-heavy and the
   AD-heavy simulation composition (same true sequence A->B->C, one abnormality model).

Inputs:
  figures/inputs/simulation/toy_target_shift.txt   output of simulation/scripts/dev/toy_target_shift.py
  simulation/scripts/dev/toy_target_shift.py       read for the two compositions of d
  figures/inputs/realdata/reference_counts_weights_ESS.csv   optional; if present, the Table 2 counts
                                                   used in c are checked against it

Only stored aggregate files are read; nothing is re-fitted. Everything outside the data axes of
c and d is drawn on one full-figure canvas in mm, measured from the top-left corner.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

ROOT = S.ROOT
COUNTS = S.REAL_INPUTS / "reference_counts_weights_ESS.csv"   # optional cross-check of panel c
TOY = S.SIM_INPUTS / "toy_target_shift.txt"
TOY_PY = ROOT / "simulation/scripts/dev/toy_target_shift.py"
S.require_inputs([TOY, TOY_PY], "Fig. 1")
S.SRC.mkdir(parents=True, exist_ok=True)

W_MM, H_MM = 180.0, 140.0
INK = S.PAL["ink"]
BLUE = S.PAL["blue_main"]
BLUE_BG = S.PAL["blue_bg"]
BLUE_EDGE = "#B7CBE6"
GREY_DARK = S.PAL["grey_dark"]
GREY_MID = S.PAL["grey_mid"]
GREY_LINE = "#9A9A9A"
LANE_BG = "#EFEFEF"
GOLD = S.PAL["highlight"]
G1, G2 = S.GROUP[1], S.GROUP[2]
DIAGS = ["CN", "MCI", "AD"]

# ----------------------------------------------------------------------------------------------
# Data: panel c
# ----------------------------------------------------------------------------------------------
# CN/MCI/AD counts of the APOE groups in the six-event ADNI and NACC analysis populations, as
# printed in Table 2. When the aggregate count table of the ADNI/NACC analysis is available, the
# counts are checked against it.
EXPECTED_A = {
    ("adni", "e2"): (50, 39, 8), ("adni", "e33"): (222, 273, 65), ("adni", "e4"): (98, 299, 149),
    ("nacc", "e2"): (106, 15, 22), ("nacc", "e33"): (480, 104, 175), ("nacc", "e4"): (281, 94, 307),
}
cnt = {(c, g, d): n for (c, g), ns in EXPECTED_A.items() for d, n in zip(DIAGS, ns)}
if COUNTS.is_file():
    stored = {}
    with open(COUNTS) as fh:
        for r in csv.DictReader(fh):
            if r["application"] != "actual_APOE":
                continue
            stored[(r["cohort"], r["group"], r["diagnosis"])] = int(r["cell_n"])
    for key, exp in EXPECTED_A.items():
        got = tuple(stored[(key[0], key[1], d)] for d in DIAGS)
        assert got == exp, (key, got, exp)
    print(f"[check] panel c counts agree with {S.source_name(COUNTS)}")
else:
    print(f"[note] {S.source_name(COUNTS)} not found; panel c uses the Table 2 counts in this script")

# simulation design constants (group sizes and composition of a published ADNI APOE comparison)
SIM = [("CN-heavy group", (57, 6, 12)), ("Intermediate group", (244, 66, 101)), ("AD-heavy group", (110, 156, 219))]

blocks = [
    ("ADNI", [(S.GENOTYPE_LABEL[g], tuple(cnt[("adni", g, d)] for d in DIAGS)) for g in ("e2", "e33", "e4")]),
    ("NACC", [(S.GENOTYPE_LABEL[g], tuple(cnt[("nacc", g, d)] for d in DIAGS)) for g in ("e2", "e33", "e4")]),
    ("Simulation", SIM),
]

# key sentence check: % CN in e2 vs e4
pcn = {b: [c[0] / sum(c) for _, c in rows] for b, rows in blocks}
for b, exp in (("ADNI", (0.515, 0.179)), ("NACC", (0.741, 0.412))):
    assert abs(pcn[b][0] - exp[0]) < 0.002 and abs(pcn[b][2] - exp[1]) < 0.002, (b, pcn[b])

with open(S.SRC / "Fig1c.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["source", "group", "diagnosis", "n", "group_n", "proportion"])
    for b, rows in blocks:
        for g, c in rows:
            for d, n in zip(DIAGS, c):
                wr.writerow([b, g, d, n, sum(c), f"{n / sum(c):.6f}"])

# ----------------------------------------------------------------------------------------------
# Data: panel d
# ----------------------------------------------------------------------------------------------
txt = TOY.read_text()
first = txt.split("===")[1]
auc = dict(zip("ABC", [float("0" + v) for v in re.findall(r"\((?:AUC )?(\.\d+)\)", first.splitlines()[0])]))
assert auc == {"A": 0.76, "B": 0.96, "C": 0.86}, auc


def parse_losses(tag):
    line = next(ln for ln in first.splitlines() if ln.strip().startswith(tag))
    return {k: float(v) for k, v in re.findall(r"'([ABC]<[ABC]<[ABC])': np\.float64\(([0-9.]+)\)", line)}


LOSS = {"CN-heavy composition": parse_losses("e2 composition"), "AD-heavy composition": parse_losses("e4 composition")}
ORDERS = ["A<B<C", "A<C<B", "B<A<C", "B<C<A", "C<A<B", "C<B<A"]
TRUE = "A<B<C"
best = {k: min(v, key=v.get) for k, v in LOSS.items()}
assert best["CN-heavy composition"] == "A<C<B" and abs(LOSS["CN-heavy composition"]["A<C<B"] - 0.2921) < 0.002
assert abs(LOSS["CN-heavy composition"][TRUE] - 0.3524) < 0.002
assert best["AD-heavy composition"] == "B<A<C" and abs(LOSS["AD-heavy composition"]["B<A<C"] - 0.3809) < 0.002
assert abs(LOSS["AD-heavy composition"][TRUE] - 0.3892) < 0.002

# the toy example's two compositions are the simulation CN-heavy and AD-heavy groups of panel c
toy_comp = {m[0]: tuple(int(v) for v in m[1:])
            for m in re.findall(r"\('(e2|e4)', \((\d+), (\d+), (\d+)\)\)", TOY_PY.read_text())}
assert toy_comp == {"e2": SIM[0][1], "e4": SIM[2][1]}, toy_comp
COMP = {"CN-heavy composition": toy_comp["e2"], "AD-heavy composition": toy_comp["e4"]}

with open(S.SRC / "Fig1d.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["composition", "ordering", "expected_ranking_loss", "best_fitting", "true_sequence",
                 "n_CN", "n_MCI", "n_AD", "AUC_A", "AUC_B", "AUC_C"])
    for comp, d in LOSS.items():
        for o in ORDERS:
            wr.writerow([comp, o.replace("<", "->"), f"{d[o]:.4f}", int(o == best[comp]), int(o == TRUE),
                         *COMP[comp], *(f"{auc[e]:.2f}" for e in "ABC")])

# panel a, box 2: illustrative compositions and the resulting weights (schematic values)
P1 = np.array([6, 3, 1]) / 10.0      # illustrative CN-heavy group
P2 = np.array([1, 3, 6]) / 10.0      # illustrative AD-heavy group
PC = (P1 + P2) / 2.0                 # common composition 3.5/3/3.5
W1, W2 = PC / P1, PC / P2
assert np.allclose(W1, [0.35 / 0.6, 1.0, 3.5]) and np.allclose(W2, [3.5, 1.0, 0.35 / 0.6])
with open(S.SRC / "Fig1a_illustrative.csv", "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["group", "diagnosis", "observed_share", "common_share", "weight"])
    for gname, P, Wt in (("Group 1", P1, W1), ("Group 2", P2, W2)):
        for d, p, pc, w in zip(DIAGS, P, PC, Wt):
            wr.writerow([gname, d, f"{p:.2f}", f"{pc:.2f}", f"{w:.3f}"])

# ==============================================================================================
# Figure and mm canvas (x right, y down, origin top-left)
# ==============================================================================================
fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))
cv = fig.add_axes([0, 0, 1, 1], zorder=0)
cv.set_xlim(0, W_MM)
cv.set_ylim(H_MM, 0)
cv.axis("off")
_R = fig.canvas.get_renderer()
PX_MM = fig.dpi / 25.4
CAP = 0.716 * 25.4 / 72          # Arial cap height per pt, in mm


def cap_mm(size):
    return CAP * size


def text_w(s, size, weight="normal"):
    t = cv.text(0, 0, s, fontsize=size, fontweight=weight)
    w = t.get_window_extent(_R).width / PX_MM
    t.remove()
    return w


def data_axes(x, y_top, w, h):
    ax = S.add_axes_mm(fig, x, y_top, w, h, zorder=1)
    ax.set_facecolor("none")
    return ax


def rbox(x, y, w, h, fc=BLUE_BG, ec=BLUE_EDGE, lw=0.5, r=1.6, z=1):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=z)
    cv.add_patch(p)
    return p


def arrow(x0, y0, x1, y1, color=BLUE, lw=0.8, ms=5, style="-|>", rad=0.0, z=4):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                        connectionstyle=f"arc3,rad={rad}", zorder=z, shrinkA=0, shrinkB=0)
    cv.add_patch(a)
    return a


def mk(x, y, g, ms=3.4, z=5, mec="white", mew=0.3):
    cv.plot(x, y, g["marker"], ms=ms * S.MARKER_SCALE[g["marker"]], mfc=g["color"], mec=mec, mew=mew, zorder=z)


def chip(xc, yc, w, h=2.5, z=1.5):
    cv.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.45",
                                fc=GOLD, ec="none", zorder=z))


def comp_bar(x, y, w, h, props, labels=False, pct=False, z=3, lw=0.3):
    left = x
    for dname, p in zip(DIAGS, props):
        cv.add_patch(Rectangle((left, y), w * p, h, fc=S.DIAG[dname], ec=INK, lw=lw, zorder=z))
        if labels:
            cv.text(left + w * p / 2, y + h / 2 + cap_mm(5) / 2, dname, fontsize=5, ha="center", va="baseline",
                    color=S.DIAG_TEXT[dname], zorder=z + 1)
        left += w * p
    if pct:
        cv.text(x + 0.6, y + h / 2 + cap_mm(5.5) / 2, S.fmt_pct(100 * props[0]), fontsize=5.5, ha="left",
                va="baseline", color=INK, zorder=z + 1)


TITLE_DX = 4.5           # title starts this far right of its letter (as in Figs 2-5)
# Two rows; every y below is an offset from its row's letter baseline (method row on top)
ROW_M = 3.6              # letter baseline, top row (a, b: the method); 3.6 keeps the 'b' ascender >= 1.5 mm from the trim
ROW_D = 79.5             # letter / title baseline, bottom row (c, d: why it is needed)
AX_BOTTOM = ROW_D + 51.5  # common x-axis baseline of c and d
XLAB_TOP = ROW_D + 55.8   # common top of the x-axis titles of c and d

# ==============================================================================================
# Panel c — diagnostic composition
# ==============================================================================================
S.panel_title_mm(fig, "Diagnostic composition", x=1.0 + TITLE_DX, y_top=ROW_D)
lx = 36.8
for d in DIAGS:
    cv.add_patch(Rectangle((lx, ROW_D - 1.75), 2.6, 1.75, fc=S.DIAG[d], ec=INK, lw=0.3, zorder=2))
    cv.text(lx + 3.3, ROW_D, d, fontsize=6, va="baseline", ha="left", color=INK)
    lx += 3.3 + text_w(d, 6) + 2.6

A_X, A_W, A_TOP = 23.5, 36.5, ROW_D + 4.0
ax_a = data_axes(A_X, A_TOP, A_W, AX_BOTTOM - A_TOP)
y = 0.0
ypos, heads = [], []
for b, rows in blocks:
    hy = 0.0 if not heads else y - 1.15      # header sits close to its own bars
    heads.append((b, hy))
    y = hy - 0.85
    for k, (g, c) in enumerate(rows):
        if k:
            y -= 1.0
        ypos.append((y, g, c))
BH = 0.66
for yy, g, c in ypos:
    tot = sum(c)
    left = 0.0
    for d, n in zip(DIAGS, c):
        wdt = 100.0 * n / tot
        ax_a.barh(yy, wdt, left=left, height=BH, color=S.DIAG[d], edgecolor=INK, linewidth=0.3, zorder=2)
        left += wdt
    ax_a.text(1.6, yy, S.fmt_pct(100.0 * c[0] / tot), va="center", ha="left", fontsize=5.5, color=INK, zorder=3)
    ax_a.text(101.5, yy, f"{tot}", va="center", ha="left", fontsize=5.5, color=INK, clip_on=False)
ax_a.set_yticks([p[0] for p in ypos])
ax_a.set_yticklabels([p[1] for p in ypos])
ax_a.tick_params(axis="y", length=0, pad=2)
ax_a.spines["left"].set_visible(False)
ax_a.set_xlim(0, 100)
ax_a.set_ylim(y - 0.45, heads[0][1] + 0.35)
ax_a.set_xticks([0, 25, 50, 75, 100])
from matplotlib import transforms  # noqa: E402

trans_a = transforms.blended_transform_factory(fig.transFigure, ax_a.transData)
for b, yy in heads:
    ax_a.text((1.0 + TITLE_DX) / W_MM, yy, b, transform=trans_a, fontsize=6.5, fontweight="bold", color=INK,
              va="center", ha="left")
    if b == "Simulation":
        ax_a.text(0.0, yy, "design from a published ADNI comparison", fontsize=5.5, color=GREY_MID,
                  va="center", ha="left")
ax_a.text(101.5, heads[0][1], "n", fontsize=5.5, style="italic", color=INK, va="center", ha="left")
cv.text(A_X + A_W / 2, XLAB_TOP, "Participants (%)", fontsize=7, color=INK, ha="center", va="top")

# ==============================================================================================
# Panel d — three-event example (two small multiples on a shared scale)
# ==============================================================================================
B_X0 = 72.5 + TITLE_DX
S.panel_title_mm(fig, "Same true sequence A→B→C, one abnormality model", x=B_X0, y_top=ROW_D)
# event separability line, the sharp event B on a gold chip (not blue: blue means CONCORD)
yb_line = ROW_D + 4.7
x = B_X0
pieces = [("Event separability (AUC):", "normal", 1.4), ("A", "normal", 0.8), (f"{auc['A']:.2f},", "normal", 1.8),
          ("B", "bold", 1.3), (f"{auc['B']:.2f},", "normal", 1.4), ("C", "normal", 0.8), (f"{auc['C']:.2f}", "normal", 0)]
for s, wt, gap in pieces:
    wpx = text_w(s, 6, wt)
    if s == "B":
        chip(x + wpx / 2, yb_line - cap_mm(6) / 2, wpx + 0.9)
        cv.text(x, yb_line, s, fontsize=6, fontweight="bold", color=INK, va="baseline", ha="left", zorder=2)
    else:
        cv.text(x, yb_line, s, fontsize=6, fontweight=wt, color=GREY_DARK, va="baseline", ha="left", zorder=2)
    x += wpx + gap

B_TOP, B_W = ROW_D + 18.5, 34.0
B_H = AX_BOTTOM - B_TOP
SUB_Y, SUBBAR_Y = ROW_D + 10.7, ROW_D + 12.4   # composition subtitle baseline / bar top
BSUB = [(B_X0, 89.0, "CN-heavy composition"), (129.5, 142.0, "AD-heavy composition")]
XL = (0.26, 0.72)
YL = (-0.6, len(ORDERS) - 0.4)
ys = {o: len(ORDERS) - 1 - i for i, o in enumerate(ORDERS)}
b_axes = []
for x_lab, x_ax, comp in BSUB:
    d = LOSS[comp]
    # composition subtitle and the composition itself (same bar as panel c)
    cv.text(x_lab, SUB_Y, comp, fontsize=6.5, fontweight="bold", color=INK, va="baseline", ha="left")
    c = np.array(COMP[comp], dtype=float)
    comp_bar(x_ax, SUBBAR_Y, B_W, 2.6, c / c.sum(), pct=True)
    cv.text(x_ax - 1.4, SUBBAR_Y + 1.3 + cap_mm(5.5) / 2, "CN/MCI/AD", fontsize=5.5, color=GREY_MID, ha="right",
            va="baseline")

    ax = data_axes(x_ax, B_TOP, B_W, B_H)
    b_axes.append(ax)
    ax.axhspan(ys[TRUE] - 0.42, ys[TRUE] + 0.42, color="#EFEFEF", zorder=0, lw=0)
    for o in ORDERS:
        ax.hlines(ys[o], *XL, color="#E2E2E2", lw=0.5, zorder=0.6)
    ax.text(XL[1] - 0.006, ys[TRUE], "true", fontsize=5.5, color=GREY_DARK, va="center", ha="right", style="italic",
            bbox=dict(boxstyle="square,pad=0.12", fc="#EFEFEF", ec="none"))
    for o in ORDERS:
        if o == best[comp]:
            ax.plot(d[o], ys[o], "o", ms=4.8, mfc=INK, mec="white", mew=0.4, zorder=4)
            ax.text(d[o] + 0.017, ys[o], f"{d[o]:.3f}  best fit", fontsize=5.5, fontweight="bold", color=INK,
                    va="center", ha="left", bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"))
        elif o == TRUE:
            ax.plot(d[o], ys[o], "o", ms=3.6, mfc="white", mec=INK, mew=0.7, zorder=4)
            ax.text(d[o] + 0.02, ys[o], f"{d[o]:.3f}", fontsize=5.5, color=GREY_DARK, va="center", ha="left",
                    bbox=dict(boxstyle="square,pad=0.12", fc="#EFEFEF", ec="none"))
        else:
            ax.plot(d[o], ys[o], "o", ms=3.0, mfc="#A9A9A9", mec="white", mew=0.3, zorder=3)
    ax.set_xlim(*XL)
    ax.set_ylim(*YL)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7])

    # ordering labels, drawn letter by letter so that event B keeps its gold chip
    lw_, aw_, gap_ = text_w("A", 6), text_w("→", 6), 0.6
    x_right = x_ax - 1.3
    for o in ORDERS:
        yc = B_TOP + (YL[1] - ys[o]) / (YL[1] - YL[0]) * B_H
        isb, ist = o == best[comp], o == TRUE
        col = INK if (isb or ist) else "#6A6A6A"
        wt = "bold" if isb else "normal"
        seq = o.split("<")
        xx = x_right - (3 * lw_ + 2 * aw_ + 4 * gap_)
        for k, e in enumerate(seq):
            if e == "B" and (isb or ist):
                chip(xx + lw_ / 2, yc, lw_ + 0.8)
            cv.text(xx + lw_ / 2, yc + cap_mm(6) / 2, e, fontsize=6, fontweight="bold" if (e == "B" or isb) else wt,
                    color=col, ha="center", va="baseline", zorder=2)
            xx += lw_ + gap_
            if k < 2:
                cv.text(xx + aw_ / 2, yc + cap_mm(6) / 2, "→", fontsize=6, color=col, ha="center", va="baseline",
                        zorder=2)
                xx += aw_ + gap_
cv.text((BSUB[0][1] + BSUB[1][1] + B_W) / 2, XLAB_TOP, "Expected ranking loss (lower = better fit)", fontsize=7,
        color=INK, ha="center", va="top")

# ==============================================================================================
# Panel a — separately fitted DEBM (grey lane) vs CONCORD (blue lane)
# ==============================================================================================
BOXES = [(3.0, 35.5), (42.0, 39.5), (85.0, 41.0)]
SEP = S.METHOD["separate"]
CON = S.METHOD["concord"]

# grey lane: separately fitted DEBM
LANE_Y, LANE_H = ROW_M + 1.9, 6.6
cv.plot(1.0 + TITLE_DX + 0.7, ROW_M - cap_mm(6.5) / 2, SEP["marker"], ms=S.ms("separate", 3.4), mfc=SEP["color"], mec="white",
        mew=0.3)
cv.text(1.0 + TITLE_DX + 2.2, ROW_M, SEP["label"], fontsize=6.5, fontweight="bold", color=S.TEXT_COLOR["separate"],
        va="baseline", ha="left")
lane_txt = ["Abnormality model fitted per group", "Observed CN/MCI/AD proportions", "Orderings differ with composition"]
for (bx, bw), t in zip(BOXES, lane_txt):
    rbox(bx, LANE_Y, bw, LANE_H, fc=LANE_BG, ec="none", lw=0, r=1.2)
    cv.text(bx + bw / 2, LANE_Y + LANE_H / 2 + cap_mm(6) / 2, t, fontsize=6, color=GREY_DARK, ha="center",
            va="baseline", zorder=3)
for i in range(2):
    xa = BOXES[i][0] + BOXES[i][1] + 0.5
    xb = BOXES[i + 1][0] - 0.5
    arrow(xa, LANE_Y + LANE_H / 2, xb, LANE_Y + LANE_H / 2, color=GREY_LINE, lw=0.8, ms=5)

# blue lane: CONCORD
CON_BASE = ROW_M + 13.7
cv.plot(1.0 + TITLE_DX + 0.7, CON_BASE - cap_mm(6.5) / 2, CON["marker"], ms=S.ms("concord", 3.4), mfc=CON["color"], mec="white",
        mew=0.3)
cv.text(1.0 + TITLE_DX + 2.2, CON_BASE, CON["label"], fontsize=6.5, fontweight="bold", color=BLUE, va="baseline", ha="left")
T, BOT = ROW_M + 15.5, ROW_M + 68.0
CAP1, CAP2 = ROW_M + 61.9, ROW_M + 65.3     # common baselines of the bottom captions in a and b
titles = ["One abnormality model fitted\nto the pooled sample",
          "Weights give every group the\nsame CN/MCI/AD proportions",
          "Each group's ordering\nestimated from weighted scores"]
for (bx, bw), t in zip(BOXES, titles):
    rbox(bx, T, bw, BOT - T)
    cv.text(bx + 1.6, T + 1.5, t, fontsize=6.5, fontweight="bold", color=BLUE, va="top", ha="left", linespacing=1.15,
            zorder=3)
for i in range(2):
    xa = BOXES[i][0] + BOXES[i][1] + 0.5
    xb = BOXES[i + 1][0] - 0.5
    arrow(xa, T + (BOT - T) / 2, xb, T + (BOT - T) / 2, lw=1.0, ms=6)

# --- box 1: one abnormality model from the pooled sample
bx, bw = BOXES[0]
cx = bx + bw / 2
y_mk = ROW_M + 30.1
g1x = [bx + 4.0 + 1.9 * k for k in range(5)]
g2x = [bx + bw - 4.0 - 1.9 * k for k in range(5)][::-1]
for xm in g1x:
    mk(xm, y_mk, G1)
for xm in g2x:
    mk(xm, y_mk, G2)
cv.text(np.mean(g1x), y_mk - 2.0, G1["label"], fontsize=5.5, color=INK, ha="center", va="baseline")
cv.text(np.mean(g2x), y_mk - 2.0, G2["label"], fontsize=5.5, color=INK, ha="center", va="baseline")
arrow(np.mean(g1x), y_mk + 1.8, cx - 2.0, y_mk + 6.2, color=GREY_DARK, lw=0.6, ms=4)
arrow(np.mean(g2x), y_mk + 1.8, cx + 2.0, y_mk + 6.2, color=GREY_DARK, lw=0.6, ms=4)
cv.text(cx, y_mk + 8.9, "Pooled sample", fontsize=5.5, color=INK, ha="center", va="baseline")
# densities (illustrative)
x0p, x1p, ybase, hmax = bx + 3.0, bx + bw - 3.0, CAP2 - 3.4, 16.5
z = np.linspace(-3.2, 5.7, 300)
fn = np.exp(-0.5 * z ** 2)
fa = 0.8 * np.exp(-0.5 * ((z - 2.6) / 1.1) ** 2)
edges = np.linspace(-3.2, 5.7, 17)
cent = 0.5 * (edges[:-1] + edges[1:])
hist = 0.55 * np.exp(-0.5 * cent ** 2) + 0.45 * 0.8 * np.exp(-0.5 * ((cent - 2.6) / 1.1) ** 2)


def zx(v):
    return x0p + (v - z[0]) / (z[-1] - z[0]) * (x1p - x0p)


binw = (x1p - x0p) / (len(edges) - 1)
for c_, h_ in zip(cent, hist):
    cv.add_patch(Rectangle((zx(c_) - binw / 2 + 0.12, ybase - h_ * hmax), binw - 0.24, h_ * hmax, fc="#D0DBEA",
                           ec="none", zorder=2))
cv.plot(zx(z), ybase - fn * hmax, color=INK, lw=0.8, zorder=3)
cv.plot(zx(z), ybase - fa * hmax, color=INK, lw=0.8, ls=(0, (2.2, 1.2)), zorder=3)
cv.plot([x0p, x1p], [ybase, ybase], color=INK, lw=0.6, zorder=3)
cv.text(zx(-0.9), ybase - 1.03 * hmax, "Normal", fontsize=5.5, color=INK, ha="right", va="baseline")
cv.text(zx(3.4), ybase - 0.84 * hmax, "Abnormal", fontsize=5.5, color=INK, ha="left", va="baseline")
cv.text(cx, CAP2, "Biomarker value", fontsize=5.5, color=INK, ha="center", va="baseline")

# --- box 2: weights to a common composition (illustrative 6/3/1 vs 1/3/6, common 3.5/3/3.5)
bx, bw = BOXES[1]
bar_o, bar_c, barh = 14.0, 15.5, 3.2
xo = bx + 4.4
xc_ = bx + bw - 1.6 - bar_c
cv.text(xo + bar_o / 2, ROW_M + 28.5, "Observed", fontsize=5.5, color=INK, ha="center", va="baseline")
cv.text(xc_ + bar_c / 2, ROW_M + 28.5, "Common", fontsize=5.5, color=INK, ha="center", va="baseline")
for g, P, Wt, yrow in ((G1, P1, W1, ROW_M + 30.8), (G2, P2, W2, ROW_M + 45.1)):
    mk(bx + 2.3, yrow + barh / 2, g, ms=3.6)
    comp_bar(xo, yrow, bar_o, barh, P)
    arrow(xo + bar_o + 0.7, yrow + barh / 2, xc_ - 0.7, yrow + barh / 2, color=BLUE, lw=0.8, ms=4.5)
    comp_bar(xc_, yrow, bar_c, barh, PC, labels=True)
    # ten participants under the observed bar; marker area proportional to the weight
    k = 0
    for di, nd in enumerate(np.rint(P * 10).astype(int)):
        for _ in range(nd):
            xm = xo + (k + 0.5) * bar_o / 10
            mk(xm, yrow + barh + 2.4, g, ms=2.35 * np.sqrt(Wt[di]), mew=0.25)
            k += 1
        # the weight itself, under the participants of that diagnosis; the MCI "1" is nudged 0.5 mm
        # away from the neighbouring "3.5"
        xw = xo + (k - nd / 2) * bar_o / 10
        if DIAGS[di] == "MCI":
            xw += -0.5 if Wt[2] > Wt[0] else 0.5
        cv.text(xw, yrow + barh + 5.9, f"{Wt[di]:.2f}".rstrip("0").rstrip("."), fontsize=5, color=GREY_DARK,
                ha="center", va="baseline")
cv.text(bx + bw / 2, CAP1, "weight = common share ÷ group share", fontsize=6, color=INK, ha="center",
        va="baseline")
cv.text(bx + bw / 2, CAP2, "per diagnosis; marker area shows weight", fontsize=5.5, color=GREY_MID, ha="center",
        va="baseline")

# --- box 3: orderings from weighted scores (neutral event boxes)
bx, bw = BOXES[2]
EV = ["A", "B", "C", "D", "E"]
ew, eh, egap = 4.4, 4.4, 1.25
row_w = 5 * ew + 4 * egap
ex0 = bx + (bw - row_w) / 2 + 1.6
seqs = [(G1, ["A", "B", "C", "D", "E"], ROW_M + 30.5), (G2, ["A", "B", "D", "C", "E"], ROW_M + 44.7)]
pos = {}
for g, seq, yrow in seqs:
    mk(ex0 - 3.0, yrow + eh / 2, g, ms=3.6)
    for i, e in enumerate(seq):
        xe = ex0 + i * (ew + egap)
        rbox(xe, yrow, ew, eh, fc="white", ec=INK, lw=0.5, r=0.7, z=3)
        cv.text(xe + ew / 2, yrow + eh / 2 + cap_mm(6) / 2, e, fontsize=6, fontweight="bold", color=INK,
                ha="center", va="baseline", zorder=4)
        pos[(g["label"], e)] = xe + ew / 2
for e in EV:
    xa_, xb_ = pos[("Group 1", e)], pos[("Group 2", e)]
    crossing = xa_ != xb_
    cv.plot([xa_, xb_], [seqs[0][2] + eh + 0.3, seqs[1][2] - 0.3], color=BLUE if crossing else "#B5B5B5",
            lw=1.0 if crossing else 0.5, zorder=2, solid_capstyle="round")
# bracket under the swapped pair
xl_, xr_ = ex0 + 2 * (ew + egap), ex0 + 3 * (ew + egap) + ew
yb_ = seqs[1][2] + eh + 1.0
cv.plot([xl_, xl_, xr_, xr_], [yb_ - 0.6, yb_, yb_, yb_ - 0.6], color=BLUE, lw=0.6, zorder=3)
cv.text((xl_ + xr_) / 2, yb_ + 1.0 + cap_mm(5.5), "1 discordant pair", fontsize=5.5, fontweight="bold", color=BLUE,
        ha="center", va="baseline")
cv.text(bx + bw / 2, CAP1, "Distance between group orderings", fontsize=5.5, color=INK, ha="center", va="baseline")
cv.text(bx + bw / 2, CAP2, f"({S.LABEL['distance'].lower()})", fontsize=5.5, color=GREY_MID, ha="center", va="baseline")

# ==============================================================================================
# Panel b — within-diagnosis permutation
# ==============================================================================================
DX0, DX1 = 129.2, 178.3
rbox(DX0, LANE_Y, DX1 - DX0, BOT - LANE_Y)
cv.text(DX0 + 1.6, LANE_Y + 1.5, "Group labels permuted\nwithin diagnosis", fontsize=6.5, fontweight="bold",
        color=BLUE, va="top", ha="left", linespacing=1.15, zorder=3)
y_leg = ROW_M + 13.9
for k, g in enumerate((G1, G2)):
    xx = DX0 + 5.0 + k * 14.0
    mk(xx, y_leg - cap_mm(5.5) / 2, g, ms=3.4)
    cv.text(xx + 1.6, y_leg, g["label"], fontsize=5.5, color=INK, ha="left", va="baseline")
# participants within each diagnosis (1 = Group 1, 2 = Group 2); one within-row swap each
rows_d = [("CN", [1, 1, 2, 1, 1, 2], (1, 2), ROW_M + 21.5),
          ("MCI", [1, 2, 1, 1, 2, 2], (3, 4), ROW_M + 29.1),
          ("AD", [2, 2, 1, 2, 2, 1], (2, 3), ROW_M + 36.7)]
xs = [DX0 + 12.5 + 5.3 * k for k in range(6)]
for dname, members, (i0, i1), yr in rows_d:
    cv.add_patch(Rectangle((xs[0] - 2.6, yr - 2.2), xs[-1] - xs[0] + 5.2, 4.4, fc="white", ec="none", zorder=2,
                           alpha=0.9))
    cv.text(xs[0] - 3.6, yr + cap_mm(6) / 2, dname, fontsize=6, color=INK, ha="right", va="baseline")
    for xm, gm in zip(xs, members):
        mk(xm, yr, S.GROUP[gm], ms=3.6)
    arrow(xs[i0], yr - 1.4, xs[i1], yr - 1.4, color=BLUE, lw=0.6, ms=3.5, style="<|-|>", rad=-0.9, z=6)
# a swap across diagnoses is not allowed: greyed, crossed-out arc between the CN and AD rows
xa_d = xs[-1] + 1.4
arc_rad = -0.32
arrow(xa_d, rows_d[0][3], xa_d, rows_d[2][3], color="#AFAFAF", lw=0.6, ms=3.5, style="<|-|>", rad=arc_rad, z=6)
xm_arc = xa_d + 0.5 * abs(arc_rad) * (rows_d[2][3] - rows_d[0][3])
ym_arc = (rows_d[0][3] + rows_d[2][3]) / 2
for sgn in (-1, 1):
    cv.plot([xm_arc - 0.8, xm_arc + 0.8], [ym_arc - 0.8 * sgn, ym_arc + 0.8 * sgn], color=GREY_DARK, lw=0.8,
            zorder=7, solid_capstyle="round")
# permutation distribution of the distance (illustrative)
hx0, hx1, hbase, hh = DX0 + 4.5, DX1 - 3.0, CAP1 - 3.4, 12.0
nb = 16
cent = np.linspace(0.5, nb - 0.5, nb)
shape = cent ** 2.2 * np.exp(-cent / 2.2)
shape = shape / shape.max()
obs_bin = 12.0
bwid = (hx1 - hx0) / nb
for c_, h_ in zip(cent, shape):
    tail = c_ > obs_bin
    cv.add_patch(Rectangle((hx0 + (c_ - 0.5) * bwid + 0.12, hbase - h_ * hh), bwid - 0.24, h_ * hh,
                           fc=S.PAL["blue_light"] if tail else S.PAL["neutral"], ec="none", zorder=2))
cv.plot([hx0, hx1], [hbase, hbase], color=INK, lw=0.6, zorder=3)
xobs = hx0 + obs_bin * bwid
# tail of the permutation distribution (no formula: E/B would clash with the events of a and d;
# the caption defines the P value)
cv.text(xobs + 0.8, hbase - 0.32 * hh, "≥ observed", fontsize=5.5, color=S.PAL["blue_secondary"], ha="left",
        va="baseline")
cv.plot([xobs, xobs], [hbase, hbase - hh - 1.0], color=BLUE, lw=1.1, zorder=4)
cv.text(xobs, hbase - hh - 1.6, "Observed", fontsize=5.5, color=BLUE, ha="center", va="baseline",
        fontweight="bold")
cv.text((hx0 + hx1) / 2, CAP1, "Distance between group orderings", fontsize=5.5, color=INK, ha="center",
        va="baseline")
cv.text((hx0 + hx1) / 2, CAP2, f"({S.LABEL['distance'].lower()})", fontsize=5.5, color=GREY_MID, ha="center",
        va="baseline")

# ----------------------------------------------------------------------------------------------
S.panel_label_mm(fig, "a", x=1.0, y_top=ROW_M)
S.panel_label_mm(fig, "b", x=127.4, y_top=ROW_M)
S.panel_label_mm(fig, "c", x=1.0, y_top=ROW_D)
S.panel_label_mm(fig, "d", x=72.5, y_top=ROW_D)

paths = S.save(fig, "Fig1_overview")
print("\n".join(str(p) for p in paths))
