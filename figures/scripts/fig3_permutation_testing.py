"""Figure 3 | Permutation testing of group differences in simulation (180 x 141 mm).

a  False-positive rate when groups differ only in composition (ADNI-like composition, one
   common ordering): unrestricted (hollow) versus within-diagnosis (filled) permutation, per
   ordering estimator, on a log10 axis so the fold reduction reads as arrow length.
b  Paired benchmark: false-positive rate with within-diagnosis permutation under three
   measurement models.
c  Detection of a true ordering change in one group, four scenarios (rows); horizontal grouped
   bars; brackets give the gain from weighting (CONCORD - pooled-score DEBM).
d  Paired gain of CONCORD over separately fitted DEBM with paired-bootstrap 95% intervals,
   sharing the scenario rows of c. The rows of c and d are grouped by the changed group
   (AD-heavy, 485 participants; CN-heavy, 75 participants).

Inputs (figures/inputs/simulation/):
  fig2_calibration.csv     a: rejections of the one-common-ordering simulation (1,000 datasets)
  rates.csv                b, c: paired benchmark rates (Supplementary Table 2)
  paired_differences.csv   c brackets, d: paired differences with paired-bootstrap intervals
Model labels in rates.csv / paired_differences.csv: repaired = separately fitted DEBM,
shared = pooled-score DEBM, invariant_min = CONCORD.

All plotted numbers are read from stored aggregate files; nothing is refitted.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path as _Path

sys.path.insert(0, str(_Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402

S.apply_style()

CAL = S.SIM_INPUTS / "fig2_calibration.csv"
RATES = S.SIM_INPUTS / "rates.csv"
PAIRED = S.SIM_INPUTS / "paired_differences.csv"
S.require_inputs([CAL, RATES, PAIRED], "Fig. 3")

W, H = 180.0, 141.0
INK = S.PAL["ink"]
GREY = S.PAL["grey_mid"]
METHODS = ["separate", "pooled", "concord"]  # drawing order: CONCORD last (on top)
ARM = {"separate": "repaired", "pooled": "shared", "concord": "invariant_min"}
CAL_NAME = {"separate": "Standard", "pooled": "Pooled only", "concord": "CONCORD"}  # method keys of CAL
TOL_RATE, TOL_GAP = 0.002, 0.01


def C(m):
    return S.METHOD[m]["color"]


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(name, header, rows):
    S.SRC.mkdir(parents=True, exist_ok=True)
    with open(S.SRC / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def close(a, b, tol, what):
    assert abs(a - b) <= tol, f"{what}: got {a}, expected {b}"


# ----------------------------------------------------------------------------- data
cal = {(r["method"], r["operator"]): r for r in read_csv(CAL)}
rates = {(r["cell"], r["arm"], r["endpoint"]): r for r in read_csv(RATES)}
paired = {(r["cell"], r["left_arm"], r["right_arm"], r["endpoint"]): r for r in read_csv(PAIRED)}

# a: calibration under one common ordering
A = {}
for m in METHODS:
    A[m] = {}
    for op, key in (("U-all", "unrestricted"), ("D-all", "within")):
        r = cal[(CAL_NAME[m], op)]
        k, n = int(r["reject"]), int(r["R"])
        lo, hi = S.wilson(k, n)
        close(lo, float(r["low"]), 1e-6, f"a Wilson low {m} {op}")
        close(hi, float(r["high"]), 1e-6, f"a Wilson high {m} {op}")
        A[m][key] = dict(k=k, n=n, rate=k / n, lo=lo, hi=hi)
exp_a = {"separate": (0.132, 0.030), "pooled": (0.777, 0.028), "concord": (0.108, 0.042)}
for m, (u, d) in exp_a.items():
    close(A[m]["unrestricted"]["rate"], u, TOL_RATE, f"a unrestricted {m}")
    close(A[m]["within"]["rate"], d, TOL_RATE, f"a within {m}")

# b: paired benchmark false positives (within-diagnosis permutation).
# T5 = Student t with 5 degrees of freedom.
B_CELLS = [("GAUSSIAN_H0", "Independent\nGaussian noise"), ("CORR_H0", "Correlated noise\nwithin modality"),
           ("T5_H0", "Heavy-tailed noise\n" + r"($\mathit{t}$, 5 d.f.)")]
B = {}
for cell, _ in B_CELLS:
    for m in METHODS:
        r = rates[(cell, ARM[m], "any_pair")]
        k, n = int(r["rejections"]), int(r["planned_R"])
        lo, hi = S.wilson(k, n)
        close(lo, float(r["wilson_low"]), 1e-6, f"b Wilson low {cell} {m}")
        close(hi, float(r["wilson_high"]), 1e-6, f"b Wilson high {cell} {m}")
        B[(cell, m)] = dict(k=k, n=n, rate=k / n, lo=lo, hi=hi)
exp_b = {"GAUSSIAN_H0": (0.034, 0.028, 0.042), "CORR_H0": (0.037, 0.033, 0.026), "T5_H0": (0.047, 0.019, 0.025)}
for cell, vals in exp_b.items():
    for m, v in zip(METHODS, vals):
        close(B[(cell, m)]["rate"], v, TOL_RATE, f"b {cell} {m}")

# c/d: detection of a truly changed contrast
SCEN = [  # rows of c and d, grouped by the changed group (headers in GROUP_HEAD)
    ("EXACT27_G3", "27 of 91 pairs\nreversed"),
    ("SINGLE13_G3", "One event moved\n13 positions"),
    ("CORR_SINGLE13_G3", "One event moved\n13 positions,\ncorrelated noise"),
    ("EXACT27_G1", "27 of 91 pairs\nreversed"),
]
SCEN_GROUP = {"EXACT27_G3": "AD-heavy", "SINGLE13_G3": "AD-heavy", "CORR_SINGLE13_G3": "AD-heavy",
              "EXACT27_G1": "CN-heavy"}
BRACKET_CELLS = ["EXACT27_G3", "SINGLE13_G3", "CORR_SINGLE13_G3"]
Cd = {}
for cell, _ in SCEN:
    for m in METHODS:
        r = rates[(cell, ARM[m], "truly_affected_pair")]
        k, n = int(r["rejections"]), int(r["planned_R"])
        assert n == 500, f"c {cell} {m}: R = {n}"
        lo, hi = S.wilson(k, n)
        close(lo, float(r["wilson_low"]), 1e-6, f"c Wilson low {cell} {m}")
        Cd[(cell, m)] = dict(k=k, n=n, rate=k / n, lo=lo, hi=hi)
exp_c = {"EXACT27_G3": (0.136, 0.350, 0.542), "SINGLE13_G3": (0.082, 0.444, 0.542),
         "CORR_SINGLE13_G3": (0.088, 0.376, 0.510), "EXACT27_G1": (0.082, 0.184, 0.176)}
for cell, vals in exp_c.items():
    for m, v in zip(METHODS, vals):
        close(Cd[(cell, m)]["rate"], v, TOL_RATE, f"c {cell} {m}")

# weighting gain (CONCORD - pooled-score DEBM), paired file, orientation checked against rates
WG = {}
for cell in BRACKET_CELLS:
    r = paired[(cell, "invariant_min", "shared", "truly_affected_pair")]
    d = float(r["difference"])
    close(d, Cd[(cell, "concord")]["rate"] - Cd[(cell, "pooled")]["rate"], 1e-9, f"c orientation {cell}")
    WG[cell] = dict(diff=d, lo=float(r["paired_bootstrap_low"]), hi=float(r["paired_bootstrap_high"]))
for cell, v in zip(BRACKET_CELLS, (0.192, 0.098, 0.134)):
    close(WG[cell]["diff"], v, TOL_GAP, f"c bracket {cell}")

# d: paired gain CONCORD - separately fitted DEBM
D = {}
for cell, _ in SCEN:
    r = paired[(cell, "invariant_min", "repaired", "truly_affected_pair")]
    d = float(r["difference"])
    close(d, Cd[(cell, "concord")]["rate"] - Cd[(cell, "separate")]["rate"], 1e-9, f"d orientation {cell}")
    D[cell] = dict(diff=d, lo=float(r["paired_bootstrap_low"]), hi=float(r["paired_bootstrap_high"]),
                   R=int(r["paired_R"]), B=int(r["bootstrap_replicates"]))
exp_d = {"EXACT27_G3": (0.406, 0.354, 0.458), "SINGLE13_G3": (0.460, 0.412, 0.508),
         "CORR_SINGLE13_G3": (0.422, 0.376, 0.470), "EXACT27_G1": (0.094, 0.054, 0.134)}
for cell, (d, lo, hi) in exp_d.items():
    close(D[cell]["diff"], d, TOL_GAP, f"d {cell}")
    close(D[cell]["lo"], lo, TOL_GAP, f"d low {cell}")
    close(D[cell]["hi"], hi, TOL_GAP, f"d high {cell}")

# ----------------------------------------------------------------------------- source data
write_csv("Fig3a.csv", ["model", "permutation", "rejections", "datasets", "false_positive_rate",
                        "wilson_low", "wilson_high", "fold_reduction_unrestricted_over_within", "source"],
          [[S.METHOD[m]["label"], {"unrestricted": "unrestricted", "within": "within diagnosis"}[op],
            A[m][op]["k"], A[m][op]["n"], f"{A[m][op]['rate']:.4f}", f"{A[m][op]['lo']:.4f}",
            f"{A[m][op]['hi']:.4f}", f"{A[m]['unrestricted']['rate'] / A[m]['within']['rate']:.2f}",
            S.source_name(CAL)]
           for m in METHODS for op in ("unrestricted", "within")])
write_csv("Fig3b.csv", ["measurement_model", "model", "rejections", "datasets", "false_positive_rate",
                        "wilson_low", "wilson_high", "cell", "source"],
          [[lab.replace(r"$\mathit{t}$", "t").replace("\n", " "), S.METHOD[m]["label"], B[(c, m)]["k"], B[(c, m)]["n"],
            f"{B[(c, m)]['rate']:.4f}", f"{B[(c, m)]['lo']:.4f}", f"{B[(c, m)]['hi']:.4f}", c,
            S.source_name(RATES)] for c, lab in B_CELLS for m in METHODS])
rows_c = []
for c, lab in SCEN:
    for m in METHODS:
        v = Cd[(c, m)]
        rows_c.append([f"{SCEN_GROUP[c]} group: " + lab.replace("\n", " "), S.METHOD[m]["label"], v["k"], v["n"], f"{v['rate']:.4f}",
                       f"{v['lo']:.4f}", f"{v['hi']:.4f}", "", "", "", c, S.source_name(RATES)])
    if c in WG:
        g = WG[c]
        rows_c.append([f"{SCEN_GROUP[c]} group: " + lab.replace("\n", " "), "Bracket: CONCORD minus Pooled-score DEBM", "", "", "", "", "",
                       f"{g['diff']:.4f}", f"{g['lo']:.4f}", f"{g['hi']:.4f}", c, S.source_name(PAIRED)])
write_csv("Fig3c.csv", ["scenario", "model", "rejections", "datasets", "detection_rate", "wilson_low",
                        "wilson_high", "weighting_gain", "paired_bootstrap_low", "paired_bootstrap_high",
                        "cell", "source"], rows_c)
write_csv("Fig3d.csv", ["scenario", "contrast", "paired_difference", "paired_bootstrap_low",
                        "paired_bootstrap_high", "paired_datasets", "bootstrap_replicates", "cell", "source"],
          [[f"{SCEN_GROUP[c]} group: " + lab.replace("\n", " "), "CONCORD minus Separately fitted DEBM", f"{D[c]['diff']:.4f}",
            f"{D[c]['lo']:.4f}", f"{D[c]['hi']:.4f}", D[c]["R"], D[c]["B"], c, S.source_name(PAIRED)]
           for c, lab in SCEN])

# ----------------------------------------------------------------------------- figure
fig = plt.figure(figsize=(W * S.MM, H * S.MM))

LEG_TOP = 1.0                  # method legend row (top-centre)
ROW1, ROW2 = 10.5, 75.2        # letter / title baselines (mm from top)
A_TOP, A_H = 14.0, 44.0        # top row axes (a, b)
C_TOP, C_H = 78.7, 50.8        # bottom row axes (c, d)
AX_A = (25.5, A_TOP, 50.5, A_H)
axA = S.add_axes_mm(fig, *AX_A)
axB = S.add_axes_mm(fig, 92.0, A_TOP, 86.0, A_H)
axC = S.add_axes_mm(fig, 28.5, C_TOP, 80.0, C_H)
axD = S.add_axes_mm(fig, 121.5, C_TOP, 56.5, C_H)
for _ax in (axC, axD):
    _ax.patch.set_visible(False)  # let the shared scenario separators show through
LETTER_X = {"a": 1.0, "b": 82.5, "c": 1.0, "d": 114.0}
TITLE_DX = 4.5                 # title starts ~3 mm right of the letter (as in Figs 1, 2, 4 and 5)


def text_w(ax, s, fs=6, weight="normal", log=False):
    """Width of a string: in x data units (linear axis) or as a multiplicative factor (log axis)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    t = ax.text(0.5, 0.5, s, fontsize=fs, fontweight=weight, transform=ax.transAxes)
    wpx = t.get_window_extent(renderer=r).width
    t.remove()
    x0, x1 = ax.get_xlim()
    axw = ax.get_window_extent(renderer=r).width
    if log:
        return 10 ** (wpx / axw * (math.log10(x1) - math.log10(x0)))
    return wpx / axw * (x1 - x0)


def interval(ax, m, x=None, y=None, xerr=None, yerr=None, zorder=3):
    lw = 0.9 if m == "concord" else S.ERRBAR["elinewidth"]
    ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt="none", ecolor=C(m), capsize=S.ERRBAR["capsize"],
                elinewidth=lw, zorder=zorder)


def ms(m, base=4.1):
    return S.ms(m, base + (0.5 if m == "concord" else 0.0))


# ---- a: groups differ only in composition (log10 axis) -----------------------
ax = axA
YA = {"separate": 2.0, "pooled": 1.0, "concord": 0.0}
XL = (0.82, 150.0)
ax.set_xscale("log")
ax.set_xlim(*XL)
YLIM_A = (-0.56, 3.02)
ax.set_ylim(*YLIM_A)
ax.spines["bottom"].set_bounds(1, 100)
TICKS = [1, 2, 5, 10, 20, 50, 100]
ax.xaxis.set_major_locator(FixedLocator(TICKS))
ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in TICKS]))
ax.xaxis.set_minor_locator(NullLocator())
ax.set_yticks([YA[m] for m in METHODS])
ax.set_yticklabels([S.METHOD[m]["label"] for m in METHODS])
for tl, m in zip(ax.get_yticklabels(), METHODS):
    tl.set_color(INK)
    if m == "concord":
        tl.set_fontweight("bold")
        tl.set_color(C("concord"))
ax.tick_params(axis="y", length=0, pad=2.5)
ax.spines["left"].set_visible(False)
ax.set_xlabel("False-positive rate (%, log scale)")
Y5_TOP = 2.19  # 5% line stops under its "Nominal 5%" label (as in b)
ax.plot([5, 5], [YLIM_A[0], Y5_TOP], **S.REF_LINE)
ax.text(5, Y5_TOP + 0.03, "Nominal 5%", color=GREY, fontsize=6, ha="center", va="bottom")
DY_ARR = 0.25  # elbow arrow below the row: down from hollow (unrestricted), up into filled (within diagnosis)
for m in METHODS:
    y, col = YA[m], C(m)
    u, d = A[m]["unrestricted"], A[m]["within"]
    for v, hollow in ((u, True), (d, False)):
        interval(ax, m, x=100 * v["rate"], y=y,
                 xerr=[[100 * (v["rate"] - v["lo"])], [100 * (v["hi"] - v["rate"])]])
        ax.plot(100 * v["rate"], y, marker=S.METHOD[m]["marker"], ms=ms(m), mew=0.85 if hollow else 0.5,
                mfc="white" if hollow else col, mec=col, ls="none", zorder=4)
    xu, xd = 100 * u["rate"], 100 * d["rate"]
    path = Path([(xu, y - 0.11), (xu, y - DY_ARR), (xd, y - DY_ARR), (xd, y - 0.105)])
    ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=4.6, lw=0.6, color=col,
                                 shrinkA=0, shrinkB=0, zorder=2, joinstyle="miter"))
# value labels beside the markers: within diagnosis to the left, unrestricted to the right
fig.canvas.draw()
# CONCORD call-out band covers the CONCORD row label and the plot row (as Fig2a)
_r = fig.canvas.get_renderer()
_cl = [t for t in ax.get_yticklabels() if t.get_text() == "CONCORD"][0]
_axbb = ax.get_window_extent(_r)
_x0 = (_cl.get_window_extent(_r).x0 - 1.3 * fig.dpi / 25.4 - _axbb.x0) / _axbb.width
S.concord_band(ax, YA["concord"], h=0.74, xmin=_x0, xmax=1.0, clip_on=False)
for m in METHODS:
    y = YA[m]
    u, d = A[m]["unrestricted"], A[m]["within"]
    ax.text(100 * d["lo"] / 1.06, y, S.fmt_pct(100 * d["rate"]), color=S.TEXT_COLOR[m], fontsize=6,
            fontweight="bold", ha="right", va="center", zorder=5)
    ax.text(max(100 * u["hi"], 100 * u["rate"] * 1.1) * 1.05, y, S.fmt_pct(100 * u["rate"]), color=S.TEXT_COLOR[m], fontsize=6,
            ha="left", va="center", zorder=5)
# hollow / filled key (upper strip, above the first row)
kx = 1.12
for yk, hollow, lab in ((2.89, True, S.LABEL["unrestricted"]), (2.61, False, S.LABEL["within"])):
    ax.plot(kx, yk, marker="o", ms=3.4, mfc="white" if hollow else GREY, mec=GREY, mew=0.8 if hollow else 0.4,
            ls="none", clip_on=False)
    ax.text(kx * 1.16, yk, lab, fontsize=6, va="center", ha="left", color=INK)

# ---- b: within-diagnosis permutation across measurement models ---------------
ax = axB
OFF = {"separate": -0.24, "pooled": 0.0, "concord": 0.24}
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(0, 6.6)
ax.set_yticks(range(0, 7))
ax.set_xticks(range(3))
ax.set_xticklabels([lab for _, lab in B_CELLS])
ax.tick_params(axis="x", length=0, pad=3)
ax.set_ylabel("False-positive rate (%)")
ax.axhline(5, **S.REF_LINE)
ax.text(2.5, 5.1, "Nominal 5%", color=GREY, fontsize=6, ha="right", va="bottom")
for i, (cell, _) in enumerate(B_CELLS):
    for m in METHODS:
        v = B[(cell, m)]
        x = i + OFF[m]
        interval(ax, m, x=x, y=100 * v["rate"],
                 yerr=[[100 * (v["rate"] - v["lo"])], [100 * (v["hi"] - v["rate"])]])
        ax.plot(x, 100 * v["rate"], marker=S.METHOD[m]["marker"], ms=ms(m), color=C(m), mec=C(m),
                ls="none", zorder=4)
        ax.text(x + 0.05, 100 * v["rate"], S.fmt_pct(100 * v["rate"]), color=S.TEXT_COLOR[m], fontsize=6,
                ha="left", va="center", fontweight="bold" if m == "concord" else "normal", zorder=5)

# ---- c: true change in one group (horizontal grouped bars) -------------------
GY = {"EXACT27_G3": 4.2, "SINGLE13_G3": 3.2, "CORR_SINGLE13_G3": 2.2, "EXACT27_G1": 0.8}  # rows, shared with d
GROUP_HEAD = [(4.84, "Change in the AD-heavy group (485 participants)"),
              (1.40, "Change in the CN-heavy group (75 participants)")]
SEP_Y = [(3.7, "#E6E6E6"), (2.7, "#E6E6E6"), (1.64, "#BDBDBD")]  # row separators; darker between groups
BOFF = {"separate": 0.28, "pooled": 0.0, "concord": -0.28}      # bars within a scenario, CONCORD lowest
BH = 0.25
YLIM = (0.3, 5.0)
SEP = "#E1E1E1"
BRK = S.PAL["grey_dark"]  # weighting-gain brackets

ax = axC
ax.set_xlim(0, 80)
ax.set_ylim(*YLIM)
ax.spines["bottom"].set_bounds(0, 60)
ax.set_xticks([0, 20, 40, 60])
ax.set_yticks([GY[c] for c, _ in SCEN])
ax.set_yticklabels([lab for _, lab in SCEN])
ax.tick_params(axis="y", length=0, pad=3)
ax.set_xlabel("Detection rate (%)")
# the y axis line is drawn per changed group, so it does not run through the group headers
_sp = ax.spines["left"]
_sp_lw, _sp_col = _sp.get_linewidth(), _sp.get_edgecolor()
_sp.set_visible(False)
for y0, y1 in ((YLIM[0], 1.24), (1.76, 4.64)):
    ax.plot([0, 0], [y0, y1], color=_sp_col, lw=_sp_lw, solid_capstyle="butt", clip_on=False, zorder=2.5)
# scenario separators run continuously through c and d (shared rows)


def y_mm(yd):
    return C_TOP + C_H * (YLIM[1] - yd) / (YLIM[1] - YLIM[0])


for yb, colr in SEP_Y:
    fig.add_artist(Line2D([1.0 / W, 178.0 / W], [1 - y_mm(yb) / H] * 2, transform=fig.transFigure,
                          color=colr, lw=0.5, zorder=0))
# changed-group headers, left-aligned with the panel letter
for yh, txt in GROUP_HEAD:
    fig.text(1.6 / W, 1 - y_mm(yh) / H, txt, fontsize=6, fontweight="bold", color=INK, ha="left",
             va="center")
fig.canvas.draw()
lab_end = {}
for cell, _ in SCEN:
    for m in METHODS:
        v = Cd[(cell, m)]
        y = GY[cell] + BOFF[m]
        ax.barh(y, 100 * v["rate"], height=BH, color=S.FILL[m], edgecolor=C(m), linewidth=0.6, zorder=2)
        if m == "concord":  # solid bar: in-bar half of the interval in white, outside half in ink
            ax.plot([100 * v["lo"], 100 * v["rate"]], [y, y], color="white", lw=S.ERRBAR["elinewidth"],
                    solid_capstyle="butt", zorder=3)
            ax.plot([100 * v["rate"], 100 * v["hi"]], [y, y], color=INK, lw=S.ERRBAR["elinewidth"],
                    solid_capstyle="butt", zorder=3)
        else:
            ax.errorbar(100 * v["rate"], y, xerr=[[100 * (v["rate"] - v["lo"])], [100 * (v["hi"] - v["rate"])]],
                        fmt="none", ecolor=INK, zorder=3, **S.ERRBAR)
        s = S.fmt_pct(100 * v["rate"])
        bold = m == "concord"
        x0 = 100 * v["hi"] + 1.0
        ax.text(x0, y, s, color=S.TEXT_COLOR[m], fontsize=6, fontweight="bold" if bold else "normal",
                ha="left", va="center", zorder=5)
        lab_end[(cell, m)] = x0 + text_w(ax, s, 6, "bold" if bold else "normal")
# weighting-gain brackets: lines run from the pooled-score and CONCORD value labels to one common x
XB = max(lab_end[(c, m)] for c in BRACKET_CELLS for m in ("pooled", "concord")) + 2.2
for cell in BRACKET_CELLS:
    yp, yc = GY[cell] + BOFF["pooled"], GY[cell] + BOFF["concord"]
    xp, xc = lab_end[(cell, "pooled")] + 0.9, lab_end[(cell, "concord")] + 0.9
    ax.plot([xp, XB, XB, xc], [yp, yp, yc, yc], color=BRK, lw=0.5, zorder=4, solid_joinstyle="miter")
    ax.text(XB + 1.0, (yp + yc) / 2, f"+{100 * WG[cell]['diff']:.1f}", color=INK, fontsize=6,
            ha="left", va="center", zorder=5)
# bracket key on c's title baseline, right-aligned to the end of c's axes
KEY_R = 28.5 + 80.0  # mm, right edge of axC
_kt = fig.text(KEY_R / W, 1 - ROW2 / H, "Bracket: gain from weighting (CONCORD − pooled-score DEBM)", fontsize=6,
               color=INK, ha="right", va="baseline")
_dpi0 = fig.dpi
fig.set_dpi(600)  # measure at the export resolution (hinting at 150 dpi shortens text by ~2 mm)
fig.canvas.draw()
_kx1 = _kt.get_window_extent(fig.canvas.get_renderer()).x0 / fig.bbox.width * W - 1.2  # mm, end of the bracket glyph
fig.set_dpi(_dpi0)
_kyb = ROW2 - 0.76                                         # mm from top: text cap-height centre
_gx = [_kx1 - 2.2, _kx1, _kx1, _kx1 - 2.2]
_gy = [_kyb - 1.1, _kyb - 1.1, _kyb + 1.1, _kyb + 1.1]
fig.add_artist(Line2D([x / W for x in _gx], [1 - y / H for y in _gy], transform=fig.transFigure, color=BRK,
                      lw=0.5, solid_joinstyle="miter"))

# ---- d: CONCORD − separately fitted DEBM (paired gain), same rows as c --------
ax = axD
DXL = (-4.0, 56.0)
ax.set_xlim(*DXL)
ax.set_ylim(*YLIM)
ax.spines["left"].set_visible(False)
ax.set_yticks([])
ax.spines["bottom"].set_bounds(0, 50)
ax.set_xticks([0, 10, 20, 30, 40, 50])
ax.set_xlabel("Difference in detection rate\n(percentage points)", linespacing=1.1)
ax.axvline(0, **S.REF_LINE)
fig.canvas.draw()
col = C("concord")
for cell, _ in SCEN:
    y, v = GY[cell], D[cell]
    ax.plot([100 * v["lo"], 100 * v["hi"]], [y, y], color=col, lw=1.1, solid_capstyle="butt", zorder=3)
    ax.plot(100 * v["diff"], y, marker="o", ms=ms("concord"), color=col, ls="none", zorder=4)
    s = f"+{100 * v['diff']:.1f} [{100 * v['lo']:.1f}, {100 * v['hi']:.1f}]"
    hw = text_w(ax, s, 6) / 2
    xc = min(max(100 * v["diff"], DXL[0] + hw + 0.5), DXL[1] - hw)
    ax.text(xc, y + 0.17, s, color=S.TEXT_COLOR["concord"], fontsize=6, ha="center", va="bottom", zorder=5)

# ---- legend, letters, titles --------------------------------------------------
S.method_legend(fig, METHODS, y_top=LEG_TOP, bars=True)
TITLES = {"a": ("Groups differ only in composition", ROW1),
          "b": ("Within-diagnosis permutation, three noise models", ROW1),
          "c": ("True change in one group", ROW2), "d": ("Gain over separately fitted DEBM", ROW2)}
for L, (t, yb) in TITLES.items():
    S.panel_label_mm(fig, L, x=LETTER_X[L], y_top=yb)
    S.panel_title_mm(fig, t, x=LETTER_X[L] + TITLE_DX, y_top=yb)

paths = S.save(fig, "Fig3_permutation_testing")
print("\n".join(str(p) for p in paths))
