"""Fig. 5 | APOE orderings estimated with CONCORD in ADNI and NACC (180 x 184 mm).

a,b  Bump charts of the observed CONCORD orderings (CSFcog6 primary panel) for
     epsilon2, epsilon3/epsilon3 and epsilon4 carriers in ADNI (a) and NACC (b).
c,d  Event-pair precedence across 200 bootstrap refits: upper-triangular heatmaps
     of P(row event before column event), one diverging colour scale centred at 0.5.
e-g  % of bootstrap refits with (e) Abeta42 first, (f) the CSF order Abeta42 -> p-tau -> total tau,
     (g) all three CSF markers before all three cognitive scores, by genotype and cohort; Wilson 95%.

Inputs: aggregate ADNI/NACC results in figures/inputs/realdata/ (not distributed with this
repository; see figures/README.md):
  APOE_observed_orderings.csv       a,b  observed CONCORD orderings per genotype
  APOE_pair_tests.csv               a,b  adjusted P of the pairwise within-diagnosis permutation tests
  event_pair_order_summary.csv      c,d  event-pair order counts over 200 bootstrap refits
  bootstrap_orderings.csv           e-g  orderings of the 200 bootstrap refits
  apoe_cp_e4_prefix_frequencies.csv, apoe_cp_block_by_group.csv
                                    optional; if present, e-g are cross-checked against them
Genotype keys: e2 = ε2 carriers without ε4, e33 = ε3/ε3, e4 = ε4 carriers without ε2.

Reads only stored aggregate summaries; no model is fitted and no participant data are read.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import patheffects as pe  # noqa: E402
from matplotlib.cm import ScalarMappable  # noqa: E402
from matplotlib.colors import Normalize, to_rgb  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle  # noqa: E402

S.apply_style()  # includes Arial mathtext (italic P in "Adjusted P")

F_ORD = S.REAL_INPUTS / "APOE_observed_orderings.csv"
F_PAIR = S.REAL_INPUTS / "event_pair_order_summary.csv"
F_TEST = S.REAL_INPUTS / "APOE_pair_tests.csv"
BOOT = S.REAL_INPUTS / "bootstrap_orderings.csv"
XCHK_PREFIX = S.REAL_INPUTS / "apoe_cp_e4_prefix_frequencies.csv"   # optional cross-check of f
XCHK_BLOCK = S.REAL_INPUTS / "apoe_cp_block_by_group.csv"           # optional cross-check of e and g
S.require_inputs([F_ORD, F_PAIR, F_TEST, BOOT], "Fig. 5")

COHORTS = ["adni", "nacc"]
COHORT_LABEL = {"adni": "ADNI", "nacc": "NACC"}
GENOS = ["e2", "e33", "e4"]
EV = S.EVENT_ORDER  # ABETA, PTAU, TAU, MEM, EXF, LAN
CSF = ["ABETA", "PTAU", "TAU"]
COG = ["MEM", "EXF", "LAN"]

# ----------------------------------------------------------------------------- data
def read_csv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh))


orderings, ns = {}, {}
for r in read_csv(F_ORD):
    if r["panel"] == "CSFcog6" and r["analysis_role"] == "primary":
        orderings[(r["cohort"], r["genotype"])] = json.loads(r["event_order"])
        ns[(r["cohort"], r["genotype"])] = int(r["n"])

EXPECTED_ORD = {
    ("adni", "e2"): ["PTAU", "ABETA", "TAU", "LAN", "EXF", "MEM"],
    ("adni", "e33"): ["PTAU", "ABETA", "LAN", "TAU", "EXF", "MEM"],
    ("adni", "e4"): ["ABETA", "PTAU", "TAU", "LAN", "EXF", "MEM"],
    ("nacc", "e2"): ["ABETA", "LAN", "TAU", "PTAU", "MEM", "EXF"],
    ("nacc", "e33"): ["ABETA", "LAN", "PTAU", "TAU", "MEM", "EXF"],
    ("nacc", "e4"): ["ABETA", "PTAU", "TAU", "MEM", "LAN", "EXF"],
}
EXPECTED_N = {("adni", "e2"): 97, ("adni", "e33"): 560, ("adni", "e4"): 546,
              ("nacc", "e2"): 143, ("nacc", "e33"): 759, ("nacc", "e4"): 682}
for k, v in EXPECTED_ORD.items():
    assert orderings[k] == v, (k, orderings[k], v)
    assert ns[k] == EXPECTED_N[k], (k, ns[k])

# Adjusted P for epsilon3/epsilon3 vs epsilon4 (CSFcog6, primary).
adjP = {}
for r in read_csv(F_TEST):
    if r["panel"] == "CSFcog6" and r["analysis_role"] == "primary" and \
            {r["group_A"], r["group_B"]} == {"e33", "e4"}:
        v = float(r["adjusted_P"])
        assert adjP.get(r["cohort"], v) == v
        adjP[r["cohort"]] = v
for c in COHORTS:
    assert abs(adjP[c] - 0.005) < 1e-9, (c, adjP[c])

# Pairwise precedence P(row event before column event), rows/cols in fixed order.
pairs = {}
for r in read_csv(F_PAIR):
    if r["panel"] != "CSFcog6":
        continue
    a, b = r["event_A"], r["event_B"]
    na, nb, den = int(r["A_before_B_count"]), int(r["B_before_A_count"]), int(r["denominator"])
    assert na + nb == den == 200, r
    assert abs(float(r["A_before_B_fraction"]) - na / den) < 1e-9
    pairs[(r["cohort"], r["genotype"], a, b)] = (na, den)
    pairs[(r["cohort"], r["genotype"], b, a)] = (nb, den)

prec = {}   # (cohort, geno) -> 6x6 matrix (nan where undefined)
block = {}  # (cohort, geno) -> mean P(CSF before cognition) over the 9 CSF x cognition pairs
for c in COHORTS:
    for g in GENOS:
        M = np.full((6, 6), np.nan)
        for i, ei in enumerate(EV):
            for j, ej in enumerate(EV):
                if j > i:
                    k, den = pairs[(c, g, ei, ej)]
                    M[i, j] = k / den
        prec[(c, g)] = M
        block[(c, g)] = float(np.mean([M[EV.index(a), EV.index(b)] for a in CSF for b in COG]))

EXPECTED_BLOCK = {("adni", "e2"): 0.70, ("adni", "e33"): 0.75, ("adni", "e4"): 0.89,
                  ("nacc", "e2"): 0.62, ("nacc", "e33"): 0.66, ("nacc", "e4"): 0.91}
for k, v in EXPECTED_BLOCK.items():
    assert abs(block[k] - v) <= 0.005 + 1e-9, (k, block[k], v)  # expected values are 2-dp roundings

# Sanity: in epsilon4 every CSF event precedes every cognitive score in the observed ordering.
for c in COHORTS:
    o = orderings[(c, "e4")]
    assert max(o.index(e) for e in CSF) < min(o.index(e) for e in COG)
    assert o[:3] == ["ABETA", "PTAU", "TAU"]

# ----------------------------------------------------------------------------- source data
S.SRC.mkdir(parents=True, exist_ok=True)
for c, letter in (("adni", "a"), ("nacc", "b")):
    with open(S.SRC / f"Fig5{letter}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "panel", "estimator", "genotype", "n", "event", "event_label",
                    "position_1_is_earliest", "adjusted_P_e33_vs_e4"])
        for g in GENOS:
            for pos, e in enumerate(orderings[(c, g)], start=1):
                w.writerow([COHORT_LABEL[c], "CSFcog6", "CONCORD", S.GENOTYPE_LABEL[g], ns[(c, g)],
                            e, S.EVENT_LABEL[e], pos, adjP[c]])
for c, letter in (("adni", "c"), ("nacc", "d")):
    with open(S.SRC / f"Fig5{letter}.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["cohort", "panel", "genotype", "row_event", "column_event",
                    "row_before_column_count", "bootstrap_refits", "P_row_before_column",
                    "in_CSF_x_cognition_block", "block_mean_P_CSF_before_cognition"])
        for g in GENOS:
            M = prec[(c, g)]
            for i, ei in enumerate(EV):
                for j, ej in enumerate(EV):
                    if j > i:
                        k, den = pairs[(c, g, ei, ej)]
                        w.writerow([COHORT_LABEL[c], "CSFcog6", S.GENOTYPE_LABEL[g], S.EVENT_LABEL[ei],
                                    S.EVENT_LABEL[ej], k, den, f"{M[i, j]:.3f}",
                                    int(ei in CSF and ej in COG), f"{block[(c, g)]:.4f}"])

# ----------------------------------------------------------------------------- helpers
FIG_W_MM, FIG_H_MM = 180.0, 184.0
TOP_H_MM = 125.0                      # height of the a-d block
HALF = {"adni": 0.0, "nacc": 90.0}   # left edge of each cohort half (mm)
LETTER_X = 1.0                        # panel letter x inside each half (mm)
TITLE_DX = 4.5                        # panel title x offset from the letter (mm)
ROW1_BASE = 3.5                       # letter/title baseline, top row (mm from top)
ROW2_BASE = 66.5                      # letter/title baseline, bottom row (mm from top)


def mix(c, other, f):
    a, b = np.array(to_rgb(c)), np.array(to_rgb(other))
    return tuple((1 - f) * a + f * b)


CSF_TXT = S.EVENT_TEXT_COLOR["ABETA"]
COG_TXT = S.EVENT_TEXT_COLOR["MEM"]
CSF_BAND = mix(S.EVENT_COLOR["PTAU"], "white", 0.84)
HALO = [pe.Stroke(linewidth=2.9, foreground="white"), pe.Normal()]


def smooth(x0, x1, y0, y1, n=40):
    t = np.linspace(0, 1, n)
    s = t * t * (3 - 2 * t)
    return x0 + (x1 - x0) * t, y0 + (y1 - y0) * s


# ----------------------------------------------------------------------------- figure
fig = plt.figure(figsize=(S.FULL_W, FIG_H_MM * S.MM))

# ---- a, b: bump charts --------------------------------------------------------
BUMP_TOP, BUMP_H = 5.0, 55.5
BUMP_X = {"adni": 9.0, "nacc": 10.5}   # b sits 1.5 mm further right: clear gap from a's brackets
BUMP_W = 81.5
Y_LO, Y_HI = 6.45, -1.55          # y limits (inverted: position 1 at top)
for c, letter in (("adni", "a"), ("nacc", "b")):
    X0 = HALF[c]
    ax = S.add_axes_mm(fig, X0 + BUMP_X[c], BUMP_TOP, BUMP_W, BUMP_H)
    W = BUMP_W
    ax.set_xlim(0, W)
    ax.set_ylim(Y_LO, Y_HI)
    mm_per_unit = BUMP_H / abs(Y_LO - Y_HI)
    cols = {"e2": 14.0, "e4": W - 20.5}
    cols["e33"] = 0.5 * (cols["e2"] + cols["e4"])

    # y axis: positions
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_bounds(1, 6)
    ax.set_xticks([])
    ax.set_yticks(range(1, 7))
    ax.set_ylabel("Position in ordering", labelpad=2)
    ax.yaxis.set_label_coords(-0.055, 1 - (3.5 - Y_HI) / (Y_LO - Y_HI))

    # CSF band behind the three CSF nodes (and labels) of the epsilon4 column
    bx0, bx1 = cols["e4"] - 3.0, cols["e4"] + 13.2
    ax.add_patch(FancyBboxPatch((bx0, 0.56), bx1 - bx0, 2.88,
                                boxstyle="round,pad=0,rounding_size=1.6",
                                mutation_aspect=1 / mm_per_unit,
                                fc=CSF_BAND, ec="none", zorder=0.5))

    pos = {g: {e: orderings[(c, g)].index(e) + 1 for e in EV} for g in GENOS}
    for e in COG[::-1] + CSF[::-1]:  # CSF lines on top
        lc, tc = S.EVENT_LINE_COLOR[e], S.EVENT_TEXT_COLOR[e]
        xs, ys = [], []
        for g0, g1 in ((GENOS[0], GENOS[1]), (GENOS[1], GENOS[2])):
            x_, y_ = smooth(cols[g0], cols[g1], pos[g0][e], pos[g1][e])
            xs.append(x_)
            ys.append(y_)
        ax.plot(np.concatenate(xs), np.concatenate(ys), color=lc, lw=1.6, solid_capstyle="round",
                zorder=2, path_effects=HALO)
        ax.plot([cols[g] for g in GENOS], [pos[g][e] for g in GENOS], ls="none", marker="o",
                ms=5.2, mfc=lc, mec="white", mew=0.7, zorder=3)
        fw = "bold" if e in CSF else "normal"
        ax.text(cols["e2"] - 2.6, pos["e2"][e], S.EVENT_LABEL[e], ha="right", va="center",
                fontsize=6.5, color=tc, fontweight=fw)
        ax.text(cols["e4"] + 2.6, pos["e4"][e], S.EVENT_LABEL[e], ha="left", va="center",
                fontsize=6.5, color=tc, fontweight=fw)

    # column headers
    for g in GENOS:
        ax.text(cols[g], -0.38, S.GENOTYPE_LABEL[g], ha="center", va="center", fontsize=6.5,
                fontweight="bold", color=S.PAL["ink"])
        ax.text(cols[g], 0.14, rf"$\mathit{{n}}$ = {ns[(c, g)]}", ha="center", va="center", fontsize=6,
                color=S.PAL["grey_mid"])

    # significance bracket epsilon3/epsilon3 vs epsilon4 (ink)
    yb = -0.86
    ax.plot([cols["e33"], cols["e33"], cols["e4"], cols["e4"]], [yb + 0.16, yb, yb, yb + 0.16],
            color=S.PAL["ink"], lw=0.6, solid_capstyle="butt", clip_on=False)
    ax.text(0.5 * (cols["e33"] + cols["e4"]), yb - 0.1, f"Adjusted $\\it{{P}}$ = {adjP[c]:.3f}",
            ha="center", va="bottom", fontsize=6, color=S.PAL["ink"])

    # group brackets to the right of the epsilon4 labels
    xr = cols["e4"] + 14.4
    for lo, hi, lab, colr in ((1, 3, "CSF", CSF_TXT), (4, 6, "Cognition", COG_TXT)):
        ax.plot([xr - 0.7, xr, xr, xr - 0.7], [lo - 0.3, lo - 0.3, hi + 0.3, hi + 0.3],
                color=colr, lw=0.7, clip_on=False, solid_joinstyle="miter")
        ax.text(xr + 1.1, 0.5 * (lo + hi), lab, rotation=270, ha="center", va="center",
                fontsize=6, color=colr, fontweight="bold")

    S.panel_label_mm(fig, letter, x=X0 + LETTER_X, y_top=ROW1_BASE)
    S.panel_title_mm(fig, f"{COHORT_LABEL[c]}, ordering by genotype", x=X0 + LETTER_X + TITLE_DX,
                     y_top=ROW1_BASE)

# ---- c, d: precedence triangles (drawn on an mm overlay: x right, y down) ------
ov = fig.add_axes([0, 0, 1, 1])
ov.set_xlim(0, FIG_W_MM)
ov.set_ylim(FIG_H_MM, 0)
ov.set_axis_off()

CMAP, NORM = S.CMAP_DIV, Normalize(0, 1)
ROWS = EV[:5]   # Abeta ... Executive
COLS = EV[1:]   # p-tau ... Language
CELL = 5.0
GRID = 5 * CELL
GAP = 2.0                                  # between triangles of one cohort
TRI_W = 3 * GRID + 2 * GAP
GRID_L0 = {"adni": 8.5, "nacc": FIG_W_MM - 1.8 - TRI_W}  # 1.8: keeps the 1 pt outline inside 1.5 mm
HEAD_BASE = ROW2_BASE + 4.5                # genotype headers
GRID_TOP = HEAD_BASE + 11.0                # leaves room for the rotated column labels of epsilon2
MEAN_BASE = GRID_TOP + GRID + 3.6          # block-mean line
MEAN_X = 3.6 * CELL                        # block-mean value centre (under the outlined block)
LAB_PAD = 0.7                              # label to cell gap (mm)

for c, letter in (("adni", "c"), ("nacc", "d")):
    X0 = HALF[c]
    S.panel_label_mm(fig, letter, x=X0 + LETTER_X, y_top=ROW2_BASE)
    S.panel_title_mm(fig, f"{COHORT_LABEL[c]}, pairwise precedence", x=X0 + LETTER_X + TITLE_DX,
                     y_top=ROW2_BASE)
    for k, g in enumerate(GENOS):
        L, T = GRID_L0[c] + k * (GRID + GAP), GRID_TOP
        M = prec[(c, g)]
        for r, er in enumerate(ROWS):
            for q, ec in enumerate(COLS):
                i, j = EV.index(er), EV.index(ec)
                if j <= i:
                    continue
                ov.add_patch(Rectangle((L + q * CELL, T + r * CELL), CELL, CELL,
                                       fc=CMAP(NORM(M[i, j])), ec="white", lw=0.6, zorder=2))
        # staircase outline keeps near-white (undetermined) cells visible
        stair = [(L, T), (L + GRID, T), (L + GRID, T + GRID)]
        for s_ in range(5, 0, -1):
            stair += [(L + (s_ - 1) * CELL, T + s_ * CELL), (L + (s_ - 1) * CELL, T + (s_ - 1) * CELL)]
        ov.add_patch(Polygon(stair, closed=True, fc="none", ec="#B5B5B5", lw=0.4, zorder=3))
        # CSF x cognition block: rows 0-2 (CSF), columns 2-4 (Memory, Executive, Language)
        ov.add_patch(Rectangle((L + 2 * CELL, T), 3 * CELL, 3 * CELL, fc="none", ec=S.PAL["ink"],
                               lw=1.0, zorder=5))
        ov.text(L + GRID / 2, HEAD_BASE, S.GENOTYPE_LABEL[g], ha="center", va="baseline",
                fontsize=6.5, fontweight="bold", color=S.PAL["ink"])
        if g == "e2":   # row (staircase) and column (top) labels on the epsilon2 triangle only
            for r, er in enumerate(ROWS):
                ov.text(L + r * CELL - LAB_PAD, T + (r + 0.5) * CELL, S.EVENT_LABEL[er],
                        ha="right", va="center", fontsize=6, color=S.EVENT_TEXT_COLOR[er])
            for q, ec in enumerate(COLS):
                ov.text(L + (q + 0.5) * CELL, T - LAB_PAD, S.EVENT_LABEL[ec], ha="center",
                        va="bottom", rotation=90, fontsize=6, color=S.EVENT_TEXT_COLOR[ec])
            ov.text(L + MEAN_X - 3.4, MEAN_BASE, "CSF before cognition:", ha="right",
                    va="baseline", fontsize=6, color=S.PAL["grey_dark"])
        is_e4 = g == "e4"
        ov.text(L + MEAN_X, MEAN_BASE, f"{block[(c, g)]:.2f}", ha="center", va="baseline",
                fontsize=7, fontweight="bold" if is_e4 else "normal",
                color=S.PAL["blue_main"] if is_e4 else S.PAL["ink"])

# ---- shared colour bar (diverging, centred at 0.5) + key for the outlined block ----
# One legend row, centred on the page: [Column event first][colour bar][Row event first]  [key]
_r = fig.canvas.get_renderer()


def text_w_mm(txt, **kw):
    t = fig.text(0, 0, txt, fontsize=6, **kw)
    w = t.get_window_extent(_r).width / fig.dpi * 25.4
    t.remove()
    return w


CB_W, CB_H, CB_PAD, KEY_GAP, KEY_BOX = 40.0, 1.8, 1.5, 8.0, 2.8
w_col = text_w_mm("Column event first", fontweight="bold")
w_row = text_w_mm("Row event first", fontweight="bold")
w_key = text_w_mm("CSF\u2013cognition pairs")
group_w = w_col + CB_PAD + CB_W + CB_PAD + w_row + KEY_GAP + KEY_BOX + 1.2 + w_key
CB_LEFT = FIG_W_MM / 2 - group_w / 2 + w_col + CB_PAD
CB_TOP = MEAN_BASE + 8.0
cax = S.add_axes_mm(fig, CB_LEFT, CB_TOP, CB_W, CB_H)
cb = fig.colorbar(ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation="horizontal")
cb.set_ticks([0, 0.25, 0.5, 0.75, 1])
cb.set_ticklabels(["0", "0.25", "0.5", "0.75", "1"])
cb.outline.set_linewidth(0.4)
cb.outline.set_edgecolor("#9A9A9A")
cax.tick_params(length=1.6, width=0.5, pad=1.0, labelsize=6)
yc = CB_TOP + CB_H / 2
ov.text(CB_LEFT + CB_W / 2, CB_TOP - 1.3,
        "Probability that the row event precedes the column event (200 bootstrap refits)",
        ha="center", va="baseline", fontsize=6, color=S.PAL["ink"])
ov.text(CB_LEFT - CB_PAD, yc, "Column event first", ha="right", va="center", fontsize=6,
        fontweight="bold", color=S.EVENT_TEXT_COLOR["MEM"])
ov.text(CB_LEFT + CB_W + CB_PAD, yc, "Row event first", ha="left", va="center", fontsize=6,
        fontweight="bold", color=S.PAL["blue_main"])
KEY_L = CB_LEFT + CB_W + CB_PAD + w_row + KEY_GAP
ov.add_patch(Rectangle((KEY_L, yc - KEY_BOX / 2), KEY_BOX, KEY_BOX, fc="none", ec=S.PAL["ink"],
                       lw=1.0))
ov.text(KEY_L + KEY_BOX + 1.2, yc, "CSF\u2013cognition pairs", ha="left", va="center", fontsize=6,
        color=S.PAL["ink"])


# ============================================================================= e-g: gradients
# (drawn below a-d, offset by G_Y0 mm)
from collections import defaultdict  # noqa: E402
import matplotlib.text  # noqa: E402

G_Y0 = TOP_H_MM + 1.5
N_BOOT = 200

g_orders = defaultdict(list)
for r in read_csv(BOOT):
    if r["panel"] != "CSFcog6" or r["status"] != "ok":
        continue
    g_orders[(r["cohort"], r["genotype"])].append(json.loads(r["event_order"]))


def amyloid_first(o):
    return o[0] == "ABETA"


def csf_cascade(o):
    return o.index("ABETA") < o.index("PTAU") < o.index("TAU")


def csf_before_cog(o):
    return max(o.index(e) for e in CSF) < min(o.index(e) for e in COG)


G_METRICS = {"e": amyloid_first, "f": csf_cascade, "g": csf_before_cog}
g_counts = {}
for (coh, g), lst in g_orders.items():
    assert len(lst) == N_BOOT, (coh, g, len(lst))
    assert all(sorted(o) == sorted(CSF + COG) for o in lst)
    for p, f in G_METRICS.items():
        g_counts[(p, coh, g)] = (sum(f(o) for o in lst), len(lst))
G_EXPECTED = {
    "e": {"adni": [7.0, 42.5, 83.5], "nacc": [46.0, 80.5, 99.0]},
    "f": {"adni": [13.5, 66.0, 100.0], "nacc": [30.0, 54.5, 82.5]},
    "g": {"adni": [13.5, 17.0, 48.5], "nacc": [2.5, 9.0, 53.0]},
}
for p, d in G_EXPECTED.items():
    for coh, vals in d.items():
        for g, v in zip(GENOS, vals):
            k, n = g_counts[(p, coh, g)]
            assert abs(k / n - v / 100) <= 0.002, (p, coh, g, k / n, v)
# cross-checks against the separately computed refit summaries, when available
if XCHK_PREFIX.is_file():
    for r in read_csv(XCHK_PREFIX):
        k, n = g_counts[("f", r["cohort"], r["genotype"])]
        assert abs(k / n - float(r["boot_frac_csf_internal_ABETA_lt_PTAU_lt_TAU"])) <= 0.002
    print(f"[check] f agrees with {S.source_name(XCHK_PREFIX)}")
else:
    print(f"[note] {S.source_name(XCHK_PREFIX)} not found; cross-check of f skipped")
if XCHK_BLOCK.is_file():
    for r in read_csv(XCHK_BLOCK):
        if r["panel"] != "CSFcog6":
            continue
        key = {"block_all_csf_first": "g", "ABETA_rank1": "e"}.get(r["statistic"])
        if key:
            k, n = g_counts[(key, r["cohort"], r["genotype"])]
            assert abs(k / n - float(r["bootstrap_mean"])) <= 0.002, r
    print(f"[check] e and g agree with {S.source_name(XCHK_BLOCK)}")
else:
    print(f"[note] {S.source_name(XCHK_BLOCK)} not found; cross-check of e and g skipped")

G_NAME = {"adni": "ADNI", "nacc": "NACC"}
G_COL = {c: S.COHORT[G_NAME[c]]["color"] for c in COHORTS}
G_MRK = {c: S.COHORT[G_NAME[c]]["marker"] for c in COHORTS}
G_MS = {c: 4.2 * S.MARKER_SCALE[G_MRK[c]] for c in COHORTS}
GX = [0, 1, 2]
G_DODGE = {"e": {"adni": 0.0, "nacc": 0.0}, "f": {"adni": -0.07, "nacc": 0.07}, "g": {"adni": -0.07, "nacc": 0.07}}
G_TITLES = {"e": "Aβ42 placed first", "f": "CSF order Aβ42 → p-tau → total tau", "g": "All CSF before all cognition"}
G_MEASURE = {
    "e": "amyloid event placed first",
    "f": "CSF events in relative order Abeta42 < p-tau < total tau",
    "g": "all three CSF events before all three cognitive scores",
}


def pct_label(v):
    return "100%" if abs(v - 100) < 1e-9 else S.fmt_pct(v)


def spread_labels(ys, min_gap):
    ys = list(ys)
    if abs(ys[0] - ys[1]) >= min_gap:
        return ys
    m = (ys[0] + ys[1]) / 2
    if ys[0] <= ys[1]:
        return [m - min_gap / 2, m + min_gap / 2]
    return [m + min_gap / 2, m - min_gap / 2]


G_AX_TOP = G_Y0 + 6.2
G_AX_H = 42.0
G_AX_W = 37.7
G_TICK_W, G_LABEL_OUT, G_GAP = 5.0, 13.3, 2.3
G_LETTER_X, G_AX_LEFT = {"e": 1.0}, {"e": 10.8}
for prev, p in [("e", "f"), ("f", "g")]:
    G_LETTER_X[p] = G_AX_LEFT[prev] + G_AX_W + G_LABEL_OUT + G_GAP
    G_AX_LEFT[p] = G_LETTER_X[p] + G_TICK_W
G_LABEL_X = {"e": 2.22, "f": 2.17, "g": 2.17}
g_records = {p: [] for p in G_METRICS}
for p in ["e", "f", "g"]:
    ax = S.add_axes_mm(fig, G_AX_LEFT[p], G_AX_TOP, G_AX_W, G_AX_H)
    ends = {}
    for coh in ["nacc", "adni"]:
        ks = [g_counts[(p, coh, g)] for g in GENOS]
        ys = [100 * k / n for k, n in ks]
        cis = [tuple(100 * v for v in S.wilson(k, n)) for k, n in ks]
        xs = [x + G_DODGE[p][coh] for x in GX]
        yerr = [[y - c[0] for y, c in zip(ys, cis)], [c[1] - y for y, c in zip(ys, cis)]]
        ax.errorbar(xs, ys, yerr=yerr, fmt="none", ecolor=G_COL[coh], zorder=4, **S.ERRBAR)
        ax.plot(xs, ys, color=G_COL[coh], lw=1.2, zorder=5, solid_capstyle="round")
        ax.plot(xs, ys, ls="none", marker=G_MRK[coh], ms=G_MS[coh], mfc=G_COL[coh], mec="white", mew=0.5,
                zorder=6)
        ends[coh] = (xs[-1], ys[-1])
        for g, x, y, (k, n), (lo, hi) in zip(GENOS, xs, ys, ks, cis):
            g_records[p].append(dict(
                panel=f"Fig5{p}", measure=G_MEASURE[p], cohort=G_NAME[coh],
                series="Common panel (CSF + cognition, 6 events)", inventory_key="CSFcog6",
                genotype=S.GENOTYPE_LABEL[g], x_plotted=round(x, 3), n_refits=n, k_refits=k,
                percent=round(y, 2), wilson95_lo=round(lo, 2), wilson95_hi=round(hi, 2),
                source=S.source_name(BOOT)))
    ly = spread_labels([ends["adni"][1], ends["nacc"][1]], min_gap=9)
    for coh, y in zip(COHORTS, ly):
        t = ax.text(G_LABEL_X[p], y, G_NAME[coh], color=G_COL[coh], fontsize=6.5, fontweight="bold",
                    va="center", ha="left", clip_on=False)
        ax.annotate(" " + pct_label(ends[coh][1]), xy=(1, 0.5), xycoords=t, va="center", ha="left",
                    fontsize=6, color=G_COL[coh], annotation_clip=False)
    ax.set_xlim(-0.25, 2.2)
    ax.set_ylim(-4, 104)
    ax.spines["left"].set_bounds(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(GX)
    ax.set_xticklabels([S.GENOTYPE_LABEL[g] for g in GENOS])
    ax.tick_params(axis="x", length=0, pad=2.5)
    ax.set_xlabel(r"$\it{APOE}$ genotype", labelpad=2)
    if p == "e":
        ax.set_ylabel("% of bootstrap refits")
    S.panel_label_mm(fig, p, x=G_LETTER_X[p], y_top=G_Y0 + 3.5)
    S.panel_title_mm(fig, G_TITLES[p], x=G_LETTER_X[p] + TITLE_DX, y_top=G_Y0 + 3.5)

for p, recs in g_records.items():
    with open(S.SRC / f"Fig5{p}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0].keys()))
        w.writeheader()
        w.writerows(recs)

fig.canvas.draw()
_r = fig.canvas.get_renderer()
_bbs = [t.get_window_extent(_r) for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
_px = fig.dpi / 25.4
print("text extent (mm): left %.1f, right %.1f, top %.1f, bottom %.1f" % (
    min(b.x0 for b in _bbs) / _px, max(b.x1 for b in _bbs) / _px,
    FIG_H_MM - max(b.y1 for b in _bbs) / _px, FIG_H_MM - min(b.y0 for b in _bbs) / _px))

paths = S.save(fig, "Fig5_apoe")
print("\n".join(str(p) for p in paths))
for c in COHORTS:
    print(c, {g: round(block[(c, g)], 4) for g in GENOS}, "adjP", adjP[c])
for p in G_METRICS:
    print(p, {coh: [100 * g_counts[(p, coh, g)][0] / N_BOOT for g in GENOS] for coh in COHORTS})
