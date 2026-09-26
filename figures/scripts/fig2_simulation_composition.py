"""Fig. 2 | Effect of diagnostic composition on orderings estimated in simulation.

180 x 142 mm, Nature Communications print size.
  a  mean-position gap between CN-heavy and AD-heavy groups, five models (dumbbells); every printed
     number is the increase due to composition (different - same composition) at n participants;
     complete data have no same-composition 4n cell, so sample size is shown in d only
  b  per-event position shift vs event separability (AUC)
  c  systematically misordered event pairs in the three groups
  d  gap vs sample size (separately fitted DEBM and CONCORD, log10 axis); these runs have missing
     measurements, as in Fig. 3; the x axis gives the total number of participants

Simulated biomarkers are named "Biomarker 1-14" in the input order of Supplementary Table 1
(design_v2.BIOMARKERS: 1-5 modeled on CSF markers, 6-7 on cognitive scores, 8-14 on imaging
measures). The input files use the generator's code names (ABETA, ..., Precuneus).

Model labels in the input files: kde_gmm = likelihood EBM, saebm_conjugate = stage-aware EBM,
repaired = separately fitted DEBM, shared = pooled-score DEBM, invariant_min = CONCORD.
Group keys g1/g2/g3 = CN-heavy/intermediate/AD-heavy group.

Inputs (figures/inputs/simulation/): mechanism_all_cells.csv, v2_mechanism_gap_intervals.csv,
runs_random_truth_eventwise.csv, S1_event_parameters.csv (= design_v2.BIOMARKERS).
All plotted numbers are read from stored aggregate files; nothing is re-simulated.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib import colors as mcolors  # noqa: E402

CELLS = S.SIM_INPUTS / "mechanism_all_cells.csv"
INTERVALS = S.SIM_INPUTS / "v2_mechanism_gap_intervals.csv"
EVENTWISE = S.SIM_INPUTS / "runs_random_truth_eventwise.csv"
AUC_CSV = S.SIM_INPUTS / "S1_event_parameters.csv"
S.require_inputs([CELLS, INTERVALS, EVENTWISE, AUC_CSV], "Fig. 2")

STEM = "Fig2_simulation_composition"
W_MM, H_MM = 180.0, 142.0
M = S.METHOD
INK = S.PAL["ink"]
GREY = S.PAL["grey_mid"]
BLUE = S.PAL["blue_main"]
TC = S.TEXT_COLOR

TOL_GAP = 0.01
N_BASE = 971


# ----------------------------------------------------------------------------- helpers
def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def lighten(color, f):
    """Mix a colour with white; f = 0 keeps the colour, f = 1 gives white."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def fpos(x, signed=False):
    """S.fmt_pos with decimal half-up rounding (an exact .xx5 stored as .xx4999... still rounds up)."""
    from decimal import Decimal, ROUND_HALF_UP
    v = Decimal(repr(round(x, 9)))
    q = Decimal("0.001") if abs(v) < Decimal("0.1") else Decimal("0.01")
    s = f"{v.quantize(q, rounding=ROUND_HALF_UP):+}"
    return s if signed else s.lstrip("+")


def check(name, got, want, tol):
    assert abs(got - want) <= tol, f"{name}: got {got:.4f}, expected {want:.4f} (tol {tol})"


cells = read_csv(CELLS)
ivals = read_csv(INTERVALS)


def cell(run, c, engine):
    rows = [r for r in cells if r["run_id"] == run and r["cell"] == c and r["engine"] == engine]
    assert len(rows) == 1, (run, c, engine, len(rows))
    return rows[0]


def gap_ci(run, c, engine):
    rows = [r for r in ivals if r["run_id"] == run and r["cell"] == c and r["engine"] == engine]
    assert len(rows) == 1, (run, c, engine, len(rows))
    r = rows[0]
    g_cells = float(cell(run, c, engine)["mean_position_gap_g1_g3"])
    g = float(r["gap"])
    assert abs(g - g_cells) < 1e-9, (run, c, engine, g, g_cells)
    return g, float(r["low"]), float(r["high"]), int(r["R"])


# ----------------------------------------------------------------------------- data: a, c
# complete-data comparison; (run, engine, cell different composition, cell same composition), n participants
COMPLETE = {
    "likelihood": ("method_kde_mechanism_a", "kde_gmm", "REF_S_CD", "IID_S_CD"),
    "saebm": ("method_saebm_a", "saebm_conjugate", "REF_S", "IID_S"),
    "separate": ("method_saebm_reference_a", "repaired", "REF_S_CD", "IID_S_CD"),
    "pooled": ("method_saebm_reference_a", "shared", "REF_S_CD", "IID_S_CD"),
    "concord": ("method_saebm_reference_a", "invariant_min", "REF_S_CD", "IID_S_CD"),
}
EXP_A = {  # different composition, same composition
    "likelihood": (1.544, 0.906),
    "saebm": (2.391, 1.495),
    "separate": (2.010, 0.875),
    "pooled": (0.997, 0.080),
    "concord": (0.116, 0.077),
}
EXP_INC = {"likelihood": 0.64, "saebm": 0.90, "separate": 1.14, "pooled": 0.92, "concord": 0.04}
EXP_C = {
    "likelihood": (29, 21, 3),
    "saebm": (11, 0, 0),
    "separate": (23, 10, 2),
    "pooled": (13, 2, 2),
    "concord": (1, 1, 1),
}

A = {}
C = {}
for k, (run, eng, c_diff, c_same) in COMPLETE.items():
    for c_ in (c_diff, c_same):
        assert cell(run, c_, eng)["missingness_enabled"] == "False", (k, c_)
        assert int(cell(run, c_, eng)["n_subjects"]) == N_BASE, (k, c_)
    d = gap_ci(run, c_diff, eng)
    s = gap_ci(run, c_same, eng)
    A[k] = dict(diff=d, same=s, inc=d[0] - s[0])
    for name, got, want in zip(("diff", "same"), (d[0], s[0]), EXP_A[k]):
        check(f"a {k} {name}", got, want, TOL_GAP)
    check(f"a {k} increase", round(d[0] - s[0], 2), EXP_INC[k], 1e-9)
    row = cell(run, c_diff, eng)
    C[k] = tuple(int(row[f"majority_wrong_pairs_g{g}"]) for g in (1, 2, 3))
    assert C[k] == EXP_C[k], (k, C[k], EXP_C[k])
    C[k + "_err"] = tuple(float(row[f"mean_error_g{g}"]) for g in (1, 2, 3))
K_EVENTS = int(cell("method_saebm_reference_a", "REF_S_CD", "invariant_min")["K_events"])
N_PAIRS = K_EVENTS * (K_EVENTS - 1) // 2
assert N_PAIRS == 91
# the increases printed in a
EXP_INC_TEXT = {"likelihood": "+0.64", "saebm": "+0.90", "separate": "+1.14", "pooled": "+0.92", "concord": "+0.039"}
for k, want in EXP_INC_TEXT.items():
    assert fpos(A[k]["inc"], signed=True) == want, (k, fpos(A[k]["inc"], signed=True), want)

# ----------------------------------------------------------------------------- data: b
auc_rows = read_csv(AUC_CSV)
AUC = {r["event"]: float(r["component_AUC"]) for r in auc_rows}
# simulated biomarkers are numbered in the input order of Supplementary Table 1 (design_v2.BIOMARKERS)
SIM_CODES = ["ABETA", "PTAU", "TAU", "NG", "NFL", "ADAS13", "MMSE", "Hippocampus", "Entorhinal", "MidTemp",
             "Fusiform", "WholeBrain", "Ventricles", "Precuneus"]
assert [r["event"] for r in auc_rows] == SIM_CODES, [r["event"] for r in auc_rows]
EVENT_NAME = {c: f"Biomarker {i + 1}" for i, c in enumerate(SIM_CODES)}
MODALITY = {c: ("CSF" if i < 5 else "cognitive" if i < 7 else "imaging") for i, c in enumerate(SIM_CODES)}
EW_SRC = {
    "likelihood": ("method_kde_calibration_a", "REF_H0_CD", "kde_gmm"),
    "separate": ("confirm_core_a", "REF_H0", "repaired"),
    "pooled": ("method_calibration_a", "REF_H0", "shared"),
    "concord": ("method_calibration_a", "REF_H0", "invariant_min"),
}
ew = read_csv(EVENTWISE)
B = {}
for k, (run, c_, eng) in EW_SRC.items():
    rows = [r for r in ew if r["run_id"] == run and r["cell"] == c_ and r["engine"] == eng]
    assert len(rows) == 14, (k, len(rows))
    ev = [r["biomarker"] for r in rows]
    x = np.array([AUC[e] for e in ev])
    y = np.array([float(r["mean_pos_g1_minus_g3"]) for r in rows])
    se = np.array([float(r["se"]) for r in rows])
    r_ = float(np.corrcoef(x, y)[0, 1])
    slope, icpt = np.polyfit(x, y, 1)
    B[k] = dict(events=ev, x=x, y=y, se=se, r=r_, slope=slope, icpt=icpt, R=int(rows[0]["R_ok"]))
check("b r separate", B["separate"]["r"], 0.99, 0.01)
check("b r pooled", B["pooled"]["r"], 0.96, 0.01)
check("b r likelihood", B["likelihood"]["r"], 0.66, 0.01)
assert np.max(np.abs(B["concord"]["y"])) <= 0.14 + 1e-9, np.max(np.abs(B["concord"]["y"]))

# ----------------------------------------------------------------------------- data: d
# Panel d shows separately fitted DEBM and CONCORD only.
SIZES = [("n", 1), ("4n", 4), ("16n", 16)]
D_SRC = {  # method -> list of (run, diff cell, same cell) per size
    "separate": [("method_mechanism_a", "REF_S", "IID_S"), ("method_mechanism_a", "REF_S_N4", "IID_S_N4"),
                 ("method_mechanism_16n_a", "REF_S_N16", "IID_S_N16")],
    "concord": [("method_mechanism_a", "REF_S", "IID_S"), ("method_mechanism_a", "REF_S_N4", "IID_S_N4"),
                ("method_mechanism_16n_a", "REF_S_N16", "IID_S_N16")],
}
D_ENGINE = {"separate": "repaired", "concord": "invariant_min"}
EXP_D = {
    "separate": ((2.019, 2.208, 2.221), (0.814, 0.496, 0.259)),
    "concord": ((0.119, 0.057, 0.035), (0.070, 0.036, 0.044)),
}
D = {}
for k, srcs in D_SRC.items():
    D[k] = {"diff": [], "same": [], "mult": []}
    for i, (run, cd, cs) in enumerate(srcs):
        assert cell(run, cd, D_ENGINE[k])["missingness_enabled"] == "True"
        assert int(cell(run, cd, D_ENGINE[k])["n_subjects"]) == N_BASE * SIZES[i][1]
        dd = gap_ci(run, cd, D_ENGINE[k])
        ss = gap_ci(run, cs, D_ENGINE[k])
        check(f"d {k} diff {SIZES[i][0]}", dd[0], EXP_D[k][0][i], TOL_GAP)
        check(f"d {k} same {SIZES[i][0]}", ss[0], EXP_D[k][1][i], TOL_GAP)
        D[k]["diff"].append(dd)
        D[k]["same"].append(ss)
        D[k]["mult"].append(SIZES[i][1])

# ----------------------------------------------------------------------------- source data
S.SRC.mkdir(parents=True, exist_ok=True)


def write_csv(name, header, rows):
    with open(S.SRC / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


ROW_ORDER = ["likelihood", "saebm", "separate", "pooled", "concord"]
rows_a = []
for k in ROW_ORDER:
    run, eng, c_diff, c_same = COMPLETE[k]
    for cond, c_, key, n_ in (("different composition", c_diff, "diff", N_BASE),
                              ("same composition", c_same, "same", N_BASE)):
        g, lo, hi, R = A[k][key]
        rows_a.append([M[k]["label"], cond, n_, run, c_, eng, R, f"{g:.4f}", f"{lo:.4f}", f"{hi:.4f}"])
    rows_a.append([M[k]["label"], "increase due to composition (different - same)", N_BASE, run,
                   f"{c_diff}-{c_same}", eng, "", f"{A[k]['inc']:.4f}", "", ""])
write_csv("Fig2a.csv", ["model", "condition", "total_participants", "run_id", "cell", "engine", "datasets",
                        "mean_position_gap_CNheavy_ADheavy", "ci95_low", "ci95_high"], rows_a)

rows_b = []
for k in ["likelihood", "separate", "pooled", "concord"]:
    run, c_, eng = EW_SRC[k]
    b = B[k]
    for e, x, y, se in zip(b["events"], b["x"], b["y"], b["se"]):
        rows_b.append([M[k]["label"], run, c_, eng, b["R"], EVENT_NAME[e], MODALITY[e], f"{x:.2f}", f"{y:.4f}",
                       f"{se:.4f}", f"{b['r']:.4f}", f"{b['slope']:.4f}", f"{b['icpt']:.4f}"])
write_csv("Fig2b.csv", ["model", "run_id", "cell", "engine", "datasets", "biomarker", "modality", "component_AUC",
                        "position_shift_CNheavy_minus_ADheavy", "se", "pearson_r_vs_AUC",
                        "ls_slope", "ls_intercept"], rows_b)

GROUPS = ["CN-heavy", "Intermediate", "AD-heavy"]
rows_c = []
for k in ROW_ORDER:
    run, eng, c_diff, _ = COMPLETE[k]
    for gi, g in enumerate(GROUPS):
        rows_c.append([M[k]["label"], run, c_diff, eng, f"g{gi + 1}", g, C[k][gi], N_PAIRS,
                       f"{C[k + '_err'][gi]:.4f}"])
write_csv("Fig2c.csv", ["model", "run_id", "cell", "engine", "group_key", "group", "systematically_misordered_pairs",
                        "event_pairs", "mean_error_vs_truth"], rows_c)

rows_d = []
for k in ["separate", "concord"]:
    for i, (run, cd, cs) in enumerate(D_SRC[k]):
        mult = D[k]["mult"][i]
        for cond, c_, v in (("different composition", cd, D[k]["diff"][i]), ("same composition", cs, D[k]["same"][i])):
            rows_d.append([M[k]["label"], cond, f"{mult}n", mult * N_BASE, run, c_, D_ENGINE[k], v[3],
                           f"{v[0]:.4f}", f"{v[1]:.4f}", f"{v[2]:.4f}"])
write_csv("Fig2d.csv", ["model", "condition", "sample_size", "total_participants", "run_id", "cell", "engine",
                        "datasets", "mean_position_gap_CNheavy_ADheavy", "ci95_low", "ci95_high"], rows_d)

# ----------------------------------------------------------------------------- figure
fig = plt.figure(figsize=(W_MM * S.MM, H_MM * S.MM))

# layout, mm from the top-left corner of the 180 x 142 mm page
LX_L, LX_R = 1.0, 104.5            # panel-letter x (left / right column)
TITLE_DX = 4.5                     # title starts ~3 mm right of the letter (as in Figs 1 and 3-5)
ROW1_Y, ROW2_Y = 10.5, 76.5        # letter / title baselines
TOP_Y, TOP_H = 14.0, 45.5          # top-row axes
BOT_Y, BOT_H = 80.0, 53.6          # bottom-row axes
AX_R_X, AX_R_W = 117.5, 60.8       # right-column axes (b, d)
A_AX_X0, A_AX_W = 26.5, 70.8      # a: right edge flush with c's direct labels (97.3 mm)
ax_a = S.add_axes_mm(fig, A_AX_X0, TOP_Y, A_AX_W, TOP_H)
ax_b = S.add_axes_mm(fig, AX_R_X, TOP_Y, AX_R_W, TOP_H)
ax_c = S.add_axes_mm(fig, 13.5, BOT_Y, 57.0, BOT_H)
ax_d = S.add_axes_mm(fig, AX_R_X, BOT_Y, AX_R_W, BOT_H)

for x, y, s, title in ((LX_L, ROW1_Y, "a", "Same vs different composition"),
                       (LX_R, ROW1_Y, "b", "Shift of each event"),
                       (LX_L, ROW2_Y, "c", "Misordered event pairs"),
                       (LX_R, ROW2_Y, "d", "Larger samples")):
    S.panel_label_mm(fig, s, x=x, y_top=y)
    S.panel_title_mm(fig, title, x=x + TITLE_DX, y_top=y)

S.method_legend(fig, S.METHOD_ORDER, y_top=1.0)

GAP_LABEL = "Mean-position gap, CN-heavy vs AD-heavy\n(positions per event)"
assert GAP_LABEL.replace(", CN-heavy vs AD-heavy\n", " ") == S.LABEL["gap"]

# ============================== a: dumbbells
ax = ax_a
A_XL = (0, 2.6)
ypos = {k: len(ROW_ORDER) - 1 - i for i, k in enumerate(ROW_ORDER)}  # likelihood on top, CONCORD at bottom
for k in ROW_ORDER:
    is_c = k == "concord"
    col = M[k]["color"]
    mk = M[k]["marker"]
    y = ypos[k]
    d, s = A[k]["diff"], A[k]["same"]
    msz = S.ms(k, 5.8 if is_c else 5.2)
    z0 = 10 if is_c else 0          # CONCORD drawn last
    # composition-driven increase: translucent thick segment between the two compositions
    ax.plot([s[0], d[0]], [y, y], color=col, alpha=0.30, lw=6.5 if is_c else 5.5, solid_capstyle="butt",
            zorder=2 + z0)
    # different composition (n): filled marker on its 95% interval (forest-plot style: one continuous line
    # through the point). A side that would reach less than 0.3 mm beyond the marker is not drawn, so no
    # stray spur appears next to a marker (as in d).
    hw = (msz / 2) * 25.4 / 72 / A_AX_W * (A_XL[1] - A_XL[0])
    min_len = 0.3 / A_AX_W * (A_XL[1] - A_XL[0])
    lo_ = d[1] if (d[0] - hw) - d[1] >= min_len else d[0]
    hi_ = d[2] if d[2] - (d[0] + hw) >= min_len else d[0]
    if hi_ > lo_:
        ax.plot([lo_, hi_], [y, y], color=col, lw=S.ERRBAR["elinewidth"] + 0.2, zorder=4 + z0,
                solid_capstyle="butt")
    # (CONCORD's two markers overlap, so its filled marker keeps a white rim to stay distinct)
    ax.plot([d[0]], [y], ls="none", marker=mk, ms=msz, mfc=col, mec="white" if is_c else col,
            mew=0.6 if is_c else 0.4, zorder=6 + z0)
    # same composition (n): hollow marker (reference condition)
    ax.plot([s[0]], [y], ls="none", marker=mk, ms=msz, mfc="white", mec=col, mew=1.0 if is_c else 0.9,
            zorder=5 + z0)
    inc = fpos(A[k]["inc"], signed=True)
    if is_c:
        ax.text(d[2] + 0.08, y, inc, color=BLUE, fontsize=7, fontweight="bold", ha="left", va="center",
                zorder=20)
    else:
        ax.text((s[0] + d[0]) / 2, y + 0.19, inc, color=TC[k], fontsize=6, ha="center", va="bottom")

YL_A = (-0.55, 5.45)
ax.set_ylim(*YL_A)
ax.set_yticks([ypos[k] for k in ROW_ORDER])
ax.set_yticklabels([M[k]["label"] for k in ROW_ORDER], fontsize=6)
for t, k in zip(ax.get_yticklabels(), ROW_ORDER):
    if k == "concord":
        t.set_fontweight("bold")
        t.set_color(BLUE)
ax.tick_params(axis="y", length=0, pad=3)
ax.spines["left"].set_visible(False)
ax.set_xlim(*A_XL)
ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0, 2.5])
ax.set_xticklabels(["0", "0.5", "1.0", "1.5", "2.0", "2.5"])
ax.set_xlabel(GAP_LABEL, linespacing=1.1)

# light-blue call-out band behind the CONCORD row, from its row label to the end of the axis
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_lab = [t for t in ax.get_yticklabels() if t.get_text() == M["concord"]["label"]][0]
_lx0 = ax.transAxes.inverted().transform(_lab.get_window_extent(_r).p0)[0]
_pad = 1.6 / A_AX_W
S.concord_band(ax, ypos["concord"], h=0.84, xmin=_lx0 - _pad, xmax=1.0, clip_on=False)

# condition key inside the panel, above the top row
a_handles = [  # ncol=2 fills by column: row 1 reads "Same | Different", row 2 "Increase due to composition"
    Line2D([0], [0], ls="none", marker="o", ms=4.4, mfc="white", mec=GREY, mew=0.9, label="Same composition"),
    Patch(facecolor=lighten(GREY, 0.6), edgecolor="none", label="Increase due to composition"),
    Line2D([0], [0], ls="none", marker="o", ms=4.4, mfc=GREY, mec=GREY, mew=0.4, label="Different composition"),
]
ax.legend(handles=a_handles, loc="upper left", bbox_to_anchor=(-0.02, 1.02), ncol=2, fontsize=6,
          handlelength=1.3, handleheight=0.8, labelspacing=0.3, columnspacing=1.2, borderaxespad=0.0)

# ============================== b: per-event shift vs AUC
ax = ax_b
S.concord_band(ax, 0.0, h=0.28)
ax.axhline(0, **S.REF_LINE)
xx = np.array([min(AUC.values()), max(AUC.values())])  # fitted range = observed AUC range
for k in ["likelihood", "pooled", "separate"]:
    b = B[k]
    ax.plot(xx, b["slope"] * xx + b["icpt"], color=M[k]["color"], lw=0.8, zorder=2)
B_DODGE = {"likelihood": -0.003, "pooled": 0.003, "separate": 0.0, "concord": 0.0}
for k in ["likelihood", "pooled", "separate", "concord"]:
    b = B[k]
    is_c = k == "concord"
    ax.plot(b["x"] + B_DODGE[k], b["y"], ls="none", marker=M[k]["marker"], ms=S.ms(k, 4.3 if is_c else 4.0),
            mfc=M[k]["color"], mec="white", mew=0.4, zorder=8 if is_c else 4)
# CONCORD band named in the plot
ax.text(0.932, -0.36, "CONCORD, all within ±0.14", color=BLUE, fontsize=6, fontweight="bold", ha="right",
        va="top", zorder=9)
# readable names for the extreme events (separately fitted DEBM points)
sep = B["separate"]
lab = {  # event: (text x, text y, ha, leader line)
    "ADAS13": (0.912, 3.954, "right", False),
    "ABETA": (0.872, 2.80, "right", True),
    "Precuneus": (0.708, -2.392, "left", False),
    "NG": (0.686, -3.02, "left", True),
}
for e, (tx, ty, ha, lead) in lab.items():
    i = sep["events"].index(e)
    x0, y0 = sep["x"][i], sep["y"][i]
    ax.text(tx, ty, EVENT_NAME[e], fontsize=5.5, color=INK, ha=ha, va="center")
    if lead:  # short leader from the label towards the point, stopping short of the marker
        f0, f1 = 0.18, 0.72
        ax.plot([tx + (x0 - tx) * f0, tx + (x0 - tx) * f1], [ty + (y0 - ty) * f0, ty + (y0 - ty) * f1],
                color=GREY, lw=0.4, zorder=3)
b_handles = []
for k in ["separate", "pooled", "likelihood"]:
    b_handles.append(Line2D([0], [0], color=M[k]["color"], lw=0.8, marker=M[k]["marker"], ms=S.ms(k, 3.8),
                            mec="white", mew=0.4, label=rf"$\mathit{{r}}$ = {B[k]['r']:.2f}"))
lg = ax.legend(handles=b_handles, loc="upper left", bbox_to_anchor=(0.0, 1.02), fontsize=6, handlelength=1.8,
               labelspacing=0.3, borderaxespad=0.1)
for t, k in zip(lg.get_texts(), ["separate", "pooled", "likelihood"]):
    t.set_color(TC[k])
ax.set_xlim(0.665, 0.935)
ax.set_xticks([0.70, 0.75, 0.80, 0.85, 0.90])
ax.set_ylim(-3.3, 4.4)
ax.set_yticks([-2, 0, 2, 4])
ax.set_yticklabels(["−2", "0", "2", "4"])
ax.set_xlabel("Event separability (AUC)")
ax.set_ylabel(S.LABEL["shift"].replace("shift, ", "shift,\n"), linespacing=1.1)

# ============================== c: systematically misordered pairs (groups ordered by CN share)
ax = ax_c
XC = {"AD-heavy": 0, "Intermediate": 1, "CN-heavy": 2}
xg = np.array([XC[g] for g in GROUPS], float)            # CN-heavy, Intermediate, AD-heavy -> 2, 1, 0
dodge = {"likelihood": 0.0, "saebm": 0.035, "separate": 0.035, "pooled": -0.035, "concord": 0.0}
for k in ROW_ORDER:
    is_c = k == "concord"
    col = M[k]["color"]
    y = np.array(C[k], float)
    x = xg + dodge[k]
    o = np.argsort(x)
    ax.plot(x[o], y[o], color=col, lw=1.7 if is_c else 0.9, marker=M[k]["marker"],
            ms=S.ms(k, 5.2 if is_c else 4.4), mfc=col, mec="white", mew=0.5, zorder=10 if is_c else 3)
    # direct label at the line end (CN-heavy group), outside the axes: model name and count
    ax.text(2.0 + 0.15, y[0], f"{M[k]['label']}  {int(y[0])}", color=TC[k], fontsize=6.5 if is_c else 6,
            fontweight="bold" if is_c else "normal", ha="left", va="center", clip_on=False)
ax.text(2.02, 1 - 1.95, f"1 of {N_PAIRS} pairs in every group", color=BLUE, fontsize=6.5, fontweight="bold",
        ha="right", va="top")
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(["AD-heavy group", "Intermediate group", "CN-heavy group"])
ax.set_xlim(-0.1, 2.1)
ax.set_ylim(-3.4, 31)
ax.set_yticks([0, 10, 20, 30])
ax.set_ylabel(f"Systematically misordered\nevent pairs (of {N_PAIRS})", linespacing=1.1)

# ============================== d: sample size, with missing measurements (log10 gap axis)
ax = ax_d
YL_D = (0.02, 3.3)
XL_D = (1 / 1.5, 16 * 2.7)
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=10)
ax.set_xlim(*XL_D)
ax.set_ylim(*YL_D)
mm_per_dec = BOT_H / (np.log10(YL_D[1]) - np.log10(YL_D[0]))
D_DODGE = 0.11                     # log2 units (about 1.1 mm): different right, same left
drawn_bars = []
for k in ["separate", "concord"]:
    is_c = k == "concord"
    col = M[k]["color"]
    mk = M[k]["marker"]
    z0 = 10 if is_c else 0
    for cond in ("same", "diff"):
        # dodge on the log2 axis, same for both models, so the two intervals stay separable
        x = np.array(D[k]["mult"], float) * 2 ** (D_DODGE if cond == "diff" else -D_DODGE)
        v = np.array(D[k][cond])
        filled = cond == "diff"
        ax.plot(x, v[:, 0], color=col, lw=(1.7 if is_c else 1.0) if filled else 0.7,
                ls="-" if filled else (0, (2.5, 1.5)), zorder=(6 if filled else 2) + z0)
        msz = S.ms(k, (4.8 if is_c else 4.3) if filled else 4.0)
        ax.plot(x, v[:, 0], ls="none", marker=mk, ms=msz, mfc=col if filled else "white",
                mec="white" if filled else col, mew=0.5 if filled else 0.8, zorder=(7 if filled else 5) + z0)
        # 95% intervals: only the part beyond the marker edge; omitted when narrower than the marker
        r_dec = (msz / 2) * 25.4 / 72 / mm_per_dec          # marker half-height in decades
        for xi, (g, lo, hi, _) in zip(x, D[k][cond]):
            up, dn = g * 10 ** r_dec, g / 10 ** r_dec
            if hi > up:
                ax.plot([xi, xi], [up, hi], color=col, lw=S.ERRBAR["elinewidth"], zorder=(7.5 if filled else 5.5) + z0,
                        solid_capstyle="butt")
                drawn_bars.append((k, cond, xi, "upper"))
            if lo < dn:
                ax.plot([xi, xi], [lo, dn], color=col, lw=S.ERRBAR["elinewidth"], zorder=(7.5 if filled else 5.5) + z0,
                        solid_capstyle="butt")
                drawn_bars.append((k, cond, xi, "lower"))
# direct value labels at 16n
xr = 16 * 1.4
ax.text(xr, D["separate"]["diff"][-1][0], fpos(D["separate"]["diff"][-1][0]), color=TC["separate"],
        fontsize=6, ha="left", va="center")
ax.text(xr, D["separate"]["same"][-1][0], fpos(D["separate"]["same"][-1][0]), color=TC["separate"],
        fontsize=6, ha="left", va="center")
ax.text(xr, D["concord"]["diff"][-1][0], fpos(D["concord"]["diff"][-1][0]), color=BLUE, fontsize=7,
        fontweight="bold", ha="left", va="center")
ax.set_xticks([1, 4, 16])
ax.set_xticklabels([f"{m * N_BASE:,}" for m in (1, 4, 16)])
ax.set_yticks([0.03, 0.1, 0.3, 1, 3])
ax.set_yticklabels(["0.03", "0.1", "0.3", "1", "3"])
ax.minorticks_off()
ax.set_xlabel("Participants in the three groups")
ax.set_ylabel(GAP_LABEL, linespacing=1.1)
d_handles = [
    Line2D([0], [0], color=GREY, lw=1.0, marker="o", ms=3.8, mfc=GREY, mec="white", mew=0.5,
           label="Different composition"),
    Line2D([0], [0], color=GREY, lw=0.7, ls=(0, (2.5, 1.5)), marker="o", ms=3.6, mfc="white", mec=GREY, mew=0.8,
           label="Same composition"),
]
ax.legend(handles=d_handles, loc="center right", bbox_to_anchor=(1.0, 0.745), fontsize=6, handlelength=2.4,
          labelspacing=0.35, borderaxespad=0.1)

paths = S.save(fig, STEM)
print("\n".join(str(p) for p in paths))

# ----------------------------------------------------------------------------- report
print("\nPanel a (different / same composition, n; increase):")
for k in ROW_ORDER:
    print(f"  {M[k]['label']:24s} {A[k]['diff'][0]:.3f} [{A[k]['diff'][1]:.3f}, {A[k]['diff'][2]:.3f}] / "
          f"{A[k]['same'][0]:.3f}; +{A[k]['inc']:.3f}")
print("Panel b r (AUC vs position shift):")
for k in ["separate", "pooled", "likelihood", "concord"]:
    print(f"  {M[k]['label']:24s} r = {B[k]['r']:+.3f}  range [{B[k]['y'].min():+.3f}, {B[k]['y'].max():+.3f}]"
          f"  datasets {B[k]['R']}")
print("Panel c:", {k: C[k] for k in ROW_ORDER})
print("Panel c mean error g1:", {k: round(C[k + '_err'][0], 3) for k in ROW_ORDER})
print("Panel d:")
for k in ["separate", "concord"]:
    print(f"  {M[k]['label']:24s} diff {[round(v[0], 3) for v in D[k]['diff']]} same {[round(v[0], 3) for v in D[k]['same']]}")
print("Panel d interval segments drawn:", drawn_bars)
