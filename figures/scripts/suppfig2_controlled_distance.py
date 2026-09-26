"""Supplementary Figure 2 | Distance between group orderings in the controlled composition experiment.

180 x 70 mm. Two panels sharing rows and x-scale: a ADNI, b NACC. For each ordering estimator,
the mean normalised Kendall distance between the two groups' fitted orderings, shown as
discordant event pairs (%) = 100 x distance, for the matched (hollow) and different (filled)
composition splits, with paired-replicate bootstrap 95% intervals. The right-hand column gives
the paired change (different - matched) in percentage points with its 95% interval.

Inputs: aggregate ADNI/NACC results in figures/inputs/realdata/ (not distributed with this
repository; see figures/README.md): composition_distance_summary.csv and
composition_paired_contrast_summary.csv. Model labels in the input files: likelihood_ebm =
likelihood EBM, separate = separately fitted DEBM, shared = pooled-score DEBM, concord = CONCORD.

Every plotted number is read from stored aggregate files (no participant-level data, no refits).
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

STEM = "SuppFig2_controlled_distance"
F_DIST = S.REAL_INPUTS / "composition_distance_summary.csv"
F_CON = S.REAL_INPUTS / "composition_paired_contrast_summary.csv"
S.require_inputs([F_DIST, F_CON], "Supplementary Fig. 2")

ENGINE = {"likelihood_ebm": "likelihood", "separate": "separate", "shared": "pooled", "concord": "concord"}
MODELS = ["likelihood", "separate", "pooled", "concord"]
COHORTS = [("adni", "ADNI"), ("nacc", "NACC")]
M = S.METHOD
INK = S.PAL["ink"]
GREY = S.PAL["grey_mid"]
TOL = 0.0006


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def check(name, got, want, tol=TOL):
    assert abs(got - want) <= tol, f"{name}: got {got:.4f}, expected {want:.4f} (tol {tol})"


def pp(x, signed=False):
    """Percentage points, one decimal, half-up, typographic minus."""
    v = Decimal(str(round(100 * x, 9))).quantize(Decimal("0.1"), rounding=ROUND_HALF_UP)
    s = f"{v:+f}" if signed else f"{v:f}"
    return s.replace("-", "−")


# ----------------------------------------------------------------------------------------
# Data
# ----------------------------------------------------------------------------------------
dist = {}
for r in read_csv(F_DIST):
    if r["kind"] != "absolute_distance":
        continue
    assert int(r["complete_replicates"]) == int(r["planned_replicates"]) == int(r["n"]) == 200
    dist[(r["cohort"], ENGINE[r["engine_left"]], r["design"])] = dict(
        mean=float(r["mean"]), lo=float(r["mean_ci_low"]), hi=float(r["mean_ci_high"]), n=int(r["n"]),
        ci=r["ci_method"])
con = {}
for r in read_csv(F_CON):
    if r["kind"] != "different_minus_matched":
        continue
    assert r["engine_left"] == r["engine_right"] and int(r["n"]) == 200
    con[(r["cohort"], ENGINE[r["engine_left"]])] = dict(
        mean=float(r["mean"]), lo=float(r["mean_ci_low"]), hi=float(r["mean_ci_high"]), n=int(r["n"]),
        ci=r["ci_method"])

EXP_DIST = {  # (matched, different)
    ("adni", "likelihood"): (0.346, 0.361), ("adni", "separate"): (0.326, 0.379),
    ("adni", "pooled"): (0.066, 0.128), ("adni", "concord"): (0.068, 0.081),
    ("nacc", "likelihood"): (0.512, 0.538), ("nacc", "separate"): (0.396, 0.517),
    ("nacc", "pooled"): (0.114, 0.249), ("nacc", "concord"): (0.115, 0.137),
}
EXP_CON = {  # different - matched: (mean, low, high)
    ("adni", "likelihood"): (0.015, -0.017, 0.047), ("adni", "separate"): (0.053, 0.024, 0.082),
    ("adni", "pooled"): (0.062, 0.051, 0.074), ("adni", "concord"): (0.013, 0.003, 0.023),
    ("nacc", "likelihood"): (0.026, -0.013, 0.064), ("nacc", "separate"): (0.121, 0.092, 0.149),
    ("nacc", "pooled"): (0.135, 0.119, 0.151), ("nacc", "concord"): (0.023, 0.010, 0.035),
}
for (c, m), (em, ed) in EXP_DIST.items():
    check(f"{c} {m} matched", dist[(c, m, "matched")]["mean"], em)
    check(f"{c} {m} different", dist[(c, m, "different")]["mean"], ed)
    for dsg in ("matched", "different"):
        d = dist[(c, m, dsg)]
        assert d["lo"] <= d["mean"] <= d["hi"]
for (c, m), (e, lo, hi) in EXP_CON.items():
    k = con[(c, m)]
    check(f"{c} {m} change", k["mean"], e)
    check(f"{c} {m} change low", k["lo"], lo)
    check(f"{c} {m} change high", k["hi"], hi)
    # the stored paired change must equal different - matched of the stored means
    check(f"{c} {m} change = different - matched", k["mean"],
          dist[(c, m, "different")]["mean"] - dist[(c, m, "matched")]["mean"], 1e-9)
LABELS = {k: f"{pp(v['mean'], True)} [{pp(v['lo'])}, {pp(v['hi'])}]" for k, v in con.items()}
EXP_LAB = {
    ("adni", "likelihood"): "+1.5 [−1.7, 4.7]", ("adni", "separate"): "+5.3 [2.4, 8.2]",
    ("adni", "pooled"): "+6.2 [5.1, 7.4]", ("adni", "concord"): "+1.3 [0.3, 2.3]",
    ("nacc", "likelihood"): "+2.6 [−1.3, 6.4]", ("nacc", "separate"): "+12.1 [9.2, 14.9]",
    ("nacc", "pooled"): "+13.5 [11.9, 15.1]", ("nacc", "concord"): "+2.3 [1.0, 3.5]",
}
for k, want in EXP_LAB.items():
    assert LABELS[k] == want, (k, LABELS[k], want)

# ----------------------------------------------------------------------------------------
# Source data
# ----------------------------------------------------------------------------------------
S.SRC.mkdir(parents=True, exist_ok=True)
with open(S.SRC / f"{STEM}.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["panel", "cohort", "model", "quantity", "n_draws", "estimate", "ci95_low", "ci95_high", "unit",
                "interval_method", "source_file"])
    for (c, cl), panel in zip(COHORTS, "ab"):
        for m in MODELS:
            for dsg in ("matched", "different"):
                d = dist[(c, m, dsg)]
                w.writerow([panel, cl, M[m]["label"], f"discordant event pairs, {dsg} composition", d["n"],
                            round(100 * d["mean"], 6), round(100 * d["lo"], 6), round(100 * d["hi"], 6), "%",
                            d["ci"], S.source_name(F_DIST)])
            k = con[(c, m)]
            w.writerow([panel, cl, M[m]["label"], "paired change, different minus matched composition", k["n"],
                        round(100 * k["mean"], 6), round(100 * k["lo"], 6), round(100 * k["hi"], 6),
                        "percentage points", k["ci"], S.source_name(F_CON)])

# ----------------------------------------------------------------------------------------
# Figure (positions in mm from the top-left corner)
# ----------------------------------------------------------------------------------------
FW, FH = 180.0, 70.0
fig = plt.figure(figsize=(FW * S.MM, FH * S.MM))

S.method_legend(fig, MODELS, y_top=1.0)
# composition key: hollow = matched (reference), filled = different
ckey = [Line2D([0], [0], ls="none", marker="o", ms=3.8, mfc="white", mec=GREY, mew=0.8),
        Line2D([0], [0], ls="none", marker="o", ms=3.8, mfc=GREY, mec="white", mew=0.4)]
_, ky = (0, 1 - 5.4 / FH)
fig.legend(ckey, ["Matched composition", "Different composition"], loc="upper center",
           bbox_to_anchor=(0.5, ky), ncol=2, fontsize=6, frameon=False, handlelength=1.0,
           handletextpad=0.4, columnspacing=2.0, borderaxespad=0, borderpad=0)

L1 = 13.0         # letter / title baseline
TITLE_DX = 4.5
AX_TOP, AX_H, AX_W = 15.0, 44.5, 50.0
PANELS = [  # (cohort key, cohort label, letter x, axes x, column x, band x0)
    ("adni", "ADNI", 1.0, 27.0, 78.5, 1.0),
    ("nacc", "NACC", 99.0, 104.0, 155.5, 104.0),
]
COL_W = 18.5
ROW = {m: i for i, m in enumerate(MODELS)}
DODGE = 0.14      # matched above, different below the row centre
Y_LIM = (3.55, -0.95)
X_MAX = 60

for c, cl, lx, ax_x, col_x, band_x0 in PANELS:
    ax = S.add_axes_mm(fig, ax_x, AX_TOP, AX_W, AX_H)
    S.panel_label_mm(fig, "a" if c == "adni" else "b", x=lx, y_top=L1)
    S.panel_title_mm(fig, cl, x=lx + TITLE_DX, y_top=L1)
    ax.set_ylim(*Y_LIM)
    ax.set_xlim(0, X_MAX)
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60])
    ax.set_xticklabels(["0", "", "20", "", "40", "", "60"])
    ax.set_xlabel("Discordant event pairs between groups (%)")
    ax.tick_params(axis="y", length=0, pad=4.5)
    ax.spines["left"].set_visible(False)
    to_ax = lambda x_mm: (x_mm - ax_x) / AX_W  # noqa: E731
    S.concord_band(ax, ROW["concord"], 0.84, xmin=to_ax(band_x0), xmax=to_ax(col_x + COL_W), clip_on=False)
    for m in MODELS:
        y = ROW[m]
        col, mk = M[m]["color"], M[m]["marker"]
        is_c = m == "concord"
        msz = S.ms(m, 4.3 if is_c else 3.8)
        z = 10 if is_c else 3
        lw = 0.9 if is_c else 0.6
        d0, d1 = dist[(c, m, "matched")], dist[(c, m, "different")]
        y0, y1 = y - DODGE, y + DODGE
        # thin connector between the two point estimates
        ax.plot([100 * d0["mean"], 100 * d1["mean"]], [y0, y1], color=col, alpha=0.45, lw=0.6, zorder=z,
                solid_capstyle="butt")
        for dd, yy in ((d0, y0), (d1, y1)):
            ax.plot([100 * dd["lo"], 100 * dd["hi"]], [yy, yy], color=col, lw=lw, zorder=z + 1,
                    solid_capstyle="butt")
        ax.plot(100 * d1["mean"], y1, marker=mk, ms=msz, mfc=col, mec="white", mew=0.4, ls="none", zorder=z + 2)
        ax.plot(100 * d0["mean"], y0, marker=mk, ms=msz, mfc="white", mec=col, mew=0.9 if is_c else 0.8,
                ls="none", zorder=z + 3)
        ax.text(to_ax(col_x), y, LABELS[(c, m)], transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=6, color=col if is_c else S.TEXT_COLOR[m], fontweight="bold" if is_c else "normal")
    ax.text(to_ax(col_x), -0.62, "Different − matched\n(percentage points)", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=6, color=GREY, linespacing=1.15)
    ax.set_yticks([ROW[m] for m in MODELS])
    if c == "adni":
        ax.set_yticklabels([M[m]["label"] for m in MODELS])
        for tl, m in zip(ax.get_yticklabels(), MODELS):
            if m == "concord":
                tl.set_color(M["concord"]["color"])
                tl.set_fontweight("bold")
    else:
        ax.set_yticklabels([])

paths = S.save(fig, STEM)
print("\n".join(str(p) for p in paths))
for (c, m) in EXP_LAB:
    d0, d1 = dist[(c, m, "matched")], dist[(c, m, "different")]
    print(f"{c:4s} {m:10s} matched {100*d0['mean']:5.1f}% [{100*d0['lo']:.1f}, {100*d0['hi']:.1f}]  "
          f"different {100*d1['mean']:5.1f}% [{100*d1['lo']:.1f}, {100*d1['hi']:.1f}]  change {LABELS[(c, m)]}")
