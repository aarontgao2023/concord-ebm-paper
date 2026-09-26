"""Supplementary Fig. 3 | Amyloid placed first in two earlier biomarker panels (88 x 62 mm).

One panel, same encoding as main Fig. 5e. x = APOE genotype (e2, e3/e3, e4); y = % of the 500
resamples in which the amyloid event is placed first. ADNI = filled circles (S.COHORT, blue_main;
12-event panel with CSF Abeta42), NACC = filled squares (blue_secondary; 5-event panel with amyloid
PET). Thin error bars are 95% Wilson intervals over the 500 resamples per cohort x genotype.

Rows used from the input file (cohort_panel / arm_or_reference keys):
  ADNI_K12 / concord_min   and   NACC5_PET / min      (both CONCORD orderings)

Input: aggregate ADNI/NACC result in figures/inputs/realdata/ (not distributed with this
repository; see figures/README.md): apoe_hist_amyloid_position.csv
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

S.apply_style()
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.text  # noqa: E402

HIST = S.REAL_INPUTS / "apoe_hist_amyloid_position.csv"
STEM = "SuppFig3_earlier_panels"
S.require_inputs([HIST], "Supplementary Fig. 3")

GENOS = ["e2", "e33", "e4"]
COHORTS = ["adni", "nacc"]
NAME = {"adni": "ADNI", "nacc": "NACC"}

# ---------------------------------------------------------------------------- data
HIST_ROWS = {"adni": ("ADNI_K12", "concord_min"), "nacc": ("NACC5_PET", "min")}
HIST_EXPECTED = {"adni": [0.8, 25.6, 100.0], "nacc": [0.2, 46.8, 99.8]}
EXPECTED_META = {"adni": dict(K=12, event="ABETA"), "nacc": dict(K=5, event="CENTILOID")}
N_EXPECTED = 500

with open(HIST) as fh:
    rows = list(csv.DictReader(fh))
hist = {}
for coh, (cp, ref) in HIST_ROWS.items():
    for g in GENOS:
        sel = [r for r in rows if r["cohort_panel"] == cp and r["arm_or_reference"] == ref and r["genotype"] == g]
        assert len(sel) == 1, (coh, g, len(sel))
        r = sel[0]
        n = int(r["resamples"])
        assert n == N_EXPECTED, (coh, g, n)
        frac = float(r["frac_amyloid_first"])
        k = round(frac * n)
        assert abs(frac * n - k) < 1e-6, (coh, g, frac, n)
        assert int(r["K"]) == EXPECTED_META[coh]["K"] and r["amyloid_event"] == EXPECTED_META[coh]["event"], r
        hist[(coh, g)] = dict(frac=frac, n=n, k=k, K=int(r["K"]), event=r["amyloid_event"], cohort_panel=cp, ref=ref)
    for g, v in zip(GENOS, HIST_EXPECTED[coh]):
        assert abs(100 * hist[(coh, g)]["frac"] - v) < 1e-6, (coh, g, hist[(coh, g)]["frac"], v)

# ---------------------------------------------------------------------------- style (as Fig. 5e-g)
COL = {c: S.COHORT[NAME[c]]["color"] for c in COHORTS}
MRK = {c: S.COHORT[NAME[c]]["marker"] for c in COHORTS}
MS = {c: 4.2 * S.MARKER_SCALE[MRK[c]] for c in COHORTS}          # equal visual area
CUR_LW = 1.2
X = [0, 1, 2]
DODGE = {"adni": -0.07, "nacc": 0.07}      # Wilson intervals overlap at e2 and e4
AMY_MEASURE = {"ABETA": "CSF Aβ42", "CENTILOID": "amyloid PET"}
DESC = {c: f"{AMY_MEASURE[EXPECTED_META[c]['event']]} ({EXPECTED_META[c]['K']} events)" for c in COHORTS}


def pct_label(v):
    """One decimal; an exact 100% (500/500) stays '100%'."""
    return "100%" if abs(v - 100) < 1e-9 else S.fmt_pct(v)


# ---------------------------------------------------------------------------- layout (mm, from top-left)
W, H = 88, 62
fig = plt.figure(figsize=(S.HALF_W, H * S.MM))
AX_TOP, AX_H = 6.2, 48.0
AX_LEFT = 10.8
LABEL_OUT = 23.8       # end-label column beyond the axes' right edge
AX_W = W - AX_LEFT - LABEL_OUT - 1.5
XLIM = (-0.25, 2.2)
YLIM = (-4, 104)
LABEL_X = 2.22

ax = S.add_axes_mm(fig, AX_LEFT, AX_TOP, AX_W, AX_H)
records = []
ends = {}
for coh in ["nacc", "adni"]:                 # ADNI drawn on top
    ks = [(hist[(coh, g)]["k"], hist[(coh, g)]["n"]) for g in GENOS]
    ys = [100 * k / n for k, n in ks]
    cis = [tuple(100 * v for v in S.wilson(k, n)) for k, n in ks]
    xs = [x + DODGE[coh] for x in X]
    yerr = [[y - c[0] for y, c in zip(ys, cis)], [c[1] - y for y, c in zip(ys, cis)]]
    ax.errorbar(xs, ys, yerr=yerr, fmt="none", ecolor=COL[coh], zorder=4, **S.ERRBAR)
    ax.plot(xs, ys, color=COL[coh], lw=CUR_LW, zorder=5, solid_capstyle="round")
    ax.plot(xs, ys, ls="none", marker=MRK[coh], ms=MS[coh], mfc=COL[coh], mec="white", mew=0.5, zorder=6)
    ends[coh] = (xs[-1], ys[-1])
    for g, x, y, (k, n), (lo, hi) in zip(GENOS, xs, ys, ks, cis):
        h = hist[(coh, g)]
        records.append(dict(
            figure="Supplementary Fig. 3", measure="amyloid event placed first", cohort=NAME[coh],
            biomarker_panel=DESC[coh], inventory_key=f"{h['cohort_panel']}/{h['ref']}",
            amyloid_measure=AMY_MEASURE[h["event"]], n_events=h["K"], genotype=S.GENOTYPE_LABEL[g],
            x_plotted=round(x, 3), n_resamples=n, k_amyloid_first=k, percent=round(y, 2),
            wilson95_lo=round(lo, 2), wilson95_hi=round(hi, 2), source=S.source_name(HIST)))

# direct end labels: cohort (bold) + e4 value, biomarker panel on the second line.
# The two e4 values are 100% and 99.8%, so NACC's label sits below ADNI's with a leader line.
LY = {"adni": ends["adni"][1], "nacc": ends["adni"][1] - 19.0}
for coh in COHORTS:
    y = LY[coh]
    t = ax.text(LABEL_X, y, NAME[coh], color=COL[coh], fontsize=6.5, fontweight="bold", va="center",
                ha="left", clip_on=False)
    ax.annotate(" " + pct_label(ends[coh][1]), xy=(1, 0.5), xycoords=t, va="center", ha="left",
                fontsize=6, color=COL[coh], annotation_clip=False)
    ax.annotate(DESC[coh], xy=(LABEL_X, y), xycoords="data", xytext=(0, -7.5), textcoords="offset points",
                va="center", ha="left", fontsize=6, color=COL[coh], annotation_clip=False)
x0, y0 = ends["nacc"]
ax.plot([x0 + 0.04, LABEL_X - 0.03], [y0 - 2.5, LY["nacc"] + 1.0], color=COL["nacc"], lw=0.5,
        solid_capstyle="butt", zorder=4, clip_on=False)

ax.set_xlim(*XLIM)
ax.set_ylim(*YLIM)
ax.spines["left"].set_bounds(0, 100)
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xticks(X)
ax.set_xticklabels([S.GENOTYPE_LABEL[g] for g in GENOS])
ax.tick_params(axis="x", length=0, pad=2.5)
ax.set_xlabel(r"$\mathit{APOE}$ genotype", labelpad=2)
ax.set_ylabel("% of bootstrap refits")

S.panel_title_mm(fig, "Amyloid placed first in two other biomarker panels", x=1.6, y_top=3.5)

# ---------------------------------------------------------------------------- outputs
S.SRC.mkdir(parents=True, exist_ok=True)
with open(S.SRC / f"{STEM}.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
    w.writeheader()
    w.writerows(records)

fig.canvas.draw()
_r = fig.canvas.get_renderer()
_bbs = [t.get_window_extent(_r) for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
_px = fig.dpi / 25.4
print("text extent (mm): left %.1f, right %.1f, top %.1f, bottom %.1f" % (
    min(b.x0 for b in _bbs) / _px, max(b.x1 for b in _bbs) / _px,
    H - max(b.y1 for b in _bbs) / _px, H - min(b.y0 for b in _bbs) / _px))

paths = S.save(fig, STEM)
print("\n".join(str(x) for x in paths))
for coh in COHORTS:
    print(NAME[coh], [f"{100 * hist[(coh, g)]['frac']:.1f}%" for g in GENOS],
          [f"{hist[(coh, g)]['k']}/{hist[(coh, g)]['n']}" for g in GENOS])
