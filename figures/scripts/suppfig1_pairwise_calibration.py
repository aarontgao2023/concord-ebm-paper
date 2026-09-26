"""Supplementary Fig. 1 | Calibration of pairwise within-diagnosis permutation (180 x 62 mm).

The pairwise test behind the APOE comparisons (Table 3): for each pair of groups only the two
compared groups are permuted within diagnosis; the third group stays in the pooled model.
Simulated groups: CN-heavy group (g1, 75 participants), intermediate group (g2, 411),
AD-heavy group (g3, 485). Legacy labels in the stored files: e2 = g1, e33 = g2, e4 = g3.

a  No group differs: CONCORD per-contrast rejection rate for the three contrasts
   (R = 1,000 simulated datasets), Wilson 95% intervals.
b  One group changed (27 of 91 pairs reversed): rejection rate of the contrast between the two
   unchanged groups (R = 500), Wilson 95% intervals.
Dashed line: Bonferroni per-contrast level 0.05 / 3 = 1/60 = 1.67%.

Inputs (figures/inputs/simulation/):
  rates_exact_final.txt    a: block "== method_pair_complete_a" (output of
                           simulation/scripts/dev/rates_exact_all.py); b is cross-checked against
                           block "== method_pair_partial_a"
  Fig3f_partial_null.csv   b: CONCORD rejections of the contrast between the two unchanged groups
                           when one group is changed (run method_pair_partial_a)
Model label invariant_min = CONCORD.

All plotted numbers are read from stored aggregate files; nothing is refitted.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import nc_style as S  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

S.apply_style()

RATES_TXT = S.SIM_INPUTS / "rates_exact_final.txt"
PARTIAL = S.SIM_INPUTS / "Fig3f_partial_null.csv"
S.require_inputs([RATES_TXT, PARTIAL], "Supplementary Fig. 1")

W, H = 180.0, 62.0
INK = S.PAL["ink"]
GREY = S.PAL["grey_mid"]
COL = S.METHOD["concord"]["color"]
LEVEL = 1 / 60  # Bonferroni per-contrast level, 0.05 / 3

LEGACY = {"e2": 1, "e33": 2, "e4": 3}                     # legacy genotype labels -> simulated group
GNAME = {1: "CN-heavy", 2: "intermediate", 3: "AD-heavy"}
GSIZE = {1: 75, 2: 411, 3: 485}
CONTRASTS = [(1, 2), (1, 3), (2, 3)]                      # row order, top to bottom (both panels)
CHANGED = {(2, 3): 1, (1, 3): 2, (1, 2): 3}               # b: contrast tested -> group whose ordering changed


def close(a, b, tol, what):
    assert abs(a - b) <= tol, f"{what}: got {a}, expected {b}"


def contrast_label(c):
    a, b = GNAME[c[0]], GNAME[c[1]]
    return f"{a[0].upper() + a[1:]} vs {b}"


# ----------------------------------------------------------------------------- a: parse text file
lines = RATES_TXT.read_text().splitlines()
start = lines.index("== method_pair_complete_a")
block = []
for ln in lines[start + 1:]:
    if ln.startswith("=="):
        break
    block.append(ln)
ROW_RE = re.compile(r"^(\S+)\s+(\S+)\s+diagnosis_pair_(e2|e33|e4)_(e2|e33|e4)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
                    r"([\d.]+)\s+\[([\d.]+),\s*([\d.]+)\]")
A = {}
for ln in block:
    m = ROW_RE.match(ln)
    if not m:
        continue
    cell, engine, l1, l2, rows, det, rej, rate, lo_s, hi_s = m.groups()
    if cell != "REF_H0" or engine != "invariant_min":
        continue
    c = (LEGACY[l1], LEGACY[l2])
    k, n = int(rej), int(rows)
    assert int(det) == n, f"a {c}: {det} of {n} datasets determinate"
    lo, hi = S.wilson(k, n)
    close(k / n, float(rate), 5e-4, f"a rate {c}")
    close(lo, float(lo_s), 5e-4, f"a Wilson low {c}")
    close(hi, float(hi_s), 5e-4, f"a Wilson high {c}")
    A[c] = dict(k=k, n=n, rate=k / n, lo=lo, hi=hi, legacy=f"REF_H0 diagnosis_pair_{l1}_{l2}")
assert set(A) == set(CONTRASTS), f"a: contrasts found {sorted(A)}"
for c, k in zip(CONTRASTS, (13, 18, 20)):
    assert (A[c]["k"], A[c]["n"]) == (k, 1000), f"a {c}: {A[c]['k']}/{A[c]['n']}, expected {k}/1000"

# ----------------------------------------------------------------------------- b: partial null csv
with open(PARTIAL, newline="") as fh:
    rows_b = [r for r in csv.DictReader(fh) if r["method"] == "CONCORD"]
B = {}
for r in rows_b:
    g = int(r["displaced_group"])
    c = tuple(int(x) for x in r["unchanged_pair"].split("–"))
    assert CHANGED[c] == g, f"b: changed group {g} does not match unchanged contrast {c}"
    k, n = int(r["reject"]), int(r["R"])
    lo, hi = S.wilson(k, n)
    close(k / n, float(r["rate"]), 1e-9, f"b rate {c}")
    close(lo, float(r["low"]), 1e-6, f"b Wilson low {c}")
    close(hi, float(r["high"]), 1e-6, f"b Wilson high {c}")
    B[c] = dict(k=k, n=n, rate=k / n, lo=lo, hi=hi, changed=g)
assert set(B) == set(CONTRASTS), f"b: contrasts found {sorted(B)}"
for g, k in ((1, 2), (2, 5), (3, 12)):
    c = [cc for cc, gg in CHANGED.items() if gg == g][0]
    assert (B[c]["k"], B[c]["n"]) == (k, 500), f"b changed {g}: {B[c]['k']}/{B[c]['n']}, expected {k}/500"

# cross-check b against the PWR_*_K27 rows of the text file (method_pair_partial_a block)
start = lines.index("== method_pair_partial_a")
for ln in lines[start + 1:]:
    if ln.startswith("=="):
        break
    m = ROW_RE.match(ln)
    if not m or m.group(2) != "invariant_min":
        continue
    cell, l1, l2, rej = m.group(1), m.group(3), m.group(4), int(m.group(7))
    c = (LEGACY[l1], LEGACY[l2])
    changed = {"PWR_E2_K27": 1, "PWR_E33_K27": 2, "PWR_E4_K27": 3}.get(cell)
    if changed is not None and CHANGED[c] == changed:
        assert rej == B[c]["k"], f"b cross-check {cell} {c}: {rej} vs {B[c]['k']}"

# ----------------------------------------------------------------------------- source data
S.SRC.mkdir(parents=True, exist_ok=True)
with open(S.SRC / "SuppFig1_pairwise_calibration.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["panel", "scenario", "contrast_tested", "changed_group", "group_sizes", "model", "rejections",
                "datasets", "rejection_rate", "wilson_low", "wilson_high", "per_contrast_level", "source",
                "source_row"])
    for c in CONTRASTS:
        v = A[c]
        w.writerow(["a", "No group differs", contrast_label(c), "none",
                    f"{GSIZE[c[0]]} vs {GSIZE[c[1]]}", "CONCORD", v["k"], v["n"], f"{v['rate']:.4f}",
                    f"{v['lo']:.4f}", f"{v['hi']:.4f}", f"{LEVEL:.4f}", S.source_name(RATES_TXT),
                    f"method_pair_complete_a {v['legacy']}"])
    for c in CONTRASTS:
        v = B[c]
        w.writerow(["b", "One group changed (27 of 91 pairs reversed)", contrast_label(c),
                    f"{GNAME[v['changed']]} group", f"{GSIZE[c[0]]} vs {GSIZE[c[1]]}", "CONCORD", v["k"],
                    v["n"], f"{v['rate']:.4f}", f"{v['lo']:.4f}", f"{v['hi']:.4f}", f"{LEVEL:.4f}",
                    S.source_name(PARTIAL), f"displaced_group {v['changed']}, CONCORD"])

# ----------------------------------------------------------------------------- figure
fig = plt.figure(figsize=(W * S.MM, H * S.MM))
ROW = 4.5                       # letter / title baseline (mm from top)
AX_TOP, AX_H = 7.5, 44.5
AX_W = 50.0
PANELS = {"a": (36.0, "No group differs", A), "b": (126.0, "One group changed (27 of 91 pairs reversed)", B)}
LETTER_X = {"a": 1.0, "b": 91.0}
TITLE_DX = 4.5
Y = {c: 2 - i for i, c in enumerate(CONTRASTS)}
YLIM = (-0.55, 2.95)
XLIM = (0, 4.9)
MS = S.ms("concord", 4.6)

for L, (x0, title, data) in PANELS.items():
    ax = S.add_axes_mm(fig, x0, AX_TOP, AX_W, AX_H)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.spines["bottom"].set_bounds(0, 4)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel("Rejection rate (%)")
    ax.spines["left"].set_visible(False)
    ax.set_yticks([Y[c] for c in CONTRASTS])
    if L == "a":
        ax.set_yticklabels([contrast_label(c) for c in CONTRASTS])
    else:
        ax.set_yticklabels([f"{contrast_label(c)}\n({GNAME[data[c]['changed']]} group changed)"
                            for c in CONTRASTS])
    ax.tick_params(axis="y", length=0, pad=3)
    for y in (0, 1, 2):  # light row bands (as Fig. 4e)
        ax.axhspan(y - 0.36, y + 0.36, color="#F2F2F2", lw=0, zorder=0)
    ytop = 2.52
    ax.plot([100 * LEVEL] * 2, [YLIM[0], ytop], **S.REF_LINE)
    ax.text(100 * LEVEL, ytop + 0.04, f"Per-comparison level {100 * LEVEL:.2f}%", color=GREY, fontsize=6,
            ha="center", va="bottom")
    for c in CONTRASTS:
        v, y = data[c], Y[c]
        ax.plot([100 * v["lo"], 100 * v["hi"]], [y, y], color=COL, lw=0.9, solid_capstyle="butt", zorder=3)
        ax.plot(100 * v["rate"], y, marker="o", ms=MS, color=COL, mec=COL, ls="none", zorder=4)
        ax.text(100 * v["hi"] + 0.1, y, S.fmt_pct(100 * v["rate"]), color=S.TEXT_COLOR["concord"], fontsize=6,
                fontweight="bold", ha="left", va="center", zorder=5,
                bbox=dict(boxstyle="square,pad=0.1", fc="#F2F2F2", ec="none"))
    S.panel_label_mm(fig, L, x=LETTER_X[L], y_top=ROW)
    S.panel_title_mm(fig, title, x=LETTER_X[L] + TITLE_DX, y_top=ROW)

paths = S.save(fig, "SuppFig1_pairwise_calibration")
print("\n".join(str(p) for p in paths))
for L, (_, _, data) in PANELS.items():
    for c in CONTRASTS:
        v = data[c]
        print(f"{L} {contrast_label(c):28s} {v['k']}/{v['n']} = {100 * v['rate']:.1f}% "
              f"[{100 * v['lo']:.1f}, {100 * v['hi']:.1f}]")
