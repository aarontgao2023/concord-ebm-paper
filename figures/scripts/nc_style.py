"""Shared figure style for the CONCORD Nature Communications figures.

House style follows figures4papers (ChenLiu-1996): minimalist axes (top/right spines off),
frameless legends, one strong blue for the proposed method, softer tones for comparators,
neutral greys for secondary elements, light-blue call-out boxes for the method.
Sizes follow Nature Communications print rules: 180 mm double-column width, Arial,
5-7 pt text, 8 pt bold lowercase panel letters, lines >= 0.25 pt, editable vector text.

Paths (all relative to this repository unless overridden by an environment variable):
  inputs       figures/inputs/                CONCORD_FIG_INPUTS
                 simulation/  simulation-only aggregate files (in the repository)
                 realdata/    ADNI/NACC aggregate files (not in the repository; see figures/README.md)
  figures      figures/output/                CONCORD_OUTPUT_DIR
  source data  figures/source_data/           CONCORD_SOURCE_DATA_DIR
  Arial        system font folders            CONCORD_FONT_DIR (folder with Arial.ttf, Arial Bold.ttf,
                                                                Arial Italic.ttf)
"""
from __future__ import annotations

import math
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

FIGDIR = Path(__file__).resolve().parents[1]      # figures/
ROOT = FIGDIR.parent                               # repository root


def _env_dir(name, default):
    v = os.environ.get(name)
    return Path(v).expanduser() if v else default


INPUTS = _env_dir("CONCORD_FIG_INPUTS", FIGDIR / "inputs")
SIM_INPUTS = INPUTS / "simulation"
REAL_INPUTS = INPUTS / "realdata"
OUT = _env_dir("CONCORD_OUTPUT_DIR", FIGDIR / "output")
SRC = _env_dir("CONCORD_SOURCE_DATA_DIR", FIGDIR / "source_data")


def require_inputs(paths, item):
    """Stop with a clear message if any input file is missing (real-data aggregates are not in the
    repository; they are placed in figures/inputs/realdata/ by the user)."""
    missing = [Path(p) for p in paths if not Path(p).is_file()]
    if not missing:
        return
    lines = [f"{item}: {len(missing)} input file(s) not found:"]
    lines += [f"  - {p}" for p in missing]
    if any(p.parent == REAL_INPUTS for p in missing):
        lines += [
            "These are aggregate ADNI/NACC results that are not distributed with this repository.",
            f"Place them in {REAL_INPUTS}",
            "(or set CONCORD_FIG_INPUTS to a folder that contains simulation/ and realdata/).",
            "See figures/README.md for the list of files and how they are produced.",
        ]
    sys.exit("\n".join(lines))


def source_name(p):
    """Input path as recorded in the source-data files: relative to the inputs folder."""
    p = Path(p)
    try:
        return "inputs/" + p.relative_to(INPUTS).as_posix()
    except ValueError:
        return p.name

MM = 1 / 25.4
FULL_W = 180 * MM
HALF_W = 88 * MM

# figures4papers palette
PAL = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "blue_light": "#7FA6D6",
    "blue_bg": "#E6EEF8",
    "green_1": "#DDF3DE",
    "green_2": "#AADCA9",
    "green_3": "#8BCF8B",
    "red_1": "#F6CFCB",
    "red_2": "#E9A6A1",
    "red_strong": "#B64342",
    "neutral": "#CFCECE",
    "grey_mid": "#767676",
    "grey_dark": "#4D4D4D",
    "ink": "#272727",
    "highlight": "#FFD700",
    "teal": "#42949E",
    "violet": "#9A4D8E",
}

# Method encoding, used identically in every figure (colour + marker).
METHOD = {
    "concord": dict(label="CONCORD", color=PAL["blue_main"], marker="o"),
    "pooled": dict(label="Pooled-score DEBM", color=PAL["blue_light"], marker="s"),
    "separate": dict(label="Separately fitted DEBM", color="#C8574F", marker="^"),
    "likelihood": dict(label="Likelihood EBM", color="#E0A33B", marker="D"),
    "saebm": dict(label="Stage-aware EBM", color="#B8619F", marker="v"),
}
METHOD_ORDER = ["likelihood", "saebm", "separate", "pooled", "concord"]

# Diagnosis: figures4papers green ramp (light CN -> deeper AD); no method or event uses green,
# and the ramp stays ordered in greyscale by lightness.
DIAG = {"CN": "#E3F2E1", "MCI": "#A8D5A4", "AD": "#4F9458"}
DIAG_TEXT = {"CN": "#272727", "MCI": "#272727", "AD": "white"}

# Events: CSF in a teal family, cognition in a warm brown family.
EVENT_COLOR = {
    "ABETA": "#1B6F7A",
    "PTAU": "#42949E",
    "TAU": "#8CC4CB",
    "MEM": "#8C4A2F",
    "EXF": "#C07A4E",
    "LAN": "#E2B48C",
}
EVENT_LABEL = {
    "ABETA": "Aβ42",
    "PTAU": "p-tau",
    "TAU": "Total tau",
    "MEM": "Memory",
    "EXF": "Executive",
    "LAN": "Language",
}
EVENT_ORDER = ["ABETA", "PTAU", "TAU", "MEM", "EXF", "LAN"]
GENOTYPE_LABEL = {"e2": "ε2", "e33": "ε3/ε3", "e4": "ε4"}

FONT_SIZE = 6.5


ARIAL_FILES = {"arial.ttf", "arial bold.ttf", "arial italic.ttf",   # macOS names
               "arialbd.ttf", "ariali.ttf"}                         # Windows / msttcorefonts names


def _register_arial():
    """Register Arial Regular, Bold and Italic with matplotlib.

    The fonts are looked up in CONCORD_FONT_DIR if set, otherwise in the system font folders.
    Without Arial, matplotlib falls back to Helvetica or DejaVu Sans; the figures are still drawn,
    but text widths (and so the layout checks) differ from the published figures."""
    font_dir = os.environ.get("CONCORD_FONT_DIR")
    found = [p for p in font_manager.findSystemFonts(fontpaths=[font_dir] if font_dir else None)
             if Path(p).name.lower() in ARIAL_FILES]
    for p in sorted(found):
        font_manager.fontManager.addfont(p)
    if not found:
        print("[warning] Arial not found; set CONCORD_FONT_DIR to a folder with Arial.ttf, "
              "Arial Bold.ttf and Arial Italic.ttf", file=sys.stderr)


def apply_style():
    _register_arial()
    plt.rcParams.update(
        {
            "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": FONT_SIZE,
            "axes.labelsize": 7,
            "axes.titleweight": "bold",
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "legend.frameon": False,
            "legend.handlelength": 1.4,
            "legend.handletextpad": 0.4,
            "legend.borderaxespad": 0.2,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.6,
            "axes.edgecolor": PAL["ink"],
            "axes.labelcolor": PAL["ink"],
            "axes.labelpad": 2.5,
            "xtick.color": PAL["ink"],
            "ytick.color": PAL["ink"],
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.major.pad": 1.5,
            "ytick.major.pad": 1.5,
            "lines.linewidth": 1.0,
            "lines.markersize": 3.5,
            "errorbar.capsize": 0,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.dpi": 600,
            "figure.dpi": 150,
            "mathtext.default": "regular",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.sf": "Arial",
            "mathtext.cal": "Arial",
            "mathtext.tt": "Arial",
            "axes.titlesize": 6.5,
            "axes.titlelocation": "left",
        }
    )


def panel_label(fig, ax, letter, dx=-0.05, dy=0.02, x=None, y=None):
    """Bold lowercase 8 pt panel letter placed at the top-left of an axes (figure coords)."""
    bb = ax.get_position()
    fx = bb.x0 + dx if x is None else x
    fy = bb.y1 + dy if y is None else y
    fig.text(fx, fy, letter, fontsize=8, fontweight="bold", va="bottom", ha="left", color="black")


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return c - h, c + h


def qa_text(fig, min_pt=5.0, max_pt=7.0, letter_pt=8.0, margin_mm=1.5, verbose=True):
    """Check every visible text object: size within [min_pt, max_pt] (panel letters 8 pt),
    inside the canvas by margin_mm, and report overlapping text boxes. Returns list of problems."""
    from matplotlib.text import Text
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    px_per_mm = fig.dpi / 25.4
    probs, boxes = [], []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        fs = t.get_fontsize()
        if not (min_pt - 1e-6 <= fs <= max_pt + 1e-6 or abs(fs - letter_pt) < 1e-6):
            probs.append(f"size {fs:.2f} pt: {t.get_text()[:40]!r}")
        if "_{" in t.get_text() or "^{" in t.get_text():
            if fs * 0.7 < min_pt - 1e-6:
                probs.append(f"mathtext sub/superscript {fs*0.7:.2f} pt: {t.get_text()[:40]!r}")
        bb = t.get_window_extent(r)
        m = (0.8 if abs(fs - letter_pt) < 1e-6 else margin_mm) * px_per_mm
        if bb.x0 < m - 0.5 or bb.y0 < m - 0.5 or bb.x1 > W - m + 0.5 or bb.y1 > H - m + 0.5:
            probs.append(f"within {margin_mm} mm of edge: {t.get_text()[:40]!r}")
        boxes.append((bb, t.get_text()[:30]))
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][0], boxes[j][0]
            ox = min(a.x1, b.x1) - max(a.x0, b.x0)
            oy = min(a.y1, b.y1) - max(a.y0, b.y0)
            if ox > 1.0 and oy > 1.0:
                probs.append(f"text overlap: {boxes[i][1]!r} / {boxes[j][1]!r}")
    if verbose:
        w_mm, h_mm = fig.get_size_inches() * 25.4
        print(f"[QA] canvas {w_mm:.1f} x {h_mm:.1f} mm; {len(probs)} text issue(s)")
        for p in probs:
            print("   -", p)
    return probs


def save(fig, stem, qa=True):
    if qa:
        qa_text(fig)
    OUT.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in ("pdf", "png", "svg"):
        p = OUT / f"{stem}.{ext}"
        fig.savefig(p, dpi=600 if ext == "png" else None, facecolor="white")
        paths.append(p)
    plt.close(fig)
    return paths


def method_handles(keys, ms=3.5):
    from matplotlib.lines import Line2D

    return [
        Line2D([0], [0], color=METHOD[k]["color"], marker=METHOD[k]["marker"], ms=ms, lw=1.0,
               label=METHOD[k]["label"])
        for k in keys
    ]


# ---------------------------------------------------------------- shared helpers
# Darker text versions of method colours (light colours are hard to read as 6 pt text).
TEXT_COLOR = {
    "concord": PAL["blue_main"],
    "pooled": "#4F6F97",
    "separate": "#B0433C",
    "likelihood": "#A8741C",
    "saebm": "#9A4585",
}
# Tinted fills for comparator bars; CONCORD stays solid (figures4papers: only the method is saturated).
FILL = {
    "concord": PAL["blue_main"],
    "pooled": "#C9D9EE",
    "separate": "#EDB3AE",
    "likelihood": "#F2D29B",
    "saebm": "#E6C0DA",
}
EVENT_TEXT_COLOR = dict(EVENT_COLOR, TAU="#3C8189", LAN="#9A6A3F", EXF="#A5633A")  # >= ~4:1 contrast at 6 pt
EVENT_LINE_COLOR = dict(EVENT_COLOR, TAU="#6FB0B8", LAN="#D39B69")
MARKER_SCALE = {"o": 1.0, "s": 0.9, "^": 1.12, "v": 1.12, "D": 0.82}
REF_LINE = dict(color="#9A9A9A", lw=0.6, ls=(0, (3, 2)), zorder=0.5)
ERRBAR = dict(capsize=0, elinewidth=0.6)
GROUP = {  # schematic groups (Fig. 1 only): ochre circle / slate-grey square
    1: dict(label="Group 1", color="#D9A441", marker="o"),
    2: dict(label="Group 2", color="#5B6B7A", marker="s"),
}
COHORT = {
    "ADNI": dict(color=PAL["blue_main"], marker="o"),
    "NACC": dict(color=PAL["blue_secondary"], marker="s"),
}
LABEL = {
    "gap": "Mean-position gap (positions per event)",
    "shift": "Position shift, CN-heavy − AD-heavy (positions)",
    "within": "Within-diagnosis permutation",
    "unrestricted": "Unrestricted permutation",
    "distance": "Discordant event pairs",
    "interval": "95% interval",
}


def ms(key, base=3.5):
    """Marker size for a method, equalising visual area across marker shapes."""
    return base * MARKER_SCALE.get(METHOD[key]["marker"], 1.0)


def fmt_pct(x):
    return f"{x:.1f}%"


def fmt_pos(x, signed=False):
    s = f"{x:+.3f}" if abs(x) < 0.1 else f"{x:+.2f}"
    return s if signed else s.lstrip("+")


def _mm_to_fig(fig, x_mm, y_top_mm):
    w_mm, h_mm = fig.get_size_inches() * 25.4
    return x_mm / w_mm, 1 - y_top_mm / h_mm


def add_axes_mm(fig, x, y_top, w, h, **kw):
    """Add axes by position in mm from the figure's top-left corner."""
    w_mm, h_mm = fig.get_size_inches() * 25.4
    return fig.add_axes([x / w_mm, 1 - (y_top + h) / h_mm, w / w_mm, h / h_mm], **kw)


def panel_label_mm(fig, letter, x=1.0, y_top=3.5):
    """8 pt bold panel letter with its baseline-top at (x, y_top) mm from the top-left."""
    fx, fy = _mm_to_fig(fig, x, y_top)
    return fig.text(fx, fy, letter, fontsize=8, fontweight="bold", ha="left", va="baseline")


def panel_title_mm(fig, text, x, y_top=3.5, color=None):
    """6.5 pt bold panel title on the same baseline as the panel letter."""
    fx, fy = _mm_to_fig(fig, x, y_top)
    return fig.text(fx, fy, text, fontsize=6.5, fontweight="bold", ha="left", va="baseline",
                    color=color or PAL["ink"])


def method_legend(fig, keys, y_top=3.0, bars=False, x_center=0.5, ncol=None, fontsize=6):
    """One-row method legend (top-centre). CONCORD label bold blue. bars=True adds a tinted patch."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles = []
    for k in keys:
        m = METHOD[k]
        if bars:
            handles.append((Patch(facecolor=FILL[k], edgecolor=m["color"], lw=0.6),
                            Line2D([0], [0], color=m["color"], marker=m["marker"], ms=ms(k), lw=0)))
        else:
            handles.append(Line2D([0], [0], color=m["color"], marker=m["marker"], ms=ms(k),
                                  lw=1.4 if k == "concord" else 1.0, mec="white", mew=0.4))
    labels = [METHOD[k]["label"] for k in keys]
    fx, fy = _mm_to_fig(fig, 0, y_top)
    from matplotlib.legend_handler import HandlerTuple

    leg = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(x_center, fy),
                     ncol=ncol or len(keys), frameon=False, fontsize=fontsize, columnspacing=1.6,
                     handlelength=1.8, handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
    for txt, k in zip(leg.get_texts(), keys):
        if k == "concord":
            txt.set_color(PAL["blue_main"])
            txt.set_fontweight("bold")
    return leg


def concord_band(ax, y, h=0.8, orientation="h", **kw):
    """Light-blue call-out band behind a CONCORD row (house style)."""
    style = dict(color=PAL["blue_bg"], zorder=0, lw=0)
    style.update(kw)
    if orientation == "h":
        return ax.axhspan(y - h / 2, y + h / 2, **style)
    return ax.axvspan(y - h / 2, y + h / 2, **style)


def _div_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "concord_div", ["#8C4A2F", "#D9A57F", "#FFFFFF", "#7FA6D6", PAL["blue_main"]])


def _seq_cmap():
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "concord_seq", ["#E9EFF7", "#7FA6D6", PAL["blue_secondary"], PAL["blue_main"]])


CMAP_DIV = _div_cmap()
CMAP_SEQ = _seq_cmap()
