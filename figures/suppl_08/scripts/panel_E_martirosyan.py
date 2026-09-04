"""
Supplemental Figure 8, panel E: independent replication of the subtype ordering in a second published substantia nigra cohort.

Input
    data/human/martirosyan_replication.csv
Output
    figures/suppl_08/output/panel_E_martirosyan.png
    figures/suppl_08/output/panel_E_martirosyan_transparent.png
    figures/suppl_08/output/panel_E_martirosyan.svg
    figures/suppl_08/output/panel_E_checksums.csv
Run
    python3 figures/suppl_08/scripts/panel_E_martirosyan.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, scipy

Six subpopulations. The statistic is the Spearman rank correlation over those six, recomputed at draw time and asserted against the recorded value; no subpopulation reaches q < 0.05 and the script asserts that too.

The input md5 is asserted before drawing, so a changed input fails instead of silently redrawing.

Three files are written: the white PNG is the review copy and the md5 of record, the transparent PNG and the SVG are the placement art. Transparency is verified in the written files.
"""
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedFormatter, FixedLocator
from scipy.stats import spearmanr

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
SRC = REPO_ROOT / "data" / "human" / "martirosyan_replication.csv"
FINAL = FIGURE_DIR / "output"
FINAL.mkdir(exist_ok=True)

# ------------------------------------------------------------------ input guard
INPUT_MD5 = "8a6a84d9c89a279d423f801a92a7d4ab"
assert hashlib.md5(SRC.read_bytes()).hexdigest() == INPUT_MD5, "input CSV changed"

mar = pd.read_csv(SRC)
assert len(mar) == 6, f"expected 6 subpopulations, got {len(mar)}"
assert list(mar["cluster"]) == ["neurons3", "neurons0", "neurons5",
                                "neurons1", "neurons2", "neurons4"], "row order changed"
assert mar["vulnerable"].sum() == 3, "expected 3 vulnerable / 3 resilient"
assert (mar["q"] >= 0.05).all(), "no subpopulation should reach q < 0.05"
assert (mar["lo"] < mar["OR"]).all() and (mar["OR"] < mar["hi"]).all(), "CI does not bracket OR"

rho, pval = spearmanr(mar["PD_pct"], mar["OR"])
assert abs(rho + 0.9428571428571428) < 1e-9, f"rho drifted: {rho}"
assert abs(pval - 0.004804664723032055) < 1e-9, f"p drifted: {pval}"

# -------------------------------------------------------------------- constants
C_VULN, C_RES = "#BA1028", "#88B838"       # lineage palette shared with the other panels
TEXT, GUIDE, FAINT = "#1A1A1A", "#6E6E6E", "#B3B3B3"
LEGEND_BG, LEGEND_EDGE = "#F0F0F0", "#D5D5D5"
EXPECTED_PD_PCT = 28.0                     # Martirosyan overall PD/control ratio

FS_AXLAB, FS_TICK = 8.0, 6.5               # unchanged from panel_d_v03
FS_PT = FS_LEG = FS_NOTE = 7.5             # adopted in-plot type
BORDERPAD, FRAMEPAD = 0.30, 0.12           # t3b
LABELSPACING, HANDLETEXTPAD = 0.32, 0.32   # t3b

# placement B: the two points differ by 1.2% in OR (0.8075 vs 0.7976), so they sit
# at the same height by construction and only vertical label offset separates them.
LABELS = {
    "neurons3": (5, 5, "left"),
    "neurons0": (5, -8, "left"),
    "neurons5": (0, 10, "center"),
    "neurons1": (0, -10, "center"),
    "neurons4": (5, 4, "left"),
    "neurons2": (-6, 5, "right"),
}

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": FS_AXLAB,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",                # editable text in Illustrator
})


def draw():
    fig, ax = plt.subplots(figsize=(3.4, 3.0), dpi=300)
    xs = mar["PD_pct"].to_numpy(float)
    ys = mar["OR"].to_numpy(float)

    ax.axhline(1.0, ls=(0, (1, 2)), lw=0.8, color=GUIDE, zorder=1)
    ax.axvline(EXPECTED_PD_PCT, ls=(0, (3, 2)), lw=0.8, color=FAINT, zorder=1)

    for _, r in mar.iterrows():
        col = C_VULN if r["vulnerable"] else C_RES
        ax.plot([r["PD_pct"]] * 2, [r["lo"], r["hi"]], color=col, lw=0.9,
                alpha=0.55, solid_capstyle="butt", zorder=2)

    # trend fitted in log2 space, i.e. in the space the axis shows
    sl, ic = np.polyfit(xs, np.log2(ys), 1)
    xf = np.linspace(xs.min() - 1.5, xs.max() + 1.5, 60)
    ax.plot(xf, 2 ** (sl * xf + ic), color=FAINT, lw=0.9, zorder=1)

    for _, r in mar.iterrows():
        col = C_VULN if r["vulnerable"] else C_RES
        ax.scatter(r["PD_pct"], r["OR"], s=46, color=col, edgecolor="white",
                   linewidth=0.7, zorder=3)
        dx, dy, ha = LABELS[r["cluster"]]
        ax.annotate(f"Subpop {r['cluster'].replace('neurons', '')}",
                    xy=(r["PD_pct"], r["OR"]), xytext=(dx, dy),
                    textcoords="offset points", ha=ha, va="center",
                    fontsize=FS_PT, color=TEXT, zorder=4)

    # y limits derived from what is drawn, so no CI bound can be cropped
    ylo, yhi = float(mar["lo"].min()) / 1.45, float(mar["hi"].max()) * 1.45
    ax.set_yscale("log", base=2)
    ticks = [t for t in (0.125, 0.25, 0.5, 1, 2, 4, 8) if ylo <= t <= yhi]
    ax.yaxis.set_major_locator(FixedLocator(ticks))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.yaxis.set_major_formatter(FixedFormatter(
        [("%g" % t) if t >= 1 else str(t) for t in ticks]))
    ax.set_ylim(ylo, yhi)
    ax.set_xlim(xs.min() - 4.0, xs.max() + 4.5)

    ax.set_xlabel("PD-donor nuclei (%)", color=TEXT)
    ax.set_ylabel("M11 enrichment (OR)", color=TEXT)
    ax.annotate("expected 28%", xy=(EXPECTED_PD_PCT, 1.0),
                xycoords=("data", "axes fraction"), xytext=(0, 5),
                textcoords="offset points", ha="center", va="bottom",
                fontsize=FS_PT * 0.95, color=FAINT)
    ax.annotate(f"Spearman $\\rho$ = {rho:.2f}\np = {pval:.3f}",
                xy=(0.97, 0.96), xycoords="axes fraction", ha="right",
                va="top", fontsize=FS_NOTE, color=TEXT)

    leg = ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="", mfc=C_VULN, mec="white", ms=6.5,
               label="vulnerable"),
        Line2D([0], [0], marker="o", ls="", mfc=C_RES, mec="white", ms=6.5,
               label="resilient")],
        loc="lower left", fontsize=FS_LEG, frameon=True,
        handletextpad=HANDLETEXTPAD, labelspacing=LABELSPACING,
        borderpad=BORDERPAD, title="SN subpopulation", title_fontsize=FS_LEG)
    leg.get_title().set_color(TEXT)
    for t in leg.get_texts():
        t.set_color(TEXT)
    fr = leg.get_frame()
    fr.set_facecolor(LEGEND_BG)
    fr.set_edgecolor(LEGEND_EDGE)
    fr.set_linewidth(0.6)
    fr.set_boxstyle("square", pad=FRAMEPAD)

    ax.spines[["left", "bottom"]].set_color(GUIDE)
    ax.spines[["left", "bottom"]].set_linewidth(0.8)
    ax.tick_params(colors=TEXT, length=2.5, width=0.7)

    fig.tight_layout()

    # ---- data-fidelity and specification guards, before anything is written
    lo_ax, hi_ax = ax.get_ylim()
    assert lo_ax <= ys.min() and ys.max() <= hi_ax, "a point is outside the axes"
    assert lo_ax <= mar["lo"].min(), (
        f"lower CI cropped ({mar['lo'].min():.4f} < {lo_ax:.4f})")
    assert mar["hi"].max() <= hi_ax, (
        f"upper CI cropped ({mar['hi'].max():.4f} > {hi_ax:.4f})")
    assert ax.xaxis.label.get_fontsize() == FS_AXLAB
    assert ax.yaxis.label.get_fontsize() == FS_AXLAB
    assert ax.get_xticklabels()[0].get_fontsize() == FS_TICK
    assert leg.get_title().get_text() == "SN subpopulation"
    assert leg.get_title().get_fontsize() == FS_LEG
    assert [t.get_text() for t in leg.get_texts()] == ["vulnerable", "resilient"]
    # no point label may still carry the removed TH+ suffix
    assert not any("TH" in t.get_text() for t in ax.texts), "TH+ suffix reappeared"
    return fig


fig = draw()
png = FINAL / "panel_E_martirosyan.png"
png_tr = FINAL / "panel_E_martirosyan_transparent.png"
svg = FINAL / "panel_E_martirosyan.svg"
# White PNG = review copy and md5 of record. Transparent PNG + SVG = placement art.
fig.savefig(png, bbox_inches="tight", facecolor="white")
fig.savefig(png_tr, bbox_inches="tight", transparent=True)
fig.savefig(svg, bbox_inches="tight", transparent=True)
plt.close(fig)

# ---- verify transparency in the written files, do not assume savefig honoured it
svg_text = svg.read_text()
assert 'style="fill: #ffffff"' not in svg_text, "SVG carries a white canvas fill"

import matplotlib.image as mpimg
arr = mpimg.imread(png_tr)
assert arr.shape[2] == 4, f"transparent PNG has no alpha channel: shape {arr.shape}"
corners = [arr[0, 0, 3], arr[0, -1, 3], arr[-1, 0, 3], arr[-1, -1, 3]]
assert max(corners) == 0.0, f"transparent PNG corners are not transparent: {corners}"
assert arr[..., 3].max() > 0.99, "transparent PNG has no opaque ink at all"
opaque_frac = float((arr[..., 3] > 0).mean())
# the white PNG must, by contrast, be fully opaque
arr_w = mpimg.imread(png)
assert arr_w.shape[2] == 3 or arr_w[..., 3].min() == 1.0, "white PNG is not opaque"

rows = [{"file": p.name,
         "bytes": p.stat().st_size,
         "md5": hashlib.md5(p.read_bytes()).hexdigest(),
         "check": "md5 (valid)" if p.suffix == ".png" else "byte size (SVG md5 is NOT valid here)"}
        for p in (png, png_tr, svg)]
rows.append({"file": SRC.name, "bytes": SRC.stat().st_size,
             "md5": INPUT_MD5, "check": "input, md5 (valid)"})
pd.DataFrame(rows).to_csv(FINAL / "panel_E_checksums.csv", index=False)

print(pd.DataFrame(rows).to_string(index=False))
print(f"\nSpearman rho = {rho:.4f}  p = {pval:.6f}  (recomputed, n = {len(mar)})")
print(f"CI span drawn: {mar['lo'].min():.4f} to {mar['hi'].max():.4f}")
print(f"transparent PNG: alpha present, 4 corners alpha=0, {100 * opaque_frac:.1f}% of pixels carry ink")
print(f"transparent SVG: 0 white canvas rects, {svg_text.count('<text')} editable <text> elements")
