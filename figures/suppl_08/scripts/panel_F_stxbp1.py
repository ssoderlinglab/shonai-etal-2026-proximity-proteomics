"""
Supplemental Figure 8, panel F: STXBP1 abundance in dopamine neurons surviving in disease donors, by subtype.

Input
    data/human/stxbp1_reproduction_table.csv
Output
    figures/suppl_08/output/panel_F_stxbp1.png
    figures/suppl_08/output/panel_F_stxbp1_transparent.png
    figures/suppl_08/output/panel_F_stxbp1.svg
    figures/suppl_08/output/panel_F_checksums.csv
Run
    python3 figures/suppl_08/scripts/panel_F_stxbp1.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib

Ten subtypes. Estimates are donor-level pseudobulk group coefficients with age as a covariate, corrected across the ten subtypes; the values are read from the frozen input and are not recomputed. Horizontal lines are the 95 per cent confidence interval of the group coefficient.

The input md5 is asserted before drawing. The script also asserts that the largest estimate belongs to CALB1_RBP4, which is the fact the legend turns on.

Marker size does not encode significance; the script asserts a single marker size.
"""
import hashlib
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
SRC = REPO_ROOT / "data" / "human" / "stxbp1_reproduction_table.csv"
FINAL = FIGURE_DIR / "output"
FINAL.mkdir(exist_ok=True)

# ------------------------------------------------------------------ input guard
INPUT_MD5 = "d630b4795072bdf32bdf8636040c5fc9"
assert hashlib.md5(SRC.read_bytes()).hexdigest() == INPUT_MD5, "input CSV changed"

stx = pd.read_csv(SRC)
assert len(stx) == 10, f"expected 10 subtypes, got {len(stx)}"
for col in ("subtype", "class", "n_PD_LBD", "n_control", "log2FC_ageadj",
            "ci_lo", "ci_hi", "p_ageadj", "FDR"):
    assert col in stx.columns, f"missing column {col}"
assert (stx["log2FC_ageadj"] >= 0).all(), "a subtype effect is negative"
assert (stx["ci_lo"] < stx["log2FC_ageadj"]).all(), "CI does not bracket the estimate"
assert (stx["log2FC_ageadj"] < stx["ci_hi"]).all(), "CI does not bracket the estimate"
# the caption-critical fact: the top effect is a RESISTANT subtype
assert stx.loc[stx["log2FC_ageadj"].idxmax(), "subtype"] == "CALB1_RBP4"

# -------------------------------------------------------------------- constants
SOX6, CALB1 = "#BA1028", "#88B838"     # lineage palette, shared with the other panels
TEXT, GUIDE = "#1A1A1A", "#6E6E6E"
SHADE, SHADE_EDGE = "#F0F0F0", "#D5D5D5"

FS_AXLAB, FS_TICK, FS_LEG, FS_STAR = 8.0, 6.5, 6.2, 9.0
XLIM = (-0.5, 1.05)
XLABEL = ("STXBP1 expression in surviving DA neurons\n"
          "PD+LBD vs control, log$_2$FC (95% CI)")

plt.rcParams.update({
    "font.size": 8, "font.family": "DejaVu Sans", "axes.labelsize": FS_AXLAB,
    "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none",                # editable text in Illustrator
})

df = stx.sort_values("log2FC_ageadj", ascending=True).reset_index(drop=True)
df["vuln"] = df["subtype"].str.startswith("SOX6")


def stars(q):
    """** q < 0.05, * q < 0.10. Both thresholds must appear in the caption."""
    return "**" if q < 0.05 else ("*" if q < 0.10 else "")


df["stars"] = df["FDR"].map(stars)
# the tiering must match the data, not an expectation of it
assert list(df["stars"]) == ["", "", "", "", "", "", "*", "**", "*", "**"], \
    f"star tiering changed: {list(df['stars'])}"
assert (df.loc[df["stars"] == "*", "FDR"] > 0.05).all(), \
    "a single-star row now clears q < 0.05 - the caption wording would be wrong"

SHORT = [s.replace("SOX6_", "SOX6\u00b7").replace("CALB1_", "CALB1\u00b7")
          .replace("CRYM_CCDC68", "CRYM") for s in df["subtype"]]


def draw():
    fig, ax = plt.subplots(figsize=(3.4, 3.0), dpi=300)

    ax.axvline(0, color=GUIDE, lw=0.9, zorder=1)
    for i, r in df.iterrows():
        c = SOX6 if r["vuln"] else CALB1
        ax.plot([r["ci_lo"], r["ci_hi"]], [i, i], "-", color=c, lw=1.5, zorder=2)
        ax.plot(r["log2FC_ageadj"], i, "o", color=c, ms=5.5, mec="white",
                mew=0.7, zorder=3)

    # variant s2: asterisks aligned in a column at the right edge
    for i, r in df.iterrows():
        if r["stars"]:
            ax.annotate(r["stars"], xy=(1.0, i),
                        xycoords=("axes fraction", "data"), xytext=(-2, -0.5),
                        textcoords="offset points", ha="right", va="center",
                        fontsize=FS_STAR, color=TEXT)

    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(SHORT, fontsize=FS_TICK)
    for tick, v in zip(ax.get_yticklabels(), df["vuln"]):
        tick.set_color(SOX6 if v else CALB1)
    ax.set_xlabel(XLABEL, color=TEXT)
    ax.set_xlim(*XLIM)
    ax.margins(y=0.04)

    leg = ax.legend(handles=[
        Line2D([0], [0], marker="o", ls="", mfc=SOX6, mec="white", ms=6,
               label="SOX6 lineage"),
        Line2D([0], [0], marker="o", ls="", mfc=CALB1, mec="white", ms=6,
               label="CALB1 lineage")],
        loc="lower right", fontsize=FS_LEG, frameon=True, handletextpad=0.4,
        labelspacing=0.32, borderpad=0.3)
    fr = leg.get_frame()
    fr.set_facecolor(SHADE)
    fr.set_edgecolor(SHADE_EDGE)
    fr.set_linewidth(0.6)
    fr.set_boxstyle("square", pad=0.12)
    for t in leg.get_texts():
        t.set_color(TEXT)

    ax.spines[["left", "bottom"]].set_color(GUIDE)
    ax.spines[["left", "bottom"]].set_linewidth(0.8)
    ax.tick_params(colors=TEXT, length=2.5, width=0.7)
    fig.tight_layout()

    # ---- data-fidelity and specification guards, before anything is written
    x0, x1 = ax.get_xlim()
    assert df["ci_lo"].min() >= x0, (
        f"a CI lower bound is cropped ({df['ci_lo'].min():.4f} < {x0})")
    assert df["ci_hi"].max() <= x1, (
        f"a CI upper bound is cropped ({df['ci_hi'].max():.4f} > {x1})")
    assert "expression" in ax.get_xlabel(), "the axis label lost 'expression'"
    assert "95% CI" in ax.get_xlabel(), "the axis label lost the CI definition"
    assert [t.get_text() for t in leg.get_texts()] == ["SOX6 lineage", "CALB1 lineage"], \
        "the legend is no longer the two lineages only"
    assert [t.get_text() for t in ax.get_yticklabels()] == SHORT
    # marker size must not encode significance again
    sizes = {ln.get_markersize() for ln in ax.lines if ln.get_marker() == "o"}
    assert len(sizes) == 1, f"marker size varies again: {sizes}"
    return fig


fig = draw()
png = FINAL / "panel_F_stxbp1.png"
png_tr = FINAL / "panel_F_stxbp1_transparent.png"
svg = FINAL / "panel_F_stxbp1.svg"
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
assert arr.shape[2] == 4, f"transparent PNG has no alpha channel: {arr.shape}"
corners = [arr[0, 0, 3], arr[0, -1, 3], arr[-1, 0, 3], arr[-1, -1, 3]]
assert max(corners) == 0.0, f"transparent PNG corners are not transparent: {corners}"
assert arr[..., 3].max() > 0.99, "transparent PNG has no opaque ink at all"
opaque_frac = float((arr[..., 3] > 0).mean())
arr_w = mpimg.imread(png)
assert arr_w.shape[2] == 3 or arr_w[..., 3].min() == 1.0, "white PNG is not opaque"

rows = [{"file": p.name, "bytes": p.stat().st_size,
         "md5": hashlib.md5(p.read_bytes()).hexdigest(),
         "check": "md5 (valid)" if p.suffix == ".png"
                  else "byte size (SVG md5 is NOT valid here)"}
        for p in (png, png_tr, svg)]
rows.append({"file": SRC.name, "bytes": SRC.stat().st_size, "md5": INPUT_MD5,
             "check": "input, md5 (valid)"})
pd.DataFrame(rows).to_csv(FINAL / "panel_F_checksums.csv", index=False)

print(pd.DataFrame(rows).to_string(index=False))
print("\nrows as drawn, top of plot first:")
print(df.iloc[::-1][["subtype", "class", "n_PD_LBD", "n_control",
                     "log2FC_ageadj", "ci_lo", "ci_hi", "FDR", "stars"]].to_string(index=False))
print(f"\nCI span drawn: {df['ci_lo'].min():.4f} to {df['ci_hi'].max():.4f} "
      f"(axis {XLIM[0]} to {XLIM[1]})")
print(f"transparent PNG: alpha present, 4 corners alpha=0, {100 * opaque_frac:.1f}% of pixels carry ink")
print(f"transparent SVG: 0 white canvas rects, {svg_text.count('<text')} editable <text> elements")
