"""
Supplemental Figure 8, panel B: convergent module enrichment across the ten published human dopamine neuron subtypes.

Input
    data/human/subtype_enrichment_M11.json
Output
    figures/suppl_08/output/panel_B_subtype_ranking.{png,svg}
    figures/suppl_08/output/panel_B_checksums.csv
Run
    python3 figures/suppl_08/scripts/panel_B_subtype_ranking.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib

Odds ratios and their significance are read from the frozen input; nothing is recomputed here. Benjamini-Hochberg correction is over the ten subtypes.

SOX6_AGTR1 is drawn at full strength because it is the subtype the source publication reports as selectively depleted, not because vulnerability was measured here.

The ranking is threshold-dependent: marker-set size and odds ratio are correlated (Spearman rho -0.818, p 0.0038) and the three highest subtypes hold three of the shortest marker lists, so the order is not size-independent.
"""

import hashlib
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
matplotlib.rcParams["pdf.fonttype"] = 42

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
IN = os.path.join(REPO_ROOT, "data", "human")
OUT = os.path.join(FIGURE_DIR, "output")
os.makedirs(OUT, exist_ok=True)

# Frozen input, derived from the published per-subtype marker table and
# reproducing the project's module by subtype enrichment matrix to 2 d.p.
DATA = json.load(open(os.path.join(IN, "subtype_enrichment_M11.json")))

# Lineage palette, fixed with panel A (sampled from Luppi Fig 7A, then assigned
# SOX6 = red / CALB1 = green on Daichi's instruction — note this INVERTS the
# source figure's own assignment, so no caption may claim it follows Luppi).
#
# Daichi 2026-08-12: the two lineage colours are MUTED (blended toward white)
# so the full-strength red can be reserved for SOX6_AGTR1. The saturated hue
# now carries one meaning only — "this is the subtype the source paper reports
# as selectively depleted" — instead of competing with nine other bars.
#
# A DARKER red was tried first and rejected: at bar size a very dark fill reads
# as *less* prominent than its lighter neighbours, not more.
SOX6_FULL = "#BA1028"      # reserved for SOX6_AGTR1 only
CALB1_FULL = "#88B838"
SOX6 = "#D6707E"           # lighten(#BA1028, 0.40)
CALB1 = "#B8D488"          # lighten(#88B838, 0.40)
AGTR1 = SOX6_FULL
TEXT = "#1A1A1A"
GUIDE = "#6E6E6E"

FS_TICK, FS_LAB, FS_LEG, FS_NOTE = 7.2, 8.0, 6.8, 5.6


def stars(q):
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "n.s."


def pretty(name):
    return name.replace("_", "\u00b7")


def base_axes(ax):
    ax.set_xlabel("M11 enrichment odds ratio", fontsize=FS_LAB, color=TEXT)
    ax.axvline(1.0, ls=(0, (2, 2)), lw=0.8, color=GUIDE, zorder=1)
    ax.set_xlim(0, 8.6)
    ax.invert_yaxis()
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GUIDE)
    ax.tick_params(labelsize=FS_TICK, colors=TEXT, length=2.5, width=0.7)


# "selective vulnerability" is the source paper's own phrase (Results: "the
# selective vulnerability of SOX6_AGTR1 cells"). It is an attribution to Kamath
# et al., NOT a property measured in this study — the caption must say so.
#
# Daichi 2026-08-12 (chosen layout, v5_final): drop the note from the y-axis
# gutter and carry the phrase in the LEGEND instead. This is the better place
# for it: in the gutter the phrase sits beside a tick label and reads like an
# axis category, i.e. like something this panel measured; in the legend it is
# unambiguously a key entry describing what the highlighted colour denotes.
NOTE = "selective vulnerability"
AGTR1_LEGEND = "SOX6_AGTR1\n(selective vulnerability in PD)"

# the design the figure carries; the rejected variants are not built here
VARIANTS = ["v5_final"]


def draw(ax, variant):
    """All variants: muted lineage fills, full-strength red for SOX6_AGTR1.
    They differ only in how the AGTR1 tick label is annotated."""
    y = np.arange(len(DATA))
    for i, r in enumerate(DATA):
        is_ag = r["subtype"] == "SOX6_AGTR1"
        col = AGTR1 if is_ag else (SOX6 if r["lineage"] == "SOX6" else CALB1)
        ec, lw = ("none", 0.0)
        if variant in ("v3_note_outline", "v5_final") and is_ag:
            ec, lw = TEXT, 0.9
        ax.barh(i, r["OR"], height=0.72, color=col, edgecolor=ec, lw=lw, zorder=3)
        ax.text(r["OR"] + 0.16, i, stars(r["q"]), va="center", ha="left",
                fontsize=FS_TICK, color=TEXT if r["q"] < 0.05 else GUIDE, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels([pretty(r["subtype"]) for r in DATA])
    for t, r in zip(ax.get_yticklabels(), DATA):
        if r["subtype"] == "SOX6_AGTR1":
            t.set_fontweight("bold")
            t.set_color(AGTR1)
    base_axes(ax)

    # the annotation
    idx = next(i for i, r in enumerate(DATA) if r["subtype"] == "SOX6_AGTR1")
    if variant in ("v2_note", "v3_note_outline"):
        # a second line under the tick label, in the y-tick gutter
        ax.text(-0.14, idx + 0.30, NOTE, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=FS_NOTE, style="italic",
                color=AGTR1, clip_on=False, zorder=6)
    elif variant == "v4_note_axlabel":
        # inside the plot, at the base of the bar — keeps the gutter clean
        ax.text(0.16, idx + 0.02, NOTE, ha="left", va="center",
                fontsize=FS_NOTE, style="italic", color="white", zorder=6)


def legend(ax, variant):
    lab = AGTR1_LEGEND if variant == "v5_final" else "SOX6_AGTR1 (depleted in PD)"
    h = [Patch(facecolor=AGTR1,
               edgecolor=TEXT if variant in ("v3_note_outline", "v5_final") else "none",
               lw=0.9, label=lab),
         Patch(facecolor=SOX6, label="SOX6 lineage"),
         Patch(facecolor=CALB1, label="CALB1 lineage")]
    if variant == "v5_final":
        # The AGTR1 entry is long. On one line it overran the two shortest bars
        # and their significance marks; moved below the axes it hung outside the
        # panel. Kept INSIDE the axes, wrapped onto two lines and set slightly
        # smaller, in the empty wedge the short bottom bars leave free. Verified
        # numerically against every bar patch and every text in the axes.
        ax.legend(handles=h, loc="lower right", fontsize=FS_LEG - 0.9,
                  frameon=False, handlelength=1.0, handleheight=0.9,
                  borderpad=0.2, labelspacing=0.42,
                  borderaxespad=0.35)
    else:
        ax.legend(handles=h, loc="lower right", fontsize=FS_LEG, frameon=False,
                  handlelength=1.1, handleheight=1.0, borderpad=0.3,
                  labelspacing=0.45)


def main():
    manifest = []
    for v in VARIANTS:
        fig = plt.figure(figsize=(4.2, 3.0))
        ax = fig.add_axes([0.36, 0.14, 0.60, 0.83])
        draw(ax, v)
        legend(ax, v)
        for ext in ("png", "svg"):
            p = os.path.join(OUT, f"panel_B_subtype_ranking.{ext}")
            fig.savefig(p, transparent=True, bbox_inches="tight",
                        pad_inches=0.02, dpi=600 if ext == "png" else None)
            manifest.append((
                f"panel_B_subtype_ranking.{ext}", os.path.getsize(p),
                hashlib.md5(open(p, "rb").read()).hexdigest()
                if ext == "png" else "SVG-md5-not-a-valid-check"))
        plt.close(fig)

    mp = os.path.join(OUT, "panel_B_checksums.csv")
    with open(mp, "w") as f:
        f.write("# PNG md5 is a valid identity check. SVG md5 is NOT:\n"
                "# matplotlib stamps a <metadata> date and randomises clip and\n"
                "# glyph ids, so a re-render changes SVG bytes while the artwork\n"
                "# is identical. Compare SVGs after stripping <metadata> and\n"
                "# normalising ids.\n")
        f.write("file,bytes,md5\n")
        for r in sorted(manifest):
            f.write(f"{r[0]},{r[1]},{r[2]}\n")
    print(f"{len(manifest)} files -> {OUT}")


if __name__ == "__main__":
    main()
