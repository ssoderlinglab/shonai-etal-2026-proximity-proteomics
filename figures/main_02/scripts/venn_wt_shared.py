"""
Figure 2, panel B: three-way Venn diagram of the WT proximity proteomes.

Input
    data/ibioid/proxiome_sets/{aSyn,Vps35,Lrrk2}_WT.csv
Output
    figures/main_02/output/venn_wt_shared.{png,svg}
Run
    python3 figures/main_02/scripts/venn_wt_shared.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, matplotlib-venn, shapely

Set sizes are asserted against the values reported in the manuscript (425, 361, 311; 86 in all three; 769 unique).
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path
import pandas as pd
from shapely.geometry import Point

# ── Paths ───────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid', 'proxiome_sets')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')

DPI = 300
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"

# ── Visual config (inherited verbatim) ───────────────────────────────
CONVERGENT_BORDER_LW = 3.0
CONVERGENT_BORDER_COLOR = "#333333"
LABEL_FONTSIZE = 14

# ── Expected set sizes, as reported in the manuscript ───────────────
EXPECTED = {"aSyn": 425, "Vps35": 361, "Lrrk2": 311, "3way": 86, "union": 769}


def load_set(path):
    """Load protein accessions from proxiome CSV."""
    df = pd.read_csv(path)
    return set(df["Accession"].dropna())


def make_convergent_border_patch(cx, cy, r):
    """
    Build a matplotlib PathPatch outlining the region where >=2 of 3 circles
    overlap. Uses shapely for robust polygon boolean operations. The union of
    all pairwise intersections forms a single connected region, so no internal
    border lines appear at the 3-way center.
    """
    circles = [Point(cx[i], cy[i]).buffer(r, resolution=128) for i in range(3)]
    ab = circles[0].intersection(circles[1])
    ac = circles[0].intersection(circles[2])
    bc = circles[1].intersection(circles[2])
    convergent = ab.union(ac).union(bc)

    def shapely_to_mpl_path(geom):
        verts = []
        codes = []
        polygons = list(geom.geoms) if geom.geom_type == 'MultiPolygon' else [geom]
        for poly in polygons:
            ext = np.array(poly.exterior.coords)
            verts.extend(ext.tolist())
            codes.append(Path.MOVETO)
            codes.extend([Path.LINETO] * (len(ext) - 2))
            codes.append(Path.CLOSEPOLY)
        return Path(verts, codes)

    mpl_path = shapely_to_mpl_path(convergent)
    return PathPatch(
        mpl_path,
        facecolor="none",
        edgecolor=CONVERGENT_BORDER_COLOR,
        linewidth=CONVERGENT_BORDER_LW,
        zorder=3,
    )


def draw_venn3(set_a, set_b, set_c,
               label_a, label_b, label_c,
               color_a, color_b, color_c,
               convergent_label, output_path=None):
    """3-circle Venn with standard set algebra."""
    only_a   = len(set_a - set_b - set_c)
    only_b   = len(set_b - set_a - set_c)
    only_c   = len(set_c - set_a - set_b)
    ab_only  = len((set_a & set_b) - set_c)
    ac_only  = len((set_a & set_c) - set_b)
    bc_only  = len((set_b & set_c) - set_a)
    abc      = len(set_a & set_b & set_c)
    total_a  = len(set_a)
    total_b  = len(set_b)
    total_c  = len(set_c)

    computed_union = only_a + only_b + only_c + ab_only + ac_only + bc_only + abc
    actual_union = len(set_a | set_b | set_c)
    assert computed_union == actual_union, (
        f"Region sum ({computed_union}) != union ({actual_union})"
    )

    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    # v02: opaque white background. ax.axis("off") hides the spines and ticks but does NOT
    # make the axes patch transparent, so the patch must be given a facecolor explicitly and
    # set_visible(True) re-asserted after axis("off") turned the patch off.
    fig.patch.set_facecolor("white")
    fig.patch.set_alpha(1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.patch.set_facecolor("white")
    ax.patch.set_alpha(1.0)
    ax.patch.set_visible(True)

    r = 1.55
    d = 0.90
    angle = np.radians(90)
    angles = [angle, angle + np.radians(120), angle + np.radians(240)]
    cx = [d * np.cos(a) for a in angles]
    cy = [d * np.sin(a) for a in angles]

    colors = [color_a, color_b, color_c]
    for i in range(3):
        ax.add_patch(Circle((cx[i], cy[i]), r, color=colors[i], alpha=0.30, zorder=1))

    # Bold border around convergent region (>=2 overlap)
    ax.add_patch(make_convergent_border_patch(cx, cy, r))

    # Region count labels
    s = 0.62
    ax.text(cx[0], cy[0] + r * s, str(only_a),
            ha="center", va="center", fontsize=15, fontweight="bold", zorder=5)
    ax.text(cx[1] - r * s * 0.7, cy[1] - r * s * 0.4, str(only_b),
            ha="center", va="center", fontsize=15, fontweight="bold", zorder=5)
    ax.text(cx[2] + r * s * 0.7, cy[2] - r * s * 0.4, str(only_c),
            ha="center", va="center", fontsize=15, fontweight="bold", zorder=5)
    ax.text((cx[0]+cx[1])/2 - r*0.25, (cy[0]+cy[1])/2 + r*0.05, str(ab_only),
            ha="center", va="center", fontsize=13, fontweight="bold", zorder=5)
    ax.text((cx[0]+cx[2])/2 + r*0.25, (cy[0]+cy[2])/2 + r*0.05, str(ac_only),
            ha="center", va="center", fontsize=13, fontweight="bold", zorder=5)
    ax.text((cx[1]+cx[2])/2, (cy[1]+cy[2])/2 - r*0.3, str(bc_only),
            ha="center", va="center", fontsize=13, fontweight="bold", zorder=5)
    ax.text(np.mean(cx), np.mean(cy), str(abc),
            ha="center", va="center", fontsize=15, fontweight="bold",
            color="black", zorder=5)

    # Bait labels
    offset_labels = 0.38
    ax.text(cx[0], cy[0] + r + offset_labels,
            f"{label_a}\n(n={total_a})",
            ha="center", va="bottom", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.text(cx[1] - r * 0.55, cy[1] - r - offset_labels * 0.5,
            f"{label_b}\n(n={total_b})",
            ha="center", va="top", fontsize=LABEL_FONTSIZE, fontweight="bold")
    ax.text(cx[2] + r * 0.55, cy[2] - r - offset_labels * 0.5,
            f"{label_c}\n(n={total_c})",
            ha="center", va="top", fontsize=LABEL_FONTSIZE, fontweight="bold")

    # Convergent region label
    conv_label_x = np.mean(cx) + r * 1.55
    conv_label_y = np.mean(cy) + r * 0.72
    ax.text(conv_label_x, conv_label_y, convergent_label,
            ha="center", va="center", fontsize=15, fontstyle="italic",
            fontweight="bold", color="#333333", zorder=6)

    margin = r + d + 0.9
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin * 1.05, margin)

    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        for fmt in ["png", "svg"]:
            fpath = output_path.rsplit('.', 1)[0] + f'.{fmt}'
            # v02: transparent=False + explicit facecolor/edgecolor so the white patch set on
            # the figure is actually written to the file (savefig overrides fig.patch unless
            # facecolor is passed through).
            plt.savefig(fpath, dpi=DPI, bbox_inches="tight", transparent=False,
                        facecolor="white", edgecolor="none")
        print(f"  Saved: {fpath}")
    plt.close()

    return dict(only_a=only_a, only_b=only_b, only_c=only_c,
                ab_only=ab_only, ac_only=ac_only, bc_only=bc_only,
                abc=abc, union=actual_union)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    asyn_wt  = load_set(os.path.join(DATA_DIR, 'aSyn_WT.csv'))
    vps35_wt = load_set(os.path.join(DATA_DIR, 'Vps35_WT.csv'))
    lrrk2_wt = load_set(os.path.join(DATA_DIR, 'Lrrk2_WT.csv'))

    print("=" * 60)
    print("3-way Venn — WT Convergent Proximity Proteome (Figure 2B)")
    print("=" * 60)
    print(f"  aSyn WT:  {len(asyn_wt)} proteins")
    print(f"  Vps35 WT: {len(vps35_wt)} proteins")
    print(f"  Lrrk2 WT: {len(lrrk2_wt)} proteins")

    assert len(asyn_wt) == EXPECTED["aSyn"], len(asyn_wt)
    assert len(vps35_wt) == EXPECTED["Vps35"], len(vps35_wt)
    assert len(lrrk2_wt) == EXPECTED["Lrrk2"], len(lrrk2_wt)

    regions = draw_venn3(
        asyn_wt, vps35_wt, lrrk2_wt,
        "\u03b1-Syn WT", "VPS35 WT", "LRRK2 WT",
        color_a="#4472C4", color_b="#70AD47", color_c="#ED7D31",
        convergent_label="WT Majority\nProxiome",
        output_path=os.path.join(OUTPUT_DIR, "venn_wt_shared.png"),
    )

    assert regions["abc"] == EXPECTED["3way"], regions["abc"]
    assert regions["union"] == EXPECTED["union"], regions["union"]

    print(f"\n  Regions: only aSyn {regions['only_a']} · only Vps35 {regions['only_b']} · "
          f"only Lrrk2 {regions['only_c']}")
    print(f"           aSyn&Vps35 {regions['ab_only']} · aSyn&Lrrk2 {regions['ac_only']} · "
          f"Vps35&Lrrk2 {regions['bc_only']}")
    print(f"  3-way intersection: {regions['abc']} proteins")
    print(f"  >=2/3 (majority)  : "
          f"{regions['ab_only'] + regions['ac_only'] + regions['bc_only'] + regions['abc']} proteins")
    print(f"  Union: {regions['union']} proteins")
    print("  all expected counts OK")


if __name__ == "__main__":
    main()
