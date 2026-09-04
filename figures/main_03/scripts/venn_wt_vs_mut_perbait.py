"""
Figure 3, panels C, E and G: per-bait Venn diagrams of the WT against the mutant proximity proteome.

Input
    data/ibioid/proxiome_sets/{aSyn,Vps35,Lrrk2}_{WT,mutant}.csv
Output
    figures/main_03/output/*.{png,svg}
Run
    python3 figures/main_03/scripts/venn_wt_vs_mut_perbait.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, matplotlib-venn, shapely
"""



import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pandas as pd

# ── Paths (relative to Figure_3/) ────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid', 'proxiome_sets')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')

DPI = 300
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["svg.fonttype"] = "none"


def load_set(path):
    """Load protein accessions from proxiome CSV."""
    df = pd.read_csv(path)
    return set(df["Accession"].dropna())


def draw_venn2(set_a, set_b, label_a, label_b, title,
               color_a="#4C8DC1", color_b="#D15D5D", output_path=None):
    """
    Clean 2-circle Venn diagram.

    Overlap percentage denominator is max(|A|, |B|).
    This is a standard approach; Jaccard index (using union) is an alternative.
    """
    only_a  = len(set_a - set_b)
    only_b  = len(set_b - set_a)
    both    = len(set_a & set_b)
    total_a = len(set_a)
    total_b = len(set_b)

    # Sanity check: regions must sum correctly
    assert only_a + both == total_a, f"WT check failed: {only_a}+{both} != {total_a}"
    assert only_b + both == total_b, f"MUT check failed: {only_b}+{both} != {total_b}"

    pct_overlap = both / max(total_a, total_b) * 100

    fig, ax = plt.subplots(figsize=(5, 4.2))
    ax.set_aspect("equal")
    ax.axis("off")

    # Circle parameters
    r = 1.4
    offset = 0.85
    cx_a, cx_b = -offset, offset
    cy = 0.0

    circ_a = Circle((cx_a, cy), r, color=color_a, alpha=0.35, zorder=1)
    circ_b = Circle((cx_b, cy), r, color=color_b, alpha=0.35, zorder=1)
    ax.add_patch(circ_a)
    ax.add_patch(circ_b)

    # Unique counts
    ax.text(cx_a - r * 0.52, cy, str(only_a),
            ha="center", va="center", fontsize=18, fontweight="bold", zorder=3)
    ax.text(cx_b + r * 0.52, cy, str(only_b),
            ha="center", va="center", fontsize=18, fontweight="bold", zorder=3)
    ax.text(0, cy, str(both),
            ha="center", va="center", fontsize=18, fontweight="bold", zorder=3)

    # Labels below circles
    label_y = cy - r - 0.28
    ax.text(cx_a, label_y, f"{label_a}\n(n={total_a})",
            ha="center", va="top", fontsize=11, fontweight="bold")
    ax.text(cx_b, label_y, f"{label_b}\n(n={total_b})",
            ha="center", va="top", fontsize=11, fontweight="bold")

    # Overlap annotation
    ax.text(0, cy - r - 1.05,
            f"Overlap: {both} proteins ({pct_overlap:.1f}%)",
            ha="center", va="top", fontsize=9.5,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                      edgecolor="#AAAAAA", alpha=0.9))

    ax.set_title(title, fontsize=12, fontweight="bold", pad=6)
    ax.set_xlim(-r - offset - 0.5, r + offset + 0.5)
    ax.set_ylim(cy - r - 1.4, cy + r + 0.5)

    plt.tight_layout()
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        for fmt in ["png", "svg"]:
            fpath = output_path.rsplit('.', 1)[0] + f'.{fmt}'
            plt.savefig(fpath, dpi=DPI, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {os.path.basename(output_path)}")
    plt.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("Per-bait 2-way Venn — WT vs Mutant (Figure 3 C/E/G)")
    print("=" * 60)

    configs = [
        {
            'wt': os.path.join(DATA_DIR, 'Vps35_WT.csv'),
            'mut': os.path.join(DATA_DIR, 'Vps35_D620N.csv'),
            'label_a': 'WT', 'label_b': 'D620N',
            'title': 'VPS35 Proximity Proteome',
            'output': 'Venn_Vps35.png',
            'panel': 'C',
        },
        {
            'wt': os.path.join(DATA_DIR, 'aSyn_WT.csv'),
            'mut': os.path.join(DATA_DIR, 'aSyn_A30P.csv'),
            'label_a': 'WT', 'label_b': 'A30P',
            'title': '\u03b1-Synuclein Proximity Proteome',
            'output': 'Venn_aSyn.png',
            'panel': 'E',
        },
        {
            'wt': os.path.join(DATA_DIR, 'Lrrk2_WT.csv'),
            'mut': os.path.join(DATA_DIR, 'Lrrk2_G2019S.csv'),
            'label_a': 'WT', 'label_b': 'G2019S',
            'title': 'LRRK2 Proximity Proteome',
            'output': 'Venn_Lrrk2.png',
            'panel': 'G',
        },
    ]

    for cfg in configs:
        wt_set = load_set(cfg['wt'])
        mut_set = load_set(cfg['mut'])
        print(f"\n  Panel {cfg['panel']}: {cfg['title']}")
        print(f"    WT: {len(wt_set)}, MUT: {len(mut_set)}")

        draw_venn2(
            wt_set, mut_set,
            cfg['label_a'], cfg['label_b'],
            cfg['title'],
            color_a="#4C8DC1", color_b="#D06050",
            output_path=os.path.join(OUTPUT_DIR, cfg['output']),
        )


if __name__ == "__main__":
    main()
