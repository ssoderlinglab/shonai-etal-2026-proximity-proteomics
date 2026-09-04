"""
Supplemental Figure 5, panel B: principal component analysis of the co-fractionation profiles.

Input
    data/ibiocofrac/fractionation/LRRK2_young_for_PCA.csv
Output
    figures/suppl_05/output/*.{png,svg,csv}
Run
    python3 figures/suppl_05/scripts/fractionation_pca.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib
"""

import os
import sys
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'ibiocofrac' / 'fractionation'
OUTPUT_DIR = FIGURE_DIR / 'output'
SRC_DIR = REPO_ROOT / 'src'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES
apply_style()

INPUT_CSV = DATA_DIR / 'LRRK2_young_for_PCA.csv'

# Fraction colors: ColorBrewer Set3 (qualitative)
FRACTION_COLORS = {
    4:  "#FB8072",  # salmon
    5:  "#80B1D3",  # light blue
    6:  "#B3DE69",  # light green
    7:  "#FDB462",  # light orange
    8:  "#BC80BD",  # medium purple
    9:  "#8DD3C7",  # light teal
    10: "#FFED6F",  # light gold
}
FRACTION_ORDER = [4, 5, 6, 7, 8, 9, 10]


def classify_geno(s):
    s = str(s).lower()
    if "wt" in s: return "WT"
    elif "tg" in s or "ki" in s: return "Mut"
    elif "spqc" in s: return "SPQC"
    return "Unknown"


def draw_ellipse(ax, x, y, color, alpha=0.15):
    """95% confidence ellipse (chi² = 5.991 for 2 df)."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    mean_x, mean_y = np.mean(x), np.mean(y)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
    width, height = 2 * np.sqrt(np.maximum(eigenvalues, 0) * 5.991)
    ellipse = patches.Ellipse(
        (mean_x, mean_y), width, height, angle=angle,
        edgecolor=color, facecolor=color, alpha=alpha, linewidth=1.5
    )
    ax.add_patch(ellipse)


def main():
    data = pd.read_csv(INPUT_CSV)
    data["genotype"] = data["sample"].apply(classify_geno)
    data = data[data["genotype"] == "WT"].copy()
    data["fraction"] = data["fraction"].astype(int)

    print(f"LRRK2-WT samples: {len(data)}")
    print(f"Fractions: {sorted(data['fraction'].unique())}")

    meta_cols = ["Accession", "fraction", "replicate", "sample", "genotype"]
    meta = data[meta_cols].reset_index(drop=True)
    numeric = data.drop(columns=meta_cols, errors="ignore")
    numeric = numeric.apply(pd.to_numeric, errors="coerce").fillna(0).reset_index(drop=True)

    print(f"Proteins: {numeric.shape[1]}")

    scaler = StandardScaler()
    scaled = scaler.fit_transform(numeric)
    pca = PCA(n_components=4)
    scores = pca.fit_transform(scaled)
    pct = pca.explained_variance_ratio_ * 100
    scores_df = pd.DataFrame(scores, columns=["PC1", "PC2", "PC3", "PC4"])

    print(f"PC1: {pct[0]:.1f}%  PC2: {pct[1]:.1f}%  PC3: {pct[2]:.1f}%  PC4: {pct[3]:.1f}%")

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    for frac in FRACTION_ORDER:
        mask = meta["fraction"] == frac
        if mask.sum() == 0:
            continue
        color = FRACTION_COLORS[frac]
        x = scores_df.loc[mask, "PC1"].values
        y = scores_df.loc[mask, "PC2"].values
        draw_ellipse(ax, x, y, color)
        ax.scatter(x, y, c=color, s=120, edgecolors="black",
                   linewidths=0.8, alpha=0.95, zorder=3, label=f"F{frac}")

    ax.set_xlabel(f"PC1 ({pct[0]:.1f}%)", fontsize=FONT_SIZES['label'])
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}%)", fontsize=FONT_SIZES['label'])
    ax.set_title("Differential Centrifugation Fractionation — PCA",
                 fontsize=FONT_SIZES['title'], fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, linestyle="--", linewidth=0.4, alpha=0.4)

    legend_elements = [
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=FRACTION_COLORS[f], markeredgecolor="black",
               markersize=10, label=f"F{f}")
        for f in FRACTION_ORDER
    ]
    ax.legend(handles=legend_elements, fontsize=FONT_SIZES['tick'],
              loc="best", framealpha=0.9, edgecolor="gray",
              title="Fraction", title_fontsize=FONT_SIZES['tick'])

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig.tight_layout()
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'fractionation_pca.{ext}',
                    dpi=300, bbox_inches="tight", facecolor='white')
    plt.close(fig)

    # Save variance
    var_df = pd.DataFrame({
        "Component": [f"PC{i+1}" for i in range(4)],
        "Variance_Explained_Pct": [round(p, 2) for p in pct[:4]],
    })
    var_df.to_csv(OUTPUT_DIR / 'fractionation_pca_variance.csv', index=False)

    print(f"\nSaved: fractionation_pca.png/svg + variance CSV")


if __name__ == '__main__':
    main()
