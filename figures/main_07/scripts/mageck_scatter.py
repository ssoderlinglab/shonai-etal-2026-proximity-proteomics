"""
Figure 7, panel D: MAGeCK two-condition scatter of gene-level negative selection scores, in the presence and in the absence of A53T alpha-synuclein.

Input
    data/craave_seq/mageck/{aSyn,Baseline}.gene_summary.txt
Output
    figures/main_07/output/panel_D_mageck_scatter.{png,svg}
Run
    python3 figures/main_07/scripts/mageck_scatter.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, src/pub_style.py

Non-targeting control genes (identifier prefix nt) are removed. No other filtering is applied and no statistic is computed here: the negative score, the median log fold change and the false discovery rate are read from the MAGeCK gene summaries as produced.

Marker area is proportional to the absolute median log fold change under alpha-synuclein; marker colour is minus log10 of the false discovery rate under alpha-synuclein.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'craave_seq' / 'mageck'
OUTPUT_DIR = FIGURE_DIR / 'output'
SRC_DIR = REPO_ROOT / 'src'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, COLORS, FONT_SIZES
apply_style()

ASYN_PATH = DATA_DIR / 'aSyn.gene_summary.txt'
BASELINE_PATH = DATA_DIR / 'Baseline.gene_summary.txt'


def main():
    asyn = pd.read_csv(ASYN_PATH, sep='\t')
    baseline = pd.read_csv(BASELINE_PATH, sep='\t')

    merged = asyn[['id', 'neg|score', 'neg|fdr', 'neg|rank', 'neg|lfc']].merge(
        baseline[['id', 'neg|score', 'neg|fdr', 'neg|rank', 'neg|lfc']],
        on='id', suffixes=('_asyn', '_base'))

    # Remove non-targeting controls
    merged['is_nt'] = merged['id'].str.lower().str.startswith('nt')
    targeting = merged[~merged['is_nt']].copy()

    targeting['x'] = -np.log10(targeting['neg|score_base'].clip(lower=1e-10))
    targeting['y'] = -np.log10(targeting['neg|score_asyn'].clip(lower=1e-10))
    targeting['color_val'] = -np.log10(targeting['neg|fdr_asyn'].clip(lower=1e-10))
    targeting['size_val'] = np.abs(targeting['neg|lfc_asyn']) * 800
    targeting['size_val'] = targeting['size_val'].clip(lower=15)

    print("=== Data Summary ===")
    print(f"Targeting genes: {len(targeting)}")
    stxbp1 = targeting[targeting['id'] == 'Stxbp1']
    if len(stxbp1) > 0:
        r = stxbp1.iloc[0]
        print(f"  Stxbp1: aSyn rank={r['neg|rank_asyn']:.0f} "
              f"FDR={r['neg|fdr_asyn']:.4f} LFC={r['neg|lfc_asyn']:.4f}, "
              f"Baseline rank={r['neg|rank_base']:.0f} "
              f"FDR={r['neg|fdr_base']:.4f}")

    # Plot
    fig, ax = plt.subplots(figsize=(5.5, 5))

    max_val = max(targeting['x'].max(), targeting['y'].max()) * 1.1
    ax.plot([0, max_val], [0, max_val], '--', color='#cccccc',
            linewidth=0.7, zorder=0)

    sc = ax.scatter(targeting['x'], targeting['y'],
                    s=targeting['size_val'],
                    c=targeting['color_val'],
                    cmap='viridis', alpha=0.7,
                    edgecolors='k', linewidths=0.4,
                    vmin=0, vmax=targeting['color_val'].max(),
                    zorder=2)

    for _, row in stxbp1.iterrows():
        ax.annotate('Stxbp1', xy=(row['x'], row['y']),
                    xytext=(8, 6), textcoords='offset points',
                    fontsize=FONT_SIZES['annotation'], fontweight='bold',
                    fontstyle='italic', zorder=4)

    cbar = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.8)
    cbar.set_label(r'$-\log_{10}$(FDR) aSyn', fontsize=FONT_SIZES['tick'])
    cbar.ax.tick_params(labelsize=FONT_SIZES['annotation'])

    for lfc_val, lfc_label in [(0.1, '0.1'), (0.25, '0.25'), (0.4, '0.4')]:
        ax.scatter([], [], s=lfc_val * 800, c='gray', alpha=0.5,
                   edgecolors='k', linewidths=0.4, label=f'|LFC| = {lfc_label}')

    ax.set_xlabel(r'$-\log_{10}$(neg score) — Baseline',
                  fontsize=FONT_SIZES['label'])
    ax.set_ylabel(r'$-\log_{10}$(neg score) — aSyn',
                  fontsize=FONT_SIZES['label'])
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect('equal')
    ax.legend(loc='lower right', fontsize=FONT_SIZES['annotation'],
              handletextpad=0.5, labelspacing=0.8)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.tight_layout()
    for fmt in ['png', 'svg']:
        fig.savefig(OUTPUT_DIR / f'panel_D_mageck_scatter.{fmt}',
                    format=fmt, dpi=300, bbox_inches='tight',
                    facecolor='white')
    plt.close()
    print("Saved: panel_D_mageck_scatter.png/svg")


if __name__ == '__main__':
    main()
