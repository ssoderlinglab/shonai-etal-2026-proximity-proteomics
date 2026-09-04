"""
Supplemental Figures 1, 2 and 3, panel C: replicate scatter of each WT proximity proteome.

Input
    data/ibioid/replicates/{Vps35,aSyn,Lrrk2_WT}_Protein_Normalized.csv
Output
    figures/suppl_0{1,2,3}/output/panel_C_replicate_scatter_<bait>_WT.{png,svg}
Run
    python3 figures/suppl_01/scripts/replicate_scatter.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, scipy

Called with one bait. Intensities are log2 transformed, zeros dropped; the annotation reports the Pearson R squared of replicate 1 against replicate 2 and the mean of the three pairwise values.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

SRC_DIR = Path(__file__).resolve().parent
REPO_ROOT = SRC_DIR.parent
DATA_DIR = REPO_ROOT / 'data' / 'ibioid' / 'replicates'
FIGURE_DIR = REPO_ROOT / 'figures'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES
apply_style()

# WT conditions: file, replicate columns, figure directory. Panel C of each
# per-bait supplemental figure. Colours are the ones the published panels carry.
CONDITIONS = {
    'aSyn': {
        'file': DATA_DIR / 'aSyn_Protein_Normalized.csv',
        'reps': ['WT_1', 'WT_2', 'WT_3'],
        'color': '#1F77B4',
        'panel': 'C',
        'figure': 'suppl_02',
    },
    'Vps35': {
        'file': DATA_DIR / 'Vps35_Protein_Normalized.csv',
        'reps': ['WT_1', 'WT_2', 'WT_3'],
        'color': '#2CA02C',
        'panel': 'C',
        'figure': 'suppl_01',
    },
    'Lrrk2': {
        'file': DATA_DIR / 'Lrrk2_WT_Protein_Normalized.csv',
        'reps': ['WT_1', 'WT_2', 'WT_3'],
        'color': '#D62728',
        'panel': 'C',
        'figure': 'suppl_03',
    },
}


def load_data(filepath, cols):
    """Load and log2-transform normalized protein data."""
    df = pd.read_csv(filepath)
    data = df[cols].copy()
    data = data.replace(0, np.nan)
    data = np.log2(data)
    data = data.dropna()
    return data


def calc_r_squared(x, y):
    """Calculate Pearson R²."""
    mask = ~(np.isnan(x) | np.isnan(y))
    r, _ = stats.pearsonr(x[mask], y[mask])
    return r ** 2


def plot_scatter(name, config, output_dir):
    """Generate the replicate scatter panel for one WT condition."""
    data = load_data(config['file'], config['reps'])
    reps = config['reps']

    # Rep1 vs Rep2 for scatter
    x = data[reps[0]].values
    y = data[reps[1]].values

    # All pairwise R²
    r2_12 = calc_r_squared(data[reps[0]], data[reps[1]])
    r2_23 = calc_r_squared(data[reps[1]], data[reps[2]])
    r2_13 = calc_r_squared(data[reps[0]], data[reps[2]])
    r2_mean = (r2_12 + r2_23 + r2_13) / 3

    fig, ax = plt.subplots(figsize=(4, 4))
    ax.scatter(x, y, c=config['color'], s=5, alpha=0.3, rasterized=True)

    # Diagonal
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel(f'{reps[0]} (log2)', fontsize=FONT_SIZES['label'])
    ax.set_ylabel(f'{reps[1]} (log2)', fontsize=FONT_SIZES['label'])
    ax.set_title(f'{name} WT', fontsize=FONT_SIZES['title'], fontweight='bold')
    ax.text(0.05, 0.95,
            f'R\u00b2 = {r2_12:.3f}\n(mean = {r2_mean:.3f})',
            transform=ax.transAxes, fontsize=FONT_SIZES['tick'],
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    ax.set_aspect('equal', adjustable='box')

    os.makedirs(output_dir, exist_ok=True)
    fig.tight_layout()
    for ext in ['svg', 'png']:
        fig.savefig(output_dir / f'panel_C_replicate_scatter_{name}_WT.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    return {
        'Bait': name,
        'Panel': config['panel'],
        'n_proteins': len(data),
        'R2_Rep1v2': r2_12,
        'R2_Rep2v3': r2_23,
        'R2_Rep1v3': r2_13,
        'R2_mean': r2_mean,
    }


def main(bait):
    """Build the replicate scatter panel of one bait's supplemental figure."""
    config = CONDITIONS[bait]
    output_dir = FIGURE_DIR / config['figure'] / 'output'
    r = plot_scatter(bait, config, output_dir)
    print(f"{bait} WT  panel {r['Panel']}  n={r['n_proteins']} proteins")
    print(f"  R2 rep1-rep2 {r['R2_Rep1v2']:.4f}  rep2-rep3 {r['R2_Rep2v3']:.4f}  "
          f"rep1-rep3 {r['R2_Rep1v3']:.4f}  mean {r['R2_mean']:.4f}")
    print(f"  wrote {output_dir}")
    return r
