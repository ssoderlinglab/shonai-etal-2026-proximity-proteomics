"""
Figure 5, panel D: compartment marker enrichment of the co-fractionation modules.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/markers/Organelle_Marker_Reference.csv
Output
    figures/main_05/output/panel_D_enrichment_heatmap.{png,svg,csv}
Run
    python3 figures/main_05/scripts/compartment_marker_heatmap.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib

Rows are the 12 selected modules plus M13. Significance is the raw one-sided Fisher p, in three tiers (0.05, 0.01, 0.001), without multiple-testing correction.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibiocofrac')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))
PATH_PARTITION = os.path.join(DATA_DIR, 'modules',
                              'Partition_Gene_combined_3PD_models.csv')
PATH_MARKER_REF = os.path.join(DATA_DIR, 'markers',
                               'Organelle_Marker_Reference.csv')

from contamination_filter import CONTAMINATION_GENES as EXCLUDE_GENES

SELECTED_12 = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 14, 19]
PANEL_B_MODULES = sorted(SELECTED_12 + [13])
SYNAPSE_MODULES = {7, 8, 11, 19}
COL_SYNAPSE = '#E15759'

COMPARTMENT_ORDER = [
    'mitochondria', 'ribosome', 'ER', 'Golgi/ERGIC',
    'lysosome', 'early_endosome', 'plasma membrane',
    'cytosol', 'nuclear_pore',
    'Synaptic_Vesicle', 'Presynaptic_Active_Zone', 'Presynapse',
    'Postsynaptic_Density', 'Postsynapse',
]
COMPARTMENT_DISPLAY = {
    'mitochondria': 'Mito', 'ribosome': 'Ribo', 'ER': 'ER',
    'Golgi/ERGIC': 'Golgi/ERGIC', 'lysosome': 'Lyso',
    'early_endosome': 'EE', 'plasma membrane': 'PM',
    'cytosol': 'Cyto', 'nuclear_pore': 'Nuc',
    'Synaptic_Vesicle': 'SV', 'Presynaptic_Active_Zone': 'PreAZ',
    'Presynapse': 'Pre', 'Postsynaptic_Density': 'PSD',
    'Postsynapse': 'Post',
}


def main():
    print("Loading data...")
    partition = pd.read_csv(PATH_PARTITION)
    marker_ref = pd.read_csv(PATH_MARKER_REF)

    lrrk2_part = partition[partition['Experiment'] == 'LRRK2'][
        ['UniProt', 'Symbol', 'Membership']].copy()
    lrrk2_part.columns = ['Accession', 'Gene', 'Module']
    lrrk2_part = lrrk2_part[~lrrk2_part['Gene'].isin(EXCLUDE_GENES)]

    lrrk2_markers = marker_ref[marker_ref['Model'] == 'LRRK2'].copy()
    lrrk2_markers['Compartment'] = lrrk2_markers['Compartment'].replace({
        'Golgi': 'Golgi/ERGIC', 'ERGIC': 'Golgi/ERGIC'
    })
    available_comp = [c for c in COMPARTMENT_ORDER
                      if c in lrrk2_markers['Compartment'].unique()]

    lrrk2_all_genes = set(lrrk2_part['Gene'])
    N_total = len(lrrk2_all_genes)
    n_mods = len(PANEL_B_MODULES)
    n_comp = len(available_comp)

    pval_matrix = np.ones((n_mods, n_comp))
    sig_matrix = np.full((n_mods, n_comp), '', dtype=object)

    for mi, mod in enumerate(PANEL_B_MODULES):
        mod_genes = set(lrrk2_part[lrrk2_part['Module'] == mod]['Gene'])
        n_mod = len(mod_genes)
        for ci, comp in enumerate(available_comp):
            comp_genes = (set(lrrk2_markers[
                lrrk2_markers['Compartment'] == comp]['Gene_Name'])
                & lrrk2_all_genes)
            K = len(comp_genes)
            k = len(mod_genes & comp_genes)
            a, b, c, d = k, K - k, n_mod - k, N_total - K - (n_mod - k)
            d = max(0, d)
            try:
                _, p = stats.fisher_exact([[a, b], [c, d]],
                                          alternative='greater')
            except Exception:
                p = 1.0
            pval_matrix[mi, ci] = p
            if p < 0.001:
                sig_matrix[mi, ci] = '***'
            elif p < 0.01:
                sig_matrix[mi, ci] = '**'
            elif p < 0.05:
                sig_matrix[mi, ci] = '*'

    # Tiered: 0=non-sig, 1=p<0.05, 2=p<0.01, 3=p<0.001
    tier_matrix = np.zeros((n_mods, n_comp))
    tier_matrix[pval_matrix < 0.05]  = 1
    tier_matrix[pval_matrix < 0.01]  = 2
    tier_matrix[pval_matrix < 0.001] = 3

    print(f"Total tests: {n_mods * n_comp}")
    print(f"  p<0.001: {(pval_matrix < 0.001).sum()}")
    print(f"  p<0.01:  {(pval_matrix < 0.01).sum()}")
    print(f"  p<0.05:  {(pval_matrix < 0.05).sum()}")

    tier_colors = ['#FFFFFF', '#FEB24C', '#F03B20', '#BD0026']
    cmap_tier = ListedColormap(tier_colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
    norm_tier = BoundaryNorm(bounds, cmap_tier.N)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.imshow(tier_matrix, cmap=cmap_tier, norm=norm_tier, aspect='auto')

    for mi in range(n_mods):
        for ci in range(n_comp):
            sig = sig_matrix[mi, ci]
            if sig:
                tier = tier_matrix[mi, ci]
                txt_color = 'white' if tier >= 3 else 'black'
                ax.text(ci, mi, sig, ha='center', va='center',
                        fontsize=8, fontweight='bold', color=txt_color)

    ylabels = [f'M{mod} (n={len(lrrk2_part[lrrk2_part["Module"]==mod])})'
               for mod in PANEL_B_MODULES]
    ax.set_yticks(range(n_mods))
    ax.set_yticklabels(ylabels, fontsize=9)
    for i, mod in enumerate(PANEL_B_MODULES):
        if mod in SYNAPSE_MODULES:
            ax.get_yticklabels()[i].set_color(COL_SYNAPSE)
            ax.get_yticklabels()[i].set_fontweight('bold')

    xlabels = [COMPARTMENT_DISPLAY.get(c, c) for c in available_comp]
    ax.set_xticks(range(n_comp))
    ax.set_xticklabels(xlabels, fontsize=9, rotation=45, ha='right')

    legend_elements = [
        Patch(facecolor='#FEB24C', edgecolor='#999', label='p < 0.05 (*)'),
        Patch(facecolor='#F03B20', edgecolor='#999', label='p < 0.01 (**)'),
        Patch(facecolor='#BD0026', edgecolor='#999', label='p < 0.001 (***)'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='upper right')
    ax.set_title('Compartment Marker Enrichment',
                 fontsize=12, fontweight='bold')

    plt.tight_layout()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for fmt in ['png', 'svg']:
        out = os.path.join(OUTPUT_DIR,
                           f'panel_D_enrichment_heatmap.{fmt}')
        fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"  Saved: {out}")
    plt.close()

    rows = []
    for mi, mod in enumerate(PANEL_B_MODULES):
        for ci, comp in enumerate(available_comp):
            rows.append({
                'Module': f'M{mod}',
                'Compartment': comp,
                'p_value': pval_matrix[mi, ci],
                'tier': int(tier_matrix[mi, ci]),
                'marker': sig_matrix[mi, ci],
            })
    df = pd.DataFrame(rows)
    csv_path = os.path.join(OUTPUT_DIR,
                            'panel_D_enrichment_heatmap.csv')
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")


if __name__ == '__main__':
    main()
