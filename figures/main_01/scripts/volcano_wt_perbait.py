"""
Figure 1, panels E, I and M: annotated volcano plots of the three WT proximity proteomes.

Input
    data/ibioid/contrasts/{Vps35,aSyn,Lrrk2}_WT_vs_Control.tsv
    data/ibioid/canonical_groups_rederived.json
    data/ibioid/label_positions.json
Output
    figures/main_01/output/volcano_{vps35,asyn,lrrk2}_wt.{png,svg}
Run
    python3 figures/main_01/scripts/volcano_wt_perbait.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib
"""

import argparse
import csv
import json
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -- Paths (relative to sandbox root) ---------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid')
CONTRAST_DIR = os.path.join(DATA_DIR, 'contrasts')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')

CANONICAL_JSON = os.path.join(DATA_DIR, 'canonical_groups_rederived.json')
LABEL_POS_JSON = os.path.join(DATA_DIR, 'label_positions.json')

# -- Colours (carried from v05) ---------------------------------------
C_NONSIG = '#CCCCCC'
C_PROX   = '#666666'
C_CANON1 = '#0072B2'   # Blue
C_CANON2 = '#E69F00'   # Orange
C_BAIT   = '#D55E00'   # Vermillion

# -- Thresholds -------------------------------------------------------
FC_THRESHOLD = 1.5
LOG2FC_THRESHOLD = np.log2(FC_THRESHOLD)  # 0.585
FDR_THRESHOLD = 0.05

# -- Bait configs ------
BAIT_CONFIGS = {
    'Vps35': {
        'contrast': os.path.join(CONTRAST_DIR, 'Vps35_WT_vs_Control.tsv'),
        'bait_gene': 'Vps35',
        'canon1_label': 'Retromer / WASH complex / Endosome',
        'canon2_label': 'Endocytosis / Clathrin coat',
        'title': 'Vps35 WT iBioID',
    },
    'aSyn': {
        'contrast': os.path.join(CONTRAST_DIR, 'aSyn_WT_vs_Control.tsv'),
        'bait_gene': 'Snca',
        'canon1_label': 'Presynapse',
        'canon2_label': 'RNA splicing',
        'title': 'aSyn WT iBioID',
    },
    'Lrrk2': {
        'contrast': os.path.join(CONTRAST_DIR, 'Lrrk2_WT_vs_Control.tsv'),
        'bait_gene': 'Lrrk2',
        'canon1_label': 'Rab regulation / Clathrin',
        'canon2_label': 'Actin / Microtubule',
        'title': 'Lrrk2 WT iBioID',
    },
}


def load_comparison(path):
    """Load one contrast TSV.

    Returns gene symbols, TRUE log2 fold change, and DEqMS BH-adjusted P.
    Rows without a gene symbol or without a finite statistic are skipped.
    """
    genes, log2fc, adjp = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            gene = row.get('Gene', '')
            if gene in ('', 'NA'):
                continue
            l2 = row.get('log2FC', '')
            bh = row.get('deqms_BH', '')
            if l2 in ('', 'NA') or bh in ('', 'NA'):
                continue
            l2 = float(l2)
            bh = float(bh)
            if bh <= 0:
                bh = 1e-300
            genes.append(gene)
            log2fc.append(l2)
            adjp.append(bh)
    return genes, np.array(log2fc), np.array(adjp)


def make_volcano(bait_key, canon_data, config, label_positions):
    """Generate volcano plot for one bait with fixed label positions."""
    genes, log2fc, adjp = load_comparison(config['contrast'])
    neg_log10p = -np.log10(adjp)

    canon_1st = set(canon_data['canonical_1st']['genes'])
    canon_2nd = set(canon_data['canonical_2nd']['genes'])
    labeled   = set(canon_data['labeled']['genes'])
    bait_gene = config['bait_gene']

    # Assign categories (same rules as v05, gate now = log2FC>0.585 & BH<0.05)
    colors, sizes, zorders, markers = [], [], [], []
    for i, gene in enumerate(genes):
        sig = (log2fc[i] > LOG2FC_THRESHOLD) and (adjp[i] < FDR_THRESHOLD)
        if gene == bait_gene:
            colors.append(C_BAIT); sizes.append(150); zorders.append(50); markers.append('*')
        elif gene in canon_1st and sig:
            colors.append(C_CANON1); sizes.append(30); zorders.append(30); markers.append('o')
        elif gene in canon_2nd and sig:
            colors.append(C_CANON2); sizes.append(30); zorders.append(30); markers.append('o')
        elif sig:
            colors.append(C_PROX); sizes.append(12); zorders.append(10); markers.append('o')
        else:
            colors.append(C_NONSIG); sizes.append(8); zorders.append(1); markers.append('o')

    # -- Plot ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    x_margin = 0.5
    x_min = min(log2fc) - x_margin
    x_max = max(log2fc) + x_margin
    y_max = min(max(neg_log10p) + 0.5, 12)
    plot_y = np.clip(neg_log10p, 0, y_max)

    # Plot by layer (background first)
    for z_val in sorted(set(zorders)):
        idx = [i for i in range(len(genes)) if zorders[i] == z_val]
        for i in idx:
            ax.scatter(log2fc[i], plot_y[i],
                       c=colors[i], s=sizes[i], marker=markers[i],
                       edgecolors='none' if markers[i] == 'o' else 'black',
                       linewidths=0.5 if markers[i] == '*' else 0,
                       alpha=0.85 if z_val >= 30 else 0.4,
                       zorder=z_val, rasterized=(z_val < 10))

    # Threshold lines
    ax.axvline(x=LOG2FC_THRESHOLD, color='grey', linestyle='--', linewidth=0.5, alpha=0.5)
    ax.axhline(y=-np.log10(FDR_THRESHOLD), color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    # -- Labels from JSON (fixed positions) ---------------------------
    gene_idx = {}
    for i, gene in enumerate(genes):
        if gene in labeled or gene == bait_gene:
            gene_idx[gene] = i

    for gene, pos in label_positions.items():
        if gene not in gene_idx:
            continue
        i = gene_idx[gene]

        if gene == bait_gene:
            tc = C_BAIT
        elif gene in canon_1st:
            tc = C_CANON1
        elif gene in canon_2nd:
            tc = C_CANON2
        else:
            tc = C_PROX

        x_data = log2fc[i]
        y_data = plot_y[i]
        x_label = pos['x_label']
        y_label = pos['y_label']

        dx = x_label - x_data
        dy = y_label - y_data
        dist = np.sqrt(dx**2 + dy**2)
        if dist > 0.15:
            ax.annotate('', xy=(x_data, y_data), xytext=(x_label, y_label),
                        arrowprops=dict(arrowstyle='-', color='grey',
                                        lw=0.3, alpha=0.6),
                        zorder=35)

        ax.text(x_label, y_label, gene,
                fontsize=8.5, fontweight='bold', color=tc,
                ha='center', va='bottom', zorder=40)

    # Axes
    ax.set_xlabel('log$_2$ (Fold Change)', fontsize=15)
    ax.set_ylabel(u'\u2212log$_{10}$ (adj. P)', fontsize=15)
    ax.set_title(config['title'], fontsize=18, fontweight='bold')
    ax.tick_params(labelsize=9)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(-0.2, y_max)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    gene_set = set(genes)
    leg_scale = 1.5
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_NONSIG,
               markersize=5 * leg_scale, label='Non-significant'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_PROX,
               markersize=6 * leg_scale, label='Proximity enriched'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_CANON1,
               markersize=7 * leg_scale, label=f'{config["canon1_label"]} (n={len(canon_1st & gene_set)})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=C_CANON2,
               markersize=7 * leg_scale, label=f'{config["canon2_label"]} (n={len(canon_2nd & gene_set)})'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor=C_BAIT,
               markeredgecolor='black', markeredgewidth=0.3 * leg_scale,
               markersize=10 * leg_scale, label=f'Bait ({bait_gene})'),
    ]
    leg = ax.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=7 * leg_scale,
        framealpha=0.9,
        fancybox=False,
    )
    leg.get_frame().set_edgecolor('none')
    leg.get_frame().set_linewidth(0)

    plt.tight_layout()

    # Save
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    bait_lower = bait_key.lower()
    png_path = os.path.join(OUTPUT_DIR, f'volcano_{bait_lower}_wt.png')
    svg_path = os.path.join(OUTPUT_DIR, f'volcano_{bait_lower}_wt.svg')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close()

    # Stats
    sig_genes = set(g for i, g in enumerate(genes)
                    if log2fc[i] > LOG2FC_THRESHOLD and adjp[i] < FDR_THRESHOLD)
    print(f'\n{"="*60}')
    print(f'  {bait_key} WT')
    print(f'{"="*60}')
    print(f'  Saved: {png_path}')
    print(f'  Total plotted (DEqMS-tested): {len(genes)}')
    print(f'  Proxiome enriched (sig): {len(sig_genes)}')
    print(f'  Canonical 1st in proxiome: {len(canon_1st & sig_genes)}')
    print(f'  Canonical 2nd in proxiome: {len(canon_2nd & sig_genes)}')
    print(f'  Labels placed: {len([g for g in label_positions if g in gene_idx])}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Volcano: -log10(adj.P) y-axis, labels from label_positions.json.'
    )
    parser.add_argument(
        '--only',
        choices=sorted(BAIT_CONFIGS.keys()),
        metavar='BAIT',
        help='Regenerate only this panel (Vps35, aSyn, or Lrrk2).',
    )
    args = parser.parse_args()

    with open(CANONICAL_JSON) as f:
        all_canon = json.load(f)
    with open(LABEL_POS_JSON) as f:
        all_label_pos = json.load(f)

    to_run = list(BAIT_CONFIGS.items())
    if args.only:
        to_run = [(args.only, BAIT_CONFIGS[args.only])]

    for bait_key, config in to_run:
        label_positions = all_label_pos.get(bait_key, {})
        make_volcano(bait_key, all_canon[bait_key], config, label_positions)

    print(f'\nLabel positions: {LABEL_POS_JSON}')
    print('To adjust a label, edit x_label / y_label in the JSON and re-run.')
