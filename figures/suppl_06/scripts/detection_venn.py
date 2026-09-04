"""
Supplemental Figure 6, panel A: three-way Venn of proteins detected across the three co-fractionation cohorts.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
Output
    figures/suppl_06/output/panel_A_detection_venn.{png,svg}
Run
    python3 figures/suppl_06/scripts/detection_venn.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib, matplotlib-venn

Cohort sizes are asserted against the values reported in the manuscript: 1939, 678, 852, and 409 detected in all three.
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibiocofrac')
OUTDIR = os.path.join(FIGURE_DIR, 'output')

sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))
from contamination_filter import CONTAMINATION_GENES as EXCLUDE_GENES

COL_LRRK2 = '#1B9E77'
COL_SNCA = '#D95F02'
COL_VPS35 = '#7570B3'

EXPECTED = {'LRRK2': 1939, 'SNCA': 678, 'VPS35': 852, 'shared': 409}


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    partition = pd.read_csv(os.path.join(
        DATA_DIR, 'modules', 'Partition_Gene_combined_3PD_models.csv'))

    sets = {m: set(partition[partition['Experiment'] == m]['Symbol'].dropna())
            - set(EXCLUDE_GENES)
            for m in ('LRRK2', 'SNCA', 'VPS35')}
    lrrk2, snca, vps35 = sets['LRRK2'], sets['SNCA'], sets['VPS35']
    three_way = lrrk2 & snca & vps35

    # Guard the headline numbers. If the upstream partition map or the contamination exclusion changes,
    # this fails loudly rather than shipping a panel that disagrees with the text.
    got = {'LRRK2': len(lrrk2), 'SNCA': len(snca), 'VPS35': len(vps35),
           'shared': len(three_way)}
    assert got == EXPECTED, f'counts changed: expected {EXPECTED}, got {got}'

    fig, ax = plt.subplots(figsize=(7, 6))
    v = venn3([lrrk2, snca, vps35], set_labels=('', '', ''),
              set_colors=(COL_LRRK2, COL_SNCA, COL_VPS35), alpha=0.35, ax=ax)

    for pid in ['100', '010', '001', '110', '101', '011', '111']:
        patch = v.get_patch_by_id(pid)
        if patch:
            patch.set_edgecolor('white')
            patch.set_linewidth(1.5)

    if v.get_label_by_id('111'):
        v.get_label_by_id('111').set_fontsize(14)
        v.get_label_by_id('111').set_fontweight('bold')
    for pid in ['100', '010', '001', '110', '101', '011']:
        lbl = v.get_label_by_id(pid)
        if lbl:
            lbl.set_fontsize(10)

    ax.text(0.08, 0.92, f'LRRK2 (Tg-G2019S): {len(lrrk2)} proteins',
            transform=ax.transAxes, fontsize=10, color=COL_LRRK2,
            fontweight='bold')
    ax.text(0.08, 0.86, f'SNCA (Tg-A30P*A53T): {len(snca)} proteins',
            transform=ax.transAxes, fontsize=10, color=COL_SNCA,
            fontweight='bold')
    ax.text(0.08, 0.80, f'VPS35 (KI-D620N): {len(vps35)} proteins',
            transform=ax.transAxes, fontsize=10, color=COL_VPS35,
            fontweight='bold')
    ax.text(0.5, 0.02,
            f'3-way overlap: {len(three_way)} proteins '
            f'\u2192 used for cross-experiment comparison',
            transform=ax.transAxes, fontsize=10, ha='center', va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#F0F0F0',
                      edgecolor='#999999'))
    ax.set_title('Protein detection overlap', fontsize=13, fontweight='bold')

    fig.tight_layout()
    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(OUTDIR, f'panel_A_detection_venn.{ext}'),
                    dpi=200, bbox_inches='tight', transparent=True)
    plt.close(fig)

    print('wrote panel_A_detection_venn.png/.svg')
    print(f'  per-model (the contamination exclusion applied): {got["LRRK2"]} / {got["SNCA"]} / '
          f'{got["VPS35"]}')
    print(f'  3-way overlap: {got["shared"]}')
    removed = sorted(
        (set(partition[partition['Experiment'] == 'LRRK2']['Symbol'].dropna())
         & set(partition[partition['Experiment'] == 'SNCA']['Symbol'].dropna())
         & set(partition[partition['Experiment'] == 'VPS35']['Symbol'].dropna()))
        & set(EXCLUDE_GENES))
    print(f'  the contamination exclusion removed {len(removed)} from the 3-way set: '
          f'{", ".join(removed)}')


if __name__ == '__main__':
    main()
