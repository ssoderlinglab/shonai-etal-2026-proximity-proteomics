"""
Figure 7, panel G: ipsilateral minus contralateral TH-positive cell count difference for the four experimental conditions.

Input
    data/craave_seq/histology/cell_count.csv
    data/craave_seq/histology/ilastik_left_right_summary.csv
    data/craave_seq/histology/L_R_paired_comparison.csv
Output
    figures/main_07/output/panel_G_th_count_difference.{png,svg}
    figures/main_07/output/panel_G_th_count_difference.csv
Run
    python3 figures/main_07/scripts/th_count_difference.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib, src/pub_style.py

Section counts are averaged within each animal, then one one-sample t-test of the per-animal difference against zero is run per condition and the four p-values are Benjamini-Hochberg corrected. The stars drawn on the panel are the corrected values.

Conditions and sample sizes: non-targeting gRNA with alpha-synuclein (n = 8), Stxbp1 gRNA with alpha-synuclein (n = 16), Stxbp1 gRNA with GFP (n = 6), alpha-synuclein with and without tetanus toxin light chain (n = 6).
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'craave_seq' / 'histology'
SRC_DIR = REPO_ROOT / 'src'
sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES  # noqa: E402
apply_style()

OUTPUT_DIR = FIGURE_DIR / 'output'

CELL_COUNT_PATH = DATA_DIR / 'cell_count.csv'
GFP_PATH = DATA_DIR / 'ilastik_left_right_summary.csv'
TENT_PATH = DATA_DIR / 'L_R_paired_comparison.csv'

COLORS = {'nt': '#dfc27d', 'stxbp1': '#018571',
          'gfp': '#80cdc1', 'tent': '#a6611a'}


def q_to_stars(q):
    if pd.isna(q):
        return 'n.s.'
    if q < 0.001:
        return '***'
    if q < 0.01:
        return '**'
    if q < 0.05:
        return '*'
    return 'n.s.'


def main():
    # 1 & 2
    df = pd.read_csv(CELL_COUNT_PATH)
    df_animal = (df.groupby(['Sample', 'gRNA', 'injection_side'])['count']
                 .mean().reset_index())
    df_wide = (df_animal.pivot_table(
        index=['Sample', 'gRNA'], columns='injection_side', values='count')
        .reset_index())
    nt = df_wide[df_wide['gRNA'] == 'nt'].copy()
    stxbp1 = df_wide[df_wide['gRNA'] == 'Stxbp1'].copy()
    diff_nt = (nt['ipsi'] - nt['contra']).values
    diff_stxbp1 = (stxbp1['ipsi'] - stxbp1['contra']).values

    # 3 GFP
    gfp_raw = pd.read_csv(GFP_PATH)
    gfp_raw = gfp_raw[gfp_raw['Animal'].notna()
                       & (gfp_raw['Animal'] != 'NA')].copy()
    gfp_raw['Animal'] = gfp_raw['Animal'].astype(int)
    gfp_animal = (gfp_raw.groupby(['Animal', 'HemisphereSide'])
                  ['cells_predicted'].mean().reset_index())
    gfp_wide = (gfp_animal.pivot_table(
        index='Animal', columns='HemisphereSide',
        values='cells_predicted').reset_index())
    diff_gfp = (gfp_wide['Right'] - gfp_wide['Left']).values

    # 4 TeNT
    tent_raw = pd.read_csv(TENT_PATH)
    tent_animal = (tent_raw.groupby('AnimalID')
                   [['neuron_count_L', 'neuron_count_R']].mean().reset_index())
    diff_tent = (tent_animal['neuron_count_R']
                 - tent_animal['neuron_count_L']).values

    conditions = [
        {'label': 'nt + \u03b1-syn', 'diff': diff_nt, 'color': COLORS['nt']},
        {'label': 'Stxbp1-KO + \u03b1-syn', 'diff': diff_stxbp1,
         'color': COLORS['stxbp1']},
        {'label': 'Stxbp1-KO + GFP', 'diff': diff_gfp, 'color': COLORS['gfp']},
        {'label': '\u03b1-syn \u00b1 TeNT', 'diff': diff_tent,
         'color': COLORS['tent']},
    ]

    pvals = []
    for c in conditions:
        c['n'] = len(c['diff'])
        t, p = stats.ttest_1samp(c['diff'], 0)
        c['t'] = float(t)
        c['p'] = float(p)
        c['mean'] = float(c['diff'].mean())
        c['sem'] = float(c['diff'].std(ddof=1) / np.sqrt(c['n']))
        pvals.append(p)

    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    for c, q in zip(conditions, qvals):
        c['q'] = float(q)
        c['marker'] = q_to_stars(q)

    print("Panel G, four one-sample t-tests, Benjamini-Hochberg corrected:")
    for c in conditions:
        print(f"  {c['label']} (n={c['n']}): mean={c['mean']:.1f}±{c['sem']:.1f}"
              f" p={c['p']:.4f} q={c['q']:.4f} {c['marker']}")

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Plot
    fig, ax = plt.subplots(figsize=(4, 4))
    positions = [0, 1, 2, 3]
    bp = ax.boxplot(
        [c['diff'] for c in conditions], positions=positions, widths=0.5,
        patch_artist=True, showfliers=False, zorder=2,
        medianprops=dict(color='black', linewidth=1.2),
        whiskerprops=dict(color='gray', linewidth=0.8),
        capprops=dict(color='gray', linewidth=0.8))
    for patch, c in zip(bp['boxes'], conditions):
        patch.set_facecolor(c['color'])
        patch.set_edgecolor('black')
        patch.set_linewidth(0.8)

    for i, c in enumerate(conditions):
        np.random.seed(500 + i)
        jitter = np.random.uniform(-0.12, 0.12, size=c['n'])
        ax.scatter(positions[i] + jitter, c['diff'], s=20, color='black',
                   edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)

    ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.6,
               alpha=0.5, zorder=1)

    for i, c in enumerate(conditions):
        top = max(c['diff'])
        ax.text(positions[i], top + 15, c['marker'],
                ha='center', va='bottom', fontsize=7, zorder=5)

    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{c['label']}\n(n={c['n']})" for c in conditions],
        fontsize=FONT_SIZES['tick'], rotation=45, ha='right')
    ax.set_ylabel(r'$\mathrm{TH}^{+}$ ipsi $-$ contra',
                  fontsize=FONT_SIZES['label'])
    ax.set_title(r'$\mathbf{TH}^{+}$ Count Difference',
                 fontsize=FONT_SIZES['title'], fontweight='bold', pad=12)

    fig.tight_layout()
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_G_th_count_difference.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    pd.DataFrame([{k: v for k, v in c.items() if k != 'diff'}
                  for c in conditions]).to_csv(
        OUTPUT_DIR / 'panel_G_th_count_difference.csv', index=False)
    print("Saved figure + csv")


if __name__ == '__main__':
    main()
