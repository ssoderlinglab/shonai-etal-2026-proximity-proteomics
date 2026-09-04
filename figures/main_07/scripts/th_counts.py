"""
Figure 7, panel F: TH-positive dopamine neuron counts in the ipsilateral and contralateral hemispheres, non-targeting gRNA against Stxbp1-targeting gRNA.

Input
    data/craave_seq/histology/cell_count.csv
Output
    figures/main_07/output/panel_F_th_counts.{png,svg}
    figures/main_07/output/panel_F_th_counts_mixed_model.csv
Run
    python3 figures/main_07/scripts/th_counts.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib, src/pub_style.py

Section counts are averaged within each animal before the tests drawn on the panel: a paired t-test of ipsilateral against contralateral within each gRNA group, and a Welch t-test between gRNA groups within each hemisphere.

The gRNA by side interaction reported for this panel is a linear mixed-effects model of the section-level counts, count ~ injection_side * gRNA_group with a random intercept per animal, fitted by REML over 110 sections from 24 animals. It is printed and written to panel_F_th_counts_mixed_model.csv; it is not drawn on the panel.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'craave_seq' / 'histology'
OUTPUT_DIR = FIGURE_DIR / 'output'
SRC_DIR = REPO_ROOT / 'src'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES
apply_style()

CELL_COUNT_PATH = DATA_DIR / 'cell_count.csv'

COLORS_NN = {
    'nt_contra':     '#dfc27d',
    'nt_ipsi':       '#a6611a',
    'stxbp1_contra': '#80cdc1',
    'stxbp1_ipsi':   '#018571',
}


def p_to_stars(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'n.s.'


def bracket(ax, x1, x2, y, h, text, fs=7):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color='black', lw=0.8, zorder=5, clip_on=False)
    ax.text((x1+x2)/2, y+h+1, text, ha='center', va='bottom',
            fontsize=fs, zorder=5)


def mixed_model(df):
    """Section-level linear mixed-effects model, random intercept per animal.

    count ~ injection_side * gRNA_group, fitted by REML. Reference levels are the
    alphabetically first category of each factor: contra, and the Stxbp1 gRNA group.
    """
    d = df.copy()
    d['gRNA_group'] = d['gRNA'].astype(str) + '|' + d['Strain'].astype(str)
    d['injection_side'] = d['injection_side'].astype('category')
    d['gRNA_group'] = d['gRNA_group'].astype('category')
    fit = smf.mixedlm('count ~ injection_side * gRNA_group', data=d,
                      groups=d['Sample']).fit()
    fe = fit.fe_params.index
    table = pd.DataFrame({'coef': fit.params[fe], 'std_err': fit.bse[fe],
                          'z': fit.tvalues[fe], 'p': fit.pvalues[fe]})
    table.index.name = 'term'
    # Variance components in the units of the response. MixedLMResults.params reports the
    # random-effects entry divided by the residual variance; cov_re does not.
    table.loc['animal intercept variance', 'coef'] = float(fit.cov_re.iloc[0, 0])
    table.loc['residual variance', 'coef'] = float(fit.scale)
    return fit, table


def main():
    df = pd.read_csv(CELL_COUNT_PATH)
    df_animal = (df.groupby(['Sample', 'gRNA', 'injection_side'])['count']
                 .mean().reset_index())
    df_wide = (df_animal.pivot_table(
        index=['Sample', 'gRNA'], columns='injection_side', values='count')
        .reset_index())

    nt = df_wide[df_wide['gRNA'] == 'nt'].copy()
    stxbp1 = df_wide[df_wide['gRNA'] == 'Stxbp1'].copy()

    nt['diff'] = nt['ipsi'] - nt['contra']
    stxbp1['diff'] = stxbp1['ipsi'] - stxbp1['contra']
    nt['ratio'] = nt['ipsi'] / nt['contra']
    stxbp1['ratio'] = stxbp1['ipsi'] / stxbp1['contra']

    # Statistics
    t_nt, p_nt = stats.ttest_rel(nt['ipsi'], nt['contra'])
    t_stxbp1, p_stxbp1 = stats.ttest_rel(stxbp1['ipsi'], stxbp1['contra'])
    t_contra, p_contra = stats.ttest_ind(
        nt['contra'], stxbp1['contra'], equal_var=False)
    t_ipsi, p_ipsi = stats.ttest_ind(
        nt['ipsi'], stxbp1['ipsi'], equal_var=False)

    print("=" * 60)
    print("STATISTICS")
    print("=" * 60)
    print(f"nt (n={len(nt)}):     contra={nt['contra'].mean():.1f}±"
          f"{nt['contra'].sem():.1f}, ipsi={nt['ipsi'].mean():.1f}±"
          f"{nt['ipsi'].sem():.1f}")
    print(f"Stxbp1 (n={len(stxbp1)}): contra={stxbp1['contra'].mean():.1f}±"
          f"{stxbp1['contra'].sem():.1f}, ipsi={stxbp1['ipsi'].mean():.1f}±"
          f"{stxbp1['ipsi'].sem():.1f}")
    print(f"Within: nt p={p_nt:.4f}, Stxbp1 p={p_stxbp1:.2e}")
    print(f"Between: contra p={p_contra:.4f}, ipsi p={p_ipsi:.4f}")

    fit, lme = mixed_model(df)
    terms = [t for t in lme.index if ':' in t]
    assert len(terms) == 1, terms
    term = terms[0]
    print(f"Mixed model: {int(fit.nobs)} sections, "
          f"{df['Sample'].nunique()} animals")
    print(f"  gRNA x side interaction: beta={lme.loc[term, 'coef']:+.1f}, "
          f"p={lme.loc[term, 'p']:.3f}   [{term}]")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    lme.to_csv(OUTPUT_DIR / 'panel_F_th_counts_mixed_model.csv')

    # Plot
    fig, ax = plt.subplots(figsize=(3, 4))

    bar_width = 0.35
    pair_gap = 0
    group_centers = np.array([0.5, 1.9])
    contra_x = group_centers - (bar_width / 2 + pair_gap / 2)
    ipsi_x = group_centers + (bar_width / 2 + pair_gap / 2)

    groups = [
        {'name': 'nt-gRNA', 'contra': nt['contra'], 'ipsi': nt['ipsi'],
         'n': len(nt), 'c_contra': COLORS_NN['nt_contra'],
         'c_ipsi': COLORS_NN['nt_ipsi'], 'p': p_nt},
        {'name': 'Stxbp1-KO', 'contra': stxbp1['contra'],
         'ipsi': stxbp1['ipsi'], 'n': len(stxbp1),
         'c_contra': COLORS_NN['stxbp1_contra'],
         'c_ipsi': COLORS_NN['stxbp1_ipsi'], 'p': p_stxbp1},
    ]

    max_within_bar_y = 0.0
    for i, g in enumerate(groups):
        cm, cs = g['contra'].mean(), g['contra'].sem()
        im, isem = g['ipsi'].mean(), g['ipsi'].sem()

        ax.bar(contra_x[i], cm, width=bar_width, color=g['c_contra'],
               edgecolor='none', zorder=2)
        ax.errorbar(contra_x[i], cm, yerr=cs, fmt='none', color='black',
                    capsize=4, capthick=1.0, linewidth=1.0, zorder=3)
        ax.bar(ipsi_x[i], im, width=bar_width, color=g['c_ipsi'],
               edgecolor='none', zorder=2)
        ax.errorbar(ipsi_x[i], im, yerr=isem, fmt='none', color='black',
                    capsize=4, capthick=1.0, linewidth=1.0, zorder=3)

        # Individual data points
        np.random.seed(100 + i * 2)
        jitter_c = np.random.uniform(-bar_width * 0.3, bar_width * 0.3,
                                      size=len(g['contra']))
        ax.scatter(contra_x[i] + jitter_c, g['contra'].values, s=18,
                   color='black', edgecolors='black', linewidth=0.3,
                   zorder=4, alpha=0.7)

        np.random.seed(100 + i * 2 + 1)
        jitter_i = np.random.uniform(-bar_width * 0.3, bar_width * 0.3,
                                      size=len(g['ipsi']))
        ax.scatter(ipsi_x[i] + jitter_i, g['ipsi'].values, s=18,
                   color='black', edgecolors='black', linewidth=0.3,
                   zorder=4, alpha=0.7)

        # Within-group bracket
        top = max(cm + cs, im + isem)
        bracket(ax, contra_x[i], ipsi_x[i], top + 50, 5,
                p_to_stars(g['p']))
        max_within_bar_y = max(max_within_bar_y, top + 50 + 5)

    # Between-group brackets
    y_contra_between = max_within_bar_y + 10
    bracket(ax, contra_x[0], contra_x[1], y_contra_between, 5,
            p_to_stars(p_contra))
    y_ipsi_between = y_contra_between + 5 + 14
    bracket(ax, ipsi_x[0], ipsi_x[1], y_ipsi_between, 5,
            p_to_stars(p_ipsi))

    xtick_positions = [contra_x[0], ipsi_x[0], contra_x[1], ipsi_x[1]]
    xtick_labels = ['nt | contra', 'nt | ipsi',
                    'Stxbp1-KO | contra', 'Stxbp1-KO | ipsi']
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=FONT_SIZES['tick'],
                       rotation=45, ha='right')
    ax.set_ylabel(r'$\mathrm{TH}^{+}$ Cell Count',
                  fontsize=FONT_SIZES['label'])
    ax.set_ylim(0, 380)
    ax.set_xlim(-0.1, 2.6)
    ax.set_title('DA Neuron Counts', fontsize=FONT_SIZES['title'],
                 fontweight='bold')

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.tight_layout()
    for fmt in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_F_th_counts.{fmt}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("Saved: panel_F_th_counts.png/svg")


if __name__ == '__main__':
    main()
