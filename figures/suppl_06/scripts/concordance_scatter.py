"""
Supplemental Figure 6, panel C: protein-level fold-change concordance between cohorts, for the 409 proteins detected in all three.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
Output
    figures/suppl_06/output/panel_C_concordance_scatter.{png,svg}
    figures/suppl_06/output/panel_C_concordance.csv
Run
    python3 figures/suppl_06/scripts/concordance_scatter.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Three pairwise binomial tests of concordant direction, Benjamini-Hochberg corrected across the three.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
from statsmodels.stats.multitest import multipletests
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
np.random.seed(42)

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'ibiocofrac'
SRC_DIR = REPO_ROOT / 'src'
sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES  # noqa: E402
apply_style()

OUTPUT_DIR = FIGURE_DIR / 'output'

PATH_PARTITION = DATA_DIR / 'modules' / 'Partition_Gene_combined_3PD_models.csv'
PATH_FRAC = {
    'LRRK2': DATA_DIR / 'fractionation' / 'LRRK2_young_for_PCA.csv',
    'SNCA':  DATA_DIR / 'fractionation' / 'SNCA_young_for_PCA.csv',
    'VPS35': DATA_DIR / 'fractionation' / 'Vps35_young_for_PCA.csv',
}

COL_LRRK2 = '#1B9E77'
COL_SNCA  = '#D95F02'
COL_VPS35 = '#7570B3'
EXP_COLORS = {'LRRK2': COL_LRRK2, 'SNCA': COL_SNCA, 'VPS35': COL_VPS35}
GENO_SHORT = {
    'LRRK2': 'Tg-LRRK2 G2019S',
    'SNCA':  'Tg-SNCA A30P*A53T',
    'VPS35': 'KI-VPS35 D620N',
}

FRACTIONS_STR = ['4', '5', '6', '7', '8', '9', '10']

from contamination_filter import CONTAMINATION_GENES as EXCLUDE_GENES

COL_SYNAPSE = '#E15759'
COL_RIBO = '#59A14F'
COL_MITO = '#EDC948'
COL_OTHER = '#5B7BA5'

GROUPS = [
    ('Synapse',  'Synapse\n(M7,M8,M11,M19)',  [7, 8, 11, 19], COL_SYNAPSE),
    ('Ribo',     'Ribosome\n(M1,M10)',         [1, 10],         COL_RIBO),
    ('Mito',     'Mitochondria\n(M2,M5)',      [2, 5],          COL_MITO),
    ('M4',       'M4 (ER/Golgi)',              [4],             COL_OTHER),
    ('M6',       'M6 (Clath/tRNA)',            [6],             COL_OTHER),
    ('M9',       'M9 (Proteasome)',            [9],             COL_OTHER),
    ('M14',      'M14 (RNA-bind/ND)',          [14],            COL_OTHER),
]


def compute_protein_fc(frac_df):
    df = frac_df[frac_df['fraction'].astype(str).isin(FRACTIONS_STR)].copy()
    wt = df[df['sample'].str.startswith('wt_')]
    tg = df[df['sample'].str.startswith('tg_')]
    meta_cols = [df.columns[0], 'fraction', 'replicate', 'sample']
    prot_cols = [c for c in df.columns if c not in meta_cols]
    rows = []
    for prot in prot_cols:
        wt_mean = np.nanmean(wt[prot].values.astype(float))
        tg_mean = np.nanmean(tg[prot].values.astype(float))
        fc = (np.log2(tg_mean / wt_mean)
              if wt_mean > 0 and tg_mean > 0 else np.nan)
        rows.append({'Accession': prot, 'log2FC': fc})
    return pd.DataFrame(rows)


def load_all():
    partition = pd.read_csv(PATH_PARTITION)
    frac_data = {e: pd.read_csv(p) for e, p in PATH_FRAC.items()}
    protein_fc = {}
    for exp in ['LRRK2', 'SNCA', 'VPS35']:
        protein_fc[exp] = compute_protein_fc(frac_data[exp])
        a2g = dict(zip(partition[partition['Experiment'] == exp]['UniProt'],
                       partition[partition['Experiment'] == exp]['Symbol']))
        protein_fc[exp]['Gene'] = protein_fc[exp]['Accession'].map(a2g)
    merged = None
    for exp in ['LRRK2', 'SNCA', 'VPS35']:
        df = (protein_fc[exp][['Gene', 'log2FC']].dropna(subset=['Gene'])
              .rename(columns={'log2FC': f'{exp}_FC'})
              .drop_duplicates(subset='Gene', keep='first'))
        merged = df if merged is None else pd.merge(merged, df, on='Gene', how='outer')
    merged = merged[~merged['Gene'].isin(EXCLUDE_GENES)]
    return partition, merged


def q_marker(p, q, allow_nominal=True):
    if pd.isna(q):
        return ''
    if q < 0.001:
        return '***'
    if q < 0.01:
        return '**'
    if q < 0.05:
        return '*'
    if allow_nominal and not pd.isna(p) and p < 0.05:
        return '(*)'
    return ''


# ================ Panel B: Concordance with FDR ================
def panel_c_concordance(merged):
    in3 = merged.dropna(subset=['LRRK2_FC', 'SNCA_FC', 'VPS35_FC']).copy()
    pairs = [
        ('LRRK2_FC', 'VPS35_FC', 'LRRK2', 'VPS35'),
        ('LRRK2_FC', 'SNCA_FC',  'LRRK2', 'SNCA'),
        ('SNCA_FC',  'VPS35_FC', 'SNCA',  'VPS35'),
    ]

    rows = []
    pvals = []
    for xc, yc, xl, yl in pairs:
        x, y = in3[xc].values, in3[yc].values
        conc = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
        n_conc = int(conc.sum())
        r_val, _ = stats.pearsonr(x, y)
        rho, _ = stats.spearmanr(x, y)
        btest = stats.binomtest(n_conc, len(x), 0.5, alternative='greater')
        rows.append({
            'pair': f'{xl} vs {yl}',
            'n': len(x),
            'n_concordant': n_conc,
            'percent_concordant': 100 * n_conc / len(x),
            'pearson_r': r_val,
            'spearman_rho': rho,
            'binomial_p': btest.pvalue,
        })
        pvals.append(btest.pvalue)

    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    for row, q in zip(rows, qvals):
        row['binomial_q_BH'] = q
        row['fdr_marker'] = q_marker(row['binomial_p'], q, allow_nominal=False)

    print("Panel C (3 binomial tests):")
    for r in rows:
        print(f"  {r['pair']}: p={r['binomial_p']:.3e} q={r['binomial_q_BH']:.3e} {r['fdr_marker']}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i, (xc, yc, xl, yl) in enumerate(pairs):
        ax = axes[i]
        x, y = in3[xc].values, in3[yc].values
        conc = ((x > 0) & (y > 0)) | ((x < 0) & (y < 0))
        ax.scatter(x[~conc], y[~conc], s=10, alpha=0.25, c='#BBBBBB',
                   edgecolors='none', zorder=2)
        ax.scatter(x[conc], y[conc], s=10, alpha=0.4, c='#2C3E50',
                   edgecolors='none', zorder=3)
        slope, intercept = np.polyfit(x, y, 1)
        xline = np.linspace(x.min(), x.max(), 100)
        ax.plot(xline, slope * xline + intercept, color='#E74C3C',
                linewidth=1.5, alpha=0.8, zorder=4)
        ax.axhline(0, color='gray', linewidth=0.3)
        ax.axvline(0, color='gray', linewidth=0.3)
        r = rows[i]
        ax.text(0.03, 0.97,
                f"n = {r['n']}\n"
                f"Pearson r = {r['pearson_r']:.3f}\n"
                f"Spearman \u03c1 = {r['spearman_rho']:.3f}\n"
                f"Concordant: {r['n_concordant']}/{r['n']} "
                f"({r['percent_concordant']:.1f}%)\n"
                f"Binomial p = {r['binomial_p']:.1e}\n"
                f"BH q = {r['binomial_q_BH']:.1e} {r['fdr_marker']}",
                transform=ax.transAxes, fontsize=7.5, va='top', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9,
                          edgecolor='#CCCCCC'))
        ax.set_xlabel(f'{xl} log\u2082FC', fontsize=FONT_SIZES['label'])
        ax.set_ylabel(f'{yl} log\u2082FC', fontsize=FONT_SIZES['label'])
        ax.set_title(f'{xl} vs {yl}', fontsize=FONT_SIZES['title'],
                     fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        lim = max(abs(x.min()), abs(x.max()),
                  abs(y.min()), abs(y.max())) * 1.1
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

    plt.suptitle('Protein-level log\u2082FC concordance',
                 fontsize=FONT_SIZES['title'], fontweight='bold')
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_C_concordance_scatter.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    pd.DataFrame(rows).to_csv(
        OUTPUT_DIR / 'panel_C_concordance.csv', index=False)
    print("  Saved panel C figure + csv")


# ================ Panel C: Grouped FC with FDR ================
def main():
    _, merged = load_all()
    panel_c_concordance(merged)
    print("Done.")


if __name__ == '__main__':
    main()
