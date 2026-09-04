"""
Supplemental Figure 10, panels B, D and E: tetanus toxin light chain validation by TH-positive cell count, Vamp2 western blot quantification and multi-electrode array mean firing rate.

Input
    data/craave_seq/histology/L_R_paired_comparison.csv
    data/craave_seq/western_blot/fiji_quant.csv
    data/craave_seq/mea/combined_well_data_with_div14.csv
    data/craave_seq/mea/integrated_lmm_results_div14.csv
Output
    figures/suppl_10/output/panel_B_tent_barplot.{png,svg}
    figures/suppl_10/output/panel_D_vamp2_wb.{png,svg}
    figures/suppl_10/output/panel_E_mea_mfr.{png,svg}
Run
    python3 figures/suppl_10/scripts/tent_validation.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib, src/pub_style.py

Panels A and C are microscopy and an immunoblot image, and are not produced here.

Panel B: paired t-test of the two hemispheres, n = 6 animals. Panel D: Welch t-test, n = 4 lanes per condition. Panel E: the linear mixed-effects statistics are read from integrated_lmm_results_div14.csv as fitted and are not refitted here.

The jitter of the individual wells in panel E is seeded per series with a fixed checksum of the series name, so the dot positions are the same in every run.

A single panel can be produced with --only <letter>.
"""

import os
import sys
import zlib
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats as sp_stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'craave_seq'
OUTPUT_DIR = FIGURE_DIR / 'output'
SRC_DIR = REPO_ROOT / 'src'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES
apply_style()

# Data paths
PATH_TENT = DATA_DIR / 'histology' / 'L_R_paired_comparison.csv'
PATH_VAMP2 = DATA_DIR / 'western_blot' / 'fiji_quant.csv'
PATH_MEA_DATA = DATA_DIR / 'mea' / 'combined_well_data_with_div14.csv'
PATH_MEA_STATS = DATA_DIR / 'mea' / 'integrated_lmm_results_div14.csv'


def _p_to_stars(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'n.s.'


def _bracket(ax, x1, x2, y, h, text, fs=7):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color='black', lw=0.8, zorder=5, clip_on=False)
    ax.text((x1+x2)/2, y+h+1, text,
            ha='center', va='bottom', fontsize=fs, zorder=5)


def _bracket_wb(ax, x1, x2, y, h, text, fs=7):
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
            color='black', lw=0.8, zorder=5, clip_on=False)
    ax.text((x1+x2)/2, y+h+0.01, text,
            ha='center', va='bottom', fontsize=fs, zorder=5)


# ============================================================
# PANEL B -- TeNT TH+ Cell Count
# ============================================================
def panel_b_tent():
    """Panel B: TeNT TH+ cell count barplot.

    Per-animal averaging of L (alpha-syn only) vs R (alpha-syn + TeNT).
    Paired t-test.
    """
    COLOR_L = '#dfc27d'   # light tan — alpha-syn only
    COLOR_R = '#a6611a'   # dark brown — alpha-syn + TeNT

    tent_raw = pd.read_csv(PATH_TENT)
    tent_animal = tent_raw.groupby('AnimalID')[
        ['neuron_count_L', 'neuron_count_R']].mean().reset_index()

    n = len(tent_animal)
    count_L = tent_animal['neuron_count_L'].values
    count_R = tent_animal['neuron_count_R'].values

    t_stat, p_val = sp_stats.ttest_rel(count_R, count_L)

    print(f"  alpha-syn + TeNT (n={n}): L(asyn)={count_L.mean():.1f}+/-{count_L.std(ddof=1)/np.sqrt(n):.1f}, "
          f"R(asyn+TeNT)={count_R.mean():.1f}+/-{count_R.std(ddof=1)/np.sqrt(n):.1f}")
    print(f"  paired t p={p_val:.4f}")

    fig, ax = plt.subplots(figsize=(2, 4))

    bar_width = 0.35
    pair_gap = 0.02
    center = 0.5
    x_L = center - (bar_width / 2 + pair_gap / 2)
    x_R = center + (bar_width / 2 + pair_gap / 2)

    means = [count_L.mean(), count_R.mean()]
    sems = [count_L.std(ddof=1)/np.sqrt(n), count_R.std(ddof=1)/np.sqrt(n)]

    ax.bar(x_L, means[0], width=bar_width, color=COLOR_L, edgecolor='none', zorder=2)
    ax.errorbar(x_L, means[0], yerr=sems[0], fmt='none', color='black',
                capsize=4, capthick=1.0, linewidth=1.0, zorder=3)
    ax.bar(x_R, means[1], width=bar_width, color=COLOR_R, edgecolor='none', zorder=2)
    ax.errorbar(x_R, means[1], yerr=sems[1], fmt='none', color='black',
                capsize=4, capthick=1.0, linewidth=1.0, zorder=3)

    # Individual data points
    np.random.seed(300)
    jitter_L = np.random.uniform(-bar_width * 0.3, bar_width * 0.3, size=n)
    ax.scatter(x_L + jitter_L, count_L, s=18, color='black',
               edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)
    np.random.seed(301)
    jitter_R = np.random.uniform(-bar_width * 0.3, bar_width * 0.3, size=n)
    ax.scatter(x_R + jitter_R, count_R, s=18, color='black',
               edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)

    top = max(means[0]+sems[0], means[1]+sems[1])
    _bracket(ax, x_L, x_R, top + 30, 5, _p_to_stars(p_val))

    ax.set_xticks([x_L, x_R])
    ax.set_xticklabels(['\u03b1-syn only', '\u03b1-syn + TeNT'],
                        fontsize=FONT_SIZES['tick'], rotation=45, ha='right')
    ax.set_ylabel(r'$\mathrm{TH}^{+}$ Cell Count', fontsize=FONT_SIZES['label'])
    ax.set_ylim(0, 250)
    ax.set_title('A53T-\u03b1Syn \u00b1 TeNT', fontsize=FONT_SIZES['title'], fontweight='bold')

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_B_tent_barplot.{ext}',
                    dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    print("  Saved: panel_B_tent_barplot.png/svg")


# ============================================================
# PANEL D -- Vamp2 WB Quantification
# ============================================================
def panel_d_vamp2():
    """Panel D: Vamp2 WB quantification.

    Vamp2/total protein, normalized to Cre- mean = 1.0.
    Welch t-test.
    """
    COLOR_CTRL = '#999999'
    COLOR_TENT = '#4A4A4A'

    df = pd.read_csv(PATH_VAMP2)

    pivot = df.pivot_table(
        index=["lane", "sample", "rep"],
        columns="detection",
        values="intensity"
    ).reset_index()

    pivot["normalized"] = pivot["anti-Vamp2"] / pivot["TotalProtein"]
    cre_neg_mean = pivot.loc[pivot["sample"] == "Cre-", "normalized"].mean()
    pivot["scaled"] = pivot["normalized"] / cre_neg_mean

    cre_neg = pivot.loc[pivot["sample"] == "Cre-", "scaled"].values
    cre_pos = pivot.loc[pivot["sample"] == "Cre+", "scaled"].values

    conditions = [
        {'key': 'Cre-', 'label': 'Cre\u2212', 'color': COLOR_CTRL, 'data': cre_neg},
        {'key': 'Cre+', 'label': 'Cre+\n(TeNT)', 'color': COLOR_TENT, 'data': cre_pos},
    ]

    for c in conditions:
        c['n'] = len(c['data'])
        c['mean'] = c['data'].mean()
        c['sem'] = c['data'].std(ddof=1) / np.sqrt(c['n'])
        print(f"    {c['key']} (n={c['n']}): mean={c['mean']:.3f} +/- {c['sem']:.3f}")

    t_stat, p_val = sp_stats.ttest_ind(cre_neg, cre_pos, equal_var=False)
    print(f"    Welch t: t={t_stat:.4f}, p={p_val:.4f}")

    def fmt_p(p):
        if p < 0.001: return 'p < 0.001'
        return f'p = {p:.3f}'

    fig, ax = plt.subplots(figsize=(2.5, 4))
    positions = [0, 1]
    bar_width = 0.6

    for i, c in enumerate(conditions):
        ax.bar(positions[i], c['mean'], width=bar_width, color=c['color'],
               edgecolor='none', zorder=2)
        ax.errorbar(positions[i], c['mean'], yerr=c['sem'], fmt='none',
                    color='black', capsize=4, capthick=1.0, linewidth=1.0, zorder=3)
        np.random.seed(800 + i)
        jitter = np.random.uniform(-bar_width * 0.25, bar_width * 0.25, size=c['n'])
        ax.scatter(positions[i] + jitter, c['data'], s=25, color='black',
                   edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)

    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.6, alpha=0.5, zorder=1)

    y_top = max(max(c['data']) for c in conditions)
    _bracket_wb(ax, positions[0], positions[1], y_top + 0.06, 0.03, fmt_p(p_val))

    ax.set_xticks(positions)
    ax.set_xticklabels([c['label'] for c in conditions], fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('Vamp2 / total protein', fontsize=FONT_SIZES['label'])
    ax.set_ylim(0, y_top + 0.20)
    ax.set_title('TeNT Validation', fontsize=FONT_SIZES['title'], fontweight='bold')

    fig.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_D_vamp2_wb.{ext}',
                    dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    print("  Saved: panel_D_vamp2_wb.png/svg")


# ============================================================
# PANEL E -- MEA Mean Firing Rate
# ============================================================
def panel_e_mea():
    """Panel E: MEA mean firing rate timecourse.

    Control vs TeNT across DIV7, DIV11, DIV14.
    Stats from pre-computed LMM (FDR-corrected).
    """
    COLOR_CTRL = '#999999'
    COLOR_TENT = '#C0392B'

    df = pd.read_csv(PATH_MEA_DATA)
    df = df[df['Treatment'].isin(['control', 'TeNT'])].copy()

    res = pd.read_csv(PATH_MEA_STATS)
    mfr_stats = res[res['metric'] == 'Mean Firing Rate (Hz)'].iloc[0]

    metric = 'Mean Firing Rate (Hz)'
    div_order = ['DIV7', 'DIV11', 'DIV14']

    print(f"  Global: treat q={mfr_stats['q_treatment']:.2e}, "
          f"interact q={mfr_stats['q_interaction']:.2e}")

    for d in div_order:
        ctrl = df[(df['DIV'] == d) & (df['Treatment'] == 'control')][metric].values
        tent = df[(df['DIV'] == d) & (df['Treatment'] == 'TeNT')][metric].values
        q_val = mfr_stats[f'q_{d}']
        print(f"    {d}: control={np.mean(ctrl):.2f}+/-{sp_stats.sem(ctrl):.2f}, "
              f"TeNT={np.mean(tent):.2f}+/-{sp_stats.sem(tent):.2f}, q={q_val:.4f}")

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    x_map = {'DIV7': 0, 'DIV11': 1, 'DIV14': 2}

    for treatment, color, label in [('control', COLOR_CTRL, 'Control'),
                                     ('TeNT', COLOR_TENT, 'TeNT')]:
        sub = df[df['Treatment'] == treatment]

        means = []
        sems = []
        for d in div_order:
            vals = sub[sub['DIV'] == d][metric].values
            means.append(np.mean(vals))
            sems.append(sp_stats.sem(vals))

        x_pos = [x_map[d] for d in div_order]

        # Individual wells as faint dots
        for d in div_order:
            vals = sub[sub['DIV'] == d][metric].values
            np.random.seed(zlib.crc32((treatment + d).encode()) % 10000)
            jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
            ax.scatter(x_map[d] + jitter, vals, s=8, color=color, alpha=0.25,
                       edgecolors='none', zorder=2)

        # Plate means as faint lines
        for plate in sub['Plate'].unique():
            plate_means = []
            for d in div_order:
                plate_vals = sub[(sub['DIV'] == d) & (sub['Plate'] == plate)][metric].values
                plate_means.append(np.mean(plate_vals))
            ax.plot(x_pos, plate_means, color=color, alpha=0.2, lw=0.8, zorder=3)

        # Group mean +/- SEM
        ax.errorbar(x_pos, means, yerr=sems, color=color, lw=2.0,
                    capsize=3, capthick=1.0, marker='o', markersize=5,
                    markerfacecolor=color, markeredgecolor='white', markeredgewidth=0.5,
                    label=label, zorder=5)

    # Per-DIV significance
    def q_to_stars(q):
        if q < 0.001: return '***'
        elif q < 0.01: return '**'
        elif q < 0.05: return '*'
        return ''

    for d in div_order:
        q_val = mfr_stats[f'q_{d}']
        stars = q_to_stars(q_val)
        if stars:
            ctrl_vals = df[(df['DIV'] == d) & (df['Treatment'] == 'control')][metric].values
            y_max = np.max(ctrl_vals)
            ax.text(x_map[d], y_max - 1.0, stars, ha='center', va='bottom',
                    fontsize=12, fontweight='bold', color='#C0392B')

    # Global model annotation
    ax.text(0.02, 0.90,
            f'treat q = {mfr_stats["q_treatment"]:.2e}\n'
            f'interact q = {mfr_stats["q_interaction"]:.2e}',
            transform=ax.transAxes, va='top', ha='left', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='grey', alpha=0.8, linewidth=0.5))

    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(div_order, fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('Mean Firing Rate (Hz)', fontsize=FONT_SIZES['label'])
    ax.set_ylim(bottom=0)
    ax.set_title('TeNT Neuronal Silencing', fontsize=FONT_SIZES['title'], fontweight='bold')
    ax.legend(fontsize=8, frameon=False, loc='upper left',
              bbox_to_anchor=(0.0, 0.80))

    fig.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_E_mea_mfr.{ext}',
                    dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    print("  Saved: panel_E_mea_mfr.png/svg")


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Supplemental Figure 10 -- TeNT validation')
    parser.add_argument('--only', choices=['B', 'D', 'E'],
                        help='Generate single panel (default: all)')
    args = parser.parse_args()

    print("Supplemental Figure 10 -- TeNT validation")

    panels = args.only or 'ALL'

    if panels in ('ALL', 'B'):
        print("\nPanel B -- TeNT TH+ Cell Count:")
        panel_b_tent()

    if panels in ('ALL', 'D'):
        print("\nPanel D -- Vamp2 WB Quantification:")
        panel_d_vamp2()

    if panels in ('ALL', 'E'):
        print("\nPanel E -- MEA Mean Firing Rate:")
        panel_e_mea()

    print("\nDone.")


if __name__ == '__main__':
    main()
