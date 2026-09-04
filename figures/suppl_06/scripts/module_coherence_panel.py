"""
Supplemental Figure 6, panel B: module coherence against size-matched random protein sets, per cohort.

Input
    figures/suppl_06/output/module_coherence_density_BH.csv
Output
    figures/suppl_06/output/panel_B_module_coherence.{png,svg}
    figures/suppl_06/output/panel_B_module_coherence_source_data.csv
Run
    python3 figures/suppl_06/scripts/module_coherence_panel.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Third of three steps; run the permutation and the correction first.

Significance markers are drawn from q, not from Z: the null distributions are right-skewed.
"""

import csv
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
OUT = os.path.join(FIGURE_DIR, 'output')
SRC = os.path.join(OUT, 'module_coherence_density_BH.csv')
STEM = os.path.join(OUT, 'panel_B_module_coherence')

# COLOR_SYSTEM.md section 2, PD model identity (Dark2)
COL = {'LRRK2': '#1B9E77', 'SNCA': '#D95F02', 'VPS35': '#7570B3'}
NULL_GREY = '#BFBFBF'
LINK_GREY = '#8C8C8C'

COHORTS = ('LRRK2', 'SNCA', 'VPS35')
TITLE = {'LRRK2': 'Tg-LRRK2', 'SNCA': 'Tg-SNCA', 'VPS35': 'KI-Vps35'}
SUBTITLE = {'LRRK2': 'reference', 'SNCA': 'test', 'VPS35': 'test'}
ORDER = ['M1', 'M2', 'M4', 'M5', 'M6', 'M7', 'M8', 'M9', 'M10', 'M11', 'M14',
         'M19']
LABEL = {'M1': 'M1  ribosome', 'M2': 'M2  mitochondria', 'M4': 'M4  ER/Golgi',
         'M5': 'M5  mitochondria', 'M6': 'M6  clathrin/tRNA',
         'M7': 'M7  SV cycling', 'M8': 'M8  postsyn/cytoskeleton',
         'M9': 'M9  proteasome', 'M10': 'M10  ribosome',
         'M11': 'M11  presynaptic/SV', 'M14': 'M14  RNA-binding',
         'M19': 'M19  endocytosis'}

FS_TICK, FS_LAB, FS_TITLE, FS_SMALL = 6.6, 7.0, 7.6, 6.0


def star(q):
    return ('***' if q < 0.001 else '**' if q < 0.01
            else '*' if q < 0.05 else 'ns')


def load():
    d = {}
    for r in csv.DictReader(open(SRC)):
        if r['in_functional_12'] != 'True' or r['genotype'] != 'wt':
            continue
        q = r['q_BH_primary_24']
        d[(r['cohort'], r['module'])] = {
            'obs': float(r['meanCor']), 'null': float(r['meanCor_null_mean']),
            'sd': float(r['meanCor_null_sd']), 'z': float(r['Z_meanCor']),
            'n': int(r['n_tested_in_cohort']),
            'q': float(q) if q not in ('', None) else None}
    return d


def main():
    d = load()
    assert len(d) == 36, 'expected 36 wt cells, got %d' % len(d)
    for c in ('SNCA', 'VPS35'):
        for m in ORDER:
            assert d[(c, m)]['q'] is not None, 'missing q for %s %s' % (c, m)

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.05), sharey=True,
                             sharex=True)
    y = np.arange(len(ORDER))

    for ax, cohort in zip(axes, COHORTS):
        is_ref = cohort == 'LRRK2'
        for i, m in enumerate(ORDER):
            v = d[(cohort, m)]
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, color='#F2F2F2', lw=0, zorder=0)
            ax.plot([v['null'] - 2 * v['sd'], v['null'] + 2 * v['sd']],
                    [i, i], '-', color=NULL_GREY, lw=3.2,
                    solid_capstyle='butt', zorder=1)
            ax.plot([v['null'], v['obs']], [i, i], '-', color=LINK_GREY,
                    lw=0.6, zorder=2)
            ax.plot(v['obs'], i, 'o', ms=4.6, color=COL[cohort],
                    mec='white', mew=0.5, zorder=3)
            if not is_ref:
                ax.text(1.055, i, star(v['q']), fontsize=FS_SMALL,
                        va='center', ha='center', color='#333333')
            ax.text(1.30, i, '%d' % v['n'], fontsize=FS_SMALL, va='center',
                    ha='right', color='#666666')

        ax.set_xlim(-0.12, 1.32)
        ax.set_ylim(len(ORDER) - 0.5, -0.5)
        ax.set_xticks([0.0, 0.5, 1.0])
        ax.tick_params(labelsize=FS_TICK, length=2.5, pad=1.5)
        ax.set_title('%s\n%s' % (TITLE[cohort], SUBTITLE[cohort]),
                     fontsize=FS_TITLE, color=COL[cohort], linespacing=1.35,
                     pad=4)
        ax.text(1.30, -0.95, 'n', fontsize=FS_SMALL, ha='right',
                color='#666666', style='italic')
        if not is_ref:
            ax.text(1.055, -0.95, 'q', fontsize=FS_SMALL, ha='center',
                    color='#333333', style='italic')
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        for sp in ('left', 'bottom'):
            ax.spines[sp].set_linewidth(0.6)

    axes[0].set_yticks(y)
    axes[0].set_yticklabels([LABEL[m] for m in ORDER], fontsize=FS_TICK)

    fig.supxlabel('mean pairwise correlation across sedimentation fractions',
                  fontsize=FS_LAB, y=0.105)

    h_null = plt.Line2D([], [], color=NULL_GREY, lw=3.2,
                        label='size-matched random protein sets '
                              '(mean $\\pm$ 2 SD, 2,000 draws)')
    h_obs = plt.Line2D([], [], marker='o', ls='none', ms=4.6, color='#555555',
                       mec='white', mew=0.5, label='module, as measured')
    fig.legend(handles=[h_obs, h_null], loc='lower center', ncol=2,
               fontsize=FS_SMALL, frameon=False,
               bbox_to_anchor=(0.5, 0.005), handletextpad=0.6,
               columnspacing=1.8)

    fig.subplots_adjust(left=0.235, right=0.995, top=0.855, bottom=0.20,
                        wspace=0.10)
    # ⛔ SVG IS THE DELIVERABLE AND IT IS FULLY TRANSPARENT. Illustrator
    # placement is Daichi's step, so no opaque rectangle may survive: matplotlib
    # emits BOTH a figure patch and a per-axes patch, and passing
    # facecolor='none' to savefig only clears the FIGURE one. The axes patches
    # are cleared here. The PNG keeps a white background because it is the
    # on-screen proof, not the artwork.
    for ax in axes:
        ax.patch.set_alpha(0.0)
    fig.patch.set_alpha(0.0)
    fig.savefig(STEM + '.svg', dpi=300, facecolor='none',
                edgecolor='none', transparent=True)
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)
    fig.savefig(STEM + '.png', dpi=300, facecolor='white')
    plt.close(fig)
    print('wrote %s .svg (transparent) and .png (white proof)'
          % os.path.basename(STEM))
    for c in ('SNCA', 'VPS35'):
        marks = ' '.join('%s=%s' % (m, star(d[(c, m)]['q'])) for m in ORDER)
        print('%-6s %s' % (c, marks))
    return 0


if __name__ == '__main__':
    sys.exit(main())
