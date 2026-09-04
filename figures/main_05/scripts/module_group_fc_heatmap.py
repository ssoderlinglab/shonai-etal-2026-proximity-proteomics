"""
Figure 5, panel E: module-level protein fold change across the three PD models, by combined functional group.

Input
    figures/main_05/output/permutation_all21.csv
Output
    figures/main_05/output/panel_E_module_FC_heatmap.{png,svg}
Run
    python3 figures/main_05/scripts/module_group_fc_heatmap.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Run figures/main_05/scripts/module_group_permutation.py first.

Each cell is the group mean minus that model's proteome baseline. Significance is the size-matched permutation q value.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
OUT_FIG = os.path.join(FIGURE_DIR, 'output')
TAB = os.path.join(OUT_FIG, 'permutation_all21.csv')

sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))
from permutation_core import GROUPS, MODELS

GENO = {'LRRK2': 'Tg-LRRK2$^{G2019S}$',
        'SNCA':  'Tg-SNCA$^{A30P{\\times}A53T}$',
        'VPS35': 'KI-VPS35$^{D620N}$'}
FOCAL = {'Synapse', 'M14'}      # negative and significant in all three models


def marker(q):
    return ('***' if q < 0.001 else '**' if q < 0.01
            else '*' if q < 0.05 else '')


def main():
    os.makedirs(OUT_FIG, exist_ok=True)
    P = pd.read_csv(TAB)
    assert len(P) == 21, f"expected 21 rows in {TAB}, got {len(P)}"

    order = [g[0] for g in GROUPS]
    labels = {k: lab for k, lab, _ in GROUPS}

    # ---- the two claims this figure makes, asserted before drawing ----
    for gk in FOCAL:
        sub = P[P['group'] == gk]
        assert (sub['delta'] < 0).all(), \
            f"{gk} is not negative in all three models -- do not report it as convergent"
        assert (sub['q_perm_BH'] < 0.05).all(), \
            f"{gk} is not q<0.05 in all three models -- the title claim would be false"

    M = P.pivot(index='group', columns='model',
                values='delta').reindex(order)[MODELS].values

    fig, ax = plt.subplots(figsize=(5.4, 3.3))
    vmax = float(np.nanmax(np.abs(M)))
    im = ax.imshow(M, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

    for i, gk in enumerate(order):
        for j, e in enumerate(MODELS):
            r = P[(P['group'] == gk) & (P['model'] == e)].iloc[0]
            ax.text(j, i, f"{r['delta']:+.3f}{marker(r['q_perm_BH'])}",
                    ha='center', va='center', fontsize=6.8, fontweight='bold',
                    color='white' if abs(r['delta']) > vmax * 0.62 else 'black')

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([labels[g] for g in order], fontsize=6.6)
    for t, g in zip(ax.get_yticklabels(), order):
        if g in FOCAL:
            t.set_color('#C0392B')
            t.set_fontweight('bold')

    ax.set_xticks(range(3))
    ax.set_xticklabels(
        [f"{GENO[e]}\nbaseline "
         f"{P[P['model'] == e]['pool_mean'].iloc[0]:+.3f}" for e in MODELS],
        fontsize=6.4, linespacing=1.5)
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(length=0)
    ax.axhline(2.5, color='white', lw=2.0)

    # horizontal colorbar under the matrix: a vertical bar on the right forced
    # its label into either the cells or the top-right column header
    cb = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.62,
                      pad=0.06, aspect=28)
    cb.ax.tick_params(labelsize=5.8)
    cb.set_label('module mean $-$ own-proteome baseline (log$_2$FC)',
                 fontsize=6.2, labelpad=2)

    fig.suptitle('Synaptic and RNA-binding groups fall below baseline in all '
                 'three models', fontsize=7.6, x=0.01, y=0.99, ha='left')
    fig.text(0.012, 0.185,
             'Each cell is a functional group\u2019s mean protein log$_2$FC minus the mean of '
             'ALL proteins detected in\nthat dataset. Stars are size-matched permutation '
             '(200,000 draws, BH across 21 tests), so a group\ngains nothing from containing '
             'more proteins.     * q<0.05   ** q<0.01   *** q<0.001',
             fontsize=6.0, va='top', linespacing=1.5)
    fig.subplots_adjust(top=0.755, left=0.255, right=0.985, bottom=0.42)

    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(OUT_FIG,
                                 f'panel_E_module_FC_heatmap.{ext}'),
                    dpi=300, facecolor='white')
    print("wrote panel_E_module_FC_heatmap.png/.svg")
    print(f"  q<0.05 cells: {int((P['q_perm_BH'] < 0.05).sum())} of 21")
    for gk in ['Synapse', 'M14']:
        s = P[P['group'] == gk]
        vals = ', '.join(f"{r['model']} {r['delta']:+.4f} (q={r['q_perm_BH']:.3g})"
                         for _, r in s.iterrows())
        print(f"  {gk}: {vals}")

    # ---- self-checks (figure-style 9.1 + canvas containment) ----
    r = fig.canvas.get_renderer()
    ticks = set()
    for a in fig.axes:
        ticks |= set(a.get_xticklabels() + a.get_yticklabels())
    texts = [(t, t.get_window_extent(r))
             for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:20], b.get_text()[:20])
          for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
          if ba.overlaps(bb) and not (a in ticks and b in ticks)]
    fb = fig.bbox
    outside = [t.get_text()[:30] for t, b in texts
               if not (b.x0 >= fb.x0 - 0.5 and b.x1 <= fb.x1 + 0.5
                       and b.y0 >= fb.y0 - 0.5 and b.y1 <= fb.y1 + 0.5)]
    print("  non-tick text overlaps:", ov if ov else "none")
    print("  text outside the canvas:", outside if outside else "none")
    assert not ov and not outside, f"layout defect: overlaps={ov} outside={outside}"
    plt.close(fig)


if __name__ == '__main__':
    main()
