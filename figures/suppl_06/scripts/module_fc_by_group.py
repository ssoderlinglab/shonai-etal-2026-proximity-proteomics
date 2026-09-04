"""
Supplemental Figure 6, panel D: protein-level fold-change distribution by functional group, for the 409 proteins detected in all three cohorts.

Input
    figures/main_05/output/permutation_all21.csv
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
Output
    figures/suppl_06/output/panel_D_module_fc.{png,svg}
Run
    python3 figures/suppl_06/scripts/module_fc_by_group.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Run figures/main_05/scripts/module_group_permutation.py first.

The dashed line is each cohort's proteome baseline. Significance is the size-matched permutation q value.

Display range is +/-0.45, with **56 of 2021 clipped** (2.77%) drawn at the edge. The script re-computes this at run time and fails if it disagrees with this line.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
OUT_FIG = os.path.join(FIGURE_DIR, 'output')
TAB = os.path.join(REPO_ROOT, 'figures', 'main_05', 'output',
                   'permutation_all21.csv')

sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))
from permutation_core import load_all, GROUPS, MODELS, SEED

COL = {'LRRK2': '#1B9E77', 'SNCA': '#D95F02', 'VPS35': '#7570B3'}
GENO = {'LRRK2': 'Tg-LRRK2$^{G2019S}$',
        'SNCA':  'Tg-SNCA$^{A30P{\\times}A53T}$',
        'VPS35': 'KI-VPS35$^{D620N}$'}
FOCAL = {'Synapse', 'M14'}          # negative and significant in all three
LIM = 0.45                          # 97.23% of module dots inside; see docstring


def marker(q):
    return ('***' if q < 0.001 else '**' if q < 0.01
            else '*' if q < 0.05 else '')


def main():
    os.makedirs(OUT_FIG, exist_ok=True)
    merged, gsets = load_all()
    P = pd.read_csv(TAB)
    assert len(P) == 21, f"expected 21 rows, got {len(P)}"
    assert (P[P['group'] == 'Synapse']['delta'] < 0).all(), \
        "Synapse not negative in all three models -- re-read before reporting"

    rng = np.random.default_rng(SEED)
    order = [g[0] for g in GROUPS][::-1]        # Synapse ends up at the top
    labels = {k: lab for k, lab, _ in GROUPS}

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 3.9), sharey=True)
    n_clip = 0

    for j, e in enumerate(MODELS):
        ax = axes[j]
        pool = merged[['Gene', f'{e}_FC']].dropna()
        base = float(P[P['model'] == e]['pool_mean'].iloc[0])

        # alternating row bands make the horizontal scan easy
        for i, gk in enumerate(order):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='#F5F5F5',
                           edgecolor='none', zorder=0)

        ax.axvline(0, color='#DDDDDD', lw=0.6, zorder=1)
        ax.axvline(base, color='#333333', lw=1.1, ls=(0, (4, 2)), zorder=4)

        for i, gk in enumerate(order):
            row = P[(P['group'] == gk) & (P['model'] == e)].iloc[0]
            v = pool[pool['Gene'].isin(gsets[gk])][f'{e}_FC'].values
            n_clip += int((np.abs(v) > LIM).sum())

            ax.scatter(np.clip(v, -LIM, LIM),
                       np.full(len(v), i) + rng.normal(0, 0.115, len(v)),
                       s=2.6, alpha=0.20, c=COL[e], edgecolors='none',
                       zorder=3)
            # module mean: a bar tall enough to read against the dot cloud
            ax.plot([row['mean_raw'], row['mean_raw']], [i - 0.33, i + 0.33],
                    color=COL[e], lw=3.2, zorder=6, solid_capstyle='butt')
            # thin white halo so the bar reads even inside a dense cloud
            ax.plot([row['mean_raw'], row['mean_raw']], [i - 0.33, i + 0.33],
                    color='white', lw=5.0, zorder=5, solid_capstyle='butt')

            m = marker(row['q_perm_BH'])
            if m:
                ax.annotate(m, xy=(LIM * 0.97, i), ha='right', va='center',
                            fontsize=8.5, fontweight='bold', color=COL[e],
                            zorder=7)

        ax.set_xlim(-LIM, LIM)
        ax.set_ylim(-0.55, len(order) - 0.45)
        ax.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
        ax.set_xlabel('protein log$_2$FC (Mut / WT)', fontsize=6.8)
        ax.set_title(GENO[e], fontsize=7.4, color=COL[e], fontweight='bold',
                     loc='left', pad=5)
        ax.tick_params(labelsize=6.2, length=2)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.spines['left'].set_visible(False)

        if j > 0:
            ax.tick_params(axis='y', length=0, labelleft=False)
        if j == 0:
            ax.set_yticks(range(len(order)))
            ax.set_yticklabels([labels[g] for g in order], fontsize=6.6)
            for t, g in zip(ax.get_yticklabels(), order):
                if g in FOCAL:
                    t.set_color('#C0392B')
                    t.set_fontweight('bold')
            ax.tick_params(axis='y', length=0)

    # ---- one key line under the panels, anchored in FIGURE coords so it
    # cannot collide with the x-axis labels regardless of the margin

    fig.suptitle('Protein-level log$_2$FC by functional group',
                 fontsize=7.8, x=0.5, y=0.985, ha='center')
    fig.subplots_adjust(top=0.845, left=0.185, right=0.995, bottom=0.155,
                        wspace=0.06)

    for ext in ('png', 'svg'):
        fig.savefig(os.path.join(OUT_FIG, f'panel_D_module_fc.{ext}'),
                    dpi=300, facecolor='white')
    print("wrote panel_D_module_fc.png/.svg")
    # Total module-dot count is itself derived, not hardcoded.
    n_total = sum(
        len(merged[merged['Gene'].isin(gsets[gk])][f'{e}_FC'].dropna())
        for gk, _, _ in GROUPS for e in MODELS)
    pct_clip = 100 * n_clip / n_total
    print(f"  display range +/-{LIM}; dots clipped: {n_clip} of {n_total} "
          f"({pct_clip:.2f}%), inside {100 - pct_clip:.2f}% "
          f"-- USE THESE NUMBERS IN THE CAPTION")
    # Guard against the docstring drifting away from the executed result: the
    # clipped count quoted in the module docstring must match what was just
    # computed.  A LIM change therefore fails the run until the text is updated.
    import re as _re
    claimed = _re.search(r'\*\*(\d+) of (\d+) clipped', __doc__)
    assert claimed, "docstring no longer states the clipped count in the expected form"
    assert (int(claimed.group(1)), int(claimed.group(2))) == (n_clip, n_total), (
        f"docstring says {claimed.group(1)} of {claimed.group(2)} clipped but this "
        f"run gives {n_clip} of {n_total} -- update the docstring")

    # geometric self-check
    r = fig.canvas.get_renderer()
    ticks = set()
    for a in fig.axes:
        ticks |= set(a.get_xticklabels() + a.get_yticklabels())
    texts = [(t, t.get_window_extent(r))
             for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:22], b.get_text()[:22])
          for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
          if ba.overlaps(bb) and not (a in ticks and b in ticks)]
    print("  non-tick text overlaps:", ov if ov else "none")
    fb = fig.bbox
    outside = [t.get_text()[:30] for t, b in texts
               if not (b.x0 >= fb.x0 - 0.5 and b.x1 <= fb.x1 + 0.5
                       and b.y0 >= fb.y0 - 0.5 and b.y1 <= fb.y1 + 0.5)]
    print("  text outside the canvas:", outside if outside else "none")
    assert not outside, f"text clipped by the canvas: {outside}"
    plt.close(fig)


if __name__ == '__main__':
    main()
