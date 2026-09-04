"""
Figure 3, panels D, F and H: diverging GO enrichment shift between the WT and mutant proximity proteomes, per bait.

Input
    data/ibioid/go_enrichment/{bait}_{WT,mutant}.Component.tsv
Output
    figures/main_03/output/*.{png,svg}
Run
    python3 figures/main_03/scripts/go_shift_perbait.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib
"""

import argparse
import csv
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
from matplotlib.colorbar import ColorbarBase
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform
import textwrap

# ── Paths ───────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid', 'go_enrichment')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')

# ── Filter parameters (project-wide standard, inherited) ─────────────
BG_MAX = 500
FDR_THRESHOLD = 0.05
TOP_N = 10
GROUP_DISTANCE_THRESHOLD = 0.85

CMAP_WT = LinearSegmentedColormap.from_list(
    'fdr_wt', ['#004D40', '#009688', '#80CBC4', '#E0F2F1'], N=256)
CMAP_MUT = LinearSegmentedColormap.from_list(
    'fdr_mut', ['#BF360C', '#FF6F61', '#FFAB91', '#FBE9E7'], N=256)

FDR_VMIN = 1e-10
FDR_VMAX = 1e-1
BAND_COLORS = ['#F5F5F5', '#FFFFFF']
LOLLIPOP_STEM_CAPSTYLE = 'butt'
DOT_EDGE_DARKEN = 0.55
DOT_EDGE_LINEWIDTH = 0.6

DATASETS = {
    'Vps35_CC': {
        'wt_tsv': os.path.join(DATA_DIR, 'Vps35_WT.Component.tsv'),
        'mut_tsv': os.path.join(DATA_DIR, 'Vps35_D620N.Component.tsv'),
        'title': 'Vps35 D620N \u2014 Cellular Component',
        'mut_label': 'D620N',
    },
    'aSyn_CC': {
        'wt_tsv': os.path.join(DATA_DIR, 'aSyn_WT.Component.tsv'),
        'mut_tsv': os.path.join(DATA_DIR, 'aSyn_A30P.Component.tsv'),
        'title': '\u03b1Syn A30P \u2014 Cellular Component',
        'mut_label': 'A30P',
    },
    'Lrrk2_CC': {
        'wt_tsv': os.path.join(DATA_DIR, 'Lrrk2_WT.Component.tsv'),
        'mut_tsv': os.path.join(DATA_DIR, 'Lrrk2_G2019S.Component.tsv'),
        'title': 'Lrrk2 G2019S \u2014 Cellular Component',
        'mut_label': 'G2019S',
    },
}


def load_enrichment_tsv(filepath):
    terms = {}
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            tid = row['#term ID']
            terms[tid] = {
                'term_id': tid,
                'description': row['term description'],
                'gene_count': int(row['observed gene count']),
                'bg_count': int(row['background gene count']),
                'strength': float(row['strength']),
                'signal': float(row['signal']),
                'fdr': float(row['false discovery rate']),
                'proteins': set(row['matching proteins in your network (labels)'].split(','))
            }
    return terms


def compute_signal_diff(wt_terms, mut_terms):
    # Iterate a sorted list rather than a set, so the row order of `results` does not depend
    # on PYTHONHASHSEED. Term selection is data-driven and unaffected; this only fixes how
    # ties at identical signal_diff are resolved.
    all_ids = sorted(set(wt_terms.keys()) | set(mut_terms.keys()))
    results = []
    for tid in all_ids:
        wt = wt_terms.get(tid)
        mut = mut_terms.get(tid)
        bg = wt['bg_count'] if wt else mut['bg_count']
        if bg > BG_MAX:
            continue
        wt_sig = wt is not None and wt['fdr'] < FDR_THRESHOLD
        mut_sig = mut is not None and mut['fdr'] < FDR_THRESHOLD
        if not (wt_sig or mut_sig):
            continue
        wt_signal = wt['signal'] if wt else 0.0
        mut_signal = mut['signal'] if mut else 0.0
        results.append({
            'term_id': tid,
            'description': wt['description'] if wt else mut['description'],
            'signal_diff': wt_signal - mut_signal,
            'direction': 'WT-enriched' if wt_signal - mut_signal > 0 else 'MUT-enriched',
            'wt_fdr': wt['fdr'] if wt else 1.0,
            'mut_fdr': mut['fdr'] if mut else 1.0,
            'wt_gene_count': wt['gene_count'] if wt else 0,
            'mut_gene_count': mut['gene_count'] if mut else 0,
            'proteins': (wt['proteins'] if wt else set()) | (mut['proteins'] if mut else set()),
        })
    return results


def select_top_terms(all_terms):
    wt_terms = [t for t in all_terms if t['signal_diff'] > 0]
    mut_terms = [t for t in all_terms if t['signal_diff'] < 0]
    # Secondary key term_id breaks signal_diff ties deterministically.
    wt_terms.sort(key=lambda x: (-x['signal_diff'], x['term_id']))
    mut_terms.sort(key=lambda x: (x['signal_diff'], x['term_id']))
    return wt_terms[:TOP_N], mut_terms[:TOP_N]


def group_by_similarity(terms):
    n = len(terms)
    if n <= 1:
        return terms, [0] * n
    if n == 2:
        inter = len(terms[0]['proteins'] & terms[1]['proteins'])
        union = len(terms[0]['proteins'] | terms[1]['proteins'])
        jacc_dist = 1.0 - (inter / union if union > 0 else 0.0)
        gid = [0, 0] if jacc_dist < GROUP_DISTANCE_THRESHOLD else [0, 1]
        return terms, gid

    dist_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            inter = len(terms[i]['proteins'] & terms[j]['proteins'])
            union = len(terms[i]['proteins'] | terms[j]['proteins'])
            jacc = inter / union if union > 0 else 0.0
            dist_matrix[i, j] = 1.0 - jacc
            dist_matrix[j, i] = 1.0 - jacc

    condensed = squareform(dist_matrix)
    if np.all(condensed == condensed[0]):
        return terms, list(range(n))

    Z = linkage(condensed, method='average')
    order = leaves_list(Z)
    clusters = fcluster(Z, t=GROUP_DISTANCE_THRESHOLD, criterion='distance')
    ordered_terms = [terms[i] for i in order]
    ordered_clusters = [clusters[i] for i in order]

    seen, group_ids, counter = {}, [], 0
    for c in ordered_clusters:
        if c not in seen:
            seen[c] = counter
            counter += 1
        group_ids.append(seen[c])
    return ordered_terms, group_ids


def fdr_to_color(fdr_val, cmap):
    norm = LogNorm(vmin=FDR_VMIN, vmax=FDR_VMAX, clip=True)
    return cmap(norm(max(fdr_val, FDR_VMIN)))


def darken_face_color(face, factor=DOT_EDGE_DARKEN):
    r, g, b = to_rgb(face)
    return (r * factor, g * factor, b * factor)


def draw_group_bands(ax, group_ids, y_min, y_max):
    if not group_ids:
        return
    runs, run_start = [], 0
    for i in range(1, len(group_ids)):
        if group_ids[i] != group_ids[i - 1]:
            runs.append((run_start, i - 1, group_ids[run_start]))
            run_start = i
    runs.append((run_start, len(group_ids) - 1, group_ids[run_start]))
    for band_idx, (start, end, gid) in enumerate(runs):
        ax.axvspan(start - 0.5, end + 0.5, color=BAND_COLORS[band_idx % 2],
                   alpha=1.0, zorder=0)


def plot_horizontal_shift(wt_terms, mut_terms, title, mut_label, output_prefix):
    wt_grouped, wt_gids = group_by_similarity(wt_terms)
    mut_grouped, mut_gids = group_by_similarity(mut_terms)
    mut_max = max(mut_gids) + 1 if mut_gids else 0
    wt_gids_offset = [g + mut_max for g in wt_gids]

    combined = list(mut_grouped) + list(wt_grouped)
    combined_gids = list(mut_gids) + list(wt_gids_offset)
    if not combined:
        print(f"  No terms to plot for {output_prefix}")
        return

    n = len(combined)
    x_pos = np.arange(n)
    diffs = [t['signal_diff'] for t in combined]
    labels = [t['description'] for t in combined]

    colors = []
    for t in combined:
        if t['signal_diff'] > 0:
            colors.append(fdr_to_color(t['wt_fdr'], CMAP_WT))
        else:
            colors.append(fdr_to_color(t['mut_fdr'], CMAP_MUT))

    gene_counts = [max(t['wt_gene_count'], t['mut_gene_count']) for t in combined]

    def wrap_label(text, max_width=25):
        return text if len(text) <= max_width else textwrap.fill(text, width=max_width)

    labels = [wrap_label(lbl) for lbl in labels]

    fig_width = max(10, n * 0.8)
    fig, ax = plt.subplots(figsize=(fig_width, 3.0))
    fig.subplots_adjust(right=0.82)

    y_max = max(abs(d) for d in diffs) * 1.15
    draw_group_bands(ax, combined_gids, -y_max, y_max)

    for x, d, c in zip(x_pos, diffs, colors):
        ax.plot([x, x], [0, d], color=c, linewidth=6,
                solid_capstyle=LOLLIPOP_STEM_CAPSTYLE, zorder=3)

    sizes = [max(100, gc * 10) for gc in gene_counts]
    dot_edges = [darken_face_color(c) for c in colors]
    ax.scatter(x_pos, diffs, c=colors, s=sizes, zorder=5,
               edgecolors=dot_edges, linewidth=DOT_EDGE_LINEWIDTH)

    ax.axhline(y=0, color='#333333', linewidth=0.8, linestyle='-', zorder=2)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=8, fontweight='bold')
    ax.set_ylabel('Signal difference (WT \u2212 MUT)', fontsize=11, fontweight='bold')
    ax.set_ylim(-y_max, y_max)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_title(title, fontsize=24, fontweight='bold', pad=12)

    ax.text(0.30, 0.97, '\u2191 WT-enriched', transform=ax.transAxes,
            fontsize=18, color='#004D40', fontweight='bold', va='top', zorder=10)
    ax.text(0.70, 0.03, '\u2193 MUT-enriched', transform=ax.transAxes,
            fontsize=18, color='#BF360C', fontweight='bold', ha='right',
            va='bottom', zorder=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    norm = LogNorm(vmin=FDR_VMIN, vmax=FDR_VMAX)
    cbar_wt_ax = fig.add_axes([0.84, 0.60, 0.01, 0.25])
    cb_wt = ColorbarBase(cbar_wt_ax, cmap=CMAP_WT, norm=norm, orientation='vertical')
    cb_wt.set_label('FDR (WT)', fontsize=9, fontweight='bold')
    cb_wt.ax.tick_params(labelsize=7)

    cbar_mut_ax = fig.add_axes([0.84, 0.28, 0.01, 0.25])
    cb_mut = ColorbarBase(cbar_mut_ax, cmap=CMAP_MUT, norm=norm, orientation='vertical')
    cb_mut.set_label(f'FDR ({mut_label})', fontsize=9, fontweight='bold')
    cb_mut.ax.tick_params(labelsize=7)

    legend_face = '#888888'
    legend_edge = darken_face_color(legend_face)
    legend_dots = []
    for s in [5, 20, 50]:
        legend_dots.append(Line2D([0], [0], linestyle='none', marker='o',
                                  markerfacecolor=legend_face,
                                  markeredgecolor=legend_edge,
                                  markeredgewidth=DOT_EDGE_LINEWIDTH,
                                  markersize=np.sqrt(max(100, s * 10) / np.pi),
                                  label=f'{s} genes'))
    leg_ax = fig.add_axes([0.83, 0.05, 0.06, 0.08])
    leg_ax.set_axis_off()
    leg_ax.legend(handles=legend_dots, loc='center', fontsize=8,
                  title='Gene count', title_fontsize=9,
                  framealpha=0.9, frameon=True, edgecolor='#cccccc')

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for fmt in ['svg', 'png']:
        outpath = os.path.join(OUTPUT_DIR, f'{output_prefix}_shift.{fmt}')
        fig.savefig(outpath, format=fmt, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"  Saved: {output_prefix}_shift.svg/png")


def process_dataset(name, config):
    print(f"\nProcessing {name}...")
    wt_terms = load_enrichment_tsv(config['wt_tsv'])
    mut_terms = load_enrichment_tsv(config['mut_tsv'])
    print(f"  Loaded: {len(wt_terms)} WT terms, {len(mut_terms)} MUT terms")

    all_diffs = compute_signal_diff(wt_terms, mut_terms)
    print(f"  After filtering (bg<={BG_MAX}, FDR<{FDR_THRESHOLD}): {len(all_diffs)} terms")

    wt_top, mut_top = select_top_terms(all_diffs)
    print(f"  Selected: {len(wt_top)} WT-enriched + {len(mut_top)} MUT-enriched")

    plot_horizontal_shift(wt_top, mut_top, config['title'], config['mut_label'], name)


def main():
    parser = argparse.ArgumentParser(
        description='Diverging GO shift lollipop — CC (Figure 3)')
    parser.add_argument('--only', type=str, default=None, choices=list(DATASETS.keys()))
    args = parser.parse_args()
    if args.only:
        process_dataset(args.only, DATASETS[args.only])
    else:
        for name, config in DATASETS.items():
            process_dataset(name, config)


if __name__ == '__main__':
    main()
