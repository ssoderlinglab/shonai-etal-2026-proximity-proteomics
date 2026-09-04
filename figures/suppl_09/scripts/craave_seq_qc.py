"""
Supplemental Figure 9, panels B, E, G, H, I, K and L: CrAAVe-seq multiplicity of infection, sgRNA recovery, screen quality control, Stxbp1 knockout validation and GFP control histology.

Input
    data/craave_seq/qc/quantification_table.xlsx
    data/craave_seq/qc/gRNA_recovery_estimation.xlsx
    data/craave_seq/screen/raw_counts/{Baseline,aSyn}/*.txt
    data/craave_seq/screen/library/sgRNA_library_CRISPR_KO.csv
    data/craave_seq/processed/processed_data_CRISPR_KO.csv
    data/craave_seq/western_blot/normalized_lane_table.csv
    data/craave_seq/histology/ilastik_left_right_summary.csv
Output
    figures/suppl_09/output/panel_B_moi_evaluation.{png,svg}
    figures/suppl_09/output/panel_E_gRNA_recovery.{png,svg}
    figures/suppl_09/output/panel_G_sgRNA_violin.{png,svg}
    figures/suppl_09/output/panel_H_pca.{png,svg}
    figures/suppl_09/output/panel_I_sgRNA_heatmap.{png,svg}
    figures/suppl_09/output/panel_K_munc18_wb.{png,svg}
    figures/suppl_09/output/panel_L_gfp_barplot.{png,svg}
Run
    python3 figures/suppl_09/scripts/craave_seq_qc.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib, seaborn, scikit-learn, openpyxl, src/pub_style.py

Panels A, C, D, F and J are schematics, microscopy and an immunoblot image, and are not produced here.

Panel K: Welch t-test of each Stxbp1 gRNA against the non-targeting gRNA, n = 3 lanes per condition. Panel L: paired t-test of ipsilateral against contralateral counts, n = 6 animals. Panels B, E, G, H and I are descriptive.

A single panel can be produced with --only <letter>.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
from pub_style import apply_style, FONT_SIZES, COLORS
apply_style()

# ============================================================
# DATA PATHS
# ============================================================
PATH_MOI = DATA_DIR / 'qc' / 'quantification_table.xlsx'
PATH_GRNA = DATA_DIR / 'qc' / 'gRNA_recovery_estimation.xlsx'
RAW_COUNTS_DIR = DATA_DIR / 'screen' / 'raw_counts'
LIBRARY_FILE = DATA_DIR / 'screen' / 'library' / 'sgRNA_library_CRISPR_KO.csv'
PATH_PROCESSED = DATA_DIR / 'processed' / 'processed_data_CRISPR_KO.csv'
PATH_WB = DATA_DIR / 'western_blot' / 'normalized_lane_table.csv'
PATH_GFP = DATA_DIR / 'histology' / 'ilastik_left_right_summary.csv'

# Sample order for CRISPR screen QC panels
SAMPLE_ORDER = [f"{c}_{t}_Rep{r}"
                for c in ['Baseline', 'aSyn']
                for t in ['Pre', 'Post']
                for r in [1, 2, 3, 4]]

# Colors for CRISPR screen QC
QC_COLORS = {
    'Baseline_Pre': '#9ecae1',
    'Baseline_Post': '#2171b5',
    'aSyn_Pre': '#fdae6b',
    'aSyn_Post': '#e6550d',
}
QC_GROUP_LABELS = {
    'Baseline_Pre': 'non-DA\nBaseline',
    'Baseline_Post': 'DA neuron\nBaseline',
    'aSyn_Pre': 'non-DA\naSyn',
    'aSyn_Post': 'DA neuron\naSyn',
}


def _load_crispr_counts():
    """Load raw sgRNA count matrix (shared by panels G and H)."""
    lib = pd.read_csv(LIBRARY_FILE)
    lib['NgsID'] = lib['NgsID'].astype(str).str.lstrip('0')
    lib['Gene_NgsID_Key'] = lib['Gene'] + '-' + lib['NgsID']

    all_data = {}
    for cond in ['Baseline', 'aSyn']:
        cond_dir = RAW_COUNTS_DIR / cond
        for tp in ['Pre', 'Post']:
            for rep in [1, 2, 3, 4]:
                sn = f"{cond}_{tp}_Rep{rep}"
                df = pd.read_csv(cond_dir / f"{sn}.txt", sep=r'\s+',
                                 header=None, names=['Count', 'Gene_NgsID_Key'])
                df = df.merge(lib[['Gene_NgsID_Key', 'ssOligo ID']],
                              on='Gene_NgsID_Key', how='inner')
                all_data[sn] = df.set_index('ssOligo ID')['Count']
    return pd.DataFrame(all_data)


# ============================================================
# PANEL B -- MOI Evaluation Barplot
# ============================================================
def panel_b_moi():
    """Panel B: MOI evaluation barplot.
    """
    df = pd.read_excel(PATH_MOI, sheet_name='Sheet1')

    moi_cols = ['MOI=1', 'MOI=2', 'MOI=3']
    slice_data = df[~df['slice'].isin(['average', 'SD', 'SE'])].copy()
    avg_row = df[df['slice'] == 'average'].iloc[0]
    se_row = df[df['slice'] == 'SE'].iloc[0]

    means = [avg_row[col] for col in moi_cols]
    sems = [se_row[col] for col in moi_cols]

    print("  Data Summary:")
    for col, m, s in zip(moi_cols, means, sems):
        print(f"    {col}: {m:.1f} +/- {s:.1f}% (SEM, n=8 sections)")

    fig, ax = plt.subplots(figsize=(3.2, 2.8))
    x_pos = np.arange(len(moi_cols))
    bar_colors = [COLORS['baseline'], COLORS['asyn'], COLORS['ns']]

    bars = ax.bar(x_pos, means, yerr=sems, capsize=4,
                  color=bar_colors, edgecolor='black', linewidth=0.6,
                  alpha=0.6, width=0.6,
                  error_kw={'linewidth': 0.8, 'capthick': 0.8})

    # Individual data points
    jitter_width = 0.12
    for i, col in enumerate(moi_cols):
        y_vals = slice_data[col].values
        x_jitter = x_pos[i] + np.random.default_rng(42).uniform(
            -jitter_width, jitter_width, len(y_vals))
        ax.scatter(x_jitter, y_vals, color=bar_colors[i], edgecolor='black',
                   linewidth=0.4, s=20, zorder=3, alpha=0.9)

    ax.set_ylabel('Infected cells (%)', fontsize=FONT_SIZES['label'])
    ax.set_title('MOI Evaluation', fontsize=FONT_SIZES['title'],
                 fontweight='normal', pad=0, loc='center')
    ax.text(0.5, 0.93, r'(2.5 $\times$ 10$^{11}$ GC, retro-orbital)',
            transform=ax.transAxes, ha='center', va='bottom',
            fontsize=FONT_SIZES['annotation'], color='#555555')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['MOI = 1', 'MOI = 2', r'MOI $\geq$ 3'])
    ax.set_ylim(0, 105)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.yaxis.grid(True, linestyle='-', alpha=0.15, linewidth=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.subplots_adjust(top=0.90)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_B_moi_evaluation.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: panel_B_moi_evaluation.png/svg")


# ============================================================
# PANEL E -- gRNA Recovery Dot Plot
# ============================================================
def panel_e_grna_recovery():
    """Panel E: gRNA recovery in DA neurons.
    """
    df = pd.read_excel(PATH_GRNA, sheet_name='Sheet1')

    pre = df[df['recombination'] == 'pre'].sort_values('mouse_ID').reset_index(drop=True)
    post = df[df['recombination'] == 'post'].sort_values('mouse_ID').reset_index(drop=True)
    assert list(pre['mouse_ID']) == list(post['mouse_ID']), "Mouse IDs do not match!"
    n_mice = len(pre)

    total_library = int(df['total'].iloc[0])
    recovery_pct = (post['read_counts'].mean() / total_library) * 100

    print(f"  Non-DA (Cre-): mean={pre['read_counts'].mean():.1f} +/- {pre['read_counts'].sem():.1f}")
    print(f"  DA (Cre+):     mean={post['read_counts'].mean():.1f} +/- {post['read_counts'].sem():.1f}")
    print(f"  Recovery: {recovery_pct:.1f}% of library ({post['read_counts'].mean():.0f}/{total_library})")

    fig, ax = plt.subplots(figsize=(3, 4))

    labels = ['Non-DA cells\n(Cre$^{-}$)', 'DA neurons\n(Cre$^{+}$)']
    x_positions = [0, 0.7]

    # Individual dots
    np.random.seed(42)
    for idx, (x_pos, group_df) in enumerate([(x_positions[0], pre),
                                              (x_positions[1], post)]):
        values = group_df['read_counts'].values
        jitter = (np.random.rand(len(values)) - 0.5) * 0.12
        ax.scatter(np.full(len(values), x_pos) + jitter, values,
                   color='black', s=35, alpha=0.85, zorder=2, edgecolors='none')

    # Mean + SEM bars
    mean_length = 0.25
    lw = 1.0
    cap_length = 0.10
    for x_pos, group_df in [(x_positions[0], pre), (x_positions[1], post)]:
        vals = group_df['read_counts'].values
        mean_val = vals.mean()
        sem_val = vals.std(ddof=1) / np.sqrt(len(vals))
        ax.plot([x_pos - mean_length / 2, x_pos + mean_length / 2],
                [mean_val, mean_val], color='black', linewidth=lw, zorder=3)
        ax.plot([x_pos, x_pos], [mean_val - sem_val, mean_val + sem_val],
                color='black', linewidth=lw * 0.7, zorder=3)
        for y in [mean_val - sem_val, mean_val + sem_val]:
            ax.plot([x_pos - cap_length / 2, x_pos + cap_length / 2],
                    [y, y], color='black', linewidth=lw * 0.7, zorder=3)

    # Library reference line
    ax.axhline(y=total_library, color='#aaaaaa', linestyle='--', linewidth=0.7, zorder=0)
    ax.text(x_positions[-1] + 0.18, total_library,
            f'Library\n({total_library:,})',
            va='center', ha='left', fontsize=FONT_SIZES['annotation'], color='#888888')

    ax.set_ylabel("Counts of Unique gRNAs", fontsize=FONT_SIZES['label'])
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZES['tick'])

    y_max = max(df['read_counts'].max(), total_library)
    ax.set_ylim(0, y_max * 1.12)

    y_n_pos = y_max * 0.02
    for x_pos in x_positions:
        ax.text(x_pos, y_n_pos, f'n={n_mice}', ha='center',
                fontsize=FONT_SIZES['annotation'], color='#666666')

    ax.grid(True, alpha=0.15, linestyle='-', axis='y', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(x_positions[0] - 0.35, x_positions[-1] + 0.60)
    ax.set_title("Effective gRNA recovery\nin DA neurons",
                 fontsize=FONT_SIZES['title'], pad=0, y=0.95)

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_E_gRNA_recovery.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: panel_E_gRNA_recovery.png/svg")


# ============================================================
# PANEL F -- sgRNA Count Distribution (Violin)
# ============================================================
def panel_g_violin(count_matrix):
    """Panel G: sgRNA read count distribution violin plots.
    """
    fig, ax = plt.subplots(figsize=(7.5, 4.0), facecolor='white')

    log2_data = np.log2(count_matrix[SAMPLE_ORDER] + 1)

    group_pos_ranges = {}
    positions = []
    colors_list = []
    pos = 0.0

    for cond in ['Baseline', 'aSyn']:
        for tp in ['Pre', 'Post']:
            group_key = f"{cond}_{tp}"
            group_start_pos = pos
            for rep in [1, 2, 3, 4]:
                positions.append(pos)
                colors_list.append(QC_COLORS[group_key])
                pos += 1.0
            group_pos_ranges[group_key] = (group_start_pos, pos - 1.0)
            pos += 0.6
        pos += 0.6

    data_list = [log2_data[s].values for s in SAMPLE_ORDER]
    vp = ax.violinplot(data_list, positions=positions,
                       showmedians=False, showextrema=False, widths=0.7)

    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors_list[i])
        body.set_alpha(0.75)
        body.set_edgecolor('none')

    for i, s in enumerate(SAMPLE_ORDER):
        med = np.median(log2_data[s])
        ax.scatter(positions[i], med, color='white', s=14, zorder=5,
                   edgecolors='#333333', linewidths=0.5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"Rep.{r}" for _ in range(4) for r in [1, 2, 3, 4]],
                       fontsize=6)
    ax.tick_params(axis='x', length=0, pad=2)

    for group_key, (gstart, gend) in group_pos_ranges.items():
        mid = (gstart + gend) / 2
        inv = ax.transData + ax.transAxes.inverted()
        x_frac = inv.transform((mid, 0))[0]
        ax.text(x_frac, -0.10, QC_GROUP_LABELS[group_key],
                transform=ax.transAxes,
                ha='center', va='top', fontsize=7.5, fontweight='bold',
                color=QC_COLORS[group_key], linespacing=1.15)

    ax.set_ylabel('log$_2$(read count + 1)')

    gap_pos = (group_pos_ranges['Baseline_Post'][1] +
               group_pos_ranges['aSyn_Pre'][0]) / 2
    ax.axvline(gap_pos, color='#dddddd', linestyle='-', linewidth=0.6)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(positions[0] - 0.6, positions[-1] + 0.6)

    plt.subplots_adjust(bottom=0.22)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_G_sgRNA_violin.{ext}',
                    dpi=300, facecolor='white', bbox_inches='tight')
    plt.close()
    print("  Saved: panel_G_sgRNA_violin.png/svg")


# ============================================================
# PANEL G -- PCA of sgRNA Counts
# ============================================================
def panel_h_pca(count_matrix):
    """Panel H: PCA of CPM-normalized sgRNA counts.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    raw = count_matrix[SAMPLE_ORDER]
    cpm = raw.div(raw.sum(axis=0), axis=1) * 1e6
    log2_cpm = np.log2(cpm + 1).T

    scaler = StandardScaler()
    scaled = scaler.fit_transform(log2_cpm)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)

    x_min, x_max = pcs[:, 0].min(), pcs[:, 0].max()
    y_min, y_max = pcs[:, 1].min(), pcs[:, 1].max()
    x_pad = (x_max - x_min) * 0.12
    y_pad = (y_max - y_min) * 0.12

    fig, ax = plt.subplots(figsize=(4.5, 4), facecolor='white')

    for i, sample in enumerate(SAMPLE_ORDER):
        parts = sample.split('_')
        cond, tp = parts[0], parts[1]
        group_key = f"{cond}_{tp}"
        color = QC_COLORS[group_key]
        marker = 'o' if tp == 'Pre' else '^'
        ax.scatter(pcs[i, 0], pcs[i, 1], c=color, marker=marker, s=90,
                   edgecolors='white', linewidths=0.7, zorder=3)

    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')

    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=QC_COLORS['Baseline_Pre'],
               markersize=8, markeredgecolor='white', markeredgewidth=0.5,
               label='non-DA, Baseline'),
        Line2D([0], [0], marker='^', color='w',
               markerfacecolor=QC_COLORS['Baseline_Post'],
               markersize=8, markeredgecolor='white', markeredgewidth=0.5,
               label='DA neuron, Baseline'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=QC_COLORS['aSyn_Pre'],
               markersize=8, markeredgecolor='white', markeredgewidth=0.5,
               label='non-DA, aSyn'),
        Line2D([0], [0], marker='^', color='w',
               markerfacecolor=QC_COLORS['aSyn_Post'],
               markersize=8, markeredgecolor='white', markeredgewidth=0.5,
               label='DA neuron, aSyn'),
    ]
    ax.legend(handles=legend_elements, fontsize=7.5, loc='lower right',
              framealpha=0.9, edgecolor='#cccccc')

    ax.axhline(0, color='#eeeeee', linewidth=0.4)
    ax.axvline(0, color='#eeeeee', linewidth=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_H_pca.{ext}',
                    dpi=300, facecolor='white', bbox_inches='tight')
    plt.close()

    print(f"  PC1={pca.explained_variance_ratio_[0]*100:.1f}%, "
          f"PC2={pca.explained_variance_ratio_[1]*100:.1f}%")
    print("  Saved: panel_H_pca.png/svg")


# ============================================================
# PANEL H -- sgRNA Heatmap (Combined Layout)
# ============================================================
def panel_i_heatmap():
    """Panel I: sgRNA heatmap (Baseline | aSyn).
    """
    import seaborn as sns
    from matplotlib.gridspec import GridSpec

    df = pd.read_csv(PATH_PROCESSED)
    rep_map = {'Rep.1': 'Rep1', 'Rep.2': 'Rep2', 'Rep.3': 'Rep3', 'Rep.4': 'Rep4'}
    df['rep_clean'] = df['Replicate'].map(rep_map)

    heatmap_data = {}
    for cond in ['aSyn', 'Baseline']:
        sub = df[df['Condition'] == cond].copy()
        pre = sub.pivot(index='sgRNA_ID', columns='rep_clean', values='Proportion_Pre')
        pre.columns = [f'Pre_{c}' for c in sorted(pre.columns)]
        pre = pre[sorted(pre.columns)]
        post = sub.pivot(index='sgRNA_ID', columns='rep_clean', values='Proportion_Post')
        post.columns = [f'Post_{c}' for c in sorted(post.columns)]
        post = post[sorted(post.columns)]
        ratio = sub.pivot(index='sgRNA_ID', columns='rep_clean', values='Ratio_Post_Pre')
        ratio.columns = [f'Ratio_{c}' for c in sorted(ratio.columns)]
        ratio = ratio[sorted(ratio.columns)]
        heatmap_data[cond] = pd.concat([pre, post, ratio], axis=1)

    # Sort by aSyn Ratio mean (descending)
    ratio_cols = [c for c in heatmap_data['aSyn'].columns if c.startswith('Ratio_')]
    asyn_ratio_mean = heatmap_data['aSyn'][ratio_cols].mean(axis=1).sort_values(ascending=False)
    sort_order = asyn_ratio_mean.index
    for cond in ['aSyn', 'Baseline']:
        heatmap_data[cond] = heatmap_data[cond].reindex(sort_order)

    # Shared color scales
    all_prop = []
    for cond in ['aSyn', 'Baseline']:
        pcols = [c for c in heatmap_data[cond].columns if c.startswith(('Pre_', 'Post_'))]
        all_prop.append(heatmap_data[cond][pcols])
    all_prop_df = pd.concat(all_prop)
    prop_vmax = all_prop_df.mean().mean() + 2 * all_prop_df.std().mean()

    all_ratio = []
    for cond in ['aSyn', 'Baseline']:
        rcols = [c for c in heatmap_data[cond].columns if c.startswith('Ratio_')]
        all_ratio.append(np.log10(heatmap_data[cond][rcols].replace(0, np.nan).fillna(1e-10)))
    all_ratio_df = pd.concat(all_ratio)
    ratio_center = all_ratio_df.mean().mean()
    ratio_spread = 2 * all_ratio_df.std().mean()
    ratio_vmin = ratio_center - ratio_spread
    ratio_vmax = ratio_center + ratio_spread

    # Font sizes
    F_CONDITION = 16
    F_SECTION = 11
    F_SUBSECTION = 9
    F_REP = 8
    F_YAXIS = 12
    F_CBAR_LABEL = 10
    F_CBAR_TICK = 8

    fig = plt.figure(figsize=(8, 4))
    gs = GridSpec(1, 5, figure=fig, width_ratios=[8, 4, 1.2, 8, 4], wspace=0.03)

    panel_specs = [
        ('Baseline', 'prop', 0),
        ('Baseline', 'ratio', 1),
        ('aSyn', 'prop', 3),
        ('aSyn', 'ratio', 4),
    ]
    rep_labels = ['1', '2', '3', '4']

    axes = {}
    for cond, ptype, col in panel_specs:
        ax = fig.add_subplot(gs[0, col])
        axes[(cond, ptype)] = ax
        data = heatmap_data[cond]
        n_rows = len(data)

        if ptype == 'prop':
            pre_cols = sorted([c for c in data.columns if c.startswith('Pre_')])
            post_cols = sorted([c for c in data.columns if c.startswith('Post_')])
            plot_data = data[pre_cols + post_cols]

            sns.heatmap(plot_data, ax=ax, cmap='coolwarm',
                        vmin=0, vmax=prop_vmax,
                        cbar=False, yticklabels=False,
                        xticklabels=rep_labels + rep_labels,
                        rasterized=True)
            ax.tick_params(axis='x', labelsize=F_REP, rotation=0)

            n_pre = len(pre_cols)
            n_post = len(post_cols)

            ax.text(n_pre / 2, n_rows + 40, r'Cre$^-$',
                    ha='center', va='top', fontsize=F_SECTION, fontweight='bold')
            ax.text(n_pre / 2, n_rows + 65, '(Input)',
                    ha='center', va='top', fontsize=F_SUBSECTION, fontstyle='italic')
            ax.text(n_pre + n_post / 2, n_rows + 40, r'Cre$^+$',
                    ha='center', va='top', fontsize=F_SECTION, fontweight='bold')
            ax.text(n_pre + n_post / 2, n_rows + 65, '(Dopamine neuron)',
                    ha='center', va='top', fontsize=F_SUBSECTION, fontstyle='italic')
            ax.axvline(n_pre, color='white', linewidth=2)

        else:  # ratio
            rcols = sorted([c for c in data.columns if c.startswith('Ratio_')])
            plot_data = np.log10(data[rcols].replace(0, np.nan).fillna(1e-10))
            sns.heatmap(plot_data, ax=ax, cmap='PRGn',
                        vmin=ratio_vmin, vmax=ratio_vmax,
                        cbar=False, yticklabels=False,
                        xticklabels=rep_labels,
                        rasterized=True)
            ax.tick_params(axis='x', labelsize=F_REP, rotation=0)
            ax.text(len(rcols) / 2, n_rows + 40, 'gRNA Effect',
                    ha='center', va='top', fontsize=F_SECTION, fontweight='bold')

        ax.set_xlabel('')
        ax.set_ylabel('')

    # Gap axis
    gap_ax = fig.add_subplot(gs[0, 2])
    gap_ax.axis('off')

    axes[('Baseline', 'prop')].set_title('Baseline (no virus)', fontsize=F_CONDITION,
                                          x=0.75, fontweight='bold', pad=12)
    axes[('aSyn', 'prop')].set_title('A53T-h-aSYN', fontsize=F_CONDITION,
                                      x=0.75, fontweight='bold', pad=12)
    axes[('Baseline', 'prop')].set_ylabel('Unique sgRNA', fontsize=F_YAXIS * 1.5,
                                           fontweight='bold')

    # Colorbars
    cbar_prop_ax = fig.add_axes([0.92, 0.55, 0.013, 0.28])
    sm_prop = plt.cm.ScalarMappable(cmap='coolwarm',
                                     norm=plt.Normalize(vmin=0, vmax=prop_vmax))
    fig.colorbar(sm_prop, cax=cbar_prop_ax).ax.tick_params(labelsize=F_CBAR_TICK)

    cbar_ratio_ax = fig.add_axes([0.92, 0.14, 0.013, 0.28])
    sm_ratio = plt.cm.ScalarMappable(cmap='PRGn',
                                      norm=plt.Normalize(vmin=ratio_vmin, vmax=ratio_vmax))
    fig.colorbar(sm_ratio, cax=cbar_ratio_ax).ax.tick_params(labelsize=F_CBAR_TICK)

    CBAR_LABEL_X = 0.975
    fig.text(CBAR_LABEL_X, 0.55 + 0.28 / 2, 'Proportion within library',
             rotation=90, va='center', ha='left', fontsize=F_CBAR_LABEL)
    fig.text(CBAR_LABEL_X, 0.14 + 0.28 / 2, r'$\log_{10}$(Cre$^+$/Cre$^-$)',
             rotation=90, va='center', ha='left', fontsize=F_CBAR_LABEL)

    plt.subplots_adjust(right=0.90, bottom=0.12, top=0.92)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_I_sgRNA_heatmap.{ext}',
                    dpi=300, bbox_inches='tight', transparent=True)
    plt.close()
    print("  Saved: panel_I_sgRNA_heatmap.png/svg")


# ============================================================
# PANEL J -- Munc18 WB Quantification
# ============================================================
def panel_k_wb_quant():
    """Panel K: Munc18 WB quantification (Stxbp1-KO validation).
    """
    from scipy import stats

    df = pd.read_csv(PATH_WB)
    df = df[df['loading'] == 'sample'].dropna(subset=['condition']).copy()
    plot_col = 'normalized_to_nt_mean'

    COLOR_NT = '#999999'
    COLOR_KO = '#4A4A4A'

    conditions = [
        {'key': 'nt',   'label': 'non-target\nsgRNA',    'color': COLOR_NT},
        {'key': 'del3', 'label': 'Stxbp1-KO\n(del3)',    'color': COLOR_KO},
        {'key': 'del6', 'label': 'Stxbp1-KO\n(del6)',    'color': COLOR_KO},
    ]

    for c in conditions:
        vals = df[df['condition'] == c['key']][plot_col].values
        c['data'] = vals
        c['n'] = len(vals)
        c['mean'] = vals.mean()
        c['sem'] = vals.std(ddof=1) / np.sqrt(len(vals))
        print(f"    {c['key']} (n={c['n']}): mean={c['mean']:.3f} +/- {c['sem']:.3f}")

    # Welch t-tests vs nt
    nt_data = conditions[0]['data']
    for c in conditions[1:]:
        t_stat, p_val = stats.ttest_ind(nt_data, c['data'], equal_var=False)
        c['p_vs_nt'] = p_val
        print(f"    nt vs {c['key']}: Welch t p={p_val:.4f}")

    def fmt_p(p):
        if p < 0.001:
            return 'p < 0.001'
        return f'p = {p:.3f}'

    def bracket(ax, x1, x2, y, h, text, fs=7):
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
                color='black', lw=0.8, zorder=5, clip_on=False)
        ax.text((x1+x2)/2, y+h+0.01, text,
                ha='center', va='bottom', fontsize=fs, zorder=5)

    fig, ax = plt.subplots(figsize=(3, 4))
    positions = [0, 1, 2]
    bar_width = 0.6

    for i, c in enumerate(conditions):
        ax.bar(positions[i], c['mean'], width=bar_width, color=c['color'],
               edgecolor='none', zorder=2)
        ax.errorbar(positions[i], c['mean'], yerr=c['sem'], fmt='none',
                    color='black', capsize=4, capthick=1.0, linewidth=1.0, zorder=3)
        np.random.seed(700 + i)
        jitter = np.random.uniform(-bar_width * 0.25, bar_width * 0.25, size=c['n'])
        ax.scatter(positions[i] + jitter, c['data'], s=25, color='black',
                   edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)

    ax.axhline(y=1.0, color='gray', linestyle='--', linewidth=0.6, alpha=0.5, zorder=1)

    y_top = max(max(c['data']) for c in conditions)
    bracket(ax, positions[0], positions[1], y_top + 0.08, 0.04,
            fmt_p(conditions[1]['p_vs_nt']))
    bracket(ax, positions[0], positions[2], y_top + 0.22, 0.04,
            fmt_p(conditions[2]['p_vs_nt']))

    ax.set_xticks(positions)
    ax.set_xticklabels([c['label'] for c in conditions], fontsize=FONT_SIZES['tick'])
    ax.set_ylabel('Munc18 / total protein', fontsize=FONT_SIZES['label'])
    ax.set_ylim(0, y_top + 0.40)
    ax.set_title('Stxbp1/Munc18-1 KO Validation',
                 fontsize=FONT_SIZES['title'], fontweight='bold')

    fig.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_K_munc18_wb.{ext}',
                    dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    print("  Saved: panel_K_munc18_wb.png/svg")


# ============================================================
# PANEL K -- GFP Control Barplot
# ============================================================
def panel_l_gfp():
    """Panel L: GFP control barplot (TH+ cell count).
    """
    from scipy import stats

    COLOR_CONTRA = '#80cdc1'
    COLOR_IPSI = '#018571'

    gfp_raw = pd.read_csv(PATH_GFP)
    gfp_raw = gfp_raw[gfp_raw['Animal'].notna() & (gfp_raw['Animal'] != 'NA')].copy()
    gfp_raw['Animal'] = gfp_raw['Animal'].astype(int)
    gfp_animal = gfp_raw.groupby(['Animal', 'HemisphereSide'])['cells_predicted'].mean().reset_index()
    gfp_wide = gfp_animal.pivot_table(index='Animal', columns='HemisphereSide',
                                       values='cells_predicted').reset_index()

    n = len(gfp_wide)
    contra = gfp_wide['Left'].values
    ipsi = gfp_wide['Right'].values

    t_stat, p_val = stats.ttest_rel(ipsi, contra)

    print(f"    Stxbp1-KO + GFP (n={n}): contra={contra.mean():.1f}+/-{contra.std(ddof=1)/np.sqrt(n):.1f}, "
          f"ipsi={ipsi.mean():.1f}+/-{ipsi.std(ddof=1)/np.sqrt(n):.1f}")
    print(f"    paired t p={p_val:.4f}")

    def p_to_stars(p):
        if p < 0.001: return '***'
        elif p < 0.01: return '**'
        elif p < 0.05: return '*'
        return 'n.s.'

    def bracket(ax, x1, x2, y, h, text, fs=7):
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y],
                color='black', lw=0.8, zorder=5, clip_on=False)
        ax.text((x1+x2)/2, y+h+1, text,
                ha='center', va='bottom', fontsize=fs, zorder=5)

    fig, ax = plt.subplots(figsize=(2, 4))

    bar_width = 0.35
    pair_gap = 0.02
    center = 0.5
    x_contra = center - (bar_width / 2 + pair_gap / 2)
    x_ipsi = center + (bar_width / 2 + pair_gap / 2)

    means = [contra.mean(), ipsi.mean()]
    sems = [contra.std(ddof=1)/np.sqrt(n), ipsi.std(ddof=1)/np.sqrt(n)]

    ax.bar(x_contra, means[0], width=bar_width, color=COLOR_CONTRA, edgecolor='none', zorder=2)
    ax.errorbar(x_contra, means[0], yerr=sems[0], fmt='none', color='black',
                capsize=4, capthick=1.0, linewidth=1.0, zorder=3)
    ax.bar(x_ipsi, means[1], width=bar_width, color=COLOR_IPSI, edgecolor='none', zorder=2)
    ax.errorbar(x_ipsi, means[1], yerr=sems[1], fmt='none', color='black',
                capsize=4, capthick=1.0, linewidth=1.0, zorder=3)

    np.random.seed(200)
    jitter_c = np.random.uniform(-bar_width * 0.3, bar_width * 0.3, size=n)
    ax.scatter(x_contra + jitter_c, contra, s=18, color='black',
               edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)
    np.random.seed(201)
    jitter_i = np.random.uniform(-bar_width * 0.3, bar_width * 0.3, size=n)
    ax.scatter(x_ipsi + jitter_i, ipsi, s=18, color='black',
               edgecolors='black', linewidth=0.3, zorder=4, alpha=0.7)

    top = max(means[0]+sems[0], means[1]+sems[1])
    bracket(ax, x_contra, x_ipsi, top + 30, 5, p_to_stars(p_val))

    ax.set_xticks([x_contra, x_ipsi])
    ax.set_xticklabels(['KO+GFP | contra', 'KO+GFP | ipsi'],
                        fontsize=FONT_SIZES['tick'], rotation=45, ha='right')
    ax.set_ylabel(r'$\mathrm{TH}^{+}$ Cell Count', fontsize=FONT_SIZES['label'])
    ax.set_ylim(0, 320)
    ax.set_title('AAV2/1-CMV-GFP', fontsize=FONT_SIZES['title'], fontweight='bold')

    plt.tight_layout()

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'panel_L_gfp_barplot.{ext}',
                    dpi=300, transparent=True, bbox_inches='tight')
    plt.close()
    print("  Saved: panel_L_gfp_barplot.png/svg")


# ============================================================
# MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Supplemental Figure 9 -- CrAAVe-seq quality control and validation')
    parser.add_argument('--only', choices=['B', 'E', 'G', 'H', 'I', 'K', 'L'],
                        help='Generate single panel (default: all)')
    args = parser.parse_args()

    print("Supplemental Figure 9 -- CrAAVe-seq quality control and validation")

    panels = args.only or 'ALL'

    # Count matrix shared by panels G and H
    count_matrix = None
    if panels in ('ALL', 'G', 'H'):
        print("\nLoading CRISPR screen count matrix...")
        count_matrix = _load_crispr_counts()
        print(f"  Shape: {count_matrix.shape}")

    if panels in ('ALL', 'B'):
        print("\nPanel B -- MOI evaluation:")
        panel_b_moi()

    if panels in ('ALL', 'E'):
        print("\nPanel E -- gRNA recovery:")
        panel_e_grna_recovery()

    if panels in ('ALL', 'G'):
        print("\nPanel G -- sgRNA count distribution:")
        panel_g_violin(count_matrix)

    if panels in ('ALL', 'H'):
        print("\nPanel H -- PCA:")
        panel_h_pca(count_matrix)

    if panels in ('ALL', 'I'):
        print("\nPanel I -- sgRNA heatmap:")
        panel_i_heatmap()

    if panels in ('ALL', 'K'):
        print("\nPanel K -- Munc18 western blot quantification:")
        panel_k_wb_quant()

    if panels in ('ALL', 'L'):
        print("\nPanel L -- GFP control:")
        panel_l_gfp()

    print("\nDone.")


if __name__ == '__main__':
    main()
