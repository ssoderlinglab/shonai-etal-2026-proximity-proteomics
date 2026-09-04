"""
Supplemental Figure 5, panels D and E: dopamine neuron marker enrichment heatmap and Cochran-Mantel-Haenszel forest plot.

Input
    data/ibiocofrac/pilot/20230130_pilot.csv
    data/ibiocofrac/markers/{eLife_APEX_DA_marker,20230228_DANmarker_for_Venn}.csv
    data/external/annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv
    data/external/PanglaoDB_markers_27_Mar_2020.tsv
Output
    figures/suppl_05/output/*.{png,svg,csv}
Run
    python3 figures/suppl_05/scripts/dan_enrichment_cmh.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib

Fractions 4, 7 and 10. Enrichment threshold log2(TurboID/Control) > 1.0. Permutation test with 10,000 draws.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'ibiocofrac'
EXTERNAL_DIR = REPO_ROOT / 'data' / 'external'
OUTPUT_DIR = FIGURE_DIR / 'output'
SRC_DIR = REPO_ROOT / 'src'

sys.path.insert(0, str(SRC_DIR))
from pub_style import apply_style, FONT_SIZES
apply_style()

PILOT_PATH = DATA_DIR / 'pilot' / '20230130_pilot.csv'
HOBSON_PATH = DATA_DIR / 'markers' / 'eLife_APEX_DA_marker.csv'
VENN_PATH = DATA_DIR / 'markers' / '20230228_DANmarker_for_Venn.csv'
SAUNDERS_PATH = EXTERNAL_DIR / 'annotation.BrainCellAtlas_Saunders_version_2018.04.01.csv'
PANGLAO_PATH = EXTERNAL_DIR / 'PanglaoDB_markers_27_Mar_2020.tsv'


# ============================================================
# Data loading
# ============================================================
def load_pilot():
    pilot = pd.read_csv(PILOT_PATH)
    pilot['GN'] = pilot['description'].str.extract(r'GN=([A-Za-z0-9\.]+)')
    for f in [4, 7, 10]:
        pilot[f'C{f}_ps'] = pilot[f'C{f}'] + 1
        pilot[f'T{f}_ps'] = pilot[f'T{f}'] + 1
        pilot[f'log2_TC_F{f}'] = np.log2(pilot[f'T{f}_ps'] / pilot[f'C{f}_ps'])
    return pilot


def load_dan_markers():
    hobson = pd.read_csv(HOBSON_PATH)
    dan_hobson = set(hobson['GN'].dropna().unique())
    venn = pd.read_csv(VENN_PATH)
    dan_saunders = set(venn[venn['dataset'] == 'Saunders_2018']['gene'].unique())
    return dan_hobson | dan_saunders


def build_marker_lists_v2(detected_genes):
    """Build cell-type marker lists from Saunders + PanglaoDB (exclusive, same-granularity)."""
    dan_combined = load_dan_markers()
    detected_upper = {g.upper(): g for g in detected_genes}

    # Saunders markers from SN
    saunders = pd.read_csv(SAUNDERS_PATH, index_col=0, encoding='latin1',
                           on_bad_lines='skip')
    sn = saunders[saunders['tissue'] == 'SN']
    saunders_markers = {
        'GABAergic Neurons': set(), 'Glutamatergic Neurons': set(),
        'Astrocytes': set(), 'Oligodendrocytes': set(), 'Microglia': set()
    }
    for _, row in sn.iterrows():
        cls = row['class']
        cm = str(row['class_marker'])
        tm = str(row['type_marker'])
        if cls == 'NEURON':
            if 'Gad1' in cm or 'Gad2' in cm:
                for g in ['Gad1', 'Gad2']:
                    saunders_markers['GABAergic Neurons'].add(g)
                for t in (tm.split('-') if '-' in tm else [tm]):
                    saunders_markers['GABAergic Neurons'].add(t)
            elif cm in ['Slc17a6', 'Slc17a7']:
                saunders_markers['Glutamatergic Neurons'].add(cm)
                for t in (tm.split('-') if '-' in tm else [tm]):
                    saunders_markers['Glutamatergic Neurons'].add(t)
        elif cls == 'ASTROCYTE':
            saunders_markers['Astrocytes'].update([cm, tm])
        elif cls == 'OLIGODENDROCYTE':
            saunders_markers['Oligodendrocytes'].update([cm, tm])
        elif cls == 'MICROGLIA':
            saunders_markers['Microglia'].update([cm, tm])

    # PanglaoDB exclusive markers
    pdb = pd.read_csv(PANGLAO_PATH, sep='\t')
    brain_mm = pdb[(pdb['species'].str.contains('Mm', na=False)) &
                   (pdb['organ'].str.contains('Brain', na=False))].copy()
    target_pdb_cts = ['Dopaminergic neurons', 'GABAergic neurons',
                      'Glutaminergic neurons', 'Astrocytes',
                      'Oligodendrocytes', 'Microglia']
    gene_target_map = {}
    for ct in target_pdb_cts:
        for g in set(brain_mm[brain_mm['cell type'] == ct][
                'official gene symbol'].str.upper()):
            if g not in gene_target_map:
                gene_target_map[g] = set()
            gene_target_map[g].add(ct)

    pdb_ct_mapping = {
        'GABAergic neurons': 'GABAergic Neurons',
        'Glutaminergic neurons': 'Glutamatergic Neurons',
        'Astrocytes': 'Astrocytes',
        'Oligodendrocytes': 'Oligodendrocytes',
        'Microglia': 'Microglia',
    }
    non_dan_v2 = {}
    for pdb_ct, our_ct in pdb_ct_mapping.items():
        genes = set(saunders_markers.get(our_ct, set()))
        sub = brain_mm[brain_mm['cell type'] == pdb_ct]
        for _, r in sub.iterrows():
            g_upper = str(r['official gene symbol']).upper()
            if len(gene_target_map.get(g_upper, set()) - {pdb_ct}) == 0:
                if g_upper in detected_upper:
                    genes.add(detected_upper[g_upper])
                else:
                    genes.add(str(r['official gene symbol']))
        non_dan_v2[our_ct] = genes
    non_dan_v2['Dopaminergic Neurons'] = dan_combined
    return non_dan_v2


# ============================================================
# CMH test
# ============================================================
def cmh_full(gene_set, pilot_df, fractions=[4, 7, 10], threshold=1.0):
    tables = []
    for frac in fractions:
        col = f'log2_TC_F{frac}'
        is_enr = pilot_df[col] > threshold
        is_marker = pilot_df['GN'].isin(gene_set)
        a = float((is_enr & is_marker).sum())
        b = float((is_enr & ~is_marker).sum())
        c = float((~is_enr & is_marker).sum())
        d = float((~is_enr & ~is_marker).sum())
        tables.append((a, b, c, d))
    numerator = 0; denominator_var = 0; or_num = 0; or_den = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        r1 = a + b; r2 = c + d; c1 = a + c; c2 = b + d
        E_a = r1 * c1 / n
        V_a = r1 * r2 * c1 * c2 / (n * n * (n - 1)) if n > 1 else 0
        numerator += (a - E_a)
        denominator_var += V_a
        or_num += a * d / n
        or_den += b * c / n
    if denominator_var == 0:
        return np.nan, np.nan, np.nan, np.nan
    chi2 = (abs(numerator) - 0.5) ** 2 / denominator_var
    p_value = 1 - stats.chi2.cdf(chi2, df=1)
    common_or = or_num / or_den if or_den > 0 else np.inf
    ln_or = np.log(common_or) if common_or > 0 and np.isfinite(common_or) else np.nan
    P_num = 0; P_den = 0; Q_num = 0; Q_den = 0; PQ = 0
    for a, b, c, d in tables:
        n = a + b + c + d
        if n == 0:
            continue
        R = a * d / n; S = b * c / n
        P_num += (a + d) * a * d / (n * n)
        P_den += R
        Q_num += (b + c) * b * c / (n * n)
        Q_den += S
        PQ += ((a + d) * b * c + (b + c) * a * d) / (n * n)
    if P_den > 0 and Q_den > 0:
        se = np.sqrt(P_num / (2 * P_den**2) + PQ / (2 * P_den * Q_den)
                     + Q_num / (2 * Q_den**2))
        ci_low = np.exp(ln_or - 1.96 * se)
        ci_high = np.exp(ln_or + 1.96 * se)
    else:
        ci_low = ci_high = np.nan
    return common_or, ci_low, ci_high, p_value


def cmh_or_fast(gene_set, pilot_df, fractions=[4, 7, 10], threshold=1.0):
    or_num = 0; or_den = 0
    for frac in fractions:
        col = f'log2_TC_F{frac}'
        is_enr = pilot_df[col] > threshold
        is_marker = pilot_df['GN'].isin(gene_set)
        a = float((is_enr & is_marker).sum())
        b = float((is_enr & ~is_marker).sum())
        c = float((~is_enr & is_marker).sum())
        d = float((~is_enr & ~is_marker).sum())
        n = a + b + c + d
        if n == 0:
            continue
        or_num += a * d / n
        or_den += b * c / n
    if or_den == 0:
        return np.inf if or_num > 0 else np.nan
    return or_num / or_den


def p_to_stars(p):
    if np.isnan(p): return ''
    if p < 0.001: return '***'
    if p < 0.01: return '**'
    if p < 0.05: return '*'
    return 'n.s.'


# ============================================================
# Panel D: DAN Heatmap
# ============================================================
def panel_d_heatmap(pilot, dan_combined):
    pilot['is_DAN'] = pilot['GN'].isin(dan_combined)
    dan_detected = pilot[pilot['is_DAN'] & pilot['GN'].notna()].copy()
    for f in [4, 7, 10]:
        dan_detected[f'log2_TC_F{f}'] = np.log2(
            (dan_detected[f'T{f}'] + 1) / (dan_detected[f'C{f}'] + 1))
    dan_detected['log2_TC_mean'] = dan_detected[
        ['log2_TC_F4', 'log2_TC_F7', 'log2_TC_F10']].mean(axis=1)
    dan_detected = dan_detected.sort_values('log2_TC_mean', ascending=False)
    n_total = len(dan_detected)

    print(f"  DAN markers detected: {n_total}")

    heatmap_raw = dan_detected.set_index('GN')[
        ['C4', 'T4', 'C7', 'T7', 'C10', 'T10']].copy()
    heatmap_log = np.log10(heatmap_raw + 1)
    col_labels = ['Ctrl\nF4', 'TurboID\nF4', 'Ctrl\nF7',
                  'TurboID\nF7', 'Ctrl\nF10', 'TurboID\nF10']

    fig_height = max(9, n_total * 0.3)
    fig, ax = plt.subplots(figsize=(6.5, fig_height))
    sns.heatmap(heatmap_log, ax=ax, cmap='Reds', vmin=0, annot=False,
                linewidths=0.8, linecolor='white',
                cbar_kws={'label': 'Log$_{10}$(Intensity + 1)', 'shrink': 0.4},
                xticklabels=col_labels)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8.5,
                       fontstyle='italic', rotation=0)
    ax.set_xticklabels(col_labels, fontsize=8.5, rotation=0, ha='center')
    for x_pos in [2, 4]:
        ax.axvline(x=x_pos, color='black', linewidth=1.5, zorder=10)
    ax.set_title('DAN Marker Enrichment by TurboID',
                 fontsize=FONT_SIZES['title'], fontweight='bold', pad=10)
    ax.set_ylabel('')

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig.tight_layout()
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'dan_heatmap.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: dan_heatmap.png/svg")


# ============================================================
# Panel E: CMH Forest Plot
# ============================================================
def panel_e_forest(pilot, marker_lists):
    detected_genes = set(pilot['GN'].dropna().unique())
    detected_list = list(detected_genes)

    ct_plot = ['Dopaminergic Neurons', 'Glutamatergic Neurons',
               'GABAergic Neurons', 'Astrocytes', 'Oligodendrocytes']

    # Compute observed ORs
    obs = {}
    for ct in ct_plot:
        genes = marker_lists[ct]
        n_det = len(genes & detected_genes)
        if n_det == 0:
            obs[ct] = {'or': np.nan, 'ci_lo': np.nan, 'ci_hi': np.nan,
                        'p': np.nan, 'n_det': 0}
        else:
            or_val, ci_lo, ci_hi, pval = cmh_full(genes, pilot)
            obs[ct] = {'or': or_val, 'ci_lo': ci_lo, 'ci_hi': ci_hi,
                        'p': pval, 'n_det': n_det}

    # Permutation test for DAN
    np.random.seed(42)
    N_PERM = 10000
    dan_detected = marker_lists['Dopaminergic Neurons'] & detected_genes
    n_dan = len(dan_detected)
    perm_ors = []
    for _ in range(N_PERM):
        rg = set(np.random.choice(detected_list, size=n_dan, replace=False))
        ov = cmh_or_fast(rg, pilot)
        if np.isfinite(ov):
            perm_ors.append(ov)
    perm_ors = np.array(perm_ors)
    perm_p = np.mean(perm_ors >= obs['Dopaminergic Neurons']['or'])
    perm_p_str = "p < 0.0001" if perm_p == 0 else f"p = {perm_p:.4f}"

    print(f"  Permutation: mean={perm_ors.mean():.3f}, "
          f"95th={np.percentile(perm_ors, 95):.3f}")
    print(f"  DAN OR={obs['Dopaminergic Neurons']['or']:.3f}, "
          f"perm {perm_p_str}")

    # Print all results
    print("\n  Cell-type enrichment results:")
    for ct in ct_plot:
        d = obs[ct]
        print(f"    {ct}: OR={d['or']:.2f}, p={d['p']:.4f} "
              f"{p_to_stars(d['p'])}, n_det={d['n_det']}")

    # Plot
    ct_colors = {
        'Dopaminergic Neurons': '#C0392B',
        'Glutamatergic Neurons': '#999999',
        'GABAergic Neurons': '#999999',
        'Astrocytes': '#999999',
        'Oligodendrocytes': '#999999',
    }
    y = np.arange(len(ct_plot))
    fig, ax = plt.subplots(figsize=(7.5, 3.8))

    for i, ct in enumerate(ct_plot):
        d = obs[ct]
        if d['n_det'] == 0:
            continue
        or_val = d['or']
        if or_val == 0:
            or_val = 0.03
        ci_lo = min(d['ci_lo'], or_val) if not np.isnan(d['ci_lo']) else 0.01
        ci_hi = max(d['ci_hi'], or_val) if not np.isnan(d['ci_hi']) else 200
        ci_lo = max(ci_lo, 0.005)
        ci_hi = min(ci_hi, 300)
        color = ct_colors[ct]
        ms = 7 + d['n_det'] * 0.12

        ax.errorbar(or_val, i,
                    xerr=[[or_val - ci_lo], [ci_hi - or_val]],
                    fmt='D', color=color, markersize=ms,
                    markeredgecolor='black', markeredgewidth=0.6,
                    ecolor=color, elinewidth=1.8, capsize=3, zorder=3)

        sig = p_to_stars(d['p'])
        or_str = f"OR = {or_val:.2f} [{ci_lo:.1f}\u2013{ci_hi:.1f}]"
        p_str = f"p = {d['p']:.1e}" if not np.isnan(d['p']) else ''

        if ct == 'Dopaminergic Neurons':
            ax.text(ci_hi * 1.5, i - 0.08, f'{or_str}', va='center',
                    fontsize=8, fontweight='bold', color='#333333')
            ax.text(ci_hi * 1.5, i + 0.25,
                    f'{p_str} {sig}  (perm. {perm_p_str})',
                    va='center', fontsize=7.5, fontweight='bold',
                    color='#C0392B')
        else:
            ax.text(ci_hi * 1.5, i, f'{or_str}\n{p_str} {sig}',
                    va='center', fontsize=8,
                    fontweight='bold' if d['p'] < 0.05 else 'normal',
                    color='#333333')

    ax.axvline(x=1, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.set_xscale('log')
    ax.set_xlim(0.01, 500)
    ax.set_yticks(y)
    ax.set_yticklabels(ct_plot, fontsize=FONT_SIZES['tick'])
    ax.invert_yaxis()
    ax.set_xlabel('Odds Ratio (95% CI)', fontsize=FONT_SIZES['label'])
    ax.set_title('Cell-Type Enrichment in TurboID Proteomics\n'
                 '(Cochran-Mantel-Haenszel test)',
                 fontsize=FONT_SIZES['title'], fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    for ext in ['svg', 'png']:
        fig.savefig(OUTPUT_DIR / f'cmh_forest.{ext}',
                    dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print("  Saved: cmh_forest.png/svg")


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Suppl Figure 2 — DAN heatmap + CMH forest plot')
    parser.add_argument('--only', choices=['D', 'E'],
                        help='Generate single panel (default: both)')
    args = parser.parse_args()

    print("=" * 60)
    print("SUPPL FIGURE 2 — DAN Enrichment Analysis")
    print("=" * 60)

    pilot = load_pilot()
    dan_combined = load_dan_markers()
    detected_genes = set(pilot['GN'].dropna().unique())

    if args.only != 'E':
        print("\nPanel D — DAN Heatmap:")
        panel_d_heatmap(pilot, dan_combined)

    if args.only != 'D':
        print("\nPanel E — CMH Forest Plot:")
        marker_lists = build_marker_lists_v2(detected_genes)
        panel_e_forest(pilot, marker_lists)

    print("\nDone.")


if __name__ == '__main__':
    main()
