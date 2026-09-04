"""
Cross-model functional-group statistics: size-matched permutation.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
Output
    figures/main_05/output/permutation_all21.csv
Run
    python3 src/permutation_core.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Imported by figures/main_05/scripts/module_group_fc_heatmap.py and figures/suppl_06/scripts/module_fc_by_group.py. Run it through figures/main_05/scripts/module_group_permutation.py.

200,000 draws without replacement, seed 20260803, Benjamini-Hochberg across the 21 group x model cells.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SRC_DIR)
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibiocofrac')
OUT_TABLES = os.path.join(REPO_ROOT, 'figures', 'main_05', 'output')

sys.path.insert(0, SRC_DIR)
from contamination_filter import CONTAMINATION_GENES as EXCLUDE

SEED = 20260803
N_PERM = 200_000

MODELS = ['LRRK2', 'SNCA', 'VPS35']
FRACTIONS = ['4', '5', '6', '7', '8', '9', '10']

# Functional groups exactly as the shipped panels define them.  Module numbers
# are LRRK2-partition numbers; module numbering does NOT correspond across the
# three partitions , so the grouping is defined once, on LRRK2.
GROUPS = [
    ('Synapse', 'Synapse (M7,8,11,19)', [7, 8, 11, 19]),
    ('Ribo',    'Ribosome (M1,10)',     [1, 10]),
    ('Mito',    'Mitochondria (M2,5)',  [2, 5]),
    ('M4',      'M4 (ER/Golgi)',        [4]),
    ('M6',      'M6 (Clath/tRNA)',      [6]),
    ('M9',      'M9 (Proteasome)',      [9]),
    ('M14',     'M14 (RNA-bind/ND)',    [14]),
]

PATH_PARTITION = os.path.join(
    DATA_DIR, 'modules', 'Partition_Gene_combined_3PD_models.csv')
PATH_FRAC = {
    'LRRK2': os.path.join(DATA_DIR, 'fractionation', 'LRRK2_young_for_PCA.csv'),
    'SNCA':  os.path.join(DATA_DIR, 'fractionation', 'SNCA_young_for_PCA.csv'),
    'VPS35': os.path.join(DATA_DIR, 'fractionation', 'Vps35_young_for_PCA.csv'),
}


# ---------------------------------------------------------------- data loading

def compute_protein_fc(frac_df):
    """Byte-for-byte the shipped `compute_protein_fc`. Collapses replicate AND
    fraction. Reproduced verbatim so the comparison below is like-for-like --
    do NOT 'improve' this function."""
    df = frac_df[frac_df['fraction'].astype(str).isin(FRACTIONS)].copy()
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
    frac = {e: pd.read_csv(p) for e, p in PATH_FRAC.items()}

    pfc = {}
    for e in MODELS:
        pfc[e] = compute_protein_fc(frac[e])
        a2g = dict(zip(partition[partition['Experiment'] == e]['UniProt'],
                       partition[partition['Experiment'] == e]['Symbol']))
        pfc[e]['Gene'] = pfc[e]['Accession'].map(a2g)

    merged = None
    for e in MODELS:
        d = (pfc[e][['Gene', 'log2FC']].dropna(subset=['Gene'])
             .rename(columns={'log2FC': f'{e}_FC'})
             .drop_duplicates(subset='Gene', keep='first'))
        merged = d if merged is None else merged.merge(d, on='Gene', how='outer')
    merged = merged[~merged['Gene'].isin(EXCLUDE)]

    lp = partition[partition['Experiment'] == 'LRRK2'][
        ['UniProt', 'Symbol', 'Membership']].copy()
    lp.columns = ['Accession', 'Gene', 'Module']
    lp = lp[~lp['Gene'].isin(EXCLUDE)]

    gsets = {k: set(lp[lp['Module'].isin(mods)]['Gene'])
             for k, _, mods in GROUPS}
    return merged, gsets


# ---------------------------------------------------------------- verification

def permutation_all_cells(merged, gsets, n_perm=N_PERM, seed=SEED):
    """Size-matched permutation, WITHOUT replacement, two-sided against the
    model's own proteome mean."""
    rng = np.random.default_rng(seed)
    rows = []
    for e in MODELS:
        pool = merged[['Gene', f'{e}_FC']].dropna()
        pv = pool[f'{e}_FC'].values
        pool_mean = pv.mean()
        for gk, glab, _ in GROUPS:
            v = pool[pool['Gene'].isin(gsets[gk])][f'{e}_FC'].values
            n = len(v)
            obs = v.mean()
            null = np.empty(n_perm)
            for i in range(n_perm):
                null[i] = pv[rng.choice(len(pv), n, replace=False)].mean()
            null_sd = null.std(ddof=1)
            exceed = int((np.abs(null - pool_mean) >= abs(obs - pool_mean)).sum())
            rows.append({
                'group': gk, 'group_label': glab, 'model': e, 'n': n,
                'mean_raw': obs, 'pool_mean': pool_mean,
                'delta': obs - pool_mean,
                'null_p2_5': float(np.percentile(null, 2.5)),
                'null_median': float(np.median(null)),
                'null_p97_5': float(np.percentile(null, 97.5)),
                'null_sd': null_sd,
                'z': (obs - pool_mean) / null_sd,
                'percentile_of_obs': 100 * float((null <= obs).mean()),
                'n_exceedances': exceed,
                'p_perm': exceed / n_perm if exceed > 0 else 1.0 / n_perm,
                'p_at_resolution_bound': exceed == 0,
            })
    P = pd.DataFrame(rows)
    P['q_perm_BH'] = multipletests(P['p_perm'].values, method='fdr_bh')[1]
    return P


def main():
    os.makedirs(OUT_TABLES, exist_ok=True)
    print("loading source data")
    merged, gsets = load_all()

    print(f"running size-matched permutation, {N_PERM:,} draws x 21 cells")
    P = permutation_all_cells(merged, gsets)
    P.to_csv(os.path.join(
        OUT_TABLES, 'permutation_all21.csv'), index=False)

    # Guard the headline claim: Synapse must be negative in all three models.
    syn = P[P['group'] == 'Synapse']
    assert (syn['delta'] < 0).all(), \
        "Synapse delta is not negative in all three models -- re-read before reporting"

    print(f"  q<0.05 permutation: {int((P['q_perm_BH'] < 0.05).sum())} of 21")
    print("done")


if __name__ == '__main__':
    main()
