"""
Supplemental Figure 6, panel B: module coherence permutation test.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
Output
    figures/suppl_06/output/module_coherence_density.csv
Run
    python3 figures/suppl_06/scripts/module_coherence_permutation.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

First of three steps: permutation, then Benjamini-Hochberg, then the panel.

Mean pairwise correlation across sedimentation fractions, against size-matched random protein sets, 2,000 draws, seed 20260826.
"""
import csv
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA = os.path.join(REPO_ROOT, 'data', 'ibiocofrac')
PARTITION = os.path.join(DATA, 'modules',
                         'Partition_Gene_combined_3PD_models.csv')
FRAC = {
    'LRRK2': os.path.join(DATA, 'fractionation', 'LRRK2_young_for_PCA.csv'),
    'SNCA': os.path.join(DATA, 'fractionation', 'SNCA_young_for_PCA.csv'),
    'VPS35': os.path.join(DATA, 'fractionation', 'Vps35_young_for_PCA.csv'),
}
OUT = os.path.join(FIGURE_DIR, 'output')

REFERENCE = 'LRRK2'
COHORTS = ('LRRK2', 'SNCA', 'VPS35')
GENOTYPES = ('wt', 'tg')
N_PERM = 2000
SEED = 20260826
SOFT_POWER = 6
MIN_MODULE_N = 8   # the manuscript's own detection floor (Results 75)
M11_REFERENCE_SIZE = 77
FUNCTIONAL_12 = ['1', '2', '4', '5', '6', '7', '8', '9', '10', '11', '14', '19']
META = ('fraction', 'replicate', 'sample')


def load_reference_modules():
    rows = list(csv.DictReader(open(PARTITION)))
    ref = {r['UniProt']: r['Membership']
           for r in rows if r['Experiment'].upper() == REFERENCE}
    m11 = sum(1 for v in ref.values() if v == '11')
    assert m11 == M11_REFERENCE_SIZE, (
        'reference module 11 is %d, expected %d -- input drift'
        % (m11, M11_REFERENCE_SIZE))
    by_mod = defaultdict(list)
    for acc, mod in ref.items():
        by_mod[mod].append(acc)
    return by_mod


def load_matrix(path, genotype):
    """Return (accessions, X) with X shape (n_proteins, n_fraction_columns)."""
    with open(path) as fh:
        rdr = csv.reader(fh)
        hdr = next(rdr)
        acc_col = hdr[0]
        prot_cols = [c for c in hdr if c not in (acc_col,) + META]
        i_samp = hdr.index('sample')
        idx = {c: hdr.index(c) for c in prot_cols}
        rows = []
        for row in rdr:
            if not row[i_samp].startswith(genotype + '_'):
                continue
            rows.append([row[idx[c]] for c in prot_cols])
    X = np.array([[float(v) if v not in ('', 'NA', 'NaN') else np.nan
                   for v in r] for r in rows], dtype=float).T
    return prot_cols, X


def corr_stats(X, rows, power=SOFT_POWER):
    """meanCor and meanAdj over the upper triangle of a protein subset."""
    sub = X[rows, :]
    keep = ~np.all(np.isnan(sub), axis=1)
    sub = sub[keep, :]
    if sub.shape[0] < 3:
        return np.nan, np.nan
    with np.errstate(invalid='ignore', divide='ignore'):
        C = np.corrcoef(np.nan_to_num(sub, nan=0.0))
    n = C.shape[0]
    iu = np.triu_indices(n, k=1)
    r = C[iu]
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan, np.nan
    adj = ((1.0 + r) / 2.0) ** power
    return float(np.mean(r)), float(np.mean(adj))


def main():
    by_mod = load_reference_modules()
    rng = np.random.default_rng(SEED)
    out_rows = []

    for cohort in COHORTS:
        for genotype in GENOTYPES:
            accs, X = load_matrix(FRAC[cohort], genotype)
            pos = {a: i for i, a in enumerate(accs)}
            n_cols = X.shape[1]
            for mod, members in by_mod.items():
                rows = [pos[a] for a in members if a in pos]
                if len(rows) < MIN_MODULE_N:
                    continue
                obs_cor, obs_adj = corr_stats(X, rows)
                if not np.isfinite(obs_cor):
                    continue
                null_cor = np.empty(N_PERM)
                null_adj = np.empty(N_PERM)
                allidx = np.arange(len(accs))
                for i in range(N_PERM):
                    pick = rng.choice(allidx, size=len(rows), replace=False)
                    null_cor[i], null_adj[i] = corr_stats(X, pick.tolist())
                mc, sc = np.nanmean(null_cor), np.nanstd(null_cor)
                ma, sa = np.nanmean(null_adj), np.nanstd(null_adj)
                z_cor = (obs_cor - mc) / sc if sc > 0 else np.nan
                z_adj = (obs_adj - ma) / sa if sa > 0 else np.nan
                p_cor = (1 + int(np.sum(null_cor >= obs_cor))) / (N_PERM + 1)
                out_rows.append({
                    'cohort': cohort,
                    'is_reference_cohort': cohort == REFERENCE,
                    'genotype': genotype,
                    'module': 'M' + mod,
                    'in_functional_12': mod in FUNCTIONAL_12,
                    'n_reference_members': len(members),
                    'n_tested_in_cohort': len(rows),
                    'fraction_columns': n_cols,
                    'meanCor': round(obs_cor, 4),
                    'meanCor_null_mean': round(float(mc), 4),
                    'meanCor_null_sd': round(float(sc), 4),
                    'Z_meanCor': round(float(z_cor), 3),
                    'p_perm_meanCor': p_cor,
                    'meanAdj': round(obs_adj, 5),
                    'Z_meanAdj': round(float(z_adj), 3),
                })

    with open(os.path.join(OUT, 'module_coherence_density.csv'), 'w',
              newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    meta = {'n_permutations': N_PERM, 'seed': SEED, 'soft_power': SOFT_POWER,
            'min_module_n': MIN_MODULE_N,
            'reference_partition': REFERENCE,
            'm11_reference_size_asserted': M11_REFERENCE_SIZE,
            'correlation_axis': 'sedimentation fractions, NOT biological '
                                'samples',
            'connectivity_preservation': 'deliberately not computed, 21 '
                                         'columns per genotype',
            'n_rows': len(out_rows)}
    with open(os.path.join(OUT, 'U2_framingA_meta.json'), 'w') as fh:
        json.dump(meta, fh, indent=2, sort_keys=True)

    # console: the reference-cohort control, then the two test cohorts, M11 first
    print('rows written: %d' % len(out_rows))
    print()
    print('%-6s %-3s %-5s %5s %5s %8s %8s %7s' % (
        'cohort', 'gt', 'mod', 'nref', 'ntest', 'meanCor', 'nullmean', 'Z'))
    for r in out_rows:
        if r['module'] == 'M11':
            print('%-6s %-3s %-5s %5d %5d %8.4f %8.4f %7.2f' % (
                r['cohort'], r['genotype'], r['module'],
                r['n_reference_members'], r['n_tested_in_cohort'],
                r['meanCor'], r['meanCor_null_mean'], r['Z_meanCor']))
    print()
    zs = [r['Z_meanCor'] for r in out_rows if r['in_functional_12']
          and np.isfinite(r['Z_meanCor'])]
    print('functional-12 rows: %d, Z_meanCor median %.2f, min %.2f, max %.2f'
          % (len(zs), float(np.median(zs)), min(zs), max(zs)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
