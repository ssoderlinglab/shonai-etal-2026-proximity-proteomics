"""
Supplemental Figure 6, panel B: Benjamini-Hochberg correction of the module coherence permutation p values.

Input
    figures/suppl_06/output/module_coherence_density.csv
Output
    figures/suppl_06/output/module_coherence_density_BH.csv
Run
    python3 figures/suppl_06/scripts/module_coherence_bh.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Second of three steps.

The correction family is the 48 test-cohort rows (2 test cohorts x 2 genotypes x 12 modules); the reference cohort's modules were defined in that cohort and are an internal control, not hypotheses.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
OUT = os.path.join(FIGURE_DIR, 'output')
SRC = os.path.join(OUT, 'module_coherence_density.csv')
DST = os.path.join(OUT, 'module_coherence_density_BH.csv')

N_PERM = 2000
P_BOUND = 1.0 / (N_PERM + 1)
PRIMARY_GENOTYPE = 'wt'


def bh(pvals):
    """Benjamini-Hochberg step-up. Returns q in the input order."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        val = min(prev, pvals[i] * n / k)
        q[i] = val
        prev = val
    return q


def main():
    rows = list(csv.DictReader(open(SRC)))
    for r in rows:
        r['p_perm_meanCor'] = float(r['p_perm_meanCor'])

    families = {
        'primary_24': [r for r in rows
                       if r['in_functional_12'] == 'True'
                       and r['cohort'] != 'LRRK2'
                       and r['genotype'] == PRIMARY_GENOTYPE],
        'sensitivity_48': [r for r in rows
                           if r['in_functional_12'] == 'True'
                           and r['cohort'] != 'LRRK2'],
    }

    for name, fam in families.items():
        qs = bh([r['p_perm_meanCor'] for r in fam])
        col = 'q_BH_' + name
        for r, q in zip(fam, qs):
            r[col] = round(q, 6)

    fields = list(rows[0].keys())
    for c in ('q_BH_primary_24', 'q_BH_sensitivity_48'):
        if c not in fields:
            fields.append(c)
    fields.append('p_at_resolution_bound')
    for r in rows:
        r['p_at_resolution_bound'] = (
            'TRUE' if abs(r['p_perm_meanCor'] - P_BOUND) < 1e-9 else 'FALSE')
        for c in ('q_BH_primary_24', 'q_BH_sensitivity_48'):
            r.setdefault(c, '')

    with open(DST, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    fam = families['primary_24']
    print('PRIMARY FAMILY: %d tests (%s only, test cohorts, functional 12)'
          % (len(fam), PRIMARY_GENOTYPE))
    print('%-6s %-5s %8s %10s %10s %6s' % ('cohort', 'mod', 'Z', 'p_perm',
                                           'q_BH', 'star'))
    for r in sorted(fam, key=lambda r: r['p_perm_meanCor']):
        q = r['q_BH_primary_24']
        star = ('***' if q < 0.001 else '**' if q < 0.01
                else '*' if q < 0.05 else 'ns')
        print('%-6s %-5s %8s %10.4f %10.4f %6s%s'
              % (r['cohort'], r['module'], r['Z_meanCor'],
                 r['p_perm_meanCor'], q, star,
                 '   (p at bound)' if r['p_at_resolution_bound'] == 'TRUE'
                 else ''))

    n_sig = sum(1 for r in fam if r['q_BH_primary_24'] < 0.05)
    print('\nsignificant at q < 0.05: %d of %d' % (n_sig, len(fam)))
    fam48 = families['sensitivity_48']
    n48 = sum(1 for r in fam48 if r['q_BH_sensitivity_48'] < 0.05)
    print('sensitivity family (48, both genotypes): %d of %d significant'
          % (n48, len(fam48)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
