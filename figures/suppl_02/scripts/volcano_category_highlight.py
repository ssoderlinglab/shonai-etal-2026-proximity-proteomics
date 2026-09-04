"""
Supplemental Figure 2, panel E: four stacked volcano plots of the aSyn WT proximity proteome, one functional category coloured in each.

Input
    data/ibioid/contrasts/aSyn_WT_vs_Control.tsv
    data/ibioid/go_string/aSyn_WT_enrichment_all_categories.tsv
    data/ibioid/category_gene_lists/aSyn_WT_*.csv
Output
    figures/suppl_02/output/panel_E_volcano_aSyn.{png,svg}
    figures/suppl_02/output/panel_E_volcano_values.csv
Run
    python3 figures/suppl_02/scripts/volcano_category_highlight.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, openpyxl

Calls src/perbait_category_volcano.py.
"""
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))

from perbait_category_volcano import main

if __name__ == '__main__':
    sys.exit(0 if main('aSyn') is not None else 1)
