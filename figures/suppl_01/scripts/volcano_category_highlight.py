"""
Supplemental Figure 1, panel E: four stacked volcano plots of the Vps35 WT proximity proteome, one functional category coloured in each.

Input
    data/ibioid/contrasts/Vps35_WT_vs_Control.tsv
    data/ibioid/go_string/Vps35_WT_enrichment_all_categories.tsv
    data/ibioid/category_gene_lists/Vps35_WT_*.csv
Output
    figures/suppl_01/output/panel_E_volcano_Vps35.{png,svg}
    figures/suppl_01/output/panel_E_volcano_values.csv
Run
    python3 figures/suppl_01/scripts/volcano_category_highlight.py
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
    sys.exit(0 if main('Vps35') is not None else 1)
