"""
Supplemental Figure 3, panel C: replicate scatter of the Lrrk2 WT proximity proteome.

Input
    data/ibioid/replicates/Lrrk2_WT_Protein_Normalized.csv
Output
    figures/suppl_03/output/panel_C_replicate_scatter_Lrrk2_WT.{png,svg}
Run
    python3 figures/suppl_03/scripts/replicate_scatter.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, scipy

Calls src/perbait_replicate_scatter.py.
"""
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))

from perbait_replicate_scatter import main

if __name__ == '__main__':
    sys.exit(0 if main('Lrrk2') is not None else 1)
