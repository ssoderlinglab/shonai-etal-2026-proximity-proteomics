"""
Figure 5, panel E and Supplemental Figure 6, panel D: build the permutation table both read.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
Output
    figures/main_05/output/permutation_all21.csv
Run
    python3 figures/main_05/scripts/module_group_permutation.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels

200,000 draws without replacement, seed 20260803. Runtime is a few minutes.
"""
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(SCRIPT_DIR)))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))

from permutation_core import main

if __name__ == '__main__':
    sys.exit(main())
