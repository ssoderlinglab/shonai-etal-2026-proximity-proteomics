"""
Supplemental Figure 7: the rcParams panels B, C and D were rendered under.

Input
    none; imported by figures/suppl_07/scripts/panels_BCD.py
Output
    none
Run
    python3 figures/suppl_07/scripts/panels_BCD.py
Dependencies
    python >= 3.9, matplotlib

The 22 values are pinned here so the render does not depend on interpreter state. Panel scripts call apply() at module load.
"""
import matplotlib as mpl

RC = {
    "font.size": 8.0,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.0,
    "axes.titlelocation": "left",
    "axes.linewidth": 0.6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "xtick.labelsize": 6.0,
    "ytick.labelsize": 6.0,
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.fontsize": 7.0,
    "legend.frameon": False,
    "lines.linewidth": 1.2,
    "patch.linewidth": 0.6,
    "figure.dpi": 200.0,
    "savefig.dpi": 300.0,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,   # editable text in Illustrator, not outlines
    "ps.fonttype": 42,
}


def apply():
    mpl.rcParams.update(RC)
