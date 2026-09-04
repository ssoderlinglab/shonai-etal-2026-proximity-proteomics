"""
Shared publication style for all Section 4 figures.
Import this module to ensure consistent formatting across panels.

Usage:
    import sys; sys.path.insert(0, str(Path(__file__).parent.parent))
    from pub_style import apply_style, COLORS, FONT_SIZES
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

# === Color palette (colorblind-friendly, white background) ===
COLORS = {
    # Screen conditions
    'baseline_pre':  '#9ecae1',
    'baseline_post': '#2171b5',
    'asyn_pre':      '#fdae6b',
    'asyn_post':     '#e6550d',
    'baseline':      '#3182bd',
    'asyn':          '#e6550d',
    # Histology
    'contra':        '#4daf4a',   # green — uninjected / control side
    'ipsi':          '#e41a1c',   # red — injection side
    'nt_gRNA':       '#377eb8',   # blue — non-targeting control
    'stxbp1_gRNA':   '#e41a1c',   # red — Stxbp1 KO
    'gfp':           '#984ea3',   # purple — GFP control
    'tent':          '#ff7f00',   # orange — TeNT
    # WB / validation
    'control':       '#666666',
    'ko':            '#e41a1c',
    # General
    'ns':            '#999999',   # non-significant
    'sig':           '#e41a1c',   # significant
    'highlight':     '#e41a1c',
    'nt_highlight':  '#4daf4a',   # non-targeting in MAGeCK
}

FONT_SIZES = {
    'title':    11,
    'label':    10,
    'tick':      8,
    'legend':    8,
    'annotation': 7,
}

def apply_style():
    """Apply global publication style settings."""
    plt.rcParams.update({
        'font.family':        'sans-serif',
        'font.size':          9,
        'axes.linewidth':     0.8,
        'axes.labelsize':     FONT_SIZES['label'],
        'axes.titlesize':     FONT_SIZES['title'],
        'axes.spines.top':    False,
        'axes.spines.right':  False,
        'xtick.major.width':  0.6,
        'ytick.major.width':  0.6,
        'xtick.labelsize':    FONT_SIZES['tick'],
        'ytick.labelsize':    FONT_SIZES['tick'],
        'legend.fontsize':    FONT_SIZES['legend'],
        'legend.frameon':     False,
        'figure.dpi':         300,
        'savefig.dpi':        300,
        'savefig.bbox':       'tight',
        'savefig.pad_inches': 0.05,
        'pdf.fonttype':       42,
        'ps.fonttype':        42,
    })
