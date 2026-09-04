"""
Figure 6, panels E, F and G: GO enrichment of module M11.

Input
    data/ibiocofrac/m11/go_enrichment/enrichment.{Component,Function,Process}.tsv
Output
    figures/main_06/output/*.{png,svg}
Run
    python3 figures/main_06/scripts/m11_go_enrichment.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, matplotlib, src/go_enrichment_tool
"""

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths ────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
REPO_ROOT = FIGURE_DIR.parent.parent
DATA_DIR = REPO_ROOT / 'data' / 'ibiocofrac' / 'm11' / 'go_enrichment'
OUTPUT_DIR = FIGURE_DIR / 'output'
SHARED_DIR = REPO_ROOT / 'src' / 'go_enrichment_tool'

# GO enrichment TSVs
TSV_BY_ASPECT = {
    'CC': DATA_DIR / 'enrichment.Component.tsv',
    'MF': DATA_DIR / 'enrichment.Function.tsv',
    'BP': DATA_DIR / 'enrichment.Process.tsv',
}

# ── Import shared library ───────────────────────────────────────
sys.path.insert(0, str(SHARED_DIR))
import config  # noqa: E402

config.PLOT_DEFAULTS['margin_left'] = 0.50
config.PLOT_DEFAULTS['font_family'] = [
    'DejaVu Sans', 'Arial', 'Helvetica', 'sans-serif'
]

import go_enrichment_lollipop as gel  # noqa: E402

# ── Constants ────────────────────────────────────────────────────
BG_MAX = 500
FDR_MAX = 0.05
SIMILARITY_THRESHOLD = 0.5
FIG_WIDTH = 10
TERM_LABEL_FONT_SIZE = 12

GO_TITLES = {
    'CC': 'Module 11 - Cellular Component (GO)',
    'MF': 'Module 11 - Molecular Function (GO)',
    'BP': 'Module 11 - Biological Process (GO)',
}


def fix_layout(fig):
    """Adjust subplot positions for consistent layout."""
    if fig is None:
        return
    fig.subplots_adjust(left=0.50, right=0.75, top=0.92, bottom=0.1)
    main_ax = None
    cbar_ax = None
    for ax in fig.get_axes():
        pos = ax.get_position()
        if pos.width < 0.03 and pos.x0 > 0.7:
            cbar_ax = ax
        elif pos.width > 0.1:
            main_ax = ax
    if cbar_ax is not None:
        cbar_x = 0.80
        cbar_pos = cbar_ax.get_position()
        cbar_ax.set_position([cbar_x, cbar_pos.y0,
                              cbar_pos.width, cbar_pos.height])
        if main_ax is not None:
            ax_pos = main_ax.get_position()
            cbar_center_fig_x = cbar_x + cbar_pos.width / 2
            bbox_x = ((cbar_center_fig_x - ax_pos.x0) / ax_pos.width)
            for child in main_ax.get_children():
                if hasattr(child, 'get_title') and hasattr(child, '_loc'):
                    try:
                        title_text = child.get_title().get_text()
                        if 'Gene' in title_text or 'count' in title_text:
                            child.set_bbox_to_anchor((bbox_x, 0.0))
                            child._loc = 8
                            break
                    except Exception:
                        pass


def style_labels(fig):
    """Bold y-tick labels and axis legends."""
    if fig is None:
        return
    for ax in fig.get_axes():
        pos = ax.get_position()
        if pos.width > 0.1:
            for lbl in ax.get_yticklabels():
                lbl.set_fontweight('bold')
            ax.xaxis.label.set_fontweight('bold')
            for child in ax.get_children():
                if hasattr(child, 'get_title') and hasattr(child, '_loc'):
                    try:
                        leg_title = child.get_title().get_text()
                        if 'Gene' in leg_title and 'count' in leg_title:
                            child.get_title().set_fontweight('bold')
                            break
                    except Exception:
                        pass
        elif pos.width < 0.03 and pos.x0 > 0.7:
            if getattr(ax, 'title', None) is not None:
                ax.title.set_fontweight('bold')


def load_filtered(tsv_path):
    """Load and filter: background_count <= BG_MAX, drop category column."""
    if not tsv_path.is_file():
        raise FileNotFoundError(tsv_path)
    df = gel.load_enrichment_data(str(tsv_path))
    if 'background_count' in df.columns:
        df = df[df['background_count'] <= BG_MAX].copy()
    if 'category' in df.columns:
        df = df.drop(columns=['category'])
    return df


def generate(tsv_path, aspect, max_terms, outdir):
    """Generate a single GO enrichment lollipop plot."""
    df = load_filtered(tsv_path)
    title = GO_TITLES[aspect]
    prefix = f'M11_{aspect}_top{max_terms}_bg{BG_MAX}'

    fig = gel.plot_enrichment(
        df, category=None, title=title,
        output_prefix=prefix, output_dir=str(outdir),
        theme='light', transparent_bg=False,
        color_palette='mint - blue', fig_width=FIG_WIDTH, dpi=300,
        output_formats=[], sort_by='signal',
        fdr_threshold=FDR_MAX, max_terms=max_terms,
        similarity_threshold=SIMILARITY_THRESHOLD,
        font_sizes={'term_label': TERM_LABEL_FONT_SIZE, 'title': 18, 'axis_label': 20},
    )
    fix_layout(fig)
    style_labels(fig)

    if fig is not None:
        outdir.mkdir(parents=True, exist_ok=True)
        for fmt in ('svg', 'png'):
            fig.savefig(outdir / f'{prefix}.{fmt}',
                        format=fmt, dpi=300, bbox_inches='tight')
        print(f'  OK: {prefix}')
    plt.close('all')


def main():
    parser = argparse.ArgumentParser(
        description='M11 GO enrichment lollipop (CC/MF/BP)')
    parser.add_argument('--max-terms', type=int, default=5,
                        help='Number of top terms per ontology (default: 5)')
    args = parser.parse_args()

    print(f'Conditions: BG<={BG_MAX}, FDR<={FDR_MAX}, '
          f'sort=signal, Jaccard={SIMILARITY_THRESHOLD}')
    print(f'Top terms: {args.max_terms}\n')

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for aspect, tsv in TSV_BY_ASPECT.items():
        print(f'{aspect} ({tsv.name})...')
        generate(tsv, aspect, args.max_terms, OUTPUT_DIR)

    print('\nDone.')


if __name__ == '__main__':
    main()
