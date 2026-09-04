"""
Figure 2, panels D, E and F: GO enrichment of the shared WT proximity proteome.

Input
    data/ibioid/go_enrichment/PartConv_WT.{Component,Process,Function}.tsv
Output
    figures/main_02/output/*.{png,svg}
Run
    python3 figures/main_02/scripts/go_enrichment_wt_shared.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, src/go_enrichment_tool
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths ───────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid', 'go_enrichment')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')
GO_TOOL_DIR = os.path.join(REPO_ROOT, 'src', 'go_enrichment_tool')

# Add GO enrichment tool to path
sys.path.insert(0, GO_TOOL_DIR)

import config
# Override default margin for wider term labels
config.PLOT_DEFAULTS["margin_left"] = 0.50
# Use DejaVu Sans to avoid glyph issues with special characters
config.PLOT_DEFAULTS["font_family"] = [
    "DejaVu Sans",
    "Arial",
    "Helvetica",
    "sans-serif",
]

import go_enrichment_lollipop as gel

# ── Filter parameters (project-wide standard) ────────────────────────
BG_MAX = 500
FDR_MAX = 0.05
MAX_TERMS = 10
FIG_WIDTH = 10
TERM_LABEL_FONT_SIZE = 12
TERM_LABEL_BOLD = True
AXIS_LEGEND_TITLES_BOLD = True
TITLE_FONT_SIZE = 18

COMMON_KWARGS = dict(
    theme="light",
    transparent_bg=False,
    color_palette="mint - blue",
    fig_width=FIG_WIDTH,
    dpi=300,
    output_formats=[],  # Save manually after layout fix
    sort_by="signal",
    fdr_threshold=FDR_MAX,
    max_terms=MAX_TERMS,
    font_sizes={
        "term_label": TERM_LABEL_FONT_SIZE,
        "title": TITLE_FONT_SIZE,
        "axis_label": 20,  # enlarged x-axis title (Signal label)
    },
)


def fix_layout(fig):
    """Apply 5:2.5:2.5 aspect and align legends."""
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
        cbar_ax.set_position([cbar_x, cbar_pos.y0, cbar_pos.width, cbar_pos.height])

        if main_ax is not None:
            ax_pos = main_ax.get_position()
            cbar_center_fig_x = cbar_x + cbar_pos.width / 2
            bbox_x = (cbar_center_fig_x - ax_pos.x0) / ax_pos.width

            for child in main_ax.get_children():
                if hasattr(child, 'get_title') and hasattr(child, '_loc'):
                    try:
                        title_text = child.get_title().get_text()
                        if "Gene" in title_text or "count" in title_text:
                            child.set_bbox_to_anchor((bbox_x, 0.0))
                            child._loc = 8
                            break
                    except Exception:
                        pass


def _style_go_term_y_labels(fig):
    """Bold GO term Y-axis labels."""
    if fig is None or not TERM_LABEL_BOLD:
        return
    for ax in fig.get_axes():
        if ax.get_position().width > 0.1:
            for lbl in ax.get_yticklabels():
                lbl.set_fontweight("bold")
            break


def _style_fdr_gene_signal_bold(fig):
    """Bold FDR colorbar title, gene count legend title, and X-axis label."""
    if fig is None or not AXIS_LEGEND_TITLES_BOLD:
        return
    main_ax = None
    cbar_ax = None
    for ax in fig.get_axes():
        pos = ax.get_position()
        if pos.width < 0.03 and pos.x0 > 0.7:
            cbar_ax = ax
        elif pos.width > 0.1:
            main_ax = ax
    if main_ax is not None:
        main_ax.xaxis.label.set_fontweight("bold")
        for child in main_ax.get_children():
            if hasattr(child, "get_title") and hasattr(child, "_loc"):
                try:
                    leg_title = child.get_title().get_text()
                    if "Gene" in leg_title and "count" in leg_title:
                        child.get_title().set_fontweight("bold")
                        break
                except Exception:
                    pass
    if cbar_ax is not None and getattr(cbar_ax, "title", None) is not None:
        cbar_ax.title.set_fontweight("bold")


def generate(tsv_path, title, prefix, outdir):
    """Load enrichment TSV, filter, plot, and save."""
    df = gel.load_enrichment_data(tsv_path)
    if "background_count" in df.columns:
        df = df[df["background_count"] <= BG_MAX].copy()
    if "category" in df.columns:
        df = df.drop(columns=["category"])

    fig = gel.plot_enrichment(df, category=None, title=title,
                              output_prefix=prefix, output_dir=outdir,
                              **COMMON_KWARGS)
    fix_layout(fig)
    _style_go_term_y_labels(fig)
    _style_fdr_gene_signal_bold(fig)

    if fig is not None:
        os.makedirs(outdir, exist_ok=True)
        for fmt in ["svg", "png"]:
            outpath = os.path.join(outdir, f"{prefix}.{fmt}")
            fig.savefig(outpath, format=fmt, dpi=300, bbox_inches="tight")
        print(f"    OK: {prefix}")
    plt.close("all")


def format_convergent_title(go_aspect_key):
    """
    Title for convergent GO plots.

    v02 (2026-08-05, Daichi's instruction): `(>=2/3 baits)` REMOVED as redundant — the
    Figure 2 legend already states the >=2 of 3 basis and the 242-protein count, so the
    panel title repeated it.

    ⚠ The removal is COSMETIC ONLY. The plotted data is unchanged: still the partial
    convergent (>=2 of 3 baits, 242-protein) set from `PartConv_WT_enrichment.*.tsv`.
    `PartConv` remains in every output filename, which is where the basis is recorded.
    Do NOT read the shorter title as a change of basis to the 86-protein total-convergent set.
    """
    return f"WT convergent iBioID - GO_{go_aspect_key}"


def main():
    parser = argparse.ArgumentParser(
        description="Convergent WT GO enrichment lollipop plots (Figure 2)"
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Generate CC only (for testing)",
    )
    args = parser.parse_args()

    outdir = OUTPUT_DIR

    if args.quick:
        tsv = os.path.join(DATA_DIR, "PartConv_WT.Component.tsv")
        generate(tsv, format_convergent_title("CC"),
                 "PartConv_CC_bg500_top10", outdir)
        print(f"\nQuick run done -> {outdir}/")
        return

    categories = {
        "CC": os.path.join(DATA_DIR, "PartConv_WT.Component.tsv"),
        "BP": os.path.join(DATA_DIR, "PartConv_WT.Process.tsv"),
        "MF": os.path.join(DATA_DIR, "PartConv_WT.Function.tsv"),
    }

    print(f"Output: {outdir}/\n")
    for cat_key, tsv_path in categories.items():
        print(f"PartConv {cat_key}...")
        generate(
            tsv_path,
            format_convergent_title(cat_key),
            f"PartConv_{cat_key}_bg500_top10",
            outdir,
        )

    print(f"\nAll convergent WT GO enrichment plots generated -> {outdir}/")


if __name__ == "__main__":
    main()
