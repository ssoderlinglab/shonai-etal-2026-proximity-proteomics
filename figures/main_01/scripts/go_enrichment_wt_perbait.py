"""
Figure 1, panels F, J and N: GO enrichment of the three WT proximity proteomes (Cellular Component and Biological Process).

Input
    data/ibioid/go_enrichment/{aSyn,Vps35,Lrrk2}_WT.{Component,Process}.tsv
Output
    figures/main_01/output/{aSyn,Vps35,Lrrk2}_{BP,CC}_bg500_top10.{png,svg}
Run
    python3 figures/main_01/scripts/go_enrichment_wt_perbait.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, src/go_enrichment_tool

Term labels wrap at 30 characters, except the Lrrk2 Biological Process panel (panel N), which wraps at 33.
"""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths (relative to Figure_1/) ────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# submitted Figure_1/scripts/go_enrichment_perbait_wt.py; all plotting
# configuration (filters, palette, layout, fonts) is unchanged.
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, 'data', 'ibioid', 'go_enrichment')
OUTPUT_DIR = os.path.join(FIGURE_DIR, 'output')
# Shared GO lollipop tool still lives in the GitHub repository (unchanged).
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


def generate(tsv_path, title, prefix, outdir, max_label_width=30):
    """Load enrichment TSV, filter, plot, and save."""
    df = gel.load_enrichment_data(tsv_path)
    if "background_count" in df.columns:
        df = df[df["background_count"] <= BG_MAX].copy()
    if "category" in df.columns:
        df = df.drop(columns=["category"])

    fig = gel.plot_enrichment(df, category=None, title=title,
                              output_prefix=prefix, output_dir=outdir,
                              max_label_width=max_label_width,
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


def format_ibioid_go_title(bait_display, go_aspect_key):
    return f"WT {bait_display} iBioID - GO_{go_aspect_key}"


def main():
    outdir = OUTPUT_DIR

    bait_labels = {
        "aSyn": "alpha-Synuclein",
        "Vps35": "Vps35",
        "Lrrk2": "Lrrk2",
    }
    baits = {
        "aSyn": {
            "BP": os.path.join(DATA_DIR, "aSyn_WT.Process.tsv"),
            "CC": os.path.join(DATA_DIR, "aSyn_WT.Component.tsv"),
        },
        "Vps35": {
            "BP": os.path.join(DATA_DIR, "Vps35_WT.Process.tsv"),
            "CC": os.path.join(DATA_DIR, "Vps35_WT.Component.tsv"),
        },
        "Lrrk2": {
            "BP": os.path.join(DATA_DIR, "Lrrk2_WT.Process.tsv"),
            "CC": os.path.join(DATA_DIR, "Lrrk2_WT.Component.tsv"),
        },
    }

    print(f"Output: {outdir}/\n")
    for bait, info in baits.items():
        for cat_key in ("BP", "CC"):
            print(f"{bait} {cat_key}...")
            generate(
                info[cat_key],
                format_ibioid_go_title(bait_labels[bait], cat_key),
                f"{bait}_{cat_key}_bg500_top10",
                outdir,
                max_label_width=33 if (bait == "Lrrk2" and cat_key == "BP") else 30,
            )

    print(f"\nAll per-bait WT GO enrichment plots generated -> {outdir}/")


if __name__ == "__main__":
    main()
