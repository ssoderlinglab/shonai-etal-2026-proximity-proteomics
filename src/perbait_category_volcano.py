"""
Supplemental Figures 1, 2 and 3, panel E: four stacked volcano plots per bait, one functional category coloured in each.

Input
    data/ibioid/contrasts/{Vps35,aSyn,Lrrk2}_WT_vs_Control.tsv
    data/ibioid/go_string/<network>_enrichment_all_categories.tsv
    data/ibioid/category_gene_lists/<network>_<category>.csv
    data/external/syngo_1.3_Lrrk2_WT_ontologies_with_annotations.xlsx
Output
    figures/suppl_0{1,2,3}/output/panel_E_volcano_<bait>.{png,svg}
    figures/suppl_0{1,2,3}/output/panel_E_volcano_values.csv
Run
    python3 figures/suppl_01/scripts/volcano_category_highlight.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, openpyxl

Called with one bait. Significance gate log2FC > 0.585 (1.5-fold) and DEqMS Benjamini-Hochberg FDR < 0.05.

Each panel prints the enrichment term identifier, the Benjamini-Hochberg corrected enrichment Q value and the number of category members plotted. Q values come from STRING except Lrrk2 Synapse, which has no significant GO:0045202 row in the Lrrk2 STRING enrichment and takes its Q value from SynGO 1.3 cellular component, the resource the category was defined in; that panel prints SynGO.

Category members come from the frozen category gene lists except Lrrk2 Fiber organization, which has no list and uses the matching proteins of GO:1902903 in the Lrrk2 enrichment.

Panels are drawn at final print size, 55.0 x 41.2 mm each, four stacked with 2 mm gaps.
"""
import csv
import os
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# The SVG is the deliverable, so it must be byte-reproducible. Two things in a
# default matplotlib SVG are not: an embedded <dc:date> creation timestamp, and
# the hash-derived ids matplotlib mints for reusable <path> definitions. Both
# are suppressed here, and rcParams are frozen so an inherited user matplotlibrc
# cannot change the output either.
matplotlib.rcParams["svg.hashsalt"] = "perbait_category_volcano"
matplotlib.rcParams["svg.fonttype"] = "none"   # text stays editable in Illustrator
matplotlib.rcParams["svg.image_inline"] = True
matplotlib.rcParams["path.simplify"] = False
matplotlib.rcParams["pdf.compression"] = 0

try:                                            # matplotlib >= 3.7
    matplotlib.rcParams["svg.metadata"] = False
except KeyError:
    pass

SVG_METADATA = {"Date": None}                   # drop <dc:date> from the output

SRC_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SRC_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")
FIGURE_DIR = os.path.join(REPO_ROOT, "figures")

CONTRASTS = os.path.join(DATA_DIR, "ibioid", "contrasts")
GO_DIR = os.path.join(DATA_DIR, "ibioid", "go_string")
GLDIR = os.path.join(DATA_DIR, "ibioid", "category_gene_lists")
SYNGO = os.path.join(DATA_DIR, "external")
# each bait's figure directory in manuscript v25
FIGURE_OF = {"Vps35": "suppl_01", "aSyn": "suppl_02", "Lrrk2": "suppl_03"}


# --- layout, settled in step 7 ------------------------------------------------
N_PANELS = 4
GAP_MM = 2.0
PANEL_W_MM = 55.0
PANEL_H_MM = PANEL_W_MM * 0.75
MM = 1.0 / 25.4

# --- printed point sizes ------------------------------------------------------
TITLE_PT = 7.0
STAT_PT = 4.4
AXIS_PT = 5.2
TICK_PT = 4.4
LINEWIDTH = 0.4

LOG2FC_THRESHOLD = 0.585
FDR_THRESHOLD = 0.05
C_GREY_NONSIG = "#CCCCCC"
C_GREY_SIG = "#999999"
C_BAIT = "#000000"
CAT_S = 5.0
GREY_S = 1.6

# one colour per category, shared across baits
CATEGORY_COLOUR = {
    "Synapse": "#FFD166",
    "SV_Endocytosis": "#E63946",
    "Endosome": "#4CC9F0",
    "Vesicle_mediated_transport": "#9B5DE5",
    "RNA_splicing": "#4CC9F0",
    "Cytoskeleton": "#9B5DE5",
    "Fiber_regulation": "#2A9D8F",
}
COLOUR_NOTE = {
    "Synapse": "clustogram Synapse, unchanged",
    "SV_Endocytosis": "clustogram Presynapse hue, reassigned to break a duplicate",
    "Endosome": "clustogram Endosome, unchanged",
    "Vesicle_mediated_transport": "clustogram Cytoskeleton hue, reassigned",
    "RNA_splicing": "clustogram Endosome cyan, free in this row",
    "Cytoskeleton": "clustogram Cytoskeleton, unchanged",
    "Fiber_regulation": "clustogram SV/VMT teal, vacated and reused",
}
DISPLAY = {
    "Synapse": "Synapse",
    "SV_Endocytosis": "SV + Endocytosis",
    "Endosome": "Endosome",
    "Vesicle_mediated_transport": "Vesicle-mediated transport",
    "RNA_splicing": "RNA splicing",
    "Cytoskeleton": "Cytoskeleton",
    "Fiber_regulation": "Fiber organization",
}
TERM_ID = {
    "Synapse": "GO:0045202",
    "SV_Endocytosis": "mmu04144",
    "Endosome": "GO:0005768",
    "Vesicle_mediated_transport": "GO:0016192",
    "RNA_splicing": "GO:0008380",
    "Cytoskeleton": "GO:0005856",
    "Fiber_regulation": "GO:1902903",
}

BAITS = [
    dict(network="Vps35_WT", short="Vps35", bait_gene="Vps35",
         contrast="Vps35_WT_vs_Control.tsv", syngo_tag="Vps35_WT",
         categories=["Synapse", "SV_Endocytosis", "Endosome",
                     "Vesicle_mediated_transport"]),
    dict(network="aSyn_WT", short="aSyn", bait_gene="Snca",
         contrast="aSyn_WT_vs_Control.tsv", syngo_tag="aSyn_WT",
         categories=["Synapse", "SV_Endocytosis", "RNA_splicing", "Cytoskeleton"]),
    dict(network="Lrrk2_WT", short="Lrrk2", bait_gene="Lrrk2",
         contrast="Lrrk2_WT_vs_Control.tsv", syngo_tag="Lrrk2_WT",
         categories=["Synapse", "SV_Endocytosis", "Cytoskeleton",
                     "Fiber_regulation"]),
]

SUPERSCRIPT = {"-": "\u207b", "0": "\u2070", "1": "\u00b9", "2": "\u00b2",
               "3": "\u00b3", "4": "\u2074", "5": "\u2075", "6": "\u2076",
               "7": "\u2077", "8": "\u2078", "9": "\u2079"}


def hex_to_rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)], dtype=float)


def rgb_to_lab(rgb):
    rgb = rgb / 255.0
    rgb = np.where(rgb > 0.04045, ((rgb + 0.055) / 1.055) ** 2.4, rgb / 12.92)
    m = np.array([[0.4124, 0.3576, 0.1805],
                  [0.2126, 0.7152, 0.0722],
                  [0.0193, 0.1192, 0.9505]])
    xyz = (m @ rgb) / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, xyz ** (1 / 3), 7.787 * xyz + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def delta_e(a, b):
    return float(np.linalg.norm(rgb_to_lab(hex_to_rgb(a)) - rgb_to_lab(hex_to_rgb(b))))


def load_contrast(path):
    genes, log2fc, adjp = [], [], []
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row.get("Gene", "")
            if gene in ("", "NA"):
                continue
            l2, bh = row.get("log2FC", ""), row.get("deqms_BH", "")
            if l2 in ("", "NA") or bh in ("", "NA"):
                continue
            l2, bh = float(l2), float(bh)
            genes.append(gene)
            log2fc.append(l2)
            adjp.append(bh if bh > 0 else 1e-300)
    return genes, np.array(log2fc), np.array(adjp)


def string_row(network, term_id):
    path = os.path.join(GO_DIR, "%s_enrichment_all_categories.tsv" % network)
    with open(path, encoding="utf-8") as fh:
        hits = [r for r in csv.DictReader(fh, delimiter="\t") if r["term ID"] == term_id]
    if not hits:
        return None
    preferred = [r for r in hits if r["#category"].startswith("GO ")] or hits
    row = preferred[0]
    return row if float(row["false discovery rate"]) < FDR_THRESHOLD else None


def syngo_row(tag, term_id):
    import pandas as pd
    path = os.path.join(
        SYNGO, "syngo_1.3_{}_ontologies_with_annotations.xlsx".format(tag))
    assert os.path.exists(path), "no SynGO table for " + tag
    df = pd.read_excel(path)
    sub = df[df["GO term ID"] == term_id]
    assert not sub.empty, "{} absent from SynGO for {}".format(term_id, tag)
    return sub.iloc[0]


def statistic(network, syngo_tag, category):
    """The Q value for one panel, with its source. STRING first, SynGO as the
    fallback for a category SynGO itself defines."""
    term_id = TERM_ID[category]
    row = string_row(network, term_id)
    if row is not None:
        return (row["false discovery rate"], term_id, row["#category"],
                row["term description"], row["observed gene count"], "STRING")
    assert category in ("Synapse", "Presynapse", "Postsynapse"), \
        "no fallback resource for {} / {}".format(network, category)
    sr = syngo_row(syngo_tag, term_id)
    fdr = float(sr["GSEA FDR corrected p-value"])
    assert fdr < FDR_THRESHOLD, "SynGO {} not significant for {}: {}".format(
        term_id, network, fdr)
    return ("{:.2e}".format(fdr), term_id, "SynGO 1.3 CC",
            str(sr["GO term name"]), str(sr["GSEA count foreground/input"]), "SynGO")


def category_members(network, category, plotted):
    """Panel gene set. Frozen list when one exists, else the term's own members."""
    path = os.path.join(GLDIR, "{}_{}.csv".format(network, category))
    if os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            return {r["Gene"] for r in csv.DictReader(fh)}, "frozen gene list"
    term_id = TERM_ID[category]
    tsv = os.path.join(GO_DIR, "%s_enrichment_all_categories.tsv" % network)
    with open(tsv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["term ID"] == term_id:
                labels = r.get("matching proteins in your network (labels)", "")
                return ({g.strip() for g in labels.split(",") if g.strip()},
                        "GO term matching proteins")
    raise AssertionError("no gene set for {} / {}".format(network, category))


def format_q(fdr_string):
    value = float(fdr_string)
    assert 0 < value < 1, fdr_string
    mantissa, exponent = value, 0
    while mantissa < 1.0:
        mantissa *= 10.0
        exponent -= 1
    return "Q = {:.1f} \u00d7 10{}".format(
        mantissa, "".join(SUPERSCRIPT[c] for c in str(exponent)))


def points_under_text(ax, text_artist, log2fc, plot_y, x_lim, y_lim):
    """Count data points the statistic box actually covers.

    Measured from the RENDERED text extent, not from a guessed corner fraction.
    A fixed 42 x 30 percent corner box overstates the footprint: the real block
    occupies x 0.035..0.330 and y 0.777..0.975 in axes coordinates, and two
    strongly depleted proteins (aSyn Pcid2, Lrrk2 Erlec1) sit just below it.
    """
    ax.figure.canvas.draw()
    bb = text_artist.get_window_extent(
        renderer=ax.figure.canvas.get_renderer()).transformed(ax.transAxes.inverted())
    fx = (log2fc - x_lim[0]) / (x_lim[1] - x_lim[0])
    fy = (plot_y - y_lim[0]) / (y_lim[1] - y_lim[0])
    return int(((fx >= bb.x0) & (fx <= bb.x1) & (fy >= bb.y0) & (fy <= bb.y1)).sum())


def draw_panel(ax, genes, log2fc, adjp, cat_genes, bait_gene, colour, stat_lines,
               x_lim, y_lim, show_xlabel):
    neg_log10p = -np.log10(adjp)
    sig = (log2fc > LOG2FC_THRESHOLD) & (adjp < FDR_THRESHOLD)
    in_cat = np.array([g in cat_genes for g in genes])
    is_bait = np.array([g == bait_gene for g in genes])
    plot_y = np.clip(neg_log10p, 0, y_lim[1])

    other = ~in_cat & ~is_bait
    ax.scatter(log2fc[other & ~sig], plot_y[other & ~sig], c=C_GREY_NONSIG,
               s=GREY_S, marker="o", edgecolors="none", alpha=0.35,
               zorder=1, rasterized=True)
    ax.scatter(log2fc[other & sig], plot_y[other & sig], c=C_GREY_SIG,
               s=GREY_S, marker="o", edgecolors="none", alpha=0.5,
               zorder=5, rasterized=True)
    cat = in_cat & ~is_bait
    ax.scatter(log2fc[cat], plot_y[cat], c=colour, s=CAT_S, marker="o",
               edgecolors="none", alpha=0.9, zorder=30)
    ax.scatter(log2fc[is_bait], plot_y[is_bait], c=C_BAIT, s=22, marker="*",
               edgecolors="none", zorder=50)

    ax.axvline(LOG2FC_THRESHOLD, color="grey", linestyle="--", linewidth=0.25, alpha=0.5)
    ax.axhline(-np.log10(FDR_THRESHOLD), color="grey", linestyle="--",
               linewidth=0.25, alpha=0.5)

    if show_xlabel:
        ax.set_xlabel("log$_2$ (FC)", fontsize=AXIS_PT, labelpad=1.5)
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    ax.set_ylabel("\u2212log$_{10}$ (adj. P)", fontsize=AXIS_PT, labelpad=1.5)
    ax.tick_params(labelsize=TICK_PT, width=LINEWIDTH, length=1.8, pad=1.2)
    ax.set_xlim(*x_lim)
    ax.set_ylim(*y_lim)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(LINEWIDTH)

    stat_artist = ax.text(0.035, 0.975, "\n".join(stat_lines), transform=ax.transAxes,
                          ha="left", va="top", fontsize=STAT_PT, color="#000000",
                          linespacing=1.35, zorder=60,
                          bbox=dict(boxstyle="square,pad=0.30", facecolor="#FFFFFF",
                                    edgecolor="none", alpha=0.88))

    return int(in_cat.sum()), int(sig.sum()), points_under_text(
        ax, stat_artist, log2fc, plot_y, x_lim, y_lim)


def build_column(bait):
    genes, log2fc, adjp = load_contrast(os.path.join(CONTRASTS, bait["contrast"]))
    neg = -np.log10(adjp)
    x_lim = (log2fc.min() - 0.4, log2fc.max() + 0.4)
    y_lim = (-0.2, min(neg.max() + 0.4, 12))

    palette = [CATEGORY_COLOUR[c] for c in bait["categories"]]
    assert len(set(palette)) == 4, "{}: duplicate colour {}".format(
        bait["network"], palette)
    min_de = min(delta_e(a, b) for a, b in combinations(palette, 2))
    assert min_de > 60.0, "{}: minimum dE {:.1f} too low".format(
        bait["network"], min_de)

    fig_w = PANEL_W_MM * MM
    fig_h = (N_PANELS * PANEL_H_MM + (N_PANELS - 1) * GAP_MM) * MM
    fig, axes = plt.subplots(N_PANELS, 1, figsize=(fig_w, fig_h), dpi=600)
    fig.subplots_adjust(left=0.155, right=0.985, top=0.965, bottom=0.055,
                        hspace=GAP_MM / PANEL_H_MM * 5.0)

    records = []
    for k, (ax, category) in enumerate(zip(axes, bait["categories"])):
        fdr, term_id, database, description, observed, source = statistic(
            bait["network"], bait["syngo_tag"], category)
        cat_genes, gene_set_source = category_members(bait["network"], category, genes)
        n_plotted = sum(1 for g in genes if g in cat_genes)
        stat_lines = [format_q(fdr), term_id if source == "STRING"
                      else "{} (SynGO)".format(term_id), "n = {}".format(n_plotted)]
        n_cat, n_sig, n_corner = draw_panel(
            ax, genes, log2fc, adjp, cat_genes, bait["bait_gene"],
            CATEGORY_COLOUR[category], stat_lines, x_lim, y_lim,
            show_xlabel=(k == N_PANELS - 1))
        ax.set_title(DISPLAY[category], fontsize=TITLE_PT, fontweight="bold", pad=3.0)
        records.append(dict(
            network=bait["network"], panel_position=k + 1, category=category,
            display_name=DISPLAY[category], term_id=term_id, database=database,
            term_description=description, statistic_source=source, fdr=fdr,
            observed_gene_count=observed, panel_colour=CATEGORY_COLOUR[category],
            colour_note=COLOUR_NOTE[category], gene_set_source=gene_set_source,
            category_genes_total=len(cat_genes), category_genes_plotted=n_plotted,
            significant_points=n_sig, total_plotted_points=len(genes),
            points_under_stat_box=n_corner, printed_stat=" / ".join(stat_lines)))
    return fig, records, min_de


def main(short):
    """Build the four-panel category volcano column of one bait."""
    bait = [b for b in BAITS if b["short"] == short]
    assert bait, "unknown bait " + short
    bait = bait[0]
    dest = os.path.join(FIGURE_DIR, FIGURE_OF[short], "output")
    os.makedirs(dest, exist_ok=True)

    fig, all_records, min_de = build_column(bait)
    stem = os.path.join(dest, "panel_E_volcano_{}".format(short))
    fig.savefig(stem + ".svg", facecolor="white", metadata=SVG_METADATA)
    fig.savefig(stem + ".png", dpi=600, facecolor="white")
    plt.close(fig)

    ledger = os.path.join(dest, "panel_E_volcano_values.csv")
    cols = list(all_records[0].keys())
    with open(ledger, "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols)
        wr.writeheader()
        wr.writerows(all_records)

    print("panel {:.2f} x {:.2f} mm, four stacked = {:.1f} mm column".format(
        PANEL_W_MM, PANEL_H_MM, N_PANELS * PANEL_H_MM + (N_PANELS - 1) * GAP_MM))
    print("{:6} minimum pairwise dE {:.1f}".format(short, min_de))
    for r in all_records:
        print("  {:10} {:2} {:27} {} {:>9} {:>6}  n={:>4} corner={}".format(
            r["network"], r["panel_position"], r["display_name"], r["panel_colour"],
            r["fdr"], r["statistic_source"], r["category_genes_plotted"],
            r["points_under_stat_box"]))

    assert len(all_records) == 4, "expected 4 panels, got {}".format(len(all_records))
    for r in all_records:
        assert float(r["fdr"]) < FDR_THRESHOLD, (r["network"], r["display_name"],
                                                 r["fdr"])
        assert r["points_under_stat_box"] == 0, \
            "{} {}: statistic box covers {} points".format(
                r["network"], r["display_name"], r["points_under_stat_box"])
        assert r["category_genes_plotted"] > 0, (r["network"], r["display_name"])
    return all_records
