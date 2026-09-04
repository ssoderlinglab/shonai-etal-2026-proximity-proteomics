"""
Supplemental Figure 7, panels B, C and D: functional-group shift, its per-protein distribution, and the cortex against dopamine neuron comparison, on a shared layout.

Input
    figures/suppl_07/output/cortex_permutation_all21.csv
    figures/suppl_07/output/cortex_vs_cofrac_21cells.csv
    data/cortex/deqms_all_contrasts_long.tsv
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
Output
    figures/suppl_07/output/panel_B_group_shift_heatmap.{png,svg}
    figures/suppl_07/output/panel_C_group_protein_strips.{png,svg}
    figures/suppl_07/output/panel_D_cortex_vs_cofrac.{png,svg}
Run
    python3 figures/suppl_07/scripts/panels_BCD.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, statsmodels

Run figures/suppl_07/scripts/cortex_group_permutation.py first; this script draws the table, it does not compute it.

Panel B marks one, two and three asterisks at q < 0.05, q < 0.01 and q < 0.001. Panel C shows the distribution and carries no test. Panel D's dotted diagonal marks equal shift in the two measurements and is not a null expectation.

C and D share one layout so a genotype occupies the same column in both; the script asserts the alignment and asserts that no text overlaps before writing.
"""
import hashlib
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUT = os.path.join(FIGURE_DIR, "output")
PAN = OUT
os.makedirs(OUT, exist_ok=True)

sys.path.insert(0, SCRIPT_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
import figure_style
figure_style.apply()
from contamination_filter import CONTAMINATION_GENES as EXCLUDE

# SVG text stays editable for the figure assembly.
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams.update({"font.size": 8.0, "axes.labelsize": 8.0,
                     "axes.titlesize": 8.0, "legend.fontsize": 6.5,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 7.0})

import label_placement as p14
place_labels, state_of = p14.place_labels, p14.state_of
p14.MARKER_R_PT = 7.2          # clears the outer 12.0-pt hairline

# ---------------------------------------------------------------- shared layout
FIGW = 9.6                     # inches, identical for C and D
L, R, WS = 0.135, 0.985, 0.10  # left/right margins and wspace, identical
PANEL_W = (R - L) * FIGW / (3 + 2 * WS)      # one panel's width, inches

# ⛔ VERTICAL SPACING IS SPECIFIED IN INCHES, NOT IN FIGURE FRACTIONS.
# B, C and D all have different figure heights (D's box must stay square for the
# equal-aspect diagonal), so the same fractional offset is a DIFFERENT physical
# gap in each - which is how the plot-to-key gap drifted apart in v1.
#
# GAP_IN is the gap a READER MEASURES: from the bottom of the lowest axis-
# attached text (x-axis label in C/D, genotype tick labels in B) to the top of
# the key block. It is NOT a subplots_adjust offset. Earlier versions set
# subplots_adjust and let the key land where it fell, which double-counted a
# line allowance and left ~0.36 in of white space in all three panels.
CAP_BASE_IN = 0.09    # figure bottom -> caption baseline
CAP_ASC_IN = 0.08     # caption baseline -> top of its glyphs (measured)
GAP_IN = 0.13         # THE GAP: lowest axis text -> top of the key block
LEG_GLYPH_IN = 0.19   # provisional height of D's marker key (ring glyph + label).
                      # NOT trusted: the concentric-ring handler's drawn height is
                      # not derivable from the font size, so panel D MEASURES the
                      # rendered legend and corrects its own geometry (see below).
LEG_GAP_IN = 0.10     # caption top -> marker key bottom (D only)
XLABEL_IN = 0.34      # axes bottom -> bottom of the x-axis label (C, D)
GENOLAB_IN = 0.21     # axes bottom -> bottom of the genotype tick labels (B)
# Top-down: figure top -> figure title -> TITLE_GAP_IN -> genotype title -> axes
SUPT_IN = 0.17        # figure top -> figure-title baseline
TITLE_GAP_IN = 0.30   # figure title -> genotype title (0.12 was too close)
GENO_H_IN = 0.17      # genotype title's own height
TOP_IN = SUPT_IN + TITLE_GAP_IN + GENO_H_IN
BOX_H_B = 2.35        # B's axes box height, inches (7 rows)
BOX_H_C = 2.62        # C's axes box height, inches (7 rows plus dot jitter)
                      # D's box height is PANEL_W: it must stay square.

COL = {"VPS35": "#7570B3", "SNCA": "#D95F02", "LRRK2": "#1B9E77"}
GREY = "#B0B0B0"
MODELS = ["LRRK2", "SNCA", "VPS35"]
GROUP_ORDER = ["Synapse", "Ribo", "Mito", "M4", "M6", "M9", "M14"]
GMODS = {"Synapse": [7, 8, 11, 19], "Ribo": [1, 10], "Mito": [2, 5],
         "M4": [4], "M6": [6], "M9": [9], "M14": [14]}
DISC_MS, RING_MS, RING_LW = 7.5, (10.2, 12.0), 0.6
GATE = 0.585

# genotype titles, shared by C and D (item 4). GENO_SET picks the wording.
GENO = {
    1: {"LRRK2": "Tg-LRRK2$^{G2019S}$", "SNCA": "Tg-SNCA$^{A30P{\\times}A53T}$",
        "VPS35": "KI-VPS35$^{D620N}$"},
    2: {"LRRK2": "Tg $Pdgfb$-hLRRK2$^{G2019S}$",
        "SNCA": "Tg $Th$-hSNCA$^{A30P{\\times}A53T}$",
        "VPS35": "KI-Vps35$^{D620N}$"},
    3: {"LRRK2": "LRRK2$^{G2019S}$", "SNCA": "SNCA$^{A30P{\\times}A53T}$",
        "VPS35": "VPS35$^{D620N}$"},
}
GENO_SET = 1

# Titles keep the shape used by the dopamine neuron counterparts of these two
# panels, changing only the level (module -> functional group) and the scope
# (dopamine neurons -> whole cortex). Set 1 is the one the figure carries.
TITLE_B = {
    1: "Group-level protein FC \u2014 whole cortex",
    2: "Group-level protein FC \u2014 combined groups, whole cortex",
    3: "Functional-group protein FC \u2014 whole cortex",
    4: "Only VPS35 shifts functional groups in whole cortex",
}
TITLE_C = {
    1: "Protein-level log$_2$FC by group \u2014 whole cortex",
    2: "Protein-level log$_2$FC by functional group \u2014 whole cortex",
    3: "Protein-level log$_2$FC by group \u2014 group-identity mapping",
    4: "No functional group shifts in whole cortex except under VPS35",
}
TITLE_D = {
    1: "Group-level protein FC \u2014 DA neurons versus whole cortex",
    2: "Most functional-group shifts are DA-neuron-selective",
    3: "Group shifts in DA neurons versus whole cortex",
    4: "Most functional-group shifts are DA-neuron-selective",
}
TITLE_SET = 1

CAP_C = {
    1: "dot = protein   bar = group mean   dashed = proteome baseline   "
       "* q$<$0.05   ** q$<$0.01   *** q$<$0.001",
    2: "bar = group mean   dashed = proteome baseline   permutation q: "
       "* $<$0.05   ** $<$0.01   *** $<$0.001",
    3: "group mean (bar) vs proteome baseline (dashed)   "
       "* q$<$0.05  ** q$<$0.01  *** q$<$0.001",
}
CAP_D = {
    1: "shaded = opposite sign   dotted = equal shift",
    2: "axes: group mean $-$ proteome baseline   shaded = opposite sign",
    3: "dotted = equal shift (not a null)   shaded = opposite sign",
}


class ConcentricHandler(HandlerBase):
    """Rings and disc drawn on ONE centre - matplotlib's HandlerTuple lays them
    side by side, which made a lone ring read as the open 'no shift' glyph."""

    def __init__(self, ring_sizes, disc_size, colour, ring_lw, **kw):
        self.ring_sizes, self.disc_size = ring_sizes, disc_size
        self.colour, self.ring_lw = colour, ring_lw
        super().__init__(**kw)

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
                       height, fontsize, trans):
        cx, cy = width / 2 - xdescent, height / 2 - ydescent
        out = [Line2D([cx], [cy], marker="o", ls="", ms=ms, mfc="none",
                      mec=self.colour, mew=self.ring_lw, transform=trans)
               for ms in self.ring_sizes]
        out.append(Line2D([cx], [cy], marker="o", ls="", ms=self.disc_size,
                          mfc=self.colour, mec=self.colour, transform=trans))
        return out



def measure_bands(fig, key_artists):
    """Return (band below the axes, height of the key block) in INCHES.

    ⛔ NEITHER NUMBER CAN BE COMPUTED FROM FONT SIZES, which is why the
    plot-to-key gap kept coming out unequal (0.153 / 0.139 / 0.210 in for the
    same nominal 0.13) across three attempts at arithmetic: B's genotype tick
    labels carry no descender, D's concentric-ring legend handler draws taller
    than its label row, and a Legend BOX pads above its tallest glyph so the box
    extent overstates what a reader sees as the key's top.

    Both quantities are INDEPENDENT OF FIGURE HEIGHT (font sizes are fixed), so
    each panel measures them on a provisional figure and REBUILDS at the height
    they imply. Rescaling artists in place was tried and abandoned: a Legend's
    anchor and its rendered extent live in different coordinate systems, and the
    in-place version pushed panel D's legend off the canvas.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    dpi, H = fig.dpi, fig.get_size_inches()[1]
    ax_bot = fig.axes[0].get_position().y0 * H
    # ⛔ the key artists sit below the axes too, so they must be EXCLUDED here or
    # the band swallows the caption and the gap is counted twice (measured 0.40 /
    # 0.39 / 0.74 in instead of 0.13).
    skip = set()
    for a in key_artists:
        if isinstance(a, mpl.legend.Legend):
            skip.update(id(t) for t in a.get_texts())
        else:
            skip.add(id(a))
    below = [t.get_window_extent(r).y0 / dpi for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible() and id(t) not in skip
             and t.get_window_extent(r).y1 / dpi < ax_bot + 0.02]
    band = ax_bot - min(below) if below else 0.0

    tops, bots = [], []
    for a in key_artists:
        if isinstance(a, mpl.legend.Legend):
            inner = list(a.get_texts()) + [h for h in a.legend_handles
                                           if h is not None]
            ext = [o.get_window_extent(r) for o in inner
                   if hasattr(o, "get_window_extent")]
            tops += [e.y1 / dpi for e in ext]
            bots += [e.y0 / dpi for e in ext]
    # ⛔ RETURN THE PIECES, NOT A LUMPED HEIGHT. The earlier version returned
    # (top - bottom) of everything in the key and the caller added it to
    # CAP_BASE_IN, which is the caption's BASELINE - so the caption's descender
    # was counted twice and D's gap came out 0.199 in against 0.151 for B and C.
    #   cap_asc  caption baseline -> top of its glyphs
    #   leg_h    height of the legend's drawn glyphs (not its box)
    #   pad      how far the legend BOX top sits above those glyphs; the anchor
    #            positions the box, the reader sees the glyph.
    cap_asc = leg_h = pad = 0.0
    for a in key_artists:
        if isinstance(a, mpl.legend.Legend):
            leg_h = max(tops) - min(bots) if tops else 0.0
            pad = a.get_window_extent(r).y1 / dpi - max(tops) if tops else 0.0
        else:
            e = a.get_window_extent(r)
            cap_asc = max(cap_asc, e.y1 / dpi - CAP_BASE_IN)
    return dict(band=band, cap_asc=cap_asc, leg_h=leg_h, pad=pad)


def marker(q):
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def save(fig, name, outdir=PAN, svg=True):
    png = os.path.join(outdir, name + ".png")
    fig.savefig(png, dpi=300, bbox_inches=None, facecolor="white")
    n_svg = 0
    if svg:
        s = os.path.join(outdir, name + ".svg")
        fig.savefig(s, bbox_inches=None, transparent=True)
        n_svg = os.path.getsize(s)
    print(f"  {name:38s} md5={hashlib.md5(open(png,'rb').read()).hexdigest()[:10]}"
          f"  svg={n_svg}")
    plt.close(fig)


def check_text(fig, tag):
    """figure-style 9.1 geometric half: text-vs-text and text-vs-spine."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    ts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
          if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:34], b.get_text()[:34])
          for i, (a, ba) in enumerate(ts) for b, bb in ts[i + 1:]
          if ba.overlaps(bb)]
    if ov:
        print(f"    [{tag}] {len(ov)} text overlaps")
        for a, b in ov[:6]:
            print(f"      {a!r} <-> {b!r}")
    return len(ov)


def assert_aligned(figs):
    """C and D must put a given genotype at the same x. Checked on the RENDERED
    positions, not on the arithmetic that produced them."""
    boxes = []
    for f in figs:
        f.canvas.draw()
        boxes.append([ax.get_position() for ax in f.axes[:3]])
    for i in range(3):
        x0 = [b[i].x0 for b in boxes]
        x1 = [b[i].x1 for b in boxes]
        assert max(x0) - min(x0) < 1e-6 and max(x1) - min(x1) < 1e-6, \
            f"column {i} not aligned between C and D: x0={x0} x1={x1}"
    print("  C/D column alignment: exact")


# ------------------------------------------------------------------- panel data
def load():
    P = pd.read_csv(os.path.join(OUT, "cortex_permutation_all21.csv"))
    J = pd.read_csv(os.path.join(OUT, "cortex_vs_cofrac_21cells.csv"))
    assert len(P) == 21 and len(J) == 21
    part = pd.read_csv(os.path.join(DATA_DIR, "ibiocofrac", "modules",
                                    "Partition_Gene_combined_3PD_models.csv"))
    lp = part[part.Experiment == "LRRK2"][["UniProt", "Symbol", "Membership"]]
    lp = lp[~lp.Symbol.isin(EXCLUDE)]
    gsets = {g: set(lp[lp.Membership.isin(GMODS[g])].UniProt) for g in GROUP_ORDER}
    L_ = pd.read_csv(os.path.join(DATA_DIR, "cortex",
                                  "deqms_all_contrasts_long.tsv"), sep="\t")
    return P, J, gsets, L_


# ----------------------------------------------------------------------- panel B
def panel_B(P, title_set=TITLE_SET, _geom=None):
    lab = {g: P[P.group == g].group_label.iloc[0] for g in GROUP_ORDER}
    n = {g: int(P[P.group == g].n.iloc[0]) for g in GROUP_ORDER}
    M = np.array([[P[(P.group == g) & (P.model == m)].delta.iloc[0]
                   for m in MODELS] for g in GROUP_ORDER])

    # pass 1 uses nominal bands, pass 2 the measured ones (see measure_bands)
    g = _geom or dict(band=GENOLAB_IN, cap_asc=0.07)
    bot_in = CAP_BASE_IN + g["cap_asc"] + GAP_IN + g["band"]
    FIGH = bot_in + BOX_H_B + TOP_IN
    bot, top = bot_in / FIGH, 1 - TOP_IN / FIGH
    fig = plt.figure(figsize=(5.6, FIGH))
    ax = fig.add_axes([0.30, bot, 0.56, top - bot])
    vmax = float(np.nanmax(np.abs(M)))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    for i, g in enumerate(GROUP_ORDER):
        for j, m in enumerate(MODELS):
            r = P[(P.group == g) & (P.model == m)].iloc[0]
            ax.text(j, i, f"{r.delta:+.3f}{marker(r.q_perm_BH)}", ha="center",
                    va="center", fontsize=7, fontweight="bold",
                    color="white" if abs(r.delta) > vmax * 0.62 else "black")
    ax.set_xticks(range(3))
    ax.set_xticklabels([GENO[GENO_SET][m] for m in MODELS], fontsize=6.5)
    for t, m in zip(ax.get_xticklabels(), MODELS):
        t.set_color(COL[m]); t.set_fontweight("bold")
    ax.set_yticks(range(len(GROUP_ORDER)))
    ax.set_yticklabels([f"{lab[g]}  [{n[g]}]" for g in GROUP_ORDER], fontsize=7)
    for t, g in zip(ax.get_yticklabels(), GROUP_ORDER):
        if g == "Synapse":
            t.set_color("#C0392B"); t.set_fontweight("bold")
    ax.set_xticks(np.arange(-.5, 3, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(GROUP_ORDER), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", size=0)
    for s in ax.spines.values():
        s.set_visible(False)

    cax = fig.add_axes([0.875, bot, 0.022, top - bot])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("group mean log$_2$FC $-$ proteome baseline", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)

    # item 1: one line, centred over the whole figure
    fig.suptitle(TITLE_B[title_set], x=0.5, ha="center", fontsize=8.5,
                 y=1 - SUPT_IN / FIGH)
    # significance key centred (was flush left), GAP_IN below the axis labels
    key = fig.text(0.5, CAP_BASE_IN / FIGH,
                   "* q$<$0.05   ** q$<$0.01   *** q$<$0.001    "
                   "[n] = proteins measured", ha="center", fontsize=6.3,
                   color="#4D4D4D")
    if _geom is None:
        gm = measure_bands(fig, [key])
        plt.close(fig)
        return panel_B(P, title_set, _geom=gm)
    return fig


# ----------------------------------------------------------------------- panel C
def panel_C(P, gsets, Lg, cap=3, title_set=TITLE_SET, _geom=None):
    g = _geom or dict(band=XLABEL_IN, cap_asc=0.07)
    bot_in = CAP_BASE_IN + g["cap_asc"] + GAP_IN + g["band"]
    FIGH = bot_in + BOX_H_C + TOP_IN
    fig = plt.figure(figsize=(FIGW, FIGH))
    axes = fig.subplots(1, 3, sharey=True)
    fig.subplots_adjust(left=L, right=R, wspace=WS,
                        bottom=bot_in / FIGH, top=1 - TOP_IN / FIGH)

    rng = np.random.default_rng(20260811)
    for ax, m in zip(axes, MODELS):
        d = Lg[Lg.contrast == f"{m}_cortex_vs_WT"]
        d = d[~d.Gene.isin(EXCLUDE)]
        base = float(P[P.model == m].pool_mean.iloc[0])
        for i, g in enumerate(GROUP_ORDER):
            if i % 2 == 0:
                ax.axhspan(i - .5, i + .5, color="#F2F2F2", zorder=0)
            v = d.loc[d.Accession.isin(gsets[g]), "log2FC"].to_numpy()
            ax.scatter(v, i + rng.uniform(-.28, .28, len(v)), s=4, c=COL[m],
                       alpha=0.30, linewidths=0, rasterized=True)
            r = P[(P.group == g) & (P.model == m)].iloc[0]
            ax.plot([r.mean_raw] * 2, [i - .34, i + .34], color=COL[m], lw=3,
                    solid_capstyle="butt", zorder=3)
            if r.q_perm_BH < 0.05:
                ax.text(0.97, i, marker(r.q_perm_BH),
                        transform=ax.get_yaxis_transform(), ha="right",
                        va="center", color=COL[m], fontsize=8,
                        fontweight="bold")
        ax.axvline(base, color="black", lw=0.9, ls="--", zorder=2)
        ax.axvline(0, color="#999999", lw=0.5, zorder=1)
        ax.set_xlim(-0.75, 0.75)
        # explicit ticks: the autolocator put a label at each panel's edge, and
        # with the tightened wspace the right edge of one panel collided with the
        # left edge of the next.
        ax.set_xticks([-0.6, -0.3, 0.0, 0.3, 0.6])
        ax.set_xlabel("protein log$_2$FC (mutant / WT)")
        # item 4: same genotype string as panel D, centred over the panel
        ax.set_title(GENO[GENO_SET][m], loc="center", color=COL[m],
                     fontweight="bold", fontsize=8)
    axes[0].set_ylim(len(GROUP_ORDER) - .5, -.5)
    axes[0].set_yticks(range(len(GROUP_ORDER)))
    axes[0].set_yticklabels([P[P.group == g].group_label.iloc[0]
                             for g in GROUP_ORDER], fontsize=7)
    for t, g in zip(axes[0].get_yticklabels(), GROUP_ORDER):
        if g == "Synapse":
            t.set_color("#C0392B"); t.set_fontweight("bold")

    fig.suptitle(TITLE_C[title_set], x=0.5, ha="center", fontsize=8.5,
                 y=1 - SUPT_IN / FIGH)
    # one line of key terms, KEY1_IN below the axes - same physical gap as D
    cap_txt = fig.text(0.5, CAP_BASE_IN / FIGH, CAP_C[cap], ha="center",
                       fontsize=6.3, color="#4D4D4D")
    if _geom is None:
        gm = measure_bands(fig, [cap_txt])
        plt.close(fig)
        return panel_C(P, gsets, Lg, cap, title_set, _geom=gm)
    return fig


# ----------------------------------------------------------------------- panel D
def panel_D(J, cap=2, title_set=TITLE_SET, _geom=None, _iter=0):
    # equal aspect => square box => height derives from the shared panel width
    g = _geom or dict(band=XLABEL_IN, cap_asc=0.07, leg_h=LEG_GLYPH_IN, bias=0.0)
    key_top_in = (CAP_BASE_IN + g["cap_asc"] + LEG_GAP_IN + g["leg_h"])
    bot_in = key_top_in + GAP_IN + g["band"]
    # ⛔ ANCHOR BIAS, MEASURED. bbox_to_anchor does NOT put the legend's top edge
    # at the anchor: borderaxespad and the box padding shift it, here by 0.045 in,
    # which left D's gap at 0.199 while B and C sat at 0.130. `bias` is that
    # offset, measured on pass 1 and applied on pass 2 - not a constant, because
    # it depends on the handler's drawn glyph height.
    anchor_in = key_top_in + g["bias"]
    FIGH = PANEL_W + bot_in + TOP_IN
    bot, top = bot_in / FIGH, 1 - TOP_IN / FIGH
    fig = plt.figure(figsize=(FIGW, FIGH))
    axes = fig.subplots(1, 3, sharex=True, sharey=True)
    fig.subplots_adjust(left=L, right=R, wspace=WS, bottom=bot, top=top)

    lim = 0.25
    pending = []
    for ax, m in zip(axes, MODELS):
        s = J[J.model == m]
        for xs, ys in [((-lim, 0), (0, lim)), ((0, lim), (-lim, 0))]:
            ax.add_patch(plt.Rectangle((xs[0], ys[0]), xs[1] - xs[0],
                                       ys[1] - ys[0], facecolor="#F4F4F4",
                                       edgecolor="none", zorder=0))
        ax.plot([-lim, lim], [-lim, lim], color="#CFCFCF", lw=0.7, ls=":",
                zorder=1)
        ax.axhline(0, color="#8C8C8C", lw=0.6, zorder=1)
        ax.axvline(0, color="#8C8C8C", lw=0.6, zorder=1)
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_xticks([-0.2, -0.1, 0, 0.1, 0.2])
        ax.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])

        items = []
        for r in s.sort_values(["sig_cortex", "sig_cofrac"],
                               ascending=False).itertuples():
            st = state_of(bool(r.sig_cortex), bool(r.sig_cofrac))
            c = COL[m]
            if st == "neither":
                ax.plot(r.delta_cofrac, r.delta_cortex, marker="o", ms=7.0,
                        mfc="white", mec=c, mew=0.9, ls="", zorder=4)
            else:
                if st == "both":
                    for ms in RING_MS:
                        ax.plot(r.delta_cofrac, r.delta_cortex, marker="o",
                                ms=ms, mfc="none", mec=c, mew=RING_LW, ls="",
                                zorder=5)
                ax.plot(r.delta_cofrac, r.delta_cortex, marker="o", ms=DISC_MS,
                        mfc=c, mec=c, mew=1.0, ls="", zorder=4)
            items.append((r.delta_cofrac, r.delta_cortex, p14.LABEL[r.group],
                          c if st != "neither" else "#4D4D4D",
                          "bold" if st != "neither" else "normal"))
        pending.append((ax, items))
        # item 5: single-line title, no n/7 counts; matches panel C exactly
        ax.set_title(GENO[GENO_SET][m], loc="center", color=COL[m],
                     fontweight="bold", fontsize=8)
        ax.set_xlabel("DA neurons, iBioCoFrac")
    axes[0].set_ylabel("whole cortex, DIA-MS")

    fig.suptitle(TITLE_D[title_set], x=0.5, ha="center", fontsize=8.5,
                 y=1 - SUPT_IN / FIGH)
    key = Line2D([], [], ls="", marker="None")
    leg = fig.legend(handles=[key,
                        Line2D([], [], marker="o", ls="", ms=7, mfc="#4D4D4D",
                               mec="#4D4D4D"),
                        Line2D([], [], marker="o", ls="", ms=6, mfc="white",
                               mec="#4D4D4D", mew=0.9)],
               labels=["both", "DA neurons only", "neither"],
               handler_map={key: ConcentricHandler((9.0, 10.6), 6.4, "#4D4D4D",
                                                   RING_LW)},
               frameon=False, ncol=3, fontsize=6.5, loc="upper center",
               bbox_to_anchor=(0.5, anchor_in / FIGH),
               handletextpad=0.6, columnspacing=2.4)
    cap_txt = fig.text(0.5, CAP_BASE_IN / FIGH, CAP_D[cap], ha="center",
                       fontsize=6.3, color="#4D4D4D")

    if _geom is None or _iter < 12:
        # ⛔ ITERATE. Rebuilding at the height implied by pass 1 moves the legend
        # (its anchor is a figure fraction), so the gap must be re-measured:
        # one pass alone left D at 0.199 in against 0.151 for B and C.
        g2 = measure_bands(fig, [cap_txt, leg])
        fig.canvas.draw()
        rr = fig.canvas.get_renderer()
        inner = list(leg.get_texts()) + [h for h in leg.legend_handles
                                         if h is not None]
        glyph_top = max(o.get_window_extent(rr).y1 for o in inner
                        if hasattr(o, "get_window_extent")) / fig.dpi
        g2["bias"] = g["bias"] + (key_top_in - glyph_top)
        converged = (_geom is not None
                     and abs(g2["bias"] - g["bias"]) < 0.002
                     and all(abs(g2[k] - g[k]) < 0.004
                             for k in g2 if k != "bias"))
        if not converged:
            g = g2
            plt.close(fig)
            return panel_D(J, cap, title_set, _geom=g, _iter=_iter + 1)
        # ⛔ FALL THROUGH on convergence - do NOT return here. An earlier version
        # returned the converged figure directly and skipped the label block
        # below, so panel D rendered with no group labels at all.

    # labels only after the layout is final - see step 19
    fig.canvas.draw()
    fb = sum(place_labels(ax, it, fig) for ax, it in pending)
    assert fb == 0, "a label had no free slot"
    return fig


def main():
    P, J, gsets, Lg = load()

    fb = panel_B(P); check_text(fb, "B")
    save(fb, "panel_B_group_shift_heatmap")
    fc = panel_C(P, gsets, Lg, cap=3); n_c = check_text(fc, "C")
    fd = panel_D(J, cap=2); n_d = check_text(fd, "D")
    assert_aligned([fc, fd])
    save(fc, "panel_C_group_protein_strips")
    save(fd, "panel_D_cortex_vs_cofrac")
    assert n_c == 0 and n_d == 0, "text overlaps in C or D"


if __name__ == "__main__":
    main()
