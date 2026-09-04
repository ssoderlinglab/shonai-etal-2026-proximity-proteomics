"""
Supplemental Figure 8, panel D: overlap of each proximity proteome with the marker set of the selectively depleted human dopamine neuron subtype.

Input
    data/human/proxiome_AGTR1_overlap.csv
Output
    figures/suppl_08/output/panel_D_proxiome_AGTR1.{png,svg}
    figures/suppl_08/output/panel_D_checksums.csv
Run
    python3 figures/suppl_08/scripts/panel_D_proxiome_agtr1.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib

Odds ratios, confidence intervals and q values are read from the frozen input; nothing is recomputed here. The q values drawn are corrected within the six query sets.
"""

import hashlib
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

matplotlib.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "legend.fontsize": 6.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
})

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
IN = os.path.join(REPO_ROOT, "data", "human")
OUT = os.path.join(FIGURE_DIR, "output")
os.makedirs(OUT, exist_ok=True)

D = pd.read_csv(os.path.join(IN, "proxiome_AGTR1_overlap.csv"))

# Bait colours are the reserved COLOR_SYSTEM.md model-identity hexes, used here
# with their sanctioned meaning (which PD model), so no new entry is needed.
BAIT = {"aSyn": "#D95F02", "Vps35": "#7570B3", "LRRK2": "#1B9E77"}
NS = "#B8B8B8"
TEXT = "#1A1A1A"
GUIDE = "#6E6E6E"

LBL = {"aSyn_WT": "\u03b1Syn", "Vps35_WT": "VPS35", "LRRK2_WT": "LRRK2",
       "aSyn_A30P_gained": "\u03b1Syn A30P", "Vps35_D620N_gained": "VPS35 D620N",
       "LRRK2_G2019S_gained": "LRRK2 G2019S"}

ORDER = ["aSyn_WT", "Vps35_WT", "LRRK2_WT",
         "aSyn_A30P_gained", "Vps35_D620N_gained", "LRRK2_G2019S_gained"]


def stars(q):
    return "***" if q < 1e-3 else "**" if q < 1e-2 else "*" if q < 0.05 else "n.s."


def rows():
    d = D.set_index("query_set").loc[ORDER].reset_index()
    d["sig"] = d["q_6"] < 0.05
    return d


# --------------------------------------------------------------------------
# v1 — grouped horizontal bars, WT block over gained block
# --------------------------------------------------------------------------

def draw_v1(ax, d):
    y = []
    pos = 0.0
    for i, r in d.iterrows():
        if i == 3:
            pos += 0.9                      # gap between the two groups
        y.append(pos)
        pos += 1.0
    y = np.array(y)

    for yi, (_, r) in zip(y, d.iterrows()):
        col = BAIT[r["bait"]] if r["sig"] else NS
        ax.barh(yi, r["OR"], height=0.7, color=col, zorder=3)
        ax.text(r["OR"] + 0.07, yi, stars(r["q_6"]), va="center", ha="left",
                fontsize=7, color=TEXT if r["sig"] else GUIDE)

    ax.axvline(1.0, ls=(0, (4, 3)), lw=0.8, color=GUIDE, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([LBL[q] for q in d["query_set"]])
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio, proxiome overlap with\nSOX6_AGTR1 markers")
    ax.set_xlim(0, 3.35)
    ax.text(-0.015, y[:3].mean(), "endogenous\nproxiome (WT)",
            transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.6, color=GUIDE, rotation=90, linespacing=1.2)
    ax.text(-0.015, y[3:].mean(), "gained in\nmutant",
            transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=6.6, color=GUIDE, rotation=90, linespacing=1.2)


# --------------------------------------------------------------------------
# v2 — forest plot: conditional OR with 95% CI, log x axis
# --------------------------------------------------------------------------

def draw_v2(ax, d, group_style="header", shade="#F0F0F0"):
    """
    Forest plot. Three faults of the first pass are fixed here, all of which
    Daichi flagged 2026-08-12:

    1. The rotated group captions were drawn in the y-tick gutter and landed ON
       the tick labels. Group identity is now carried by a HEADER ROW above
       each block (group_style="header") or by a bracket outside the axes
       (group_style="bracket") — neither occupies the tick gutter.
    2. The overlap counts were pinned at x = 6.4 in DATA coordinates, so a wide
       CI ran straight through them (aSyn A30P's upper bound is 5.62). They now
       live in their own column OUTSIDE the axes, in axes coordinates, where no
       interval can reach them.
    3. The two blocks did not read as different categories. The gap is wider,
       a hairline separates them, and the header rows name them.

    The x range is set from the data (widest CI = 5.62) rather than hardcoded,
    so a later data change cannot silently clip an interval.
    """
    n_wt = int((d["group"] == "WT").sum())

    # y positions: header rows at -0.75 above each block, wide gap between
    y, pos = [], 0.0
    for i in range(len(d)):
        if i == n_wt:
            pos += 1.55
        y.append(pos)
        pos += 1.0
    y = np.array(y)
    hdr_wt, hdr_gn = y[0] - 0.78, y[n_wt] - 0.78

    lo_min, hi_max = d["ci_lo"].min(), d["ci_hi"].max()
    xlo, xhi = lo_min / 2.4, hi_max * 1.75

    for yi, (_, r) in zip(y, d.iterrows()):
        col = BAIT[r["bait"]] if r["sig"] else NS
        ax.plot([r["ci_lo"], r["ci_hi"]], [yi, yi], lw=1.3, color=col,
                solid_capstyle="butt", zorder=3)
        for e in (r["ci_lo"], r["ci_hi"]):
            ax.plot([e, e], [yi - 0.15, yi + 0.15], lw=1.3, color=col, zorder=3)
        ax.plot(r["or_cond"], yi, "o", ms=4.6, color=col, zorder=4,
                mec="white", mew=0.5)
        # counts and significance in dedicated columns OUTSIDE the axes
        ax.text(1.075, yi, f"{int(r['a'])}/{int(r['n_query'])}",
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=6.4, color=GUIDE, clip_on=False)
        ax.text(1.115, yi, stars(r["q_6"]),
                transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=6.8, color=TEXT if r["sig"] else GUIDE, clip_on=False)

    ax.axvline(1.0, ls=(0, (4, 3)), lw=0.8, color=GUIDE, zorder=2)
    ax.set_xscale("log")
    ax.set_xticks([0.25, 0.5, 1, 2, 4])
    ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.set_xlim(xlo, xhi)
    ax.set_yticks(y)
    ax.set_yticklabels([LBL[q] for q in d["query_set"]])
    ax.set_ylim(y[-1] + 0.75, hdr_wt - 0.55)
    ax.set_xlabel("Odds ratio (95% CI) for overlap with\nSOX6_AGTR1 subtype markers")

    # column headers
    ax.text(1.075, hdr_wt, "overlap", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=6.2, color=GUIDE,
            style="italic", clip_on=False)

    if group_style in ("header", "shade"):
        # A header row per block. Daichi 2026-08-12: at x = 0.0 the leading
        # glyph sat ON the spine. Inset so the text starts clear of it.
        for hy, txt in ((hdr_wt, "Endogenous proxiome (WT)"),
                        (hdr_gn, "Gained in mutant")):
            ax.text(0.022, hy, txt, transform=ax.get_yaxis_transform(),
                    ha="left", va="center", fontsize=7, fontweight="bold",
                    color=TEXT, clip_on=False)

    if group_style == "header":
        # A single rule between the blocks. Daichi 2026-08-12 rejected this:
        # to separate the blocks it has to extend past the axes (the tick
        # labels and the count column sit outside them), and a line that
        # overshoots the frame reads as an artefact, not as structure.
        ax.axhline((y[n_wt - 1] + hdr_gn) / 2, xmin=-0.42, xmax=1.10,
                   lw=0.5, color="#D0D0D0", zorder=1, clip_on=False)

    elif group_style == "shade":
        # Shaded band behind the GAINED block, confined to the plotting area.
        #
        # Daichi 2026-08-12, two corrections to the first attempt:
        #   - the band ran outside the axes (xmin=-0.42, xmax=1.16, clip_on
        #     False), so it also tinted the tick labels and the count column
        #     and no longer read as part of the plot. It is now xmin=0, xmax=1
        #     and clipped: inside the frame only.
        #   - it shaded the WT block; the gained block is the one to set apart,
        #     since WT is the panel's primary result and should sit on plain
        #     ground.
        ax.axhspan((y[n_wt - 1] + hdr_gn) / 2, y[-1] + 0.72,
                   xmin=0.0, xmax=1.0, facecolor=shade,
                   edgecolor="none", zorder=0)

    else:  # "bracket": vertical rule + rotated caption outside the tick labels
        for i0, i1, txt in ((0, n_wt - 1, "Endogenous\nproxiome (WT)"),
                            (n_wt, len(y) - 1, "Gained in\nmutant")):
            ax.plot([-0.40, -0.40], [y[i0] - 0.28, y[i1] + 0.28],
                    transform=ax.get_yaxis_transform(), lw=1.1, color=TEXT,
                    clip_on=False, solid_capstyle="butt")
            ax.text(-0.53, (y[i0] + y[i1]) / 2, txt,
                    transform=ax.get_yaxis_transform(), ha="center",
                    va="center", fontsize=6.8, color=TEXT, rotation=90,
                    linespacing=1.25, clip_on=False)


# --------------------------------------------------------------------------
# v3 — paired vertical bars: WT vs gained, one pair per bait
# --------------------------------------------------------------------------

def draw_v3(ax, d):
    baits = ["aSyn", "Vps35", "LRRK2"]
    disp = {"aSyn": "\u03b1Syn", "Vps35": "VPS35", "LRRK2": "LRRK2"}
    x = np.arange(len(baits))
    w = 0.36
    for i, b in enumerate(baits):
        wt = d[(d.bait == b) & (d.group == "WT")].iloc[0]
        gn = d[(d.bait == b) & (d.group == "gained")].iloc[0]
        for off, r, hatch in ((-w / 2, wt, None), (w / 2, gn, "///")):
            col = BAIT[b] if r["sig"] else NS
            ax.bar(i + off, r["OR"], width=w, color=col, zorder=3,
                   hatch=hatch, edgecolor="white", lw=0.0)
            ax.text(i + off, r["OR"] + 0.06, stars(r["q_6"]), ha="center",
                    va="bottom", fontsize=6.6,
                    color=TEXT if r["sig"] else GUIDE)
    ax.axhline(1.0, ls=(0, (4, 3)), lw=0.8, color=GUIDE, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels([disp[b] for b in baits])
    ax.set_ylabel("Odds ratio, proxiome overlap\nwith SOX6_AGTR1 markers")
    ax.set_ylim(0, 3.45)
    ax.legend(handles=[Patch(facecolor="#8A8A8A", label="endogenous (WT)"),
                       Patch(facecolor="#8A8A8A", hatch="///",
                             edgecolor="white", label="gained in mutant")],
              loc="upper right", frameon=False, handlelength=1.2,
              borderpad=0.2, labelspacing=0.35)


# --------------------------------------------------------------------------
# v4 — WT only, with the gained sets as a small inset-style second axis
# --------------------------------------------------------------------------

def draw_v4(ax_main, ax_sub, d):
    for ax, grp, ttl in ((ax_main, "WT", "Endogenous proxiome"),
                         (ax_sub, "gained", "Gained in mutant")):
        sub = d[d.group == grp].reset_index(drop=True)
        y = np.arange(len(sub))
        for yi, (_, r) in zip(y, sub.iterrows()):
            col = BAIT[r["bait"]] if r["sig"] else NS
            ax.barh(yi, r["OR"], height=0.66, color=col, zorder=3)
            ax.text(r["OR"] + 0.07, yi, stars(r["q_6"]), va="center", ha="left",
                    fontsize=6.8, color=TEXT if r["sig"] else GUIDE)
        ax.axvline(1.0, ls=(0, (4, 3)), lw=0.8, color=GUIDE, zorder=2)
        ax.set_yticks(y)
        ax.set_yticklabels([LBL[q] for q in sub["query_set"]])
        ax.invert_yaxis()
        ax.set_xlim(0, 3.35)
        ax.set_title(ttl, fontsize=7.2, color=TEXT, pad=3, loc="left")
    ax_main.set_xticklabels([])
    ax_sub.set_xlabel("Odds ratio, proxiome overlap with\nSOX6_AGTR1 markers")


def main():
    d = rows()
    manifest = []

    def save(fig, stem):
        for ext in ("png", "svg"):
            p = os.path.join(OUT, f"{stem}.{ext}")
            fig.savefig(p, transparent=True, bbox_inches="tight",
                        pad_inches=0.02, dpi=600 if ext == "png" else None)
            manifest.append((f"{stem}.{ext}", os.path.getsize(p),
                             hashlib.md5(open(p, "rb").read()).hexdigest()
                             if ext == "png" else "SVG-md5-not-a-valid-check"))
        plt.close(fig)

    fig = plt.figure(figsize=(3.9, 2.75))
    ax = fig.add_axes([0.26, 0.24, 0.60, 0.70])
    draw_v2(ax, d, "shade", shade="#F0F0F0")
    save(fig, "panel_D_proxiome_AGTR1")

    mp = os.path.join(OUT, "panel_D_checksums.csv")
    with open(mp, "w") as f:
        f.write("# PNG md5 is a valid identity check. SVG md5 is NOT:\n"
                "# matplotlib stamps a date and randomises clip/glyph ids.\n")
        f.write("file,bytes,md5\n")
        for r in sorted(manifest):
            f.write(f"{r[0]},{r[1]},{r[2]}\n")

    for _, r in d.iterrows():
        print(f"{r['query_set']:22s} OR {r['OR']:5.2f} "
              f"[{r['ci_lo']:.2f}-{r['ci_hi']:.2f}] "
              f"{int(r['a'])}/{int(r['n_query'])} q6 {r['q_6']:.2e} "
              f"{stars(r['q_6'])}")
    print(f"{len(manifest)} files -> {OUT}")


if __name__ == "__main__":
    main()
