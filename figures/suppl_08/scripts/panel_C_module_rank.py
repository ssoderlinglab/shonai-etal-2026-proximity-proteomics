"""
Supplemental Figure 8, panel C: ranking of all co-fractionation modules by their enrichment in the selectively depleted human dopamine neuron subtype.

Input
    data/human/module_specificity_42modules.csv
Output
    figures/suppl_08/output/panel_C_module_rank.{png,svg}
    figures/suppl_08/output/panel_C_checksums.csv
Run
    python3 figures/suppl_08/scripts/panel_C_module_rank.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib

Odds ratios and Benjamini-Hochberg q values are read from the frozen input; nothing is recomputed here.
"""

import hashlib
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams.update({
    "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
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

SPEC = pd.read_csv(os.path.join(IN, "module_specificity_42modules.csv"))

# Palette. The AGTR1 red is the value panels A and B use, so
# "this bar concerns SOX6_AGTR1" reads consistently across the figure.
# The two other significant modules keep a neutral dark grey rather than a
# second hue: they are significant, not a second category of finding.
C_M11 = "#BA1028"
C_SIG = "#6E6E6E"
C_NS = "#D2D2D2"
TEXT = "#1A1A1A"
GUIDE = "#6E6E6E"

# Y-axis wording. Final form chosen by Daichi 2026-08-12 (variant B').
#
# SINGULAR "subtype" is load-bearing and must not be pluralised. SOX6_AGTR1 is
# ONE subtype, and the whole point of the source's claim is that a SINGLE
# subpopulation is selectively depleted. A plural would reinstate exactly the
# lineage-level reading this figure set was corrected to avoid.
#
# "selectively" is OMITTED here: every candidate carrying it overflowed a line
# at this axis height. It is the core of the source's claim and therefore MUST
# appear in the caption.
#
# "-log10(q)" is kept: without it the 0-6 tick scale has no stated meaning and
# the panel cannot be read standalone.
#
# Line-length budget, measured at 8 pt against this axis height: one line fits
# ~34 characters (88% of the axis). The two lines below sit at 83% and 78%.
YLAB_B = ("Enrichment in PD-vulnerable markers\n"
          "(SOX6_AGTR1 subtype), -log10(q)")

# Gene counts are printed on all three named bars (Daichi 2026-08-12). They are
# what makes M35's OR of 18.4 legible as a small-module artefact (3 of 6) rather
# than a stronger finding than M11 (16 of 62).
VARIANTS = ["vB"]


def prepare():
    sp = SPEC.sort_values("q_BH").reset_index(drop=True)
    sp["nlq"] = -np.log10(sp["q_BH"].clip(lower=1e-300))
    return sp


def draw(ax, sp, variant):
    cols = [C_M11 if r["module"] == 11 else (C_SIG if r["q_BH"] < 0.05 else C_NS)
            for _, r in sp.iterrows()]
    x = np.arange(len(sp))
    ax.bar(x, sp["nlq"], color=cols, width=0.78, zorder=3)

    ax.axhline(-np.log10(0.05), ls=(0, (5, 4)), lw=0.8, color=GUIDE, zorder=2)
    ax.text(len(sp) * 0.985, -np.log10(0.05) + 0.10, "q = 0.05",
            ha="right", va="bottom", fontsize=6, color=GUIDE)

    # index of each named module in the ranked order
    idx = {int(r["module"]): i for i, r in sp.iterrows()}

    def lab(mod):
        r = sp[sp["module"] == mod].iloc[0]
        return int(r["a"]), int(r["n_query"]), r["OR"]

    # --- M11: moved OFF the y-axis (it previously collided with it) ---
    i11 = idx[11]
    a11, n11, or11 = lab(11)
    m11_txt = f"M11\n{a11}/{n11} genes"
    ax.annotate(m11_txt, xy=(i11, sp["nlq"][i11]),
                xytext=(14, -2), textcoords="offset points",
                ha="left", va="top", fontsize=7, fontweight="bold",
                color=C_M11, linespacing=1.25,
                arrowprops=dict(arrowstyle="-", lw=0.7, color=C_M11,
                                shrinkA=1, shrinkB=2))

    # --- M19 and M35: named individually, no longer "2 other modules" ---
    i19, i35 = idx[19], idx[35]
    a19, n19, or19 = lab(19)
    a35, n35, or35 = lab(35)
    t19, t35 = f"M19  {a19}/{n19}", f"M35  {a35}/{n35}"

    ax.annotate(t19, xy=(i19, sp["nlq"][i19]), xytext=(9, 3),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=6.4, color=C_SIG,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C_SIG,
                                shrinkA=1, shrinkB=2))
    ax.annotate(t35, xy=(i35, sp["nlq"][i35]), xytext=(9, 3),
                textcoords="offset points", ha="left", va="bottom",
                fontsize=6.4, color=C_SIG,
                arrowprops=dict(arrowstyle="-", lw=0.6, color=C_SIG,
                                shrinkA=1, shrinkB=2))

    ax.set_ylabel(YLAB_B)
    ax.set_xlabel("iBioCoFrac modules (n = 42), ranked")
    ax.set_xticks([])
    ax.set_xlim(-0.9, len(sp) - 0.1)
    ax.set_ylim(0, sp["nlq"].max() * 1.16)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=2.5, width=0.7)
    return dict(M11=(a11, n11, or11), M19=(a19, n19, or19), M35=(a35, n35, or35))


def main():
    sp = prepare()
    manifest, rep = [], None
    for v in VARIANTS:
        fig = plt.figure(figsize=(3.5, 2.75))
        ax = fig.add_axes([0.30, 0.16, 0.66, 0.80])
        rep = draw(ax, sp, v)
        for ext in ("png", "svg"):
            p = os.path.join(OUT, f"panel_C_module_rank.{ext}")
            fig.savefig(p, transparent=True, bbox_inches="tight",
                        pad_inches=0.02, dpi=600 if ext == "png" else None)
            manifest.append((f"panel_C_module_rank.{ext}", os.path.getsize(p),
                             hashlib.md5(open(p, "rb").read()).hexdigest()
                             if ext == "png" else "SVG-md5-not-a-valid-check"))
        plt.close(fig)

    mp = os.path.join(OUT, "panel_C_checksums.csv")
    with open(mp, "w") as f:
        f.write("# PNG md5 is a valid identity check. SVG md5 is NOT:\n"
                "# matplotlib stamps a <metadata> date and randomises clip and\n"
                "# glyph ids, so a re-render changes SVG bytes while the\n"
                "# artwork is identical.\n")
        f.write("file,bytes,md5\n")
        for r in sorted(manifest):
            f.write(f"{r[0]},{r[1]},{r[2]}\n")

    for k, (a, n, o) in rep.items():
        print(f"{k}: overlap {a}/{n}  OR {o:.2f}")
    print(f"significant (q<0.05): {int((sp['q_BH'] < 0.05).sum())} of {len(sp)}")
    print(f"{len(manifest)} files -> {OUT}")


if __name__ == "__main__":
    main()
