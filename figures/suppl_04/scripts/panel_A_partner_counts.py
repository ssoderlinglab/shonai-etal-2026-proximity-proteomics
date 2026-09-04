"""
Supplemental Figure 4, panel A: absolute number of previously reported partners recovered in each WT proximity proteome, per database and for the union of the three.

Input
    data/ibioid/orthogonal/db_matrix_values_3DB.csv
    data/ibioid/orthogonal/recovered_partners_3DB.csv
Output
    figures/suppl_04/output/panel_A_partner_counts_{outlined,editable}.svg
    figures/suppl_04/output/panel_A_partner_counts.png
Run
    python3 figures/suppl_04/scripts/panel_A_partner_counts.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, pillow

Counts only; no statistical test is computed or shown.

The value axis is square-root scaled and keeps a true zero. The union bar is split into partners reported in one database and in two or more; the three single-database bars are not split.

Every input value is asserted on load, so a changed input fails instead of producing a plausible wrong figure.

Two SVGs are written and are not interchangeable: _outlined converts text to outlines, _editable keeps live text and needs the named fonts to render. The PNG is opaque and is for screen checking only.
"""

import hashlib
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.transforms import Bbox

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
EVIDENCE = os.path.join(REPO_ROOT, "data", "ibioid", "orthogonal")
OUT = os.path.join(FIGURE_DIR, "output")
NAME = "panel_A_partner_counts"

# COLOR_SYSTEM.md section 2 — PD model colours, not negotiable
MODEL_COLOR = {"VPS35": "#7570B3", "SNCA": "#D95F02", "LRRK2": "#1B9E77"}
BAITS = ["VPS35", "SNCA", "LRRK2"]                 # bait order in both panels
DISPLAY = {"VPS35": "Vps35", "SNCA": "aSyn", "LRRK2": "Lrrk2"}
DBS = ["BioGRID", "IntAct", "STRING_physical"]

TITLE = "Independent reports of iBioID proximity partners"      # T4
VLABEL = "Previously reported proximity partners (count)"        # X1
TICKS_X = ["any database", "BioGRID", "IntAct", "STRING\n(physical)"]

SLOT_W_PT, SLOT_H_PT = 707.0, 224.5
FIGSIZE = (SLOT_W_PT / 72.0, SLOT_H_PT / 72.0)

BARW, GAP, YMAX = 0.88, 0.018, 150
LEFT = 0.0575          # renderer-probed floor for 10 pt y ticks. Under an
                       # explicit page box a margin below this CLIPS the
                       # rotated y-axis label instead of overflowing.
BOTTOM, TOP = 0.185, 0.795
POS = [0, 1.42, 2.42, 3.42]
YTICKS = [0, 10, 25, 50, 100, 150]
FS_TICK, FS_COUNT, FS_BAIT, FS_TITLE, FS_VLAB, FS_LEG = 10, 8.5, 10, 11, 9, 6.8

BASE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.0,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "patch.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "ytick.major.size": 3.0,
    "pdf.fonttype": 42,          # TrueType, not Type 3 (Cell Press)
    "ps.fonttype": 42,
    "savefig.dpi": 300.0,
    "savefig.bbox": "standard",  # never global tight; see trap 1
    "savefig.transparent": True,
}


def tint(c, f):
    r, g, b = matplotlib.colors.to_rgb(c)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def load():
    m = pd.read_csv(os.path.join(EVIDENCE, "db_matrix_values_3DB.csv"))
    p = pd.read_csv(os.path.join(EVIDENCE, "recovered_partners_3DB.csv"))
    assert len(m) == 9 and len(p) == 236, (len(m), len(p))
    rec = m.pivot(index="db", columns="bait", values="recovered").loc[DBS, BAITS]
    assert rec.values.tolist() == [[29, 140, 22], [18, 20, 45], [16, 11, 10]]
    uniq = p.groupby("bait")["gene"].nunique().reindex(BAITS)
    assert uniq.tolist() == [41, 141, 54]
    tier = (p.groupby(["bait", "n_databases"]).size()
            .unstack(fill_value=0).reindex(BAITS)[[1, 2, 3]])
    assert tier.values.tolist() == [[25, 10, 6], [116, 20, 5], [38, 9, 7]]
    return rec, uniq, tier


def build(rec, uniq, tier):
    fig = plt.figure(figsize=FIGSIZE)
    span = (0.995 - LEFT - 2 * GAP) / 3.0
    axes = [fig.add_axes([LEFT + j * (span + GAP), BOTTOM, span, TOP - BOTTOM])
            for j in range(3)]
    for ax, b in zip(axes, BAITS):
        col, pale = MODEL_COLOR[b], tint(MODEL_COLOR[b], 0.55)
        tot = int(uniq[b])
        one = int(tier.loc[b, 1])
        multi = int(tier.loc[b, 2] + tier.loc[b, 3])
        assert one + multi == tot, (b, one, multi, tot)
        ax.bar(POS[0], one, BARW, color=pale, edgecolor="none")
        ax.bar(POS[0], multi, BARW, bottom=one, color=col, edgecolor="white",
               lw=0.5)
        ax.text(POS[0], tot + YMAX * 0.018, tot, ha="center",
                fontsize=FS_COUNT, fontweight="bold")
        for i, db in enumerate(DBS, start=1):
            v = int(rec.loc[db, b])
            ax.bar(POS[i], v, BARW, color=pale, edgecolor="none")
            ax.text(POS[i], v + YMAX * 0.018, v, ha="center", fontsize=FS_COUNT)
        ax.axvline((POS[0] + POS[1]) / 2.0, color="#999999", lw=0.7,
                   ls=(0, (1.6, 1.8)), ymax=0.97)
        # square-root value axis: compresses the aSyn/BioGRID outlier while
        # keeping a REAL zero, so bar length stays a monotone function of the
        # count with no arbitrary lower limit to justify (a log axis has none).
        ax.set_yscale("function",
                      functions=(lambda v: np.sqrt(np.abs(v)), lambda v: v ** 2))
        ax.set_ylim(0, YMAX)
        ax.set_yticks(YTICKS)
        ax.tick_params(axis="y", labelsize=FS_TICK)
        ax.set_xticks(POS)
        ax.set_xticklabels(TICKS_X, fontsize=FS_TICK)
        ax.tick_params(axis="x", length=0)
        for k, t in enumerate(ax.get_xticklabels()):
            t.set_color(col if k == 0 else "#4A4A4A")
        half = BARW / 2 + 0.10
        ax.set_xlim(POS[0] - half, POS[-1] + half)
        ax.text(0.5, 1.045, DISPLAY[b], transform=ax.transAxes, ha="center",
                va="bottom", fontsize=FS_BAIT, fontweight="bold", color=col)
        if ax is not axes[0]:
            ax.set_yticklabels([])
            ax.spines["left"].set_visible(False)
            ax.tick_params(axis="y", length=0)
    axes[0].set_ylabel(VLABEL, fontsize=FS_VLAB)
    fig.text(0.5, 0.945, TITLE, ha="center", va="bottom", fontsize=FS_TITLE)
    axes[0].legend(
        handles=[Patch(facecolor=MODEL_COLOR["VPS35"],
                       label="reported in $\\geq$2 databases"),
                 Patch(facecolor=tint(MODEL_COLOR["VPS35"], 0.55),
                       label="reported in 1 database")],
        fontsize=FS_LEG, frameon=False, loc="upper left",
        bbox_to_anchor=(0.0, 1.005), handlelength=1.1, handletextpad=0.45,
        labelspacing=0.35, borderpad=0.0)
    return fig, axes


def audit(fig, path_svg, label):
    """Verify the delivered file rather than trusting the call that wrote it."""
    import re
    from PIL import Image
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    tx = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
          if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:18], b.get_text()[:18])
          for i, (a, ba) in enumerate(tx)
          for b, bb in tx[i + 1:] if ba.overlaps(bb)]
    head = open(path_svg).read(2000)
    w = float(re.search(r'width="([\d.]+)pt"', head).group(1))
    h = float(re.search(r'height="([\d.]+)pt"', head).group(1))
    body = open(path_svg).read()
    white_bg = "#ffffff" in body.split("<g")[0].replace(" ", "").lower()
    editable = "<text" in body
    png = path_svg.replace(".svg", ".png").replace("_editable", "").replace(
        "_outlined", "")
    clip = "n/a"
    if os.path.exists(png):
        a = np.array(Image.open(png).convert("RGB")).astype(int).mean(axis=2)
        ink = a < 250
        clip = f"col0={ink[:, 0].sum()} row0={ink[0, :].sum()}"
    print(f"  {label:9s} page {w:6.1f} x {h:5.1f} pt "
          f"(target {SLOT_W_PT} x {SLOT_H_PT}) | white bg {white_bg} | "
          f"live text {editable} | collisions {len(ov)} | edge ink {clip}")
    if ov:
        print("      COLLISIONS:", ov[:5])
    assert abs(w - SLOT_W_PT) < 0.5 and abs(h - SLOT_H_PT) < 0.5, (w, h)
    assert not white_bg, "SVG carries a white background rectangle"
    return len(ov)


def main():
    os.makedirs(OUT, exist_ok=True)
    rec, uniq, tier = load()
    print("inputs asserted OK; every value is read, none recomputed")
    print(f"target page {SLOT_W_PT} x {SLOT_H_PT} pt "
          f"(aspect {SLOT_W_PT/SLOT_H_PT:.2f}:1), from the .ai ink extent\n")

    bad = 0
    for kind, fonttype in (("outlined", "path"), ("editable", "none")):
        plt.rcParams.update(BASE_RC)
        plt.rcParams["svg.fonttype"] = fonttype
        fig, _ = build(rec, uniq, tier)
        box = Bbox([[0, 0], fig.get_size_inches()])
        svg = os.path.join(OUT, f"{NAME}_{kind}.svg")
        fig.savefig(svg, bbox_inches=box, transparent=True, facecolor="none",
                    edgecolor="none")
        if kind == "outlined":       # one opaque PNG, for screen checking only
            fig.savefig(os.path.join(OUT, f"{NAME}.png"), dpi=300,
                        bbox_inches=box, transparent=False, facecolor="white",
                        edgecolor="none")
        bad += audit(fig, svg, kind)
        plt.close(fig)

    # A whole-file SVG md5 is not an identity check: matplotlib stamps a
    # <dc:date> into <metadata> and mints random clip and glyph ids on every
    # build, so the raw md5 drifts while the artwork does not. Byte size plus
    # the canonical md5 below, computed with the metadata lines deleted and the
    # generated ids normalised, are stable. Delete WHOLE LINES from <metadata>
    # to </metadata>: an inline substitution of the tags leaves the surrounding
    # newlines and yields a different hash.
    print("\nbyte size + CANONICAL md5 (metadata lines dropped, generated ids "
          "normalised).\nThe raw whole-file md5 is NOT a valid check; neither is "
          "the PNG's.")
    import re as _re
    for f in sorted(x for x in os.listdir(OUT)
                    if x.startswith(NAME) and x.endswith(".svg")):
        p = os.path.join(OUT, f)
        keep, drop = [], False
        for line in open(p):
            if "<metadata>" in line:
                drop = True
            if not drop:
                keep.append(line)
            if "</metadata>" in line:
                drop = False
        body = "".join(keep)
        body = _re.sub(r'(id|href)="#?[pm][0-9a-f]{10}"', r'\1="#GENID"', body)
        body = _re.sub(r"url\(#[pm][0-9a-f]{10}\)", "url(#GENID)", body)
        print(f"  {f:56s} {os.path.getsize(p):8d} B  "
              f"{hashlib.md5(body.encode()).hexdigest()}")
    print("\ntext collisions across both files:", bad)
    assert bad == 0, "text collision in the adopted panel"


if __name__ == "__main__":
    main()
