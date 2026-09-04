"""
Supplemental Figure 4, panel B: database support for each recovered partner reported in two or more databases.

Input
    data/ibioid/orthogonal/recovered_partners_3DB.csv
Output
    figures/suppl_04/output/panel_B_database_support_{outlined,editable}.svg
    figures/suppl_04/output/panel_B_database_support.png
Run
    python3 figures/suppl_04/scripts/panel_B_database_support.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, pillow

Rows are proteins, ordered by descending log2FC within each bait. A filled cell means the database reports that protein as a partner of the bait; a gray cell means it does not. No statistic is computed or shown.

16, 25 and 16 rows and 132 filled cells of 171 are asserted. The database columns hold the string Y or an empty value, so the test is against Y explicitly: a truth test on the raw value fills every cell and the panel loses its information while still looking plausible.

Drawn at final print size; no scaling is applied downstream.
"""

import hashlib
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
EVIDENCE = os.path.join(REPO_ROOT, "data", "ibioid", "orthogonal")
OUT = os.path.join(FIGURE_DIR, "output")
NAME = "panel_B_database_support"

MODEL_COLOR = {"VPS35": "#7570B3", "SNCA": "#D95F02", "LRRK2": "#1B9E77"}
BAITS = ["VPS35", "SNCA", "LRRK2"]
DISPLAY = {"VPS35": "Vps35", "SNCA": "aSyn", "LRRK2": "Lrrk2"}
DBS = ["BioGRID", "IntAct", "STRING_physical"]
DBLAB = ["BioGRID", "IntAct", "STRING"]     # 'physical' is defined in the caption
EMPTY = "#EDEDED"                            # 'not in this database'
EXPECT_MULTI = {"VPS35": 16, "SNCA": 25, "LRRK2": 16}   # partners in >=2 DBs
TITLE = "Proximity partners reported in multiple databases"

# panel A's grid, measured from the assembled ver3
BLOCK_W_PT, BLOCK_GAP_PT = 224.5, 12.7
PAGE_W_PT = 3 * BLOCK_W_PT + 2 * BLOCK_GAP_PT               # 698.8

# B6 parameters
LABEL_W_PT, CELL_GAP_PT = 52.0, 6.0
ROW_PITCH_PT = 18.53          # one gene row, as placed in the artwork
NMAX = 25                     # aSyn
FS_TITLE, FS_BAIT, FS_TICK, FS_GENE = 11, 10, 10, 9.5
TITLE_PAD_PT, BAIT_PAD_PT, TICK_BAND_PT = 16.0, 10.0, 18.0

BASE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.0,
    "pdf.fonttype": 42,          # TrueType, not Type 3 (Cell Press)
    "ps.fonttype": 42,
    "savefig.dpi": 300.0,
    "savefig.bbox": "standard",  # never global tight
    "savefig.transparent": True,
}


def mouse_symbol(s):
    """Source tables carry UPPERCASE symbols; the figure uses mouse nomenclature."""
    return s[0].upper() + s[1:].lower() if len(s) > 1 else s.upper()


def load():
    d = pd.read_csv(os.path.join(EVIDENCE, "recovered_partners_3DB.csv"))
    assert len(d) == 236, len(d)
    m = d[d.n_databases >= 2].copy()
    m["label"] = m.gene.map(mouse_symbol)
    out = {}
    for b in BAITS:
        s = (m[m.bait == b].sort_values("log2FC", ascending=False)
             .reset_index(drop=True))
        assert len(s) == EXPECT_MULTI[b], f"{b}: {len(s)} != {EXPECT_MULTI[b]}"
        out[b] = s
    filled = sum(int(str(getattr(r, db)).strip().upper() == "Y")
                 for b in BAITS for r in out[b].itertuples() for db in DBS)
    total = sum(len(out[b]) for b in BAITS) * len(DBS)
    assert (filled, total) == (132, 171), f"fills {filled}/{total} != 132/171"
    return out


def build(data):
    grid_w = BLOCK_W_PT - LABEL_W_PT
    cell_w = (grid_w - 2 * CELL_GAP_PT) / 3.0
    fig_w = PAGE_W_PT
    fig_h = (TITLE_PAD_PT + FS_BAIT + BAIT_PAD_PT
             + NMAX * ROW_PITCH_PT + TICK_BAND_PT)
    fig = plt.figure(figsize=(fig_w / 72.0, fig_h / 72.0))
    top_of_rows = fig_h - (TITLE_PAD_PT + FS_BAIT + BAIT_PAD_PT)

    for j, bait in enumerate(BAITS):
        s = data[bait]
        n = len(s)
        col = MODEL_COLOR[bait]
        x0 = j * (BLOCK_W_PT + BLOCK_GAP_PT) + LABEL_W_PT
        h = n * ROW_PITCH_PT
        ax = fig.add_axes([x0 / fig_w, (top_of_rows - h) / fig_h,
                           grid_w / fig_w, h / fig_h])
        ax.set_facecolor("none")          # transparent background
        for i, r in enumerate(s.itertuples()):
            y = n - 1 - i
            for k, db in enumerate(DBS):
                # 'Y' / NaN STRINGS — compare explicitly; bool(nan) is True.
                on = str(getattr(r, db)).strip().upper() == "Y"
                ax.add_patch(Rectangle(
                    (k * (cell_w + CELL_GAP_PT), y + 0.12), cell_w, 0.76,
                    facecolor=col if on else EMPTY, edgecolor="none"))
        ax.set_xlim(0, grid_w)
        ax.set_ylim(0, n)
        ax.set_xticks([k * (cell_w + CELL_GAP_PT) + cell_w / 2 for k in range(3)])
        ax.set_xticklabels(DBLAB, fontsize=FS_TICK)
        ax.set_yticks([n - 1 - i + 0.5 for i in range(n)])
        ax.set_yticklabels(s.label, fontsize=FS_GENE)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{DISPLAY[bait]}  n = {n}", fontsize=FS_BAIT, color=col,
                     fontweight="bold", loc="center", pad=BAIT_PAD_PT)
    fig.text(0.5, 1 - (TITLE_PAD_PT * 0.55) / fig_h, TITLE, ha="center",
             va="center", fontsize=FS_TITLE)
    return fig, fig_w, fig_h


def audit(fig, path_svg, label, fig_w, fig_h):
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
    body = open(path_svg).read()
    w = float(re.search(r'width="([\d.]+)pt"', body[:2000]).group(1))
    h = float(re.search(r'height="([\d.]+)pt"', body[:2000]).group(1))
    white_bg = "#ffffff" in body.split("<g")[0].replace(" ", "").lower()
    live_text = "<text" in body
    png = os.path.join(OUT, NAME + ".png")
    clip = "n/a"
    if os.path.exists(png):
        a = np.array(Image.open(png).convert("RGB")).astype(int).mean(axis=2)
        ink = a < 250
        clip = f"col0={ink[:, 0].sum()} row0={ink[0, :].sum()}"
    print(f"  {label:9s} page {w:6.1f} x {h:5.1f} pt (target {fig_w:.1f} x "
          f"{fig_h:.1f}) | white bg {white_bg} | live text {live_text} | "
          f"collisions {len(ov)} | edge ink {clip}")
    if ov:
        print("      COLLISIONS:", ov[:5])
    assert abs(w - fig_w) < 0.5 and abs(h - fig_h) < 0.5, (w, h)
    assert not white_bg, "SVG carries a white background rectangle"
    return len(ov)


def main():
    os.makedirs(OUT, exist_ok=True)
    data = load()
    print("inputs asserted OK: 16/25/16 rows, 132/171 filled cells")
    print(f"page {PAGE_W_PT:.1f} pt wide, matching panel A's ink span "
          "to within 0.2 pt\n")

    bad = 0
    for kind, fonttype in (("outlined", "path"), ("editable", "none")):
        plt.rcParams.update(BASE_RC)
        plt.rcParams["svg.fonttype"] = fonttype
        fig, fw, fh = build(data)
        box = Bbox([[0, 0], fig.get_size_inches()])
        svg = os.path.join(OUT, f"{NAME}_{kind}.svg")
        fig.savefig(svg, bbox_inches=box, transparent=True, facecolor="none",
                    edgecolor="none")
        if kind == "outlined":       # one opaque PNG, for screen checking only
            fig.savefig(os.path.join(OUT, f"{NAME}.png"), dpi=300,
                        bbox_inches=box, transparent=False, facecolor="white",
                        edgecolor="none")
        bad += audit(fig, svg, kind, fw, fh)
        plt.close(fig)

    print("\n⛔ SVG md5 IS NOT A VALID IDENTITY CHECK — matplotlib stamps a "
          "<dc:date>\n   into <metadata> AND mints random clip/glyph ids per "
          "build. Use byte size\n   plus the canonical md5 below (metadata "
          "stripped, generated ids normalised):")
    # ⛔ This MUST reproduce the shell recipe in README_ADOPTED.md byte for byte,
    # or the two disagree and neither can be trusted. sed DELETES WHOLE LINES from
    # <metadata> to </metadata>; an inline re.sub of the tags leaves the
    # surrounding newlines behind and yields a DIFFERENT hash (measured: 1f3baefc
    # inline vs d2b1af3b line-wise for the same file). Line-wise deletion below.
    #   sed '/<metadata>/,/<\/metadata>/d' f.svg \
    #     | sed -E 's/(id|href)="#?[pm][0-9a-f]{10}"/\1="#GENID"/g;
    #               s/url\(#[pm][0-9a-f]{10}\)/url(#GENID)/g' | md5
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
