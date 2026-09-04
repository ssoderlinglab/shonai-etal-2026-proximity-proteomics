"""
Supplemental Figure 7, panel D: collision-free label placement, the marker key handler and the marker state rule.

Input
    none; imported by figures/suppl_07/scripts/panels_BCD.py
Output
    none
Run
    python3 figures/suppl_07/scripts/panels_BCD.py
Dependencies
    python >= 3.9, numpy, matplotlib

place_labels walks a fixed ring by angle grid of candidate positions and takes the first that overlaps no placed label, no marker disc and lies inside the axes. It is deterministic; the only ordering dependence is the caller's own sort. It returns the number of labels that found no free slot, so the caller fails rather than shipping an unreadable panel.
"""
import matplotlib as mpl
import numpy as np
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D

LABEL = {"Synapse": "Synapse", "Ribo": "Ribosome", "Mito": "Mitochondria",
         "M4": "M4 ER/Golgi", "M6": "M6 Clathrin", "M9": "M9 Proteasome",
         "M14": "M14 RNA-binding"}
# Label search geometry: rings of increasing radius (points) x 16 angles.
# Small radii first, so a label stays adjacent to its marker when it can.
RADII = [12, 17, 23, 30, 38, 47, 57, 68]
N_ANGLES = 20
# EVERY label gets a leader. A distance threshold left the crowded right-hand
# cluster of the LRRK2 panel with three unleadered labels sitting between two
# overlapping markers - geometrically clean, but a reader cannot tell which
# label names which point. figure-style 6.9 requires the endpoint to land within
# one marker radius of the point it names, so the leader is drawn always.
MARKER_R_PT = 7.2         # half the outer 12.0-pt hairline, plus clearance


def place_labels(ax, items, fig, fontsize=6):
    """Collision-free radial label placement with leader lines.

    `items` is [(x, y, text, colour, weight), ...] in DATA coordinates. For each
    point, walks RADII x N_ANGLES candidate positions and keeps the first whose
    rendered box overlaps no already-placed label box, no plotted marker, and
    lies inside the axes. Deterministic - no randomness; the only ordering
    dependence is `items` itself, which the caller sorts.

    Returns the number of labels that had to fall back (no free slot found), so
    the caller can fail loudly rather than ship an unreadable panel.
    """
    r = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(r)
    pts = [ax.transData.transform((x, y)) for x, y, *_ in items]
    # Marker DISCS, not centres. A centre-only test passed the 9.1 bbox check
    # while two labels visibly ran across markers in the rendered PNG - the
    # markers are 7.5 pt across, so a label can clear the centre and still sit
    # on the glyph. Each marker becomes a box of half-width MARKER_R_PT.
    rpx = MARKER_R_PT * fig.dpi / 72.0
    marker_boxes = [mpl.transforms.Bbox.from_bounds(px - rpx, py - rpx,
                                                    2 * rpx, 2 * rpx)
                    for px, py in pts]
    angles = [2 * np.pi * k / N_ANGLES for k in range(N_ANGLES)]
    placed, n_fallback = [], 0

    for (x, y, text, colour, weight), pt in zip(items, pts):
        best = None
        for rad in RADII:
            for a in angles:
                dx, dy = rad * np.cos(a), rad * np.sin(a)
                t = ax.annotate(text, (x, y), textcoords="offset points",
                                xytext=(dx, dy), fontsize=fontsize, color=colour,
                                fontweight=weight, zorder=6,
                                ha="left" if dx > 1 else "right" if dx < -1 else "center",
                                va="bottom" if dy > 1 else "top" if dy < -1 else "center")
                bb = t.get_window_extent(r).expanded(1.10, 1.30)
                inside = (bb.x0 >= ax_box.x0 and bb.x1 <= ax_box.x1
                          and bb.y0 >= ax_box.y0 and bb.y1 <= ax_box.y1)
                if inside and not any(bb.overlaps(b) for b in placed) \
                        and not any(bb.overlaps(mb) for mb in marker_boxes):
                    best = (bb, rad, dx, dy)
                    break
                t.remove()
            if best:
                break
        if best is None:
            n_fallback += 1
            t = ax.annotate(text, (x, y), textcoords="offset points",
                            xytext=(RADII[-1], 0), fontsize=fontsize,
                            color=colour, fontweight=weight, zorder=6,
                            ha="left", va="center")
            best = (t.get_window_extent(r).expanded(1.10, 1.30), RADII[-1],
                    RADII[-1], 0)
        bb, rad, dx, dy = best
        if True:
            # leader from just outside the marker toward the label anchor
            norm = np.hypot(dx, dy)
            ux, uy = dx / norm, dy / norm
            trans = ax.transData.transform((x, y))
            p0 = trans + np.array([ux * MARKER_R_PT, uy * MARKER_R_PT]) * fig.dpi / 72.0
            p1 = trans + np.array([dx, dy]) * fig.dpi / 72.0 * 0.88
            (a0, b0) = ax.transData.inverted().transform(p0)
            (a1, b1) = ax.transData.inverted().transform(p1)
            ax.plot([a0, a1], [b0, b1], color=colour, lw=0.5, alpha=0.8,
                    zorder=3, solid_capstyle="butt")
        placed.append(bb)
    return n_fallback


class ConcentricHandler(HandlerBase):
    """Legend handler drawing rings and a disc CONCENTRICALLY at one position.

    ⛔ matplotlib's HandlerTuple cannot do this: with ndivide=None it still lays
    the artists out side by side, which rendered the "both" key as a ring next
    to a disc - and a lone ring is visually identical to the open "no shift"
    glyph. This handler stacks every marker size on the same centre, which is
    what the data markers actually look like.
    """

    def __init__(self, ring_sizes, disc_size, colour, ring_lw, **kw):
        self.ring_sizes, self.disc_size = ring_sizes, disc_size
        self.colour, self.ring_lw = colour, ring_lw
        super().__init__(**kw)

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width,
                       height, fontsize, trans):
        cx, cy = width / 2 - xdescent, height / 2 - ydescent
        out = []
        for ms in self.ring_sizes:
            out.append(Line2D([cx], [cy], marker="o", ls="", ms=ms, mfc="none",
                              mec=self.colour, mew=self.ring_lw,
                              transform=trans))
        out.append(Line2D([cx], [cy], marker="o", ls="", ms=self.disc_size,
                          mfc=self.colour, mec=self.colour, transform=trans))
        return out


def tint(hex_colour, frac=0.42):
    """Blend a hue toward white. Keeps the model's identity while lowering
    visual weight, so the three states form one saturation ramp per panel."""
    r, g, b = mpl.colors.to_rgb(hex_colour)
    return (r + (1 - r) * (1 - frac), g + (1 - g) * (1 - frac),
            b + (1 - b) * (1 - frac))


def state_of(sig_cortex, sig_cofrac):
    """One of the three states that actually occur in these data.

    ⛔ cortex-only is NOT a drawable state: it is empty in the current table
    (0 of 21 cells) and therefore has no key. If the statistics change and a
    cortex-only cell appears, this raises instead of silently rendering it as
    something else."""
    if sig_cortex and sig_cofrac:
        return "both"
    if sig_cofrac:
        return "da_only"
    if sig_cortex:
        raise AssertionError(
            "a cortex-only cell appeared; the 3-state key no longer covers the "
            "data - add the state and its legend entry before rendering")
    return "neither"


DISC_MS = 7.5          # the base disc, shared by all three states
RING_MS = (10.2, 12.0)  # concentric hairlines marking "both"; clear of the disc
RING_LW = 0.6
