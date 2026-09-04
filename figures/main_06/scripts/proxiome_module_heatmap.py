"""
Figure 6, panel A: overlap between the proximity proteomes and the co-fractionation modules.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibioid/proxiome_sets/*.csv
Output
    figures/main_06/output/panel_A_proxiome_module_heatmap.{png,svg}
    figures/main_06/output/panel_A_proxiome_module_enrichment.csv
Run
    python3 figures/main_06/scripts/proxiome_module_heatmap.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib
"""

import os

import numpy as np
import pandas as pd
from scipy import stats as sp_stats
from scipy.stats import false_discovery_control

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -- Paths ------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
OUT = os.path.join(FIGURE_DIR, "output")

PATH_PARTITION = os.path.join(
    REPO_ROOT, "data/ibiocofrac/modules",
    "Partition_Gene_combined_3PD_models.csv")
PATH_PROXIOME_DIR = os.path.join(REPO_ROOT, "data/ibioid/proxiome_sets")

PROXIOME_FILES = {
    "aSyn_WT": "aSyn_WT.csv",
    "aSyn_A30P": "aSyn_A30P.csv",
    "Vps35_WT": "Vps35_WT.csv",
    "Vps35_DN": "Vps35_D620N.csv",
    "Lrrk2_WT": "Lrrk2_WT.csv",
    "Lrrk2_GS": "Lrrk2_G2019S.csv",
}

# -- Display constants (unchanged from the published generator) -------
HEATMAP_ORDER = [7, 8, 11, 19, 1, 10, 2, 5, 4, 6, 9, 14]
SYNAPSE_MODULES = {7, 8, 11, 19}
MODULE_LABELS = {
    1: "Ribo (cyto)", 2: "Mito (matrix)", 4: "ER/Golgi", 5: "Mito (memb)",
    6: "Clath/tRNA", 7: "SV cycling", 8: "Post/cytosk", 9: "Proteasome",
    10: "Ribo (ER)", 11: "Pre/SV", 14: "RNA-bind/ND", 19: "CME",
}
COL_SYNAPSE = "#E15759"
PROXIOME_ORDER = ["aSyn_WT", "aSyn_A30P", "Vps35_WT", "Vps35_DN",
                  "Lrrk2_WT", "Lrrk2_GS", "Total_Conv_WT", "Total_Conv_MUT"]
PROXIOME_DISPLAY = {
    "aSyn_WT": "WT\naSyn", "aSyn_A30P": "A30P\naSyn",
    "Vps35_WT": "WT\nVps35", "Vps35_DN": "D620N\nVps35",
    "Lrrk2_WT": "WT\nLrrk2", "Lrrk2_GS": "G2019S\nLrrk2",
    "Total_Conv_WT": "Total Conv\nWT", "Total_Conv_MUT": "Total Conv\nMUT",
}
DISPLAY_CAP = 20      # ceil(max -log10(q)) == 20; nothing saturates
EXPECTED_N_TOTAL = 1957
EXPECTED_M11 = 77


def build_sets():
    """Proxiome sets on ACCESSION, plus the two 3-way convergent sets."""
    sets = {}
    for key, fname in PROXIOME_FILES.items():
        df = pd.read_csv(os.path.join(PATH_PROXIOME_DIR, fname))
        sets[key] = set(df["Accession"])
    sets["Total_Conv_WT"] = (sets["aSyn_WT"] & sets["Vps35_WT"]
                             & sets["Lrrk2_WT"])
    sets["Total_Conv_MUT"] = (sets["aSyn_A30P"] & sets["Vps35_DN"]
                              & sets["Lrrk2_GS"])
    # export values -- the accession-vs-symbol trap shows up right here
    assert len(sets["Total_Conv_WT"]) == 86, (
        f"3-way WT is {len(sets['Total_Conv_WT'])}, expected 86 "
        "(87 means the intersection ran on gene symbol)")
    assert len(sets["Total_Conv_MUT"]) == 72, (
        f"3-way MUT is {len(sets['Total_Conv_MUT'])}, expected 72")
    return sets


def build_background():
    """LRRK2 partition, the contamination exclusion NOT applied -> universe 1957, M11 = 77."""
    part = pd.read_csv(PATH_PARTITION)
    lp = part[part["Experiment"] == "LRRK2"][
        ["UniProt", "Symbol", "Membership"]].copy()
    lp.columns = ["Accession", "Gene", "Module"]

    assert lp["Accession"].is_unique, "partition accessions are not unique"
    n_total = lp["Accession"].nunique()
    assert n_total == EXPECTED_N_TOTAL, (
        f"universe is {n_total}, expected {EXPECTED_N_TOTAL} "
        "(1939 means the contamination exclusion was applied -- this panel must NOT apply it)")

    m11 = lp[lp["Module"] == 11]
    assert len(m11) == EXPECTED_M11, f"M11 is {len(m11)}, expected 77"
    # the two genes this whole question was about must BOTH be present
    for g in ("COX2", "Cldn11"):
        assert g in set(m11["Gene"]), f"{g} missing from M11 -- basis is wrong"
    return lp


def compute(lp, sets):
    n_mods, n_prox = len(HEATMAP_ORDER), len(PROXIOME_ORDER)
    background = set(lp["Accession"])
    n_total = len(background)

    k_mat = np.zeros((n_mods, n_prox), dtype=int)
    p_mat = np.ones((n_mods, n_prox))
    for mi, mod in enumerate(HEATMAP_ORDER):
        mod_acc = set(lp[lp["Module"] == mod]["Accession"])
        n_mod = len(mod_acc)
        for pi, pname in enumerate(PROXIOME_ORDER):
            in_bg = sets[pname] & background
            k = len(mod_acc & in_bg)
            k_mat[mi, pi] = k
            if k and in_bg and n_mod:
                p_mat[mi, pi] = sp_stats.hypergeom.sf(
                    k - 1, n_total, len(in_bg), n_mod)

    q_mat = false_discovery_control(
        p_mat.ravel(), method="bh").reshape(p_mat.shape)
    nlq = -np.log10(np.maximum(q_mat, 1e-300))

    assert nlq.max() <= DISPLAY_CAP, (
        f"max -log10(q) is {nlq.max():.2f}, above the cap {DISPLAY_CAP} -- "
        "the cap decision must be revisited before this panel ships")
    top = [HEATMAP_ORDER[int(np.argmax(nlq[:, i]))] for i in range(n_prox)]
    assert all(t == 11 for t in top), (
        f"M11 is no longer the top module in every column: {top}")
    return k_mat, p_mat, q_mat, nlq


def star_tier(qval):
    """3-tier significance star on the BH q-value. Adopted 2026-08-07 (Daichi).

    A single-star variant (* for q<0.05 only) was rendered alongside this and
    NOT adopted: it collapses M11's q ~ 1e-20 and M8 x Lrrk2_GS's q = 0.015 into
    the same mark. Tiers give M8 and M14 a readable bait dependence.
    """
    if qval < 0.001:
        return "***"
    if qval < 0.01:
        return "**"
    if qval < 0.05:
        return "*"
    return ""


def draw(k_mat, q_mat, nlq, sizes):
    fig, ax = plt.subplots(figsize=(9, 6.5))
    im = ax.imshow(nlq, cmap="Reds", aspect="auto", vmin=0, vmax=DISPLAY_CAP)

    # Cell text is the overlap count k; the star underneath it is the BH q tier.
    # They are DIFFERENT quantities -- k is not a significance measure (M6 x
    # aSyn_A30P has k=16 and is n.s.; M8 x Total_Conv_WT has k=7 and is
    # significant), which is exactly why the stars earn their place here.
    for mi in range(len(HEATMAP_ORDER)):
        for pi in range(len(PROXIOME_ORDER)):
            colour = "white" if nlq[mi, pi] > 12 else "black"
            star = star_tier(q_mat[mi, pi])
            if star:
                ax.text(pi, mi - 0.14, str(int(k_mat[mi, pi])),
                        ha="center", va="center", fontsize=7.5,
                        fontweight="bold", color=colour)
                ax.text(pi, mi + 0.27, star, ha="center", va="center",
                        fontsize=9, fontweight="bold", color=colour)
            else:
                ax.text(pi, mi, str(int(k_mat[mi, pi])),
                        ha="center", va="center", fontsize=7.5,
                        fontweight="bold", color=colour)

    ax.axhline(3.5, color="white", linewidth=2.5)
    ax.axvline(5.5, color="white", linewidth=2.5)

    ax.set_yticks(range(len(HEATMAP_ORDER)))
    ax.set_yticklabels(
        [f"M{m} {MODULE_LABELS[m]} (n={sizes[m]})" for m in HEATMAP_ORDER],
        fontsize=8.5)
    for i, mod in enumerate(HEATMAP_ORDER):
        if mod in SYNAPSE_MODULES:
            ax.get_yticklabels()[i].set_color(COL_SYNAPSE)
            ax.get_yticklabels()[i].set_fontweight("bold")

    ax.set_xticks(range(len(PROXIOME_ORDER)))
    ax.set_xticklabels([PROXIOME_DISPLAY[p] for p in PROXIOME_ORDER],
                       fontsize=8, ha="center")

    cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("-log10(q) hypergeometric, BH-corrected", fontsize=8)

    ax.set_title("Proxiome enrichment across co-fractionation modules",
                 fontsize=11, fontweight="bold")
    # Star key on the panel. The legend must still carry the full statement
    # (test, background, BH scope); this is the minimum so the panel is not
    # self-contradictory when read alone.
    ax.set_xlabel("*** q < 0.001    ** q < 0.01    * q < 0.05    "
                  "(number = overlap count)", fontsize=7.5, labelpad=8)
    fig.tight_layout()
    return fig


def main():
    os.makedirs(OUT, exist_ok=True)
    sets = build_sets()
    lp = build_background()
    sizes = {m: int((lp["Module"] == m).sum()) for m in HEATMAP_ORDER}
    k_mat, p_mat, q_mat, nlq = compute(lp, sets)

    fig = draw(k_mat, q_mat, nlq, sizes)
    for fmt in ("png", "svg"):
        fig.savefig(
            os.path.join(OUT, f"panel_A_proxiome_module_heatmap.{fmt}"),
            dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    rows = []
    for mi, mod in enumerate(HEATMAP_ORDER):
        for pi, pname in enumerate(PROXIOME_ORDER):
            rows.append(dict(
                module=f"M{mod}", label=MODULE_LABELS[mod],
                module_size=sizes[mod], proxiome=pname,
                overlap_k=int(k_mat[mi, pi]),
                p_value=float(p_mat[mi, pi]),
                q_value_BH=float(q_mat[mi, pi]),
                neg_log10_q=float(nlq[mi, pi]),
                sig_raw_p=bool(p_mat[mi, pi] < 0.05),
                sig_BH_q=bool(q_mat[mi, pi] < 0.05),
                star=star_tier(q_mat[mi, pi])))
    pd.DataFrame(rows).to_csv(
        os.path.join(OUT, "panel_A_proxiome_module_enrichment.csv"),
        index=False)

    i11 = HEATMAP_ORDER.index(11)
    print(f"universe N              : {lp['Accession'].nunique()}")
    print(f"M11 size                : {sizes[11]}  (COX2 + Cldn11 retained)")
    print(f"module sizes            : {[sizes[m] for m in HEATMAP_ORDER]}")
    print(f"cells sig raw p / BH q  : {int((p_mat < 0.05).sum())} / "
          f"{int((q_mat < 0.05).sum())}")
    print(f"M11 row k               : {k_mat[i11].tolist()}")
    print(f"M11 row -log10(q)       : {np.round(nlq[i11], 2).tolist()}")
    print(f"max -log10(q)           : {nlq.max():.2f}  (cap {DISPLAY_CAP}, "
          f"{int((nlq > DISPLAY_CAP).sum())} cells over)")

    n3 = int((q_mat < 0.001).sum())
    n2 = int(((q_mat >= 0.001) & (q_mat < 0.01)).sum())
    n1 = int(((q_mat >= 0.01) & (q_mat < 0.05)).sum())
    print(f"star tiers ***/**/*     : {n3}/{n2}/{n1}   n.s. {96 - n3 - n2 - n1}")
    # M11 must carry the top tier in every column, or the panel no longer says
    # what the figure claims
    i11 = HEATMAP_ORDER.index(11)
    assert all(q_mat[i11, pi] < 0.001 for pi in range(len(PROXIOME_ORDER))), \
        "M11 does not reach *** in every column -- headline claim weakened"


if __name__ == "__main__":
    main()
