"""
Figure 6, panel D: M11 protein abundance across sedimentation fractions, wild-type against transgenic.

Input
    data/ibiocofrac/fractions/LOPIT_LRRK2_young_Transformed_Rel-int_mixture.csv
    data/ibiocofrac/modules/LRRK2-TMT-Module-Results_relintensity.xlsx
    figures/main_06/output/permutation_bootstrap.json
Output
    figures/main_06/output/panel_D_m11_fraction_profile.{png,svg}
Run
    python3 figures/main_06/scripts/m11_fraction_profile.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib, openpyxl

Run module_profile_prepare.py, module_lmm.R and module_permutation_bootstrap.py first.

M11 membership is asserted at 77 proteins.
"""

import json
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, "data", "ibiocofrac")
# The mixture-normalised rel-int CSV lives in Box, NOT under the project root,
# so PD_PROJECT_ROOT cannot reach it. Override with PD_BOX_FRACTION_DIR.
OUT = os.path.join(FIGURE_DIR, "output")

GENE = "LRRK2"
MODULE = 11
WT_LABEL = "WildType"
MUT_LABEL = "Transgenic"
MIXTURES = ["M1", "M2", "M3"]
FRACTIONS = ["F4", "F5", "F6", "F7", "F8", "F9", "F10"]

PATH_TRANSFORMED = os.path.join(
    DATA_DIR, "fractions",
    f"LOPIT_{GENE}_young_Transformed_Rel-int_mixture.csv")
PATH_MODULE_RESULTS = os.path.join(
    DATA_DIR, "modules", f"{GENE}-TMT-Module-Results_relintensity.xlsx")
PATH_PERM = os.path.join(OUT, "permutation_bootstrap.json")

EXPECTED_M11 = 77

# the notebook's own plot_config, verbatim
PLOT_CONFIG = {
    "xlabel_size": 14, "ylabel_size": 14,
    "xticklabel_size": 16, "yticklabel_size": 16,
    "title_size": 24, "legend_size": 16,
    "wt_color": "green", "mut_color": "red",
    "alpha_individual": 0.05, "alpha_average": 1.0,
}


def load():
    # fail with the variable name to set, not with a bare FileNotFoundError
    if not os.path.exists(PATH_TRANSFORMED):
        raise FileNotFoundError(
            f"mixture-normalised rel-int CSV not found:\n  {PATH_TRANSFORMED}\n"
            "This file lives in Box, outside the project root. Set "
            "PD_BOX_FRACTION_DIR to the '20250227 Fraction plot mixture-norm' "
            "directory on this machine.")
    if not os.path.exists(PATH_MODULE_RESULTS):
        raise FileNotFoundError(
            f"module results workbook not found:\n  {PATH_MODULE_RESULTS}\n"
            "Set PD_PROJECT_ROOT to the Parkinson_Proteomics_Manuscript root.")
    tr = pd.read_csv(PATH_TRANSFORMED)
    part = pd.read_excel(PATH_MODULE_RESULTS, sheet_name="Partition",
                         engine="openpyxl")
    mres = pd.read_excel(PATH_MODULE_RESULTS, sheet_name="Module Results",
                         engine="openpyxl")
    return tr, part, mres


def load_permutation_bound():
    """Read the permutation p-bound from the permutation output. Never recompute it here.

    ⚠ `permutation/p_value` in that file is null and that is NOT missing data:
    0 of 200,000 size-matched draws exceeded the observed mean, so no point
    estimate exists and the reportable quantity is `p_upper_bound`.
    Reading `p_value` alone looks like a gap in the data.
    """
    with open(PATH_PERM) as fh:
        perm = json.load(fh)["permutation"]
    assert perm["exceedances"] == 0, (
        f"now reports {perm['exceedances']} exceedances; the panel's "
        "'p < bound' wording assumes zero and must be revisited")
    assert perm["n_proteins"] == EXPECTED_M11, (
        f"permutation ran on {perm['n_proteins']} proteins, not "
        f"{EXPECTED_M11} -- it does not describe this panel's module")
    return float(perm["p_upper_bound"])


def module_stats(mres):
    row = mres[mres["Module"] == f"M{MODULE}"]
    assert len(row) == 1, f"expected one M{MODULE} row, got {len(row)}"
    n_prots = int(row["nProts"].values[0])
    padjust = float(row["Padjust"].values[0])
    # *** THE POINT OF THIS SCRIPT ***
    assert n_prots == 77, (
        f"module results says nProts={n_prots}; the panel's 76 is NOT in the "
        "data and must not be reintroduced")
    assert abs(padjust - 5.126402e-16) / 5.126402e-16 < 1e-6, (
        f"Padjust {padjust:.4e} does not match the published 5.13e-16")
    return n_prots, padjust


def build_curves(tr, part):
    """Merge partition onto the transformed table, keep module MODULE only.

    The notebook merges on the protein id, then averages `rel_int_mixture` per
    protein across the three mixtures within each genotype -- so an individual
    faint line is ONE PROTEIN's mixture-averaged profile, not one sample.
    """
    # *** TRAP: BioFraction carries TRAILING SPACES in the source CSV. ***
    # The raw values are 'F4 ', 'F5 ', ... 'F9 ' but 'F10' (no space), so a
    # reindex on clean labels silently yields all-NaN for six of seven
    # fractions and a plot with one point. Strip before anything else.
    tr = tr.copy()
    tr["BioFraction"] = tr["BioFraction"].astype(str).str.strip()
    assert set(FRACTIONS) <= set(tr["BioFraction"].unique()), \
        f"fraction labels do not match after stripping: {sorted(tr['BioFraction'].unique())}"

    merged = tr.merge(part[["UniProt", "Symbol", "Membership"]],
                      left_on="Protein", right_on="UniProt", how="inner")
    mod = merged[merged["Membership"].astype(str).str.replace("M", "")
                 == str(MODULE)].copy()
    assert mod["Protein"].nunique() == 77, (
        f"pulled {mod['Protein'].nunique()} proteins for M{MODULE}, expected 77")

    curves = {}
    for label in (WT_LABEL, MUT_LABEL):
        sub = mod[mod["GenoType"] == label]
        # protein x fraction, averaged over the three mixtures
        wide = sub.pivot_table(index="Protein", columns="BioFraction",
                               values="rel_int_mixture", aggfunc="mean")
        curves[label] = wide.reindex(columns=FRACTIONS)
    return mod, curves


def draw(curves, n_prots, padjust, label_lines=None, tag="v01"):
    """Reproduce the notebook figure. `label_lines` overrides the annotation."""
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, colour in ((WT_LABEL, PLOT_CONFIG["wt_color"]),
                          (MUT_LABEL, PLOT_CONFIG["mut_color"])):
        w = curves[label]
        for _, row in w.iterrows():
            ax.plot(FRACTIONS, row.values, color=colour,
                    alpha=PLOT_CONFIG["alpha_individual"])
        ax.plot(FRACTIONS, w.mean(), color=colour,
                alpha=PLOT_CONFIG["alpha_average"], linewidth=2)

    ax.set_title(f"{GENE} Module M{MODULE}",
                 fontsize=PLOT_CONFIG["title_size"])
    ax.tick_params(axis="x", labelsize=PLOT_CONFIG["xticklabel_size"])
    ax.tick_params(axis="y", labelsize=PLOT_CONFIG["yticklabel_size"])
    ax.spines[["top", "right"]].set_visible(False)

    if label_lines is None:
        label_lines = [(f"Number of Proteins: {n_prots}", "black"),
                       (f"Padjust = {padjust:.2e}", "blue")]
    y = 0.95
    for txt, colour in label_lines:
        ax.text(0.95, y, txt, ha="right", va="top", transform=ax.transAxes,
                fontsize=PLOT_CONFIG["legend_size"], color=colour)
        y -= 0.05
    fig.tight_layout()
    return fig


def main():
    os.makedirs(OUT, exist_ok=True)
    tr, part, mres = load()
    n_prots, padjust = module_stats(mres)
    mod, curves = build_curves(tr, part)

    # ---- the published panel ---------------------------------------
    # Three annotation lines: n, the mixed-model Padjust and the size-matched
    # permutation. Both values are read from their sources and are not recomputed here.
    perm_bound = load_permutation_bound()
    fig3 = draw(curves, n_prots, padjust, label_lines=[
        (f"Number of Proteins: {n_prots}", "black"),
        (f"Padjust = {padjust:.2e}", "blue"),
        (f"Permutation p < {perm_bound:g}", "blue")])
    for fmt in ("png", "svg"):
        fig3.savefig(os.path.join(OUT, f"panel_D_m11_fraction_profile.{fmt}"),
                     dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig3)


    wt, mt = curves[WT_LABEL].mean(), curves[MUT_LABEL].mean()
    pd.DataFrame({"fraction": FRACTIONS, "WT_mean": wt.values,
                  "Tg_mean": mt.values,
                  "delta": (mt - wt).values}).to_csv(
        os.path.join(OUT, "panel_D_curve_values.csv"), index=False)

    print(f"module results : nProts={n_prots}  Padjust={padjust:.4e}")
    print(f"proteins pulled: {mod['Protein'].nunique()}")
    print(f"WT mean  : {np.round(wt.values, 4).tolist()}")
    print(f"Tg mean  : {np.round(mt.values, 4).tolist()}")
    print(f"WT peak  : {FRACTIONS[int(np.argmax(wt.values))]}  "
          f"Tg peak: {FRACTIONS[int(np.argmax(mt.values))]}")
    print(f"Tg below WT in all fractions: {bool((mt.values < wt.values).all())}")


if __name__ == "__main__":
    main()
