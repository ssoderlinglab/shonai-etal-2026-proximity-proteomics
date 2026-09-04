"""
Figure 6, panel D: prepare the module-annotated fraction tables the mixed model and the permutation read.

Input
    data/ibiocofrac/tidy/iBioCoFrac_Tidy_{LRRK2,SNCA,VPS35}.csv
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
Output
    figures/main_06/output/tidy_{LRRK2,SNCA,VPS35}_module_annotated.csv
Run
    python3 figures/main_06/scripts/module_profile_prepare.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

First of three statistics steps for panel D: prepare, then the mixed model (R), then the permutation and bootstrap.
"""
import os
import numpy as np
import pandas as pd

# --- paths (relative to the manuscript root) -------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, "data", "ibiocofrac")

TIDY_DIR = os.path.join(DATA_DIR, "tidy")
PARTITION = os.path.join(
    DATA_DIR, "modules", "Partition_Gene_combined_3PD_models.csv")
OUT = os.path.join(FIGURE_DIR, "output")

# tidy filename stem -> Experiment label used in the partition file
MODELS = {"LRRK2": "LRRK2", "SNCA": "SNCA", "VPS35": "VPS35"}


def main():
    os.makedirs(OUT, exist_ok=True)
    part = pd.read_csv(PARTITION)

    for stem, exp in MODELS.items():
        tidy = pd.read_csv(os.path.join(TIDY_DIR, f"iBioCoFrac_Tidy_{stem}.csv"))
        memb = part.loc[part.Experiment == exp, ["UniProt", "Membership"]]
        d = tidy.merge(memb, left_on="Accession", right_on="UniProt", how="inner")

        # drop non-positive intensities (cannot be log-transformed) and
        # module 0 = unassigned proteins
        d = d[(d.Norm_Intensity > 0) & (d.Membership != 0)].copy()

        # Courtland's rel_Intensity: each protein scaled to sum 1 over its own samples
        d["rel_Intensity"] = d.groupby("Accession").Norm_Intensity.transform(
            lambda x: x / x.sum()
        )
        cols = ["Accession", "Membership", "Genotype", "Fraction",
                "Replicate", "Norm_Intensity", "rel_Intensity"]
        d[cols].to_csv(os.path.join(OUT, f"tidy_{stem}_module_annotated.csv"), index=False)

        # protein-level total abundance per sample (sum over fractions), log2
        tot = (d.groupby(["Accession", "Membership", "Genotype", "Replicate"])
                 .Norm_Intensity.sum().reset_index())
        wide = tot.pivot_table(index=["Accession", "Membership"],
                               columns=["Genotype", "Replicate"],
                               values="Norm_Intensity")
        wide.columns = [f"{g}_{r}" for g, r in wide.columns]
        wide = wide.dropna()
        wide = wide[(wide > 0).all(axis=1)]
        np.log2(wide).reset_index().to_csv(
            os.path.join(OUT, f"protlevel_{stem}.csv"), index=False)

        n_mod = d.Membership.nunique()
        print(f"{stem}: {len(d)} rows, {d.Accession.nunique()} proteins, {n_mod} modules")

    # assert the cross-model non-correspondence that motivates the warning above
    s = {e: set(part.loc[(part.Experiment == e) & (part.Membership == 11), "UniProt"])
         for e in MODELS.values()}
    assert len(s["LRRK2"] & s["SNCA"]) == 0, "unexpected overlap LRRK2/SNCA module 11"
    assert len(s["LRRK2"] & s["VPS35"]) == 0, "unexpected overlap LRRK2/VPS35 module 11"
    print("\nverified: module 11 of each partition shares 0 proteins across models")
    print(f"  LRRK2 module 11 = {len(s['LRRK2'])} proteins  <- the manuscript's M11")


if __name__ == "__main__":
    main()
