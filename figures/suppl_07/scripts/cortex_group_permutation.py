"""
Supplemental Figure 7, panels B and D: build the cortical functional-group permutation table both panels read.

Input
    data/cortex/deqms_all_contrasts_long.tsv
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/permutation/permutation_all21.csv
Output
    figures/suppl_07/output/cortex_permutation_all21.csv
    figures/suppl_07/output/cortex_vs_cofrac_21cells.csv
    figures/suppl_07/output/cortex_permutation_summary.json
Run
    python3 figures/suppl_07/scripts/cortex_group_permutation.py
Dependencies
    python >= 3.9, numpy, pandas, matplotlib, statsmodels

Size-matched permutation WITHOUT replacement against the dataset's own proteome mean, 200,000 draws, seed 20260811, two-sided, Benjamini-Hochberg across the 21 group by model cells. Runtime is a few minutes.

The seven functional groups are defined on the LRRK2 partition: Synapse (M7, M8, M11, M19), Ribosome (M1, M10), Mitochondria (M2, M5), M4, M6, M9, M14. Contaminant genes are removed from both the groups and the null pool.

Run before figures/suppl_07/scripts/panels_BCD.py.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(os.path.dirname(FIGURE_DIR))
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUT = os.path.join(FIGURE_DIR, "output")
os.makedirs(OUT, exist_ok=True)

sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
from contamination_filter import CONTAMINATION_GENES as EXCLUDE

SEED = 20260811
N_PERM = 200_000
MODELS = ["LRRK2", "VPS35", "SNCA"]

GROUPS = [
    ("Synapse", "Synapse (M7,8,11,19)", [7, 8, 11, 19]),
    ("Ribo",    "Ribosome (M1,10)",     [1, 10]),
    ("Mito",    "Mitochondria (M2,5)",  [2, 5]),
    ("M4",      "M4 (ER/Golgi)",        [4]),
    ("M6",      "M6 (Clath/tRNA)",      [6]),
    ("M9",      "M9 (Proteasome)",      [9]),
    ("M14",     "M14 (RNA-bind/ND)",    [14]),
]

PATH_PARTITION = os.path.join(DATA_DIR, "ibiocofrac", "modules",
                              "Partition_Gene_combined_3PD_models.csv")


def load():
    part = pd.read_csv(PATH_PARTITION)
    lp = part[part.Experiment == "LRRK2"][["UniProt", "Symbol", "Membership"]].copy()
    assert len(lp) == 1957, f"LRRK2 partition is {len(lp)} rows, expected 1957"
    assert lp.Membership.nunique() == 42, "expected 42 LRRK2 modules"
    assert int((lp.Membership == 11).sum()) == 77, "M11 is not 77 in the partition"
    lp = lp[~lp.Symbol.isin(EXCLUDE)]

    gsets = {k: set(lp[lp.Membership.isin(mods)].UniProt) for k, _, mods in GROUPS}

    L = pd.read_csv(os.path.join(DATA_DIR, "cortex",
                                 "deqms_all_contrasts_long.tsv"), sep="\t")
    frames = {}
    for m in MODELS:
        d = L[L.contrast == f"{m}_cortex_vs_WT"][["Accession", "Gene", "log2FC"]].copy()
        d = d[~d.Gene.isin(EXCLUDE)].dropna(subset=["log2FC"])
        frames[m] = d
    return frames, gsets


def permutation_all_cells(frames, gsets, n_perm=N_PERM, seed=SEED):
    """Size-matched null, WITHOUT replacement, two-sided against the dataset's
    own proteome mean. The same null the dopamine neuron panels use."""
    rng = np.random.default_rng(seed)
    rows = []
    for m in MODELS:
        pv = frames[m].log2FC.to_numpy()
        pool_mean = pv.mean()
        accs = frames[m].Accession
        for gk, glab, _ in GROUPS:
            v = frames[m].loc[accs.isin(gsets[gk]), "log2FC"].to_numpy()
            n = len(v)
            obs = v.mean()
            null = np.empty(n_perm)
            for i in range(n_perm):
                null[i] = pv[rng.choice(len(pv), n, replace=False)].mean()
            null_sd = null.std(ddof=1)
            exceed = int((np.abs(null - pool_mean) >= abs(obs - pool_mean)).sum())
            rows.append({
                "group": gk, "group_label": glab, "model": m, "n": n,
                "n_defined": len(gsets[gk]),
                "mean_raw": obs, "pool_mean": pool_mean, "delta": obs - pool_mean,
                "null_p2_5": float(np.percentile(null, 2.5)),
                "null_median": float(np.median(null)),
                "null_p97_5": float(np.percentile(null, 97.5)),
                "null_sd": null_sd, "z": (obs - pool_mean) / null_sd,
                "percentile_of_obs": 100 * float((null <= obs).mean()),
                "n_exceedances": exceed,
                "p_perm": exceed / n_perm if exceed > 0 else 1.0 / n_perm,
                "p_at_resolution_bound": exceed == 0,
            })
    P = pd.DataFrame(rows)
    P["q_perm_BH"] = multipletests(P.p_perm.values, method="fdr_bh")[1]
    return P


def main():
    frames, gsets = load()
    print("cortex pool per model, contaminants removed:",
          {m: len(frames[m]) for m in MODELS})
    print(f"running size-matched permutation, {N_PERM:,} draws x 21 cells")
    P = permutation_all_cells(frames, gsets)
    P.to_csv(os.path.join(OUT, "cortex_permutation_all21.csv"), index=False)

    # join to the DA-neuron iBioCoFrac result. Suffixes name the PREPARATION,
    # because the two are different preparations and different mouse models.
    C = pd.read_csv(os.path.join(DATA_DIR, "ibiocofrac", "permutation",
                                 "permutation_all21.csv"))
    keep = ["group", "group_label", "model", "n", "mean_raw", "pool_mean",
            "delta", "z", "p_perm", "q_perm_BH"]
    j = (P[keep].rename(columns={c: f"{c}_cortex" for c in keep[3:]})
         .merge(C[keep].rename(columns={c: f"{c}_cofrac" for c in keep[3:]}),
                on=["group", "group_label", "model"], how="outer"))
    assert len(j) == 21, f"join produced {len(j)} cells, expected 21"
    j["same_sign"] = np.sign(j.delta_cortex) == np.sign(j.delta_cofrac)
    j["sig_cortex"] = j.q_perm_BH_cortex < 0.05
    j["sig_cofrac"] = j.q_perm_BH_cofrac < 0.05
    j = j.sort_values(["model", "group"])
    j.to_csv(os.path.join(OUT, "cortex_vs_cofrac_21cells.csv"), index=False)

    summary = dict(
        seed=SEED, n_perm=N_PERM,
        method=("size-matched permutation WITHOUT replacement, two-sided against "
                "the dataset's own proteome mean, BH across 21 cells"),
        cortex_log2FC="DEqMS contrast coefficient, 3 vs 3 mice, no fractions",
        cortex_pool={m: int(len(frames[m])) for m in MODELS},
        n_cells_sig_cortex=int(j.sig_cortex.sum()),
        n_cells_sig_cofrac=int(j.sig_cofrac.sum()),
        n_cells_same_sign=int(j.same_sign.sum()),
        n_cells_sig_both=int((j.sig_cortex & j.sig_cofrac).sum()),
        synapse_cortex={r.model: dict(delta=float(r.delta_cortex),
                                      q=float(r.q_perm_BH_cortex), n=int(r.n_cortex))
                        for r in j[j.group == "Synapse"].itertuples()},
        synapse_cofrac={r.model: dict(delta=float(r.delta_cofrac),
                                      q=float(r.q_perm_BH_cofrac), n=int(r.n_cofrac))
                        for r in j[j.group == "Synapse"].itertuples()},
    )
    with open(os.path.join(OUT, "cortex_permutation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print("\n== CORTEX, 21 cells ==")
    print(P[["model", "group", "n", "mean_raw", "pool_mean", "delta", "z",
             "q_perm_BH"]].round(4).to_string(index=False))
    print("\n== CORTEX vs DA-neuron iBioCoFrac ==")
    print(j[["model", "group", "n_cortex", "delta_cortex", "q_perm_BH_cortex",
             "n_cofrac", "delta_cofrac", "q_perm_BH_cofrac", "same_sign"]]
          .round(4).to_string(index=False))
    print(f"\nsig cortex {summary['n_cells_sig_cortex']}/21 | "
          f"sig cofrac {summary['n_cells_sig_cofrac']}/21 | "
          f"same sign {summary['n_cells_same_sign']}/21 | "
          f"sig in both {summary['n_cells_sig_both']}/21")


if __name__ == "__main__":
    main()
