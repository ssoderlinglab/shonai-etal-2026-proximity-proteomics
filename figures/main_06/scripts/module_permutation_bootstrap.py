"""
Figure 6, panel D: size-matched permutation and bootstrap for module M11.

Input
    figures/main_06/output/tidy_{LRRK2,SNCA,VPS35}_module_annotated.csv
Output
    figures/main_06/output/permutation_bootstrap.json
    figures/main_06/output/permutation_result.csv
Run
    python3 figures/main_06/scripts/module_permutation_bootstrap.py
Dependencies
    python >= 3.9, numpy, pandas, scipy, statsmodels, matplotlib

Third of three statistics steps for panel D.

200,000 permutation draws and 20,000 bootstrap resamples, seed 20260731.
"""
import os
import json
import numpy as np
import pandas as pd

SEED = 20260731
N_PERM = 200_000
N_BOOT = 20_000

HERE = os.path.dirname(os.path.abspath(__file__))
FIGURE_DIR = os.path.dirname(HERE)
DATA = os.path.join(FIGURE_DIR, "output")
OUT = os.path.join(FIGURE_DIR, "output")

# The manuscript's M11 lives in the LRRK2 partition. Module numbers do not
# correspond across models, so the focal test is defined on LRRK2 only.
FOCAL_MODEL, FOCAL_MODULE = "LRRK2", 11


def paired_log2fc(stem):
    """Mut-vs-WT log2FC per protein x fraction x replicate.

    The design is PAIRED: within each replicate the WT and Mut channels contain
    the identical protein set, so a within-replicate log2 ratio is defined and
    fraction effects largely cancel.
    """
    d = pd.read_csv(os.path.join(DATA, f"tidy_{stem}_module_annotated.csv"))
    w = d.pivot_table(index=["Accession", "Membership", "Fraction", "Replicate"],
                      columns="Genotype", values="Norm_Intensity").reset_index()
    w = w.dropna(subset=["WT", "Mut"])
    w = w[(w.WT > 0) & (w.Mut > 0)].copy()
    w["log2FC"] = np.log2(w.Mut / w.WT)
    return w


def variance_tables(stem):
    """Two variance decompositions.

    raw    - variance of a single protein's log2 intensity, partitioned across
             Fraction / Replicate / Genotype. This is the same quantity shown in
             Courtland et al. 2021 Figure 3-figure supplement 4D, which justified
             omitting the Mixture term there.
    effect - variance of the genotype EFFECT (log2FC) within the focal module,
             partitioned across Replicate / Fraction / Protein. Fraction cancels
             in the ratio, so this is the decomposition relevant to the contrast.

    NOTE on animal-level variance: each replicate is a POOL OF 5 ANIMALS, so
    between-animal variance is not recoverable from these data by construction.
    The recoverable levels are replicate (pool), protein and fraction.
    """
    d = pd.read_csv(os.path.join(DATA, f"tidy_{stem}_module_annotated.csv"))
    d = d[d.Norm_Intensity > 0].copy()
    d["logI"] = np.log2(d.Norm_Intensity)

    rows = []
    for _, g in d.groupby("Accession"):
        gm = g.logI.mean()
        ss_t = ((g.logI - gm) ** 2).sum()
        if ss_t <= 0:
            continue
        rows.append([
            sum(len(s) * (s.logI.mean() - gm) ** 2 for _, s in g.groupby(c)) / ss_t
            for c in ("Fraction", "Replicate", "Genotype")])
    raw = pd.DataFrame(np.array(rows) * 100,
                       columns=["Fraction", "Replicate", "Genotype"])

    w = paired_log2fc(stem)
    sub = w[w.Membership == FOCAL_MODULE] if stem == FOCAL_MODEL else w
    gm = sub.log2FC.mean()
    ss_t = ((sub.log2FC - gm) ** 2).sum()
    effect = {c: 100 * sum(len(s) * (s.log2FC.mean() - gm) ** 2
                           for _, s in sub.groupby(c)) / ss_t
              for c in ("Replicate", "Fraction", "Accession")}
    return raw, effect


def main():
    rng = np.random.default_rng(SEED)
    report = {"seed": SEED, "n_perm": N_PERM, "n_boot": N_BOOT}

    w = paired_log2fc(FOCAL_MODEL)
    prot = w.groupby(["Accession", "Membership"]).log2FC.mean().reset_index()
    focal = prot[prot.Membership == FOCAL_MODULE]
    n, obs = len(focal), focal.log2FC.mean()
    pool = prot.log2FC.values
    centre = pool.mean()

    # --- permutation: size-matched random protein sets ---------------------
    null = np.array([rng.choice(pool, n, replace=False).mean()
                     for _ in range(N_PERM)])
    exceed = int((np.abs(null - centre) >= abs(obs - centre)).sum())
    p_perm = exceed / N_PERM
    report["permutation"] = dict(
        n_proteins=int(n), observed_mean_log2FC=float(obs),
        proteome_mean_log2FC=float(centre), null_sd=float(null.std()),
        exceedances=exceed,
        p_value=float(p_perm) if exceed else None,
        p_upper_bound=None if exceed else float(1 / N_PERM),
        z=float((obs - null.mean()) / null.std()))

    # --- bootstrap: resample proteins -------------------------------------
    vals = focal.log2FC.values
    boot = np.array([rng.choice(vals, n, replace=True).mean()
                     for _ in range(N_BOOT)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    report["bootstrap"] = dict(mean_log2FC=float(obs), ci95_low=float(lo),
                               ci95_high=float(hi), excludes_zero=bool(hi < 0 or lo > 0))

    # --- protein-level sign test (distribution-free) -----------------------
    from scipy import stats as st
    k = int((vals < 0).sum())
    report["sign_test"] = dict(
        n=int(n), n_negative=k, pct_negative=float(100 * k / n),
        binomial_p=float(st.binomtest(k, n, 0.5, alternative="two-sided").pvalue))

    # --- competitive rank among the model's modules ------------------------
    mm = w.groupby("Membership").log2FC.mean().sort_values()
    order = list(mm.index)
    report["competitive_rank"] = dict(
        module_mean_log2FC=float(mm.loc[FOCAL_MODULE]),
        rank_most_negative=order.index(FOCAL_MODULE) + 1,
        k_modules=len(order),
        most_negative_modules={f"M{i}": round(float(v), 4)
                               for i, v in mm.head(5).items()})

    # --- per-replicate breakdown (reported, not interpreted here) ----------
    per_rep = {}
    for r, s in w[w.Membership == FOCAL_MODULE].groupby("Replicate"):
        pm = s.groupby("Accession").log2FC.mean()
        per_rep[f"Rep{r}"] = dict(mean_log2FC=float(pm.mean()),
                                  n_negative=int((pm < 0).sum()), n=int(len(pm)))
    report["per_replicate"] = per_rep

    # --- variance decompositions ------------------------------------------
    var_out = {}
    for stem in ("LRRK2", "SNCA", "VPS35"):
        raw, effect = variance_tables(stem)
        var_out[stem] = dict(
            raw_intensity_median_pct={c: round(float(raw[c].median()), 2)
                                      for c in raw.columns},
            genotype_effect_pct={k: round(float(v), 2) for k, v in effect.items()})
    report["variance_decomposition"] = var_out
    report["design_note"] = ("Each replicate is a pool of 5 animals; "
                             "animal-level variance is not recoverable by design.")

    with open(os.path.join(OUT, "permutation_bootstrap.json"), "w") as fh:
        json.dump(report, fh, indent=2)

    pd.DataFrame([report["permutation"]]).to_csv(
        os.path.join(OUT, "permutation_result.csv"), index=False)

    print(json.dumps({k: report[k] for k in
                      ("permutation", "bootstrap", "sign_test", "competitive_rank")},
                     indent=2))
    print("\nwrote: output/permutation_bootstrap.json, permutation_result.csv")


if __name__ == "__main__":
    main()
