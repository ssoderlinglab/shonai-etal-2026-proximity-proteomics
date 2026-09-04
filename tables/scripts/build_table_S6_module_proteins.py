"""
Supplemental Table S6: per-protein fold changes of the selected co-fractionation modules.

Input
    data/ibiocofrac/modules/Partition_Gene_combined_3PD_models.csv
    data/ibiocofrac/fractionation/{LRRK2,SNCA,Vps35}_young_for_PCA.csv
    data/ibiocofrac/permutation/permutation_12modules.csv
    data/ibiocofrac/m11/M11_77_gene_list.csv
Output
    tables/output/Table_S6_iBioCoFrac_Module_Proteins.xlsx
Run
    python3 tables/scripts/build_table_S6_module_proteins.py
Dependencies
    python >= 3.9, pandas, openpyxl

Modules are selected by a consistent annotation and eight or more proteins per cohort. Fold changes are recomputed through the same path the permutation uses, so the per-protein values and the module-level statistic agree.
"""
import csv
import os
import sys

import numpy as np
import openpyxl
import pandas as pd
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(TABLES_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUT = os.path.join(TABLES_DIR, "output")
os.makedirs(OUT, exist_ok=True)

PARTITION = os.path.join(DATA_DIR, "ibiocofrac", "modules",
                         "Partition_Gene_combined_3PD_models.csv")
FRACTION_DIR = os.path.join(DATA_DIR, "ibiocofrac", "fractionation")
FRACTION_PATHS = {
    "LRRK2": os.path.join(FRACTION_DIR, "LRRK2_young_for_PCA.csv"),
    "SNCA": os.path.join(FRACTION_DIR, "SNCA_young_for_PCA.csv"),
    "VPS35": os.path.join(FRACTION_DIR, "Vps35_young_for_PCA.csv"),
}
MODULE_PERM = os.path.join(DATA_DIR, "ibiocofrac", "permutation",
                           "permutation_12modules.csv")
M11_LIST = os.path.join(DATA_DIR, "ibiocofrac", "m11", "M11_77_gene_list.csv")

DEST_NAME = "Table_S6_iBioCoFrac_Module_Proteins.xlsx"
SHEET_NAME = "Module_Proteins"

FRACTIONS = ["4", "5", "6", "7", "8", "9", "10"]
REFERENCE_COHORT = "LRRK2"          # the partition the module numbers come from

# model token -> printed column. Cohort names as the manuscript writes them.
MODELS = [("LRRK2", "Log2FC_Tg_LRRK2_G2019S"),
          ("SNCA", "Log2FC_Tg_SNCA_A30P_A53T"),
          ("VPS35", "Log2FC_KI_Vps35_D620N")]

# the twelve modules the paper analyses, in the order Figure 5E prints them, with the
# labels that figure uses. Taken from the permutation output, not retyped.
MODULE_ORDER = ["M1", "M2", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M14",
                "M19"]

HEADER = (["Module", "Description", "Gene", "UniProt_Accession"]
          + [col for _t, col in MODELS])

NOT_DETECTED = "Not detected"
M11_SIZE = 77                       # settled module size


def compute_protein_fc(frac_df):
    """The permutation's compute_protein_fc, reproduced verbatim (which reproduces
    the shipped
    Figure 5 function). Collapses replicate AND fraction. Do not "improve" it: any
    change here breaks the correspondence with the reported module statistics."""
    df = frac_df[frac_df["fraction"].astype(str).isin(FRACTIONS)].copy()
    wt = df[df["sample"].str.startswith("wt_")]
    tg = df[df["sample"].str.startswith("tg_")]
    meta_cols = [df.columns[0], "fraction", "replicate", "sample"]
    prot_cols = [c for c in df.columns if c not in meta_cols]
    rows = []
    for prot in prot_cols:
        wt_mean = np.nanmean(wt[prot].values.astype(float))
        tg_mean = np.nanmean(tg[prot].values.astype(float))
        fc = (np.log2(tg_mean / wt_mean)
              if wt_mean > 0 and tg_mean > 0 else np.nan)
        rows.append({"Accession": prot, "log2FC": fc})
    return pd.DataFrame(rows)


def per_gene_log2fc(partition):
    """({model: {gene: log2FC}}, contaminant_set) through the permutation's own path.

    The contaminant set is returned rather than applied: this table carries every
    measurement, and the statistic it supports filters. See the module docstring.
    """
    sys.path.insert(0, os.path.join(REPO_ROOT, "src"))
    from contamination_filter import CONTAMINATION_GENES as EXCLUDE
    out = {}
    for token, _col in MODELS:
        pfc = compute_protein_fc(pd.read_csv(FRACTION_PATHS[token]))
        rows = partition[partition["Experiment"] == token]
        a2g = dict(zip(rows["UniProt"], rows["Symbol"]))
        pfc["Gene"] = pfc["Accession"].map(a2g)
        d = (pfc[["Gene", "log2FC"]].dropna(subset=["Gene"])
             .drop_duplicates(subset="Gene", keep="first"))
        out[token] = {g: (None if pd.isna(v) else float(v))
                      for g, v in zip(d["Gene"], d["log2FC"])}
    return out, set(EXCLUDE)


def main():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, DEST_NAME)
    if os.path.exists(dest):
        os.remove(dest)

    partition = pd.read_csv(PARTITION)
    ref = partition[partition["Experiment"] == REFERENCE_COHORT]
    perm = {(r["module"], r["model"]): r
            for r in csv.DictReader(open(MODULE_PERM))}
    labels = {r["module"]: r["label"] for r in csv.DictReader(open(MODULE_PERM))}
    fc, contaminants = per_gene_log2fc(partition)

    assert set(labels) == set(MODULE_ORDER), (
        "the permutation covers {}, this script prints {}".format(
            sorted(labels), MODULE_ORDER))

    # membership, from the reference partition. One row per protein per module.
    members = {}
    for m in MODULE_ORDER:
        num = m[1:]
        rows = ref[ref["Membership"].astype(str).str.strip() == num]
        members[m] = list(zip(rows["Symbol"], rows["UniProt"]))
        assert len({a for _s, a in members[m]}) == len(members[m]), (
            "{}: duplicate accession in the partition".format(m))

    # ⛔ M11 IS 77. Asserted against the distributed gene list, not just the partition.
    m11_list = {r["Gene"] for r in csv.DictReader(open(M11_LIST))}
    assert len(members["M11"]) == M11_SIZE, "M11 has {} rows, expected {}".format(
        len(members["M11"]), M11_SIZE)
    assert {g for g, _a in members["M11"]} == m11_list, (
        "the partition's M11 differs from the distributed 77-gene list")

    out = openpyxl.Workbook()
    ws = out.active
    ws.title = SHEET_NAME
    ws.append(HEADER)

    per_module = {}
    for m in MODULE_ORDER:
        counts = {col: 0 for _t, col in MODELS}
        for gene, acc in sorted(members[m], key=lambda p: p[0].lower()):
            row = [m, labels[m], gene, acc]
            for token, col in MODELS:
                v = fc[token].get(gene)
                if v is None:
                    row.append(NOT_DETECTED)
                else:
                    row.append(v)
                    counts[col] += 1
            ws.append(row)
        per_module[m] = counts

    # ⛔ EVERY ONE OF THE 36 MODULE DELTAS MUST FALL OUT OF THESE COLUMNS. The
    # statistic filters contaminants and this table does not, so the filter is applied
    # HERE, at the point of comparison. This is what proves the columns are the reported
    # vintage of the pipeline and not one of the 2026-03 files.
    for m in MODULE_ORDER:
        for token, col in MODELS:
            p = perm[(m, token)]
            vals = [fc[token][g] for g, _a in members[m]
                    if g not in contaminants and fc[token].get(g) is not None]
            assert len(vals) == int(p["n"]), (
                "{} {}: {} values after the contamination filter, the "
                "permutation used {}".format(m, token, len(vals), p["n"]))
            delta = sum(vals) / len(vals) - float(p["pool_mean"])
            assert abs(delta - float(p["delta"])) < 5e-4, (
                "{} {}: recomputed delta {:.6f}, permutation {:.6f}".format(
                    m, token, delta, float(p["delta"])))

    ws.freeze_panes = "E2"
    for token, col in MODELS:
        j = HEADER.index(col) + 1
        for row in ws.iter_rows(min_row=2, min_col=j, max_col=j):
            for c in row:
                c.number_format = "0.0000"
    for j, name in enumerate(HEADER, start=1):
        c = ws.cell(row=1, column=j)
        c.font = Font(bold=True)
        c.alignment = Alignment(horizontal="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = max(11, len(name) + 2)
    ws.column_dimensions["A"].width = 9
    ws.column_dimensions["B"].width = 16
    ws.column_dimensions["C"].width = 12
    ws.column_dimensions["D"].width = 18
    ws.row_dimensions[1].height = 28

    out.save(dest)

    total = sum(len(members[m]) for m in MODULE_ORDER)
    print("wrote " + os.path.relpath(dest, REPO_ROOT))
    print("1 sheet, {} rows x {} columns\n".format(total, len(HEADER)))
    print("  {:6} {:16} {:>8}  {:>10} {:>10} {:>10}".format(
        "module", "description", "proteins", "LRRK2", "SNCA", "VPS35"))
    for m in MODULE_ORDER:
        c = per_module[m]
        print("  {:6} {:16} {:>8}  {:>10} {:>10} {:>10}".format(
            m, labels[m][:16], len(members[m]),
            *[c[col] for _t, col in MODELS]))
    print("  {:6} {:16} {:>8}  {:>10} {:>10} {:>10}".format(
        "TOTAL", "", total,
        *[sum(per_module[m][col] for m in MODULE_ORDER) for _t, col in MODELS]))
    print('\n  remaining cells read "{}"'.format(NOT_DETECTED))
    print("\n  all {} module deltas reproduced from these columns "
          "(contaminants excluded at comparison time)".format(
              len(MODULE_ORDER) * len(MODELS)))
    inmod = sorted({g for m in MODULE_ORDER for g, _a in members[m]}
                   & contaminants)
    print("  contaminants inside these twelve modules: {}".format(
        ", ".join(inmod) if inmod else "none"))


if __name__ == "__main__":
    main()
