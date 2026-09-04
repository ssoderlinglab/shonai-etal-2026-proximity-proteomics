"""
Supplemental Table S1: the Vps35 proximity proteome, WT and D620N.

Input
    data/ibioid/tables/Table_S1_iBioID_Proximity_Proteomes.xlsx
Output
    tables/output/Table_S1_iBioID_Vps35_Proximity_Proteome.xlsx
Run
    python3 tables/scripts/build_table_S1_vps35.py
Dependencies
    python >= 3.9, openpyxl

Thirteen columns per sheet. The proximity proteome is the rows called Enriched plus the bait row itself; the counts written, 361 for WT and 396 for D620N, are asserted against the manuscript.

Depleted rows are kept: they are a measured outcome of the same contrast.

Only the Benjamini-Hochberg FDR is shipped. The source also carries a second correction that gives a different count, and the manuscript reports Benjamini-Hochberg.
"""
import os
import shutil

import openpyxl
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(TABLES_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUT = os.path.join(TABLES_DIR, "output")
os.makedirs(OUT, exist_ok=True)

SOURCE = os.path.join(DATA_DIR, "ibioid", "tables",
                      "Table_S1_iBioID_Proximity_Proteomes.xlsx")

DEST_NAME = "Table_S1_iBioID_Vps35_Proximity_Proteome.xlsx"

# source sheet -> (new sheet name, expected In_Proximity_Proteome TRUE count,
#                  manuscript site)
SHEETS = [
    ("Vps35_WT", "Vps35_WT", 361, "Results 30"),
    ("Vps35_D620N", "Vps35_D620N", 396, "Results 49"),
]

# The 14 columns, in this order. U1_column_decisions.csv, proposal == KEEP.
KEEP = [
    "UniProt_Accession", "Gene_Symbol", "Protein_Name", "Peptide_Count",
    "Bait_Rep1", "Bait_Rep2", "Bait_Rep3",
    "SolubleTurboID_Rep1", "SolubleTurboID_Rep2", "SolubleTurboID_Rep3",
    "log2FC", "DEqMS_BH_FDR", "Significance",
]

DROP = ["Bait_CV_percent", "SolubleTurboID_CV_percent", "fold_change", "log2FC_SE",
        "log2FC_CI_lower", "log2FC_CI_upper", "DEqMS_p_value", "DEqMS_Storey_q",
        "In_Proximity_Proteome"]

# Significance values that make up the proximity proteome. Measured: In_Proximity_
# Proteome was a PURE FUNCTION of Significance (Enriched and Bait -> True, NS,
# Depleted and Not tested -> False), so the flag column carried no information the
# call column does not. It is dropped, and the count is asserted against these two
# values instead.
#
# ⛔ THE PROTEOME IS Enriched PLUS Bait, NOT Enriched ALONE. Filtering on Enriched
# gives 360 in Vps35_WT, not the manuscript's 361; the missing row is Vps35 itself.
# Keeping Significance rather than the flag is what makes that visible: a reader sees
# 360 + 1 and can see why the total is 361.
IN_PROTEOME = ("Enriched", "Bait")


def main():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, DEST_NAME)
    if os.path.exists(dest):
        os.remove(dest)

    src = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    out = openpyxl.Workbook()
    out.remove(out.active)

    report = []
    for src_name, new_name, expect, site in SHEETS:
        ws_in = src[src_name]
        rows = list(ws_in.values)
        header = list(rows[0])
        missing = [c for c in KEEP if c not in header]
        assert not missing, "{}: columns absent from source: {}".format(src_name, missing)
        idx = [header.index(c) for c in KEEP]

        ws = out.create_sheet(new_name)
        ws.append(KEEP)
        n_true = 0
        accessions = set()
        i_sig = KEEP.index("Significance")
        i_acc = KEEP.index("UniProt_Accession")
        for r in rows[1:]:
            if r[idx[i_acc]] is None:
                continue
            vals = [r[j] for j in idx]
            ws.append(vals)
            if str(vals[i_sig]).strip() in IN_PROTEOME:
                n_true += 1
                accessions.add(vals[i_acc])

        ws.freeze_panes = "A2"
        for j, name in enumerate(KEEP, start=1):
            ws.cell(row=1, column=j).font = Font(bold=True)
            ws.cell(row=1, column=j).alignment = Alignment(horizontal="center")
            width = max(11, min(38, len(name) + 3))
            if name == "Protein_Name":
                width = 46
            ws.column_dimensions[get_column_letter(j)].width = width

        assert n_true == expect, (
            "{}: Significance in {} gives {}, manuscript {} states {}".format(
                new_name, IN_PROTEOME, n_true, site, expect))
        assert len(accessions) == expect, (
            "{}: {} TRUE rows but {} distinct accessions".format(
                new_name, n_true, len(accessions)))
        report.append((new_name, ws.max_row - 1, len(KEEP), n_true, site))

    src.close()
    out.save(dest)

    print("wrote " + os.path.relpath(dest, REPO_ROOT))
    print("{} sheets, {} columns each (was 22)\n".format(len(SHEETS), len(KEEP)))
    print("  {:14} {:>7} {:>6}  {:>10}  {}".format(
        "sheet", "rows", "cols", "in proteome", "manuscript"))
    for name, nrow, ncol, ntrue, site in report:
        print("  {:14} {:>7} {:>6}  {:>11}  {}".format(name, nrow, ncol, ntrue, site))
    print("\n  dropped: " + ", ".join(DROP))


if __name__ == "__main__":
    main()
