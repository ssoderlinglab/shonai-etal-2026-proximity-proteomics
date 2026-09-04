"""
Supplemental Table S4: proteins shared across the three proximity proteomes.

Input
    data/ibioid/tables/Table_S1_iBioID_Proximity_Proteomes.xlsx
Output
    tables/output/Table_S4_iBioID_Convergence.xlsx
Run
    python3 tables/scripts/build_table_S4_shared_proxiome.py
Dependencies
    python >= 3.9, openpyxl

The four convergence counts written are asserted against the manuscript. Membership is counted by accession, not by gene symbol: a symbol-level count gives different totals.
"""
import os

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
SOURCE_SHEET = "Convergence_Summary"

DEST_NAME = "Table_S4_iBioID_Convergence.xlsx"
SHEET_NAME = "Convergence"

# 10 columns. WT block then mutant block, each ending in its category, so the flags a
# category summarises are the three cells to its left.
KEEP = [
    "UniProt_Accession", "Gene_Symbol",
    "In_Vps35_WT", "In_aSyn_WT", "In_Lrrk2_WT", "WT_Category",
    "In_Vps35_D620N", "In_aSyn_A30P", "In_Lrrk2_G2019S", "MUT_Category",
]

DROP = ["N_Baits_WT", "N_Baits_MUT", "WT_to_MUT_Shift"]

# Results paragraph 44. category -> (expected rows, what the manuscript calls it)
EXPECT_WT = {
    "Total_Convergent": (86, "detected by all three baits"),
    "Partial_Convergent": (156, "additional, detected by at least two baits"),
}
EXPECT_SHARED = 242   # the two above
EXPECT_UNIQUE = 769   # anything detected in at least one WT bait


def main():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, DEST_NAME)
    if os.path.exists(dest):
        os.remove(dest)

    src = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    rows = list(src[SOURCE_SHEET].values)
    src.close()
    header = list(rows[0])
    missing = [c for c in KEEP if c not in header]
    assert not missing, "columns absent from source: {}".format(missing)
    idx = [header.index(c) for c in KEEP]

    out = openpyxl.Workbook()
    ws = out.active
    ws.title = SHEET_NAME
    ws.append(KEEP)

    i_acc = KEEP.index("UniProt_Accession")
    i_wt = KEEP.index("WT_Category")
    i_mut = KEEP.index("MUT_Category")
    by_wt = {}
    by_mut = {}
    for r in rows[1:]:
        if r[idx[i_acc]] is None:
            continue
        vals = [r[j] for j in idx]
        ws.append(vals)
        by_wt.setdefault(str(vals[i_wt]), set()).add(vals[i_acc])
        by_mut.setdefault(str(vals[i_mut]), set()).add(vals[i_acc])

    ws.freeze_panes = "C2"
    for j, name in enumerate(KEEP, start=1):
        ws.cell(row=1, column=j).font = Font(bold=True)
        ws.cell(row=1, column=j).alignment = Alignment(horizontal="center")
        ws.column_dimensions[get_column_letter(j)].width = max(13, len(name) + 3)

    n_rows = ws.max_row - 1
    for cat, (expect, _label) in EXPECT_WT.items():
        got = len(by_wt.get(cat, ()))
        assert got == expect, "WT {} is {}, Results 44 states {}".format(
            cat, got, expect)
    shared = len(by_wt.get("Total_Convergent", set())
                 | by_wt.get("Partial_Convergent", set()))
    assert shared == EXPECT_SHARED, "shared is {}, Results 44 states {}".format(
        shared, EXPECT_SHARED)
    unique = sum(len(v) for k, v in by_wt.items() if k != "Not_Detected")
    assert unique == EXPECT_UNIQUE, "unique is {}, Results 44 states {}".format(
        unique, EXPECT_UNIQUE)

    out.save(dest)

    print("wrote " + os.path.relpath(dest, REPO_ROOT))
    print("1 sheet, {} rows x {} columns (source had 13)\n".format(n_rows, len(KEEP)))
    print("  Results 44's four numbers, as row counts in WT_Category:")
    for cat, (expect, label) in EXPECT_WT.items():
        print("    {:20} {:>5}   {}".format(cat, len(by_wt[cat]), label))
    print("    {:20} {:>5}   total shared".format("(the two summed)", shared))
    print("    {:20} {:>5}   unique across the three WT datasets".format(
        "(not Not_Detected)", unique))
    print("\n  MUT_Category, not asserted because the manuscript does not state it:")
    for cat in ("Total_Convergent", "Partial_Convergent", "Bait_Specific",
                "Not_Detected"):
        print("    {:20} {:>5}".format(cat, len(by_mut.get(cat, ()))))
    print("\n  dropped: " + ", ".join(DROP))


if __name__ == "__main__":
    main()
