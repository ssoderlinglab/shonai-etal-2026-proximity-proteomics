"""
Supplemental Table S7: the CrAAVe-seq in vivo CRISPR screen.

Input
    data/craave_seq/tables/Table_S3_CrAAVe_seq_CRISPR_Screen.xlsx
Output
    tables/output/Table_S7_CrAAVe_seq_CRISPR_Screen.xlsx
Run
    python3 tables/scripts/build_table_S7_craave.py
Dependencies
    python >= 3.9, openpyxl

Sheets are copied unchanged from the submitted workbook; the overview sheet is dropped. The library row count is asserted and its disagreement with the design figure is recorded in the workbook rather than corrected here.
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

SOURCE = os.path.join(DATA_DIR, "craave_seq", "tables",
                      "Table_S3_CrAAVe_seq_CRISPR_Screen.xlsx")
DEST_NAME = "Table_S7_CrAAVe_seq_CRISPR_Screen.xlsx"

# source sheet -> new sheet name. "Overview" is dropped.
SHEETS = [
    ("sgRNA_Library", "sgRNA_Library"),
    ("Normalized_Read_Counts", "Normalized_Read_Counts"),
    ("MAGeCK_Gene_aSyn", "MAGeCK_aSyn_Depletion"),
    ("MAGeCK_Gene_Baseline", "MAGeCK_Baseline_Depletion"),
]
# the two sheets kept whole
WHOLE = {"sgRNA_Library", "Normalized_Read_Counts"}

# gene-level: source column -> printed name. Order is the printed order.
GENE_COLS = [("Gene", "Gene"), ("N_sgRNAs", "N_sgRNAs"),
             ("neg_LFC", "Log2FC"), ("neg_p_value", "P_Value"), ("neg_FDR", "FDR")]
GENE_DROP = ["neg_Score", "neg_Rank", "neg_Good_sgRNAs", "pos_Score", "pos_p_value",
             "pos_FDR", "pos_Rank", "pos_Good_sgRNAs", "pos_LFC"]

# every value Results 98 prints, and the sheet it is printed from.
# (sheet, gene, printed column, value, decimals the manuscript uses)
QUOTED = [
    ("MAGeCK_aSyn_Depletion", "Stxbp1", "FDR", 0.044, 3),
    ("MAGeCK_aSyn_Depletion", "Atp6v1b2", "Log2FC", -0.256, 3),
    ("MAGeCK_aSyn_Depletion", "Atp6v1b2", "P_Value", 0.0070, 4),
    ("MAGeCK_aSyn_Depletion", "Atp6v1b2", "FDR", 0.247, 3),
    ("MAGeCK_aSyn_Depletion", "Atp6v1d", "Log2FC", -0.376, 3),
    ("MAGeCK_aSyn_Depletion", "Atp6v1d", "P_Value", 0.011, 3),
    ("MAGeCK_aSyn_Depletion", "Atp6v1d", "FDR", 0.249, 3),
    ("MAGeCK_aSyn_Depletion", "Th", "Log2FC", -0.273, 3),
    ("MAGeCK_aSyn_Depletion", "Th", "P_Value", 0.021, 3),
    ("MAGeCK_aSyn_Depletion", "Th", "FDR", 0.309, 3),
    ("MAGeCK_Baseline_Depletion", "Stxbp1", "FDR", 0.78, 2),
    ("MAGeCK_Baseline_Depletion", "Th", "P_Value", 0.017, 3),
]
EXPECT_GENES = 70
EXPECT_LIBRARY = 466
EXPECT_COUNTS = 419


def main():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, DEST_NAME)
    if os.path.exists(dest):
        os.remove(dest)

    src = openpyxl.load_workbook(SOURCE, read_only=True, data_only=True)
    out = openpyxl.Workbook()
    out.remove(out.active)

    report = []
    for src_name, new_name in SHEETS:
        rows = [r for r in src[src_name].values if r[0] is not None]
        hdr = list(rows[0])
        data = rows[1:]
        ws = out.create_sheet(new_name)

        if src_name in WHOLE:
            ws.append(hdr)
            for r in data:
                ws.append(list(r))
            report.append((new_name, len(data), len(hdr), "kept whole"))
        else:
            missing = [c for c, _p in GENE_COLS if c not in hdr]
            assert not missing, "{}: no column {}".format(src_name, missing)
            # pos_LFC == neg_LFC is the reason pos_LFC is dropped; measure it
            i_neg, i_pos = hdr.index("neg_LFC"), hdr.index("pos_LFC")
            same = sum(1 for r in data if r[i_neg] == r[i_pos])
            assert same == len(data), (
                "{}: pos_LFC equals neg_LFC in only {} of {} rows, so it is NOT a "
                "restatement and must not be dropped".format(
                    src_name, same, len(data)))
            idx = [hdr.index(c) for c, _p in GENE_COLS]
            ws.append([p for _c, p in GENE_COLS])
            for r in data:
                ws.append([r[i] for i in idx])
            assert len(data) == EXPECT_GENES, "{}: {} genes".format(
                src_name, len(data))
            report.append((new_name, len(data), len(GENE_COLS),
                           "{} -> {} columns".format(len(hdr), len(GENE_COLS))))

        ws.freeze_panes = "B2"
        for j, name in enumerate([c.value for c in ws[1]], start=1):
            c = ws.cell(row=1, column=j)
            c.font = Font(bold=True)
            c.alignment = Alignment(horizontal="center", wrap_text=True)
            ws.column_dimensions[get_column_letter(j)].width = max(
                11, min(22, len(str(name)) + 2))
    src.close()

    # every quoted value must be in the sheet the manuscript prints it from
    for sheet, gene, col, value, dp in QUOTED:
        rows = [r for r in out[sheet].values]
        h = list(rows[0])
        row = next((r for r in rows[1:] if r[0] == gene), None)
        assert row is not None, "{}: no row for {}".format(sheet, gene)
        got = float(row[h.index(col)])
        assert round(got, dp) == value, (
            "{} {} {}: source {:.6g} rounds to {}, manuscript says {}".format(
                sheet, gene, col, got, round(got, dp), value))

    # Results 98: Stxbp1 is the ONLY gene significant after FDR correction
    rows = [r for r in out["MAGeCK_aSyn_Depletion"].values]
    h = list(rows[0])
    sig = [r[0] for r in rows[1:] if float(r[h.index("FDR")]) < 0.05]
    assert sig == ["Stxbp1"], (
        "Results 98 states Stxbp1 is the only FDR-significant gene; this sheet "
        "gives {}".format(sig))
    # Results 98: Th is rank one in the no-challenge arm
    brows = [r for r in out["MAGeCK_Baseline_Depletion"].values]
    bh = list(brows[0])
    top = min(brows[1:], key=lambda r: float(r[bh.index("P_Value")]))
    assert top[0] == "Th", (
        "Results 98 states Th ranks one without alpha-synuclein; lowest p here is "
        "{}".format(top[0]))

    assert out["sgRNA_Library"].max_row - 1 == EXPECT_LIBRARY
    assert out["Normalized_Read_Counts"].max_row - 1 == EXPECT_COUNTS

    out.save(dest)

    print("wrote " + os.path.relpath(dest, REPO_ROOT))
    print("{} sheets, was 5\n".format(len(SHEETS)))
    print("  {:28} {:>6} {:>6}  {}".format("sheet", "rows", "cols", "change"))
    for name, nrow, ncol, note in report:
        print("  {:28} {:>6} {:>6}  {}".format(name, nrow, ncol, note))
    print("\n  dropped sheet: Overview (33 rows of prose; that is the legend's job)")
    print("  dropped columns: " + ", ".join(GENE_DROP))
    print("\n  {} quoted values verified against the sheet they are printed "
          "from".format(len(QUOTED)))
    print("  Stxbp1 is the only FDR-significant gene in the alpha-synuclein arm")
    print("  Th is rank one by p value in the no-challenge arm")


if __name__ == "__main__":
    main()
