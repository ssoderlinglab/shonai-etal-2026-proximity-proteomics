"""
Supplemental Table S5: size-matched permutation for the individual co-fractionation modules and for the module groups, in the three models.

Input
    data/ibiocofrac/permutation/permutation_12modules.csv
    data/ibiocofrac/permutation/permutation_all21.csv
Output
    tables/output/Table_S5_iBioCoFrac_Module_Permutation.xlsx
Run
    python3 tables/scripts/build_table_S5_ibiocofrac.py
Dependencies
    python >= 3.9, openpyxl

Two sheets, 36 and 21 rows. The Benjamini-Hochberg corrections in the two sheets were applied within different family sizes, 36 and 21, and the sheets carry the family size so the two q columns are not read as one.

Every row is named by its member modules; the descriptive word goes in the description column.

data/ibiocofrac/permutation/permutation_all21.csv is the same table that figures/main_05/scripts/module_group_permutation.py regenerates. It is shipped as data so this builder does not depend on a figure having been run first.
"""
import csv
import os

import openpyxl
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TABLES_DIR = os.path.dirname(SCRIPT_DIR)
REPO_ROOT = os.path.dirname(TABLES_DIR)
DATA_DIR = os.path.join(REPO_ROOT, "data")
OUT = os.path.join(TABLES_DIR, "output")
os.makedirs(OUT, exist_ok=True)

PERM = os.path.join(DATA_DIR, "ibiocofrac", "permutation")

DEST_NAME = "Table_S5_iBioCoFrac_Module_Permutation.xlsx"

# sheet name, source csv, identifier column in the source, new identifier header,
# family size, expected rows
# sheet name, source csv, identifier column, new identifier header, label column,
# BH family size, expected source rows
SHEETS = [
    ("Individual_Modules", "permutation_12modules.csv",
     "module", "Modules", "label", 36, 36),
    ("Module_Groups", "permutation_all21.csv",
     "group", "Modules", "group_label", 21, 21),
]

# ⛔ THE IDENTIFIER COLUMN IS UNIFORM: EVERY ROW IS NAMED BY MODULE NUMBER.
# The source mixes two vocabularies in one column: word labels for the multi-module
# groups (Synapse, Ribo, Mito) and bare numbers for the single-module ones (M4, M6,
# M9, M14). A reader scanning that column cannot tell what kind of thing each row is.
# Here every row is named the same way, by its member modules, and the word goes in
# Description where the other sheet's words also live.
GROUP_ID = {"Synapse": "M7+M8+M11+M19", "Ribo": "M1+M10", "Mito": "M2+M5",
            "M4": "M4", "M6": "M6", "M9": "M9", "M14": "M14"}
GROUP_DESC = {"Synapse": "Synapse", "Ribo": "Ribosome", "Mito": "Mitochondria",
              "M4": "ER/Golgi", "M6": "Clathrin/tRNA", "M9": "Proteasome",
              "M14": "RNA-binding/ND"}

# source token -> the name the manuscript and the figures use
MODEL_DISPLAY = {"LRRK2": "Tg-LRRK2 G2019S", "SNCA": "Tg-SNCA A30P x A53T",
                 "VPS35": "KI-Vps35 D620N"}
MODELS_N = list(MODEL_DISPLAY)

# One hue per model, header and a paler banded pair for the data rows. Chosen light
# enough that black text stays readable and the sheet prints on a monochrome printer.
MODEL_FILL = ["DCE6F1", "E4DFEC", "EAF1DD"]
MODEL_FILL_ALT = ["EAF1F8", "F1EEF7", "F4F8EC"]
GROUP_ORDER = ["Synapse", "Ribo", "Mito"]

HEADER_TAIL = ["Model", "N_Proteins_Defined", "N_Proteins_Quantified",
               "Proteome_Baseline_Mean", "Delta_From_Baseline", "Q_Permutation_BH"]

# ⛔ TWO SIZE COLUMNS, AND THEY ARE NOT A DUPLICATE. The two columns answer
# different questions, so the table states both.
#
#   N_Proteins_Defined      how many proteins the module CONTAINS. Counted from the
#                           Tg-LRRK2 partition, which is where every module is defined.
#                           This is the number Results 77 prints in prose: M11 = 77.
#   N_Proteins_Quantified   how many of those were quantified in THIS cohort, i.e. the
#                           permutation n. delta is a mean over these proteins,
#                           so a reader
#                           needs it to recheck delta.
#
# For the reference cohort the two are close but not always equal (M11 77 vs 76, M7 88
# vs 86, M2 149 vs 148, M14 66 vs 64): a protein can be in the partition and still be
# missing from the matrix the permutation ran on. For the two test cohorts they diverge
# by design, and that divergence IS the mapping the manuscript describes.
#
# ⛔ NEITHER IS THE COHORT'S OWN PARTITION INDEX. The old S2 reported that (SNCA M11 =
# 30, VPS35 M11 = 40) and it is a DIFFERENT PROTEIN SET carrying the same number in a
# different partition. It appears nowhere in this table.
PARTITION = os.path.join(DATA_DIR, "ibiocofrac", "modules",
                         "Partition_Gene_combined_3PD_models.csv")
REFERENCE_COHORT = "LRRK2"

DROP = ["mean_raw", "null_p2_5", "null_median", "null_p97_5",
        "null_sd", "z", "percentile_of_obs", "n_exceedances", "p_perm",
        "p_at_resolution_bound"]

MODEL_ORDER = {"LRRK2": 0, "SNCA": 1, "VPS35": 2}

# The functional groups of the 21-test family, and which modules each spans. Taken from
# the group_label column of AN10_permutation_all21_v01.csv, which names them: a group's
# definition size is the sum over its member modules.
GROUP_MEMBERS = {
    "Synapse": ["M7", "M8", "M11", "M19"],
    "Ribo": ["M1", "M10"],
    "Mito": ["M2", "M5"],
}

# Every value Results 76, 77 and 78 print, with the sheet it must come from.
# (sheet, identifier, model, delta to 3 dp, q as the manuscript prints it, paragraph)
QUOTED = [
    ("Module_Groups", "Synapse", "LRRK2", -0.124, 3.5e-5, "Results 76"),
    ("Module_Groups", "Synapse", "SNCA", -0.051, 0.042, "Results 76"),
    ("Module_Groups", "Synapse", "VPS35", -0.052, 3.2e-4, "Results 76"),
    ("Individual_Modules", "M8", "LRRK2", -0.062, 0.005, "Results 77"),
    ("Individual_Modules", "M8", "SNCA", -0.197, 0.002, "Results 77"),
    ("Individual_Modules", "M19", "LRRK2", -0.154, 1.4e-4, "Results 77"),
    ("Individual_Modules", "M19", "SNCA", -0.195, 0.043, "Results 77"),
    ("Individual_Modules", "M7", "LRRK2", -0.079, 4.1e-4, "Results 77"),
    ("Individual_Modules", "M7", "SNCA", 0.115, 0.013, "Results 77"),
    ("Individual_Modules", "M11", "LRRK2", -0.229, 4.5e-5, "Results 77"),
    ("Individual_Modules", "M11", "SNCA", -0.094, 0.045, "Results 77"),
    ("Individual_Modules", "M11", "VPS35", -0.149, 4.5e-5, "Results 77"),
    ("Module_Groups", "M14", "LRRK2", -0.065, 0.012, "Results 78"),
    ("Module_Groups", "M14", "SNCA", -0.125, 0.042, "Results 78"),
    ("Module_Groups", "M14", "VPS35", -0.079, 0.015, "Results 78"),
    ("Module_Groups", "Mito", "VPS35", 0.119, 3.5e-5, "Results 78"),
    ("Module_Groups", "Ribo", "LRRK2", 0.055, 3.5e-5, "Results 78"),
    ("Module_Groups", "M9", "LRRK2", 0.054, 0.015, "Results 78"),
]

# Module sizes Results 77 states in prose, verbatim: "M8 (postsynaptic/cytoskeletal,
# 85 proteins...)", "M19 (clathrin-mediated endocytosis, 38 proteins)", "M7 (SV
# cycling, 86 proteins)", "M11 (presynaptic/SV, 77 proteins)".
#
# ⛔ THE PROSE MIXES THE TWO BASES, AND THAT IS WHY THE TABLE CARRIES BOTH COLUMNS.
# Measured against the Tg-LRRK2 partition (definition) and the mapped count the
# permutation used:
#
#     module   prose   definition   mapped n   the prose value is
#     M8          85           85         85   both, they coincide
#     M19         38           38         38   both, they coincide
#     M7          86           88         86   the MAPPED count
#     M11         77           77         76   the DEFINITION size
#
# So no single column reproduces all four, and asserting all four against one of them
# would be asserting something the manuscript does not say. Each prose value is
# checked against whichever column it is, and the assertion below requires that EVERY
# prose value be recoverable from SOME column of the row, which is the property that
# actually matters to a reader checking the paragraph.
QUOTED_N = [("M8", 85), ("M19", 38), ("M7", 86), ("M11", 77)]


def main():
    os.makedirs(OUT, exist_ok=True)
    dest = os.path.join(OUT, DEST_NAME)
    if os.path.exists(dest):
        os.remove(dest)

    # definition sizes, counted from the reference cohort's partition
    defined = {}
    with open(PARTITION) as fh:
        for r in csv.DictReader(fh):
            if r["Experiment"].strip() == REFERENCE_COHORT:
                key = "M" + str(r["Membership"]).strip()
                defined[key] = defined.get(key, 0) + 1
    assert defined, "no {} rows in {}".format(REFERENCE_COHORT, PARTITION)

    out = openpyxl.Workbook()
    out.remove(out.active)
    tables = {}

    for sheet_name, csv_name, id_col, id_header, label_col, family, expect_rows in \
            SHEETS:
        rows = list(csv.DictReader(open(os.path.join(PERM, csv_name))))
        assert len(rows) == expect_rows, "{}: {} rows, expected {}".format(
            csv_name, len(rows), expect_rows)

        # pivot: one entry per identifier, the three models inside it
        wide = {}
        for r in rows:
            e = wide.setdefault(r[id_col], {"label": r[label_col], "models": {}})
            e["models"][r["model"]] = r
        assert len(wide) * 3 == expect_rows, (
            "{}: {} identifiers x 3 != {} rows".format(
                csv_name, len(wide), expect_rows))

        is_groups = sheet_name == "Module_Groups"
        ws = out.create_sheet(sheet_name)

        # two header rows: the model spans above, the measures below
        # ⛔ pool_mean IS NOT A COLUMN. Measured: it takes exactly ONE value per
        # model (LRRK2 -0.028290, SNCA 0.095988, VPS35 -0.083827) across all 36 and
        # all 21 source rows, because it is the cohort's whole-proteome baseline, not
        # a per-module quantity. As a column it would repeat one number twelve times
        # under each model. It goes in the model's header instead, where a reader sees
        # it once and can still recover the raw module mean as baseline + delta.
        baselines = {}
        for r in rows:
            baselines.setdefault(r["model"], set()).add(float(r["pool_mean"]))
        assert all(len(v) == 1 for v in baselines.values()), (
            "pool_mean is not constant within a model: {}".format(baselines))
        baseline = {m: v.pop() for m, v in baselines.items()}

        top = [id_header, "Description", "N_Proteins_Defined"]
        bottom = ["", "", ""]
        for model in MODEL_DISPLAY:
            top += ["{}\nproteome baseline {:+.4f}".format(
                MODEL_DISPLAY[model], baseline[model]), "", ""]
            bottom += ["N quantified", "Delta", "Q (BH)"]
        ws.append(top)
        ws.append(bottom)
        for j in (1, 2, 3):
            ws.merge_cells(start_row=1, start_column=j, end_row=2, end_column=j)
        for k in range(len(MODEL_DISPLAY)):
            c0 = 4 + 3 * k
            ws.merge_cells(start_row=1, start_column=c0, end_row=1,
                           end_column=c0 + 2)

        data = {}
        order = sorted(wide, key=lambda i: (
            0 if i in GROUP_MEMBERS else 1,
            GROUP_ORDER.index(i) if i in GROUP_ORDER else 0,
            int(i[1:]) if i.startswith("M") and i[1:].isdigit() else 999))
        for ident in order:
            e = wide[ident]
            display = GROUP_ID.get(ident, ident) if is_groups else ident
            description = GROUP_DESC.get(ident, e["label"]) if is_groups \
                else e["label"]
            n_def = defined.get(ident)
            if n_def is None:
                members = GROUP_MEMBERS.get(ident)
                n_def = sum(defined[m] for m in members) if members else ""
            row = [display, description, n_def]
            for model in MODEL_DISPLAY:
                r = e["models"][model]
                row += [int(r["n"]), float(r["delta"]), float(r["q_perm_BH"])]
                data[(ident, model)] = (n_def, int(r["n"]), float(r["delta"]),
                                        float(r["q_perm_BH"]))
            ws.append(row)
        tables[sheet_name] = data

        ws.freeze_panes = "D3"
        # A light fill per model block, so the eye stays inside one model while
        # reading across, and banded data rows so it stays on one module.
        for k in range(len(MODELS_N)):
            c0 = 4 + 3 * k
            fill = PatternFill("solid", fgColor=MODEL_FILL[k])
            for j in range(c0, c0 + 3):
                for i in range(1, ws.max_row + 1):
                    cell = ws.cell(row=i, column=j)
                    if i > 2 and i % 2 == 1:
                        cell.fill = PatternFill("solid", fgColor=MODEL_FILL_ALT[k])
                    else:
                        cell.fill = fill
        for j in (1, 2, 3):
            for i in range(1, ws.max_row + 1):
                cell = ws.cell(row=i, column=j)
                if i <= 2:
                    cell.fill = PatternFill("solid", fgColor="D9D9D9")
                elif i % 2 == 1:
                    cell.fill = PatternFill("solid", fgColor="F2F2F2")
        thin = Side(style="thin", color="AAAAAA")
        for i in range(1, ws.max_row + 1):
            for j in range(1, ws.max_column + 1):
                ws.cell(row=i, column=j).border = Border(
                    bottom=Side(style="thin", color="D0D0D0")
                    if i == 2 else None)
        for row in (1, 2):
            for c in ws[row]:
                c.font = Font(bold=True)
                c.alignment = Alignment(horizontal="center", vertical="center",
                                        wrap_text=True)
        for i in range(3, ws.max_row + 1):
            ws.cell(row=i, column=1).font = Font(bold=True)
            ws.cell(row=i, column=1).alignment = Alignment(horizontal="left")
        ws.column_dimensions["A"].width = 13
        ws.column_dimensions["B"].width = 22
        ws.column_dimensions["C"].width = 18
        ws.row_dimensions[1].height = 30
        for k in range(len(MODEL_DISPLAY)):
            c0 = 4 + 3 * k
            for off, w in enumerate((13, 10, 11)):
                ws.column_dimensions[get_column_letter(c0 + off)].width = w
            for row in ws.iter_rows(min_row=3, min_col=c0 + 1, max_col=c0 + 1):
                for c in row:
                    c.number_format = "0.0000"
            for row in ws.iter_rows(min_row=3, min_col=c0 + 2, max_col=c0 + 2):
                for c in row:
                    c.number_format = "0.00E+00"
            # a rule at the left edge of each model block after the first
            if k:
                for i in range(1, ws.max_row + 1):
                    cell = ws.cell(row=i, column=c0)
                    cell.border = Border(
                        left=Side(style="medium", color="808080"),
                        bottom=cell.border.bottom)

    # every quoted number must be recoverable from the sheet it is supposed to be in
    for sheet_name, ident, model, delta, q, site in QUOTED:
        _nd, n, d, qq = tables[sheet_name][(ident, model)]
        assert round(d, 3) == delta, "{} {} {} delta {} != {} ({})".format(
            sheet_name, ident, model, round(d, 3), delta, site)
        # ⛔ q IS MATCHED AGAINST EVERY PLAUSIBLE ROUNDING, NOT ONE FIXED PRECISION.
        # The manuscript does not print q uniformly: Results 76's Synapse LRRK2 is
        # 3.5e-5 (two significant figures), Results 77's M8 LRRK2 is 0.005 (one), and
        # M8 SNCA prints 0.002 from a source value of 0.0015075, one figure rounding
        # UP and 33 percent away in relative terms. Neither a fixed-digit test nor a
        # relative tolerance works, so the source value must round to the printed one
        # under some consistent convention.
        forms = {"1 sig fig": float("{:.1g}".format(qq)),
                 "2 sig fig": float("{:.2g}".format(qq)),
                 "3 dp": round(qq, 3),
                 "4 dp": round(qq, 4)}
        assert any(v == q for v in forms.values()), (
            "{} {} {} q {:.6g} does not round to {} ({}); tried {}".format(
                sheet_name, ident, model, qq, q, site,
                ", ".join("{}={}".format(k, v) for k, v in forms.items())))

    prose_basis = {}
    for ident, size in QUOTED_N:
        n_def, n_quant, _d, _q = tables["Individual_Modules"][(ident, "LRRK2")]
        basis = [name for name, val in (("defined", n_def),
                                        ("quantified", n_quant)) if val == size]
        assert basis, (
            "Results 77 states {} for {}, which is neither the definition size ({}) "
            "nor the quantified count ({})".format(size, ident, n_def, n_quant))
        prose_basis[ident] = basis

    out.save(dest)

    print("wrote " + os.path.relpath(dest, REPO_ROOT))
    print("2 sheets, WIDE: one row per module, three models side by side\n")
    for sheet_name, csv_name, _i, _ih, _lc, family, expect_rows in SHEETS:
        idents = sorted({k[0] for k in tables[sheet_name]})
        print("  {:24} {:>3} source rows -> {:>2} table rows   BH over {} tests".format(
            sheet_name, expect_rows, len(idents), family))
    shared = sorted({k[0] for k in tables["Individual_Modules"]}
                    & {k[0] for k in tables["Module_Groups"]})
    print("\n  in BOTH sheets, same delta and counts, DIFFERENT q: {}".format(
        ", ".join(shared)))
    print("  {} quoted delta/q pairs verified against their own sheet".format(
        len(QUOTED)))
    print("  {} prose module sizes recovered, each from the column it is:".format(
        len(QUOTED_N)))
    for ident, size in QUOTED_N:
        print("    {:5} Results 77 says {:>3}  = {}".format(
            ident, size, " and ".join(prose_basis[ident])))
    print("\n  dropped from each source: " + ", ".join(DROP))


if __name__ == "__main__":
    main()
