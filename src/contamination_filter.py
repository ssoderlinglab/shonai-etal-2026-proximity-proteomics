"""
Contamination gene list applied before the co-fractionation module analyses.

The list holds 24 mouse gene symbols in four classes: astrocyte markers, blood and erythrocyte
proteins, epidermal and desmosome proteins, and other non-neuronal markers. Applying it takes the
LRRK2 partition background from 1957 to 1939 detected genes.

Import CONTAMINATION_GENES for membership tests, or call filter_contamination(df) to drop the rows.
Running this module prints the table.
"""

from __future__ import annotations
import pandas as pd

# ---------------------------------------------------------------------
_CONTAMINATION_TABLE = [
    # ---- Astrocyte markers ----
    ("Gfap",     "astrocyte", "GFAP, canonical astrocyte intermediate filament"),
    ("Vim",      "astrocyte", "Vimentin, astrocyte/glial marker"),
    ("Aldh1l1",  "astrocyte", "Aldh1l1, pan-astrocyte marker"),
    ("Slc1a3",   "astrocyte", "GLAST, astrocyte glutamate transporter"),
    ("S100b",    "astrocyte", "S100B, astrocyte marker"),
    ("Aqp4",     "astrocyte", "Aquaporin 4, astrocyte end-feet"),
    ("Aldoc",    "astrocyte", "Aldolase C, astrocyte-enriched"),
    ("Glul",     "astrocyte", "Glutamine synthetase, astrocyte marker"),
    # ---- Blood / erythrocyte ----
    ("Hbb-bs",   "blood",     "Hemoglobin beta, erythrocyte contamination"),
    ("Hbb-bt",   "blood",     "Hemoglobin beta-T, erythrocyte contamination"),
    ("Hba-a1",   "blood",     "Hemoglobin alpha 1, erythrocyte contamination"),
    ("Hba-a2",   "blood",     "Hemoglobin alpha 2, erythrocyte contamination"),
    ("Alb",      "blood",     "Serum albumin, plasma contamination"),
    # ---- Epidermal / desmosome ----
    ("Pkp1",     "desmosome", "Plakophilin 1, desmosome"),
    ("Dsp",      "desmosome", "Desmoplakin, desmosome"),
    ("Dsg1b",    "desmosome", "Desmoglein 1b, desmosome"),
    ("Jup",      "desmosome", "Junction plakoglobin, desmosome"),
    ("Cdsn",     "desmosome", "Corneodesmosin, cornified envelope"),
    ("Tgm1",     "desmosome", "Transglutaminase 1, keratinocyte"),
    ("Tgm3",     "desmosome", "Transglutaminase 3, keratinocyte"),
    ("Gsdma",    "desmosome", "Gasdermin A, epidermal"),
    ("Flg2",     "desmosome", "Filaggrin 2, cornified envelope"),
    # ---- Other non-neuronal ----
    ("Cldn11",   "other",     "Claudin 11, oligodendrocyte tight junction"),
    ("Lyz1",     "other",     "Lysozyme 1, myeloid"),
]

CONTAMINATION_DF: pd.DataFrame = pd.DataFrame(
    _CONTAMINATION_TABLE, columns=["Gene", "Class", "Rationale"]
)

CONTAMINATION_GENES: frozenset[str] = frozenset(CONTAMINATION_DF["Gene"])
EXCLUDE_GENES: frozenset[str] = CONTAMINATION_GENES


def filter_contamination(df: pd.DataFrame, gene_col: str = "Gene") -> pd.DataFrame:
    """
    Return a copy of `df` with contamination genes removed.

    Parameters
    ----------
    df : DataFrame containing one row per protein.
    gene_col : Name of the column holding the mouse gene symbol.

    """
    if gene_col not in df.columns:
        raise KeyError(
            f"Column {gene_col!r} not found. Available: {list(df.columns)}"
        )
    return df.loc[~df[gene_col].isin(CONTAMINATION_GENES)].copy()


def report() -> str:
    """Summary string for logging."""
    by_class = CONTAMINATION_DF.groupby("Class").size().to_dict()
    lines = [
        f"contamination filter: {len(CONTAMINATION_GENES)} genes",
    ]
    for cls, n in sorted(by_class.items()):
        lines.append(f"    {cls:<10} {n}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
    print()
    print(CONTAMINATION_DF.to_string(index=False))
