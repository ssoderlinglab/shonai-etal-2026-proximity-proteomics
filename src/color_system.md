# Color System — Parkinson Proteomics Manuscript

> **Version:** 1.0
> **Created:** 2026-03-15
> **Purpose:** Unified color definitions for all figures, scripts, and visualizations
> **Rule:** All new figures MUST reference this file. No ad-hoc color definitions.

---

## 1. Genotype (WT vs MUT)

| Label | Hex | RGB | Usage |
|-------|-----|-----|-------|
| WT | `#4895EF` | (72, 149, 239) | All WT-associated elements: nodes, bars, points, edges |
| MUT | `#E63946` | (230, 57, 70) | All MUT-associated elements |
| WT border (network) | `#1A5FAD` | (26, 95, 173) | Bait node border in network figures |
| MUT border (network) | `#9D0208` | (157, 2, 8) | Bait node border in network figures |

**Context:** Blue/Red is the most colorblind-safe pair. Maintained from existing network and volcano figures.

---

## 2. PD Model Identity (SNCA / VPS35 / LRRK2)

Source palette: **ColorBrewer Dark2**

| Model | Hex | RGB | Dark2 Position |
|-------|-----|-----|----------------|
| SNCA | `#D95F02` | (217, 95, 2) | Dark2-2 (orange) |
| VPS35 | `#7570B3` | (117, 112, 179) | Dark2-3 (purple-violet) |
| LRRK2 | `#1B9E77` | (27, 158, 119) | Dark2-1 (teal-green) |

**Usage:** PCA plots, model comparison panels, any figure where model identity (not genotype) is the variable.

**Note:** When both genotype and model identity appear in the same figure (e.g., PCA with WT and MUT from 3 models), use model color for fill and genotype for shape or border.

---

## 3. Network Functional Categories

Source palette: **ColorBrewer Set3** (selected and reassigned)

| Category | Hex | RGB | Set3 Name |
|----------|-----|-----|-----------|
| Presynaptic | `#FB8072` | (251, 128, 114) | salmon |
| Postsynaptic/PSD | `#FCCDE5` | (252, 205, 229) | light pink |
| Membrane trafficking/Clathrin/Retromer | `#B3DE69` | (179, 222, 105) | light green |
| Cytoskeleton | `#CCEBC5` | (204, 235, 197) | very light green |
| Signaling/GTPase/Kinase | `#BC80BD` | (188, 128, 189) | medium purple |
| RNA processing/Spliceosome | `#FFED6F` | (255, 237, 111) | gold |
| Protein homeostasis/Ubiquitin | `#8DD3C7` | (141, 211, 199) | light teal |
| Other/Uncharacterized | `#D9D9D9` | (217, 217, 217) | grey |

**Note:** Set3 is a pastel palette. Network categories are "secondary information" — intentionally lower saturation than Genotype/Model colors to create visual hierarchy.

**Known proximity to reserved colors:** Protein homeostasis teal ≈ LRRK2 teal, Signaling purple ≈ VPS35 purple. Accepted because network and model identity rarely appear in the same panel, and pastel vs. saturated provides sufficient distinction.

---

## 4. Cell Type Identity

| Cell Type | Hex | RGB | Notes |
|-----------|-----|-----|-------|
| Dopaminergic Neurons (DAN) | `#E74C3C` | (231, 76, 60) | Highlighted — primary focus |
| Non-DAN (all others) | `#95A5A6` | (149, 165, 166) | Grey — de-emphasized |

**Usage:** DAN enrichment figures (Step 1 CS-iLOPIT-DC), cell-type bar plots.

**If individual non-DAN types need distinction (Extended Data only):**

| Cell Type | Hex | Source |
|-----------|-----|--------|
| Glutamatergic | `#95A5A6` | grey (medium) |
| GABAergic | `#BDC3C7` | grey (light) |
| Astrocytes | `#7F8C8D` | grey (dark) |
| Oligodendrocytes | `#ABB2B9` | grey (medium-light) |
| Microglia | `#5D6D7E` | grey (darkest) |

---

## 5. Heatmaps

### Sequential (one-directional: enrichment, intensity, expression)

| Context | Colormap | Notes |
|---------|----------|-------|
| Default sequential | `Reds` | matplotlib/seaborn `Reds` |
| DAN enrichment | `Reds` | Used in Step 1 heatmap |
| Protein intensity | `Reds` | General abundance |

### Diverging (bi-directional: enrichment vs. depletion, fold-change)

| Context | Colormap | Notes |
|---------|----------|-------|
| Default diverging | `PRGn` | ColorBrewer Purple-Green. Avoids red/blue (reserved for genotype). |
| Fold-change (MUT vs WT) | `PRGn` | Purple = enriched in MUT, Green = depleted in MUT |
| CrossComparison heatmap | `PRGn` | Step 4 CS-iLOPIT-DC |

**Rationale for PRGn:** Red-Blue diverging (RdBu) would conflict with Genotype colors (WT=blue, MUT=red). Purple-Green occupies a separate color space.

---

## 6. Volcano Plots

| Category | Hex | RGB | Notes |
|----------|-----|-----|-------|
| Significantly Up | `#D62728` | (214, 39, 40) | matplotlib tab:red |
| Significantly Down | `#1F77B4` | (31, 119, 180) | matplotlib tab:blue |
| Not Significant | `#CCCCCC` | (204, 204, 204) | light grey |
| Bait | `#2CA02C` | (44, 160, 44) | matplotlib tab:green |

---

## 7. Edge Colors (Network)

| Edge Type | Hex | Notes |
|-----------|-----|-------|
| STRING database | `#AAAAAA` | grey, background |
| BioID WT | `#4895EF` | matches WT genotype color |
| BioID MUT | `#E63946` | matches MUT genotype color |

---

## 8. Fractionation Profiles (Differential Centrifugation)

Source palette: **ColorBrewer Set3** (qualitative, 7 high-visibility colors)

| Fraction | Hex | RGB | Set3 Name |
|----------|-----|-----|-----------|
| F4 (heaviest) | `#FB8072` | (251, 128, 114) | salmon |
| F5 | `#80B1D3` | (128, 177, 211) | light blue |
| F6 | `#B3DE69` | (179, 222, 105) | light green |
| F7 | `#FDB462` | (253, 180, 98) | light orange |
| F8 | `#BC80BD` | (188, 128, 189) | medium purple |
| F9 | `#8DD3C7` | (141, 211, 199) | light teal |
| F10 (lightest) | `#FFED6F` | (255, 237, 111) | light gold |

**Usage:** PCA scatter plots, fractionation profile line plots, fraction bar charts.

**Note:** These colors overlap with Network Functional Categories (same Set3 source). This is acceptable because fraction colors and network colors never appear in the same figure panel.

---

## Reserved Color Summary

These hex codes are "claimed" and must not be reused for other purposes:

| Hex | Assigned To |
|-----|-------------|
| `#4895EF` | WT genotype |
| `#E63946` | MUT genotype |
| `#D95F02` | SNCA model |
| `#7570B3` | VPS35 model |
| `#1B9E77` | LRRK2 model |
| `#E74C3C` | DAN cell type |

---

## Python Config Snippet

```python
# Copy this into any analysis script
COLOR_GENOTYPE = {"WT": "#4895EF", "MUT": "#E63946"}
COLOR_MODEL = {"SNCA": "#D95F02", "VPS35": "#7570B3", "LRRK2": "#1B9E77"}
COLOR_CELLTYPE = {"DAN": "#E74C3C", "non-DAN": "#95A5A6"}
COLOR_NETWORK = {
    "Presynaptic": "#FB8072",
    "Postsynaptic": "#FCCDE5",
    "Trafficking": "#B3DE69",
    "Cytoskeleton": "#CCEBC5",
    "Signaling": "#BC80BD",
    "RNA_processing": "#FFED6F",
    "Protein_homeostasis": "#8DD3C7",
    "Other": "#D9D9D9",
}
COLOR_VOLCANO = {"Up": "#D62728", "Down": "#1F77B4", "NS": "#CCCCCC", "Bait": "#2CA02C"}
COLOR_FRACTION = {
    "F4": "#FB8072", "F5": "#80B1D3", "F6": "#B3DE69",
    "F7": "#FDB462", "F8": "#BC80BD", "F9": "#8DD3C7", "F10": "#FFED6F",
}
CMAP_SEQUENTIAL = "Reds"
CMAP_DIVERGING = "PRGn"
```

---

## Accessibility Notes

- **Colorblind safety:** Dark2 (PD models) and Set3 (network) are both ColorBrewer certified colorblind-safe palettes. WT/MUT blue/red is the most universally distinguishable pair.
- **Print safety:** All palettes are designed to remain distinguishable in greyscale print (varying lightness values).
- **PRGn diverging:** Safe for protanopia and deuteranopia (most common forms of color vision deficiency).

---

*This file is the single source of truth for all color definitions. Update version number when changing any color.*
