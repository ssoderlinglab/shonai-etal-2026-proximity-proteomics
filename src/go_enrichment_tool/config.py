"""
Default configuration for GO enrichment lollipop plot.
All settings can be overridden via function arguments.
"""

# =============================================================================
# Color Palettes
# =============================================================================
# Each palette is defined as (low_color, high_color) for FDR colormap.
# low_color = most significant (smallest FDR), high_color = least significant.
COLOR_PALETTES = {
    "mint - blue":    ("#b2dfdb", "#01579b"),  # STRING default
    "mint - green":   ("#c8e6c9", "#1b5e20"),
    "yellow - red":   ("#fff9c4", "#b71c1c"),
    "orange - brown": ("#ffe0b2", "#3e2723"),
    "pink - purple":  ("#f8bbd0", "#4a148c"),
    "grey":           ("#e0e0e0", "#212121"),
    "viridis":        ("viridis", None),       # Use matplotlib built-in
    "plasma":         ("plasma", None),
    "cividis":        ("cividis", None),
}

# =============================================================================
# Theme Presets
# =============================================================================
THEMES = {
    "light": {
        "bg_color": "#ffffff",
        "text_color": "#333333",
        "axis_color": "#333333",
        "grid_color": "#d3d3d3",
        "stem_alpha": 0.6,
        "group_band_colors": ["#d6eaf8", "#e8f4f8"],  # Fallback; actual bands derived from FDR color
        "legend_text_color": "#333333",
    },
    "dark": {
        "bg_color": "#1e1e1e",
        "text_color": "#e0e0e0",
        "axis_color": "#e0e0e0",
        "grid_color": "#444444",
        "stem_alpha": 0.5,
        "group_band_colors": ["#2a3a4a", "#1e2e3e"],  # Fallback; actual bands derived from FDR color
        "legend_text_color": "#e0e0e0",
    },
}

# =============================================================================
# Font Sizes (all customizable)
# =============================================================================
FONT_SIZES = {
    "title": 12,
    "axis_label": 10,
    "tick_label": 9,
    "term_label": 9,
    "legend_title": 9,
    "legend_label": 8,
}

# =============================================================================
# Plot Layout
# =============================================================================
PLOT_DEFAULTS = {
    "fig_width": 10,
    "fig_height_per_term": 0.35,
    "fig_min_height": 4,
    "dpi": 150,
    "font_family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "bubble_scale": 8,         # Multiplier for bubble sizes
    "bubble_min_size": 30,     # Minimum bubble area in points^2
    "bubble_max_size": 600,    # Maximum bubble area in points^2
    "stem_linewidth": 1.5,
    "margin_left": 0.28,      # Fraction of figure width for term labels
    "margin_right": 0.15,     # Fraction for legends
}

# =============================================================================
# Filtering & Sorting
# =============================================================================
FILTER_DEFAULTS = {
    "max_terms": 10,
    "fdr_threshold": 0.05,
    "sort_by": "signal",       # "signal", "fdr", "strength", "gene_count"
    "similarity_threshold": 0.5,  # Jaccard index threshold for grouping
    "category": None,          # None = auto-detect first GO category
}

# =============================================================================
# Column Name Mapping
# =============================================================================
# Maps expected column names to possible aliases in input TSV.
COLUMN_ALIASES = {
    "category":          ["#category", "category"],
    "term_id":           ["term ID", "term_id", "GO_ID", "ID"],
    "term_description":  ["term description", "term_description", "Description", "Term"],
    "observed_count":    ["observed gene count", "observed_count", "Count", "gene_count"],
    "background_count":  ["background gene count", "background_count", "background"],
    "strength":          ["strength", "log2_enrichment", "enrichment"],
    "signal":            ["signal"],
    "fdr":               ["false discovery rate", "fdr", "FDR", "p.adjust", "padj", "q-value"],
    "matching_ids":      ["matching proteins in your network (IDs)", "matching_ids", "geneID", "Genes"],
    "matching_labels":   ["matching proteins in your network (labels)", "matching_labels", "gene_labels"],
}
