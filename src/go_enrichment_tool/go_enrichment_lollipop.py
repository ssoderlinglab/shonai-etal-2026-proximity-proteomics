#!/usr/bin/env python3
"""
GO Enrichment Lollipop Plot Generator (STRING-style)
=====================================================

Generates publication-quality lollipop plots for Gene Ontology enrichment
analysis, mimicking the STRING database visualization style.

Key features:
  - Signal score on X-axis (auto-calculated if missing)
  - Bubble size = gene count, bubble color = FDR
  - Jaccard-index-based similarity grouping with colored bands
  - Customizable themes (light/dark), color palettes, font sizes
  - SVG + PNG export with optional transparent background

Usage:
    from go_enrichment_lollipop import plot_enrichment
    plot_enrichment("enrichment.tsv", output_prefix="my_plot")

Author: Auto-generated for unified enrichment visualization pipeline
"""

import os
import textwrap
import warnings
from typing import Optional, Union

import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from config import (
    COLOR_PALETTES,
    THEMES,
    FONT_SIZES,
    PLOT_DEFAULTS,
    FILTER_DEFAULTS,
    COLUMN_ALIASES,
)


# =============================================================================
# Data Loading & Processing
# =============================================================================

def load_enrichment_data(filepath: str) -> pd.DataFrame:
    """
    Load enrichment TSV/CSV file and normalize column names.

    Parameters
    ----------
    filepath : str
        Path to TSV or CSV enrichment result file.

    Returns
    -------
    pd.DataFrame with standardized column names.
    """
    sep = "\t" if filepath.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(filepath, sep=sep)

    # Rename columns to standardized names
    rename_map = {}
    for std_name, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in df.columns and alias != std_name:
                rename_map[alias] = std_name
                break
    df = df.rename(columns=rename_map)

    # Ensure numeric types
    for col in ["observed_count", "background_count", "strength", "signal", "fdr"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def compute_signal(df: pd.DataFrame, w_strength: float = 1.0, w_fdr: float = 1.0) -> pd.Series:
    """
    Compute STRING-style 'signal' score.

    Signal is defined as the weighted harmonic mean of:
      - strength component: observed/expected ratio, normalized to [0, 1]
      - significance component: -log10(FDR), normalized to [0, 1]

    This matches STRING's definition: a balanced measure combining enrichment
    effect size and statistical significance.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain 'strength' (or 'observed_count' + 'background_count') and 'fdr'.
    w_strength : float
        Weight for the strength component (default 1.0).
    w_fdr : float
        Weight for the FDR component (default 1.0).

    Returns
    -------
    pd.Series of signal values in [0, 1].
    """
    # --- Strength component ---
    if "strength" in df.columns:
        strength_raw = df["strength"].copy()
    elif "observed_count" in df.columns and "background_count" in df.columns:
        # STRING strength = log10(observed / expected)
        # expected = background_count * (n_query / n_background_total)
        # Approximate: use observed/background as proxy
        ratio = df["observed_count"] / df["background_count"].replace(0, np.nan)
        strength_raw = np.log10(ratio.clip(lower=1e-10))
    else:
        warnings.warn("Cannot compute strength: missing columns. Setting signal to -log10(FDR) only.")
        neg_log_fdr = -np.log10(df["fdr"].clip(lower=1e-300))
        return (neg_log_fdr / neg_log_fdr.max()).clip(0, 1)

    # --- FDR component ---
    neg_log_fdr = -np.log10(df["fdr"].clip(lower=1e-300))

    # --- Normalize both to [0, 1] ---
    s_min, s_max = strength_raw.min(), strength_raw.max()
    if s_max > s_min:
        s_norm = (strength_raw - s_min) / (s_max - s_min)
    else:
        s_norm = pd.Series(0.5, index=df.index)

    f_min, f_max = neg_log_fdr.min(), neg_log_fdr.max()
    if f_max > f_min:
        f_norm = (neg_log_fdr - f_min) / (f_max - f_min)
    else:
        f_norm = pd.Series(0.5, index=df.index)

    # Shift to avoid zero in harmonic mean
    s_norm = s_norm.clip(lower=1e-6)
    f_norm = f_norm.clip(lower=1e-6)

    # --- Weighted harmonic mean ---
    total_weight = w_strength + w_fdr
    signal = total_weight / (w_strength / s_norm + w_fdr / f_norm)

    # Final normalization to [0, 1]
    sig_min, sig_max = signal.min(), signal.max()
    if sig_max > sig_min:
        signal = (signal - sig_min) / (sig_max - sig_min)
    signal = signal.clip(0, 1)

    return signal


# =============================================================================
# Jaccard Similarity Grouping
# =============================================================================

def _parse_gene_set(gene_str: str) -> set:
    """Parse comma-separated gene string into a set."""
    if pd.isna(gene_str) or gene_str == "":
        return set()
    return set(str(gene_str).split(","))


def compute_jaccard_matrix(df: pd.DataFrame, gene_col: str = "matching_ids") -> np.ndarray:
    """
    Compute pairwise Jaccard similarity matrix for gene sets.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain a column with comma-separated gene identifiers.
    gene_col : str
        Column name containing gene identifiers.

    Returns
    -------
    np.ndarray of shape (n_terms, n_terms).
    """
    # Try matching_ids first, fall back to matching_labels
    if gene_col not in df.columns:
        if "matching_labels" in df.columns:
            gene_col = "matching_labels"
        else:
            warnings.warn(f"Column '{gene_col}' not found. Skipping similarity grouping.")
            return np.zeros((len(df), len(df)))

    gene_sets = [_parse_gene_set(g) for g in df[gene_col]]
    n = len(gene_sets)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            if len(gene_sets[i]) == 0 or len(gene_sets[j]) == 0:
                sim = 0.0
            else:
                intersection = len(gene_sets[i] & gene_sets[j])
                union = len(gene_sets[i] | gene_sets[j])
                sim = intersection / union if union > 0 else 0.0
            matrix[i, j] = sim
            matrix[j, i] = sim
    return matrix


def filter_and_group_terms(
    df: pd.DataFrame,
    similarity_threshold: float = 0.5,
    max_terms: int = 10,
    sort_by: str = "signal",
) -> pd.DataFrame:
    """
    Select top terms and assign visual groups by Jaccard similarity.

    Algorithm (matching STRING's "Group terms by" behavior):
      1. Sort terms by the chosen metric (best first).
      2. Take the top N terms (no similarity-based exclusion).
      3. Assign group IDs using Jaccard similarity: terms with
         similarity >= threshold to an already-grouped term are
         placed in the same visual group.

    Parameters
    ----------
    df : pd.DataFrame
        Enrichment data (already filtered by category and FDR).
    similarity_threshold : float
        Jaccard threshold for visual grouping of similar terms.
    max_terms : int
        Maximum number of terms to display.
    sort_by : str
        Column to sort by: "signal", "fdr", "strength", "gene_count".

    Returns
    -------
    pd.DataFrame with added 'group_id' column, sorted with groups contiguous.
    """
    df = df.copy()

    # Ensure signal column exists
    if "signal" not in df.columns:
        df["signal"] = compute_signal(df)

    # Sort
    sort_map = {
        "signal": ("signal", False),
        "fdr": ("fdr", True),
        "strength": ("strength", False),
        "gene_count": ("observed_count", False),
    }
    sort_col, ascending = sort_map.get(sort_by, ("signal", False))
    df = df.sort_values(sort_col, ascending=ascending).reset_index(drop=True)

    # Take top N terms (no similarity-based exclusion)
    df_filtered = df.head(max_terms).copy().reset_index(drop=True)

    n = len(df_filtered)
    if n == 0:
        return df_filtered

    # Compute Jaccard matrix among selected terms
    jaccard_kept = compute_jaccard_matrix(df_filtered)

    # Assign group IDs: greedy grouping by similarity threshold
    group_ids = [-1] * n
    current_group = 0
    for i in range(n):
        if group_ids[i] == -1:
            group_ids[i] = current_group
            for j in range(i + 1, n):
                if group_ids[j] == -1 and jaccard_kept[i, j] >= similarity_threshold:
                    group_ids[j] = current_group
            current_group += 1

    df_filtered["group_id"] = group_ids

    # Re-sort: groups ordered by best term, within groups keep original order
    df_filtered["_group_best"] = df_filtered.groupby("group_id")[sort_col].transform(
        "min" if ascending else "max"
    )
    df_filtered = df_filtered.sort_values(
        ["_group_best", sort_col],
        ascending=[ascending, ascending],
    ).reset_index(drop=True)
    df_filtered = df_filtered.drop(columns=["_group_best"])

    return df_filtered


# =============================================================================
# Plotting
# =============================================================================

def _build_colormap(palette_name: str) -> matplotlib.colors.Colormap:
    """Create a colormap from palette name."""
    if palette_name not in COLOR_PALETTES:
        # Try as matplotlib builtin
        try:
            return plt.get_cmap(palette_name)
        except ValueError:
            warnings.warn(f"Unknown palette '{palette_name}'. Using 'mint - blue'.")
            palette_name = "mint - blue"

    low, high = COLOR_PALETTES[palette_name]
    if high is None:
        # Built-in matplotlib colormap
        return plt.get_cmap(low)
    return mcolors.LinearSegmentedColormap.from_list("custom", [low, high], N=256)


def plot_enrichment(
    data: Union[str, pd.DataFrame],
    # --- Filtering & Sorting ---
    category: Optional[str] = None,
    max_terms: int = None,
    sort_by: str = None,
    similarity_threshold: float = None,
    fdr_threshold: float = None,
    # --- Appearance ---
    theme: str = "light",
    color_palette: str = "mint - blue",
    transparent_bg: bool = False,
    font_sizes: Optional[dict] = None,
    # --- Layout ---
    fig_width: float = None,
    fig_height_per_term: float = None,
    bubble_scale: float = None,
    # --- Output ---
    output_prefix: str = "enrichment_plot",
    output_dir: str = ".",
    output_formats: list = None,
    dpi: int = None,
    # --- Advanced ---
    title: Optional[str] = None,
    x_label: str = "Signal",
    signal_weights: tuple = (1.0, 1.0),
    custom_theme: Optional[dict] = None,
    max_label_width: int = 30,
) -> plt.Figure:
    """
    Generate a STRING-style lollipop plot for GO enrichment analysis.

    Parameters
    ----------
    data : str or pd.DataFrame
        Path to enrichment TSV/CSV file, or pre-loaded DataFrame.
    category : str, optional
        Enrichment category to plot (e.g., "GO Process", "KEGG").
        If None, auto-selects the first GO category found.
    max_terms : int
        Maximum number of terms to display.
    sort_by : str
        Sort metric: "signal", "fdr", "strength", "gene_count".
    similarity_threshold : float
        Jaccard index threshold for visual grouping of similar terms.
    fdr_threshold : float
        Maximum FDR to include.
    theme : str
        "light" or "dark".
    color_palette : str
        Name of color palette for FDR coloring.
    transparent_bg : bool
        If True, output has transparent background.
    font_sizes : dict, optional
        Override specific font sizes. Keys: title, axis_label, tick_label,
        term_label, legend_title, legend_label.
    fig_width : float
        Figure width in inches.
    fig_height_per_term : float
        Figure height per term in inches.
    bubble_scale : float
        Multiplier for bubble sizes.
    output_prefix : str
        Prefix for output filenames.
    output_dir : str
        Directory for output files.
    output_formats : list
        List of formats to save (default: ["svg", "png"]).
    dpi : int
        Resolution for raster output.
    title : str, optional
        Custom plot title. Auto-generated if None.
    x_label : str
        Label for X-axis.
    signal_weights : tuple
        (w_strength, w_fdr) weights for signal calculation.
    custom_theme : dict, optional
        Override theme settings (same keys as THEMES entries).

    Returns
    -------
    matplotlib.figure.Figure
    """
    # --- Apply defaults ---
    max_terms = max_terms or FILTER_DEFAULTS["max_terms"]
    sort_by = sort_by or FILTER_DEFAULTS["sort_by"]
    similarity_threshold = similarity_threshold if similarity_threshold is not None else FILTER_DEFAULTS["similarity_threshold"]
    fdr_threshold = fdr_threshold or FILTER_DEFAULTS["fdr_threshold"]
    fig_width = fig_width or PLOT_DEFAULTS["fig_width"]
    fig_height_per_term = fig_height_per_term or PLOT_DEFAULTS["fig_height_per_term"]
    bubble_scale = bubble_scale or PLOT_DEFAULTS["bubble_scale"]
    dpi = dpi or PLOT_DEFAULTS["dpi"]
    output_formats = output_formats or ["svg", "png"]

    # Merge font sizes
    fs = FONT_SIZES.copy()
    if font_sizes:
        fs.update(font_sizes)

    # Theme
    theme_cfg = THEMES.get(theme, THEMES["light"]).copy()
    if custom_theme:
        theme_cfg.update(custom_theme)

    # --- Load data ---
    if isinstance(data, str):
        df = load_enrichment_data(data)
    else:
        df = data.copy()

    # --- Select category ---
    if "category" in df.columns:
        available_cats = df["category"].unique().tolist()
        if category is None:
            # Auto-select: prefer GO categories
            go_cats = [c for c in available_cats if "GO" in str(c)]
            category = go_cats[0] if go_cats else available_cats[0]
        df = df[df["category"] == category].copy()

    if len(df) == 0:
        warnings.warn(f"No data for category '{category}'.")
        return None

    # --- FDR filter ---
    if "fdr" in df.columns:
        df = df[df["fdr"] <= fdr_threshold].copy()

    if len(df) == 0:
        warnings.warn(f"No terms pass FDR threshold {fdr_threshold}.")
        return None

    # --- Compute signal if missing ---
    if "signal" not in df.columns:
        df["signal"] = compute_signal(df, w_strength=signal_weights[0], w_fdr=signal_weights[1])

    # --- Filter and group ---
    df_plot = filter_and_group_terms(
        df,
        similarity_threshold=similarity_threshold,
        max_terms=max_terms,
        sort_by=sort_by,
    )

    if len(df_plot) == 0:
        warnings.warn("No terms remaining after filtering.")
        return None

    n_terms = len(df_plot)

    # --- Prepare plot data ---
    terms = df_plot["term_description"].tolist()
    signals = df_plot["signal"].values
    gene_counts = df_plot["observed_count"].values
    fdrs = df_plot["fdr"].values
    group_ids = df_plot["group_id"].values

    # Reverse for bottom-to-top ordering (best at top)
    terms = terms[::-1]
    signals = signals[::-1]
    gene_counts = gene_counts[::-1]
    fdrs = fdrs[::-1]
    group_ids = group_ids[::-1]

    # Wrap long term labels (similar to STRING's two-line wrapping)
    # Characters before wrapping. Default 30 reproduces the submitted panels
    # byte-for-byte; callers may widen it per panel (see 20260805_fig1_N_BP_label_wrap).
    terms_wrapped = [textwrap.fill(t, width=max_label_width) for t in terms]

    # Compute y-positions with gaps between groups
    # Gap is relative to 1 term height (1.0 y-unit)
    TERM_HEIGHT = 1.0
    GROUP_GAP_RATIO = 0.1  # 10% of term height
    group_gap = TERM_HEIGHT * GROUP_GAP_RATIO
    y_positions = np.zeros(n_terms)
    current_y = 0.0
    for i in range(n_terms):
        if i > 0 and group_ids[i] != group_ids[i - 1]:
            current_y += group_gap  # Add gap at group boundary
        y_positions[i] = current_y
        current_y += TERM_HEIGHT

    # --- Bubble sizes ---
    gc_min, gc_max = gene_counts.min(), gene_counts.max()
    if gc_max > gc_min:
        size_norm = (gene_counts - gc_min) / (gc_max - gc_min)
    else:
        size_norm = np.full_like(gene_counts, 0.5, dtype=float)
    bubble_sizes = (
        PLOT_DEFAULTS["bubble_min_size"]
        + size_norm * (PLOT_DEFAULTS["bubble_max_size"] - PLOT_DEFAULTS["bubble_min_size"])
    ) * (bubble_scale / 8)

    # --- FDR colormap ---
    cmap = _build_colormap(color_palette)
    fdr_log = -np.log10(np.clip(fdrs, 1e-300, None))
    fdr_min, fdr_max = fdr_log.min(), fdr_log.max()
    if fdr_max > fdr_min:
        fdr_norm = (fdr_log - fdr_min) / (fdr_max - fdr_min)
    else:
        fdr_norm = np.full_like(fdr_log, 0.5)
    bubble_colors = cmap(fdr_norm)

    # --- Figure setup ---
    y_total_range = y_positions[-1] - y_positions[0] + 1 if n_terms > 0 else 1
    fig_height = max(PLOT_DEFAULTS["fig_min_height"], y_total_range * fig_height_per_term + 1.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if not transparent_bg:
        fig.patch.set_facecolor(theme_cfg["bg_color"])
        ax.set_facecolor(theme_cfg["bg_color"])
    else:
        fig.patch.set_alpha(0)
        ax.set_facecolor("none")

    font_family = PLOT_DEFAULTS["font_family"]
    if isinstance(font_family, list):
        plt.rcParams["font.sans-serif"] = font_family
        plt.rcParams["font.family"] = "sans-serif"
    else:
        plt.rcParams["font.family"] = font_family

    # --- Draw group background bands (color derived from representative term's FDR) ---
    unique_groups = []
    seen = set()
    for g in group_ids:
        if g not in seen:
            unique_groups.append(g)
            seen.add(g)

    # For each group, find the representative term (first/best in group) and
    # derive the band color from its FDR color, lightened for background use.
    group_repr_color = {}
    for gid in unique_groups:
        # Find first term in this group (representative = best ranked)
        for i, g in enumerate(group_ids):
            if g == gid:
                repr_color = bubble_colors[i][:3]  # RGB only
                # Lighten the color for background: blend toward white/dark bg
                if theme == "light":
                    band_rgb = repr_color * 0.25 + np.array([1, 1, 1]) * 0.75
                else:
                    band_rgb = repr_color * 0.3 + np.array(mcolors.to_rgb(theme_cfg["bg_color"])) * 0.7
                group_repr_color[gid] = band_rgb
                break

    # Draw bands per group as contiguous blocks (with gaps between groups)
    # Find contiguous ranges for each group in the reversed array
    group_spans = []  # [(y_start, y_end, gid), ...]
    i = 0
    while i < n_terms:
        gid = group_ids[i]
        y_start = y_positions[i] - 0.5
        y_end = y_positions[i] + 0.5
        while i + 1 < n_terms and group_ids[i + 1] == gid:
            i += 1
            y_end = y_positions[i] + 0.5
        group_spans.append((y_start, y_end, gid))
        i += 1

    for y_start, y_end, gid in group_spans:
        band_color = group_repr_color[gid]
        ax.axhspan(y_start, y_end, color=band_color, alpha=0.6, zorder=0)

    # --- Set X-axis limits first so stems extend from left edge ---
    x_data_min = signals.min()
    x_data_max = signals.max()
    x_range = x_data_max - x_data_min if x_data_max > x_data_min else 0.1
    x_left = x_data_min - x_range * 0.05
    x_right = x_data_max + x_range * 0.15
    ax.set_xlim(x_left, x_right)

    # --- Draw stems (lollipop lines from left edge to bubble) ---
    for i in range(n_terms):
        stem_color = bubble_colors[i]
        ax.plot(
            [x_left, signals[i]],
            [y_positions[i], y_positions[i]],
            color=stem_color,
            linewidth=PLOT_DEFAULTS["stem_linewidth"],
            alpha=theme_cfg["stem_alpha"],
            zorder=1,
            solid_capstyle="round",
        )

    # --- Draw bubbles ---
    scatter = ax.scatter(
        signals,
        y_positions,
        s=bubble_sizes,
        c=fdr_norm,
        cmap=cmap,
        edgecolors="white" if theme == "light" else "#444444",
        linewidths=0.5,
        zorder=2,
    )

    # --- Axis styling ---
    ax.set_yticks(y_positions)
    ax.set_yticklabels(terms_wrapped, fontsize=fs["term_label"], color=theme_cfg["text_color"])
    ax.set_xlabel(x_label, fontsize=fs["axis_label"], color=theme_cfg["text_color"])
    ax.set_ylim(y_positions[0] - 0.5, y_positions[-1] + 0.5)

    # X-axis limits already set above

    # Grid (x-axis only, subtle)
    ax.xaxis.grid(True, linestyle="--", linewidth=0.6, color=theme_cfg["grid_color"], alpha=0.7)
    ax.yaxis.grid(False)
    ax.set_axisbelow(True)

    # Tick styling
    ax.tick_params(axis="x", labelsize=fs["tick_label"], colors=theme_cfg["text_color"])
    ax.tick_params(axis="y", length=0)

    # Remove all spines (no x/y axis lines)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # --- Title ---
    if title is None:
        cat_display = category if category else "Enrichment"
        title = f"{cat_display} enrichment"
    ax.set_title(title, fontsize=fs["title"], color=theme_cfg["text_color"], pad=12, fontweight="bold")

    # --- FDR Colorbar Legend ---
    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=plt.Normalize(vmin=fdr_min, vmax=fdr_max),
    )
    sm.set_array([])

    # Position colorbar (use transform-based positioning for bbox_inches='tight')
    cbar_ax = fig.add_axes([0.86, 0.50, 0.015, 0.30])  # [left, bottom, width, height]
    cbar = fig.colorbar(sm, cax=cbar_ax)

    # FDR label ABOVE the colorbar
    cbar_ax.set_title(
        "FDR",
        fontsize=fs["legend_title"],
        color=theme_cfg["legend_text_color"],
        pad=6,
    )

    # Format colorbar ticks as scientific notation FDR values
    n_cbar_ticks = 5
    tick_positions = np.linspace(fdr_min, fdr_max, n_cbar_ticks)
    tick_labels = [f"{10**(-v):.1e}" for v in tick_positions]
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.ax.tick_params(labelsize=fs["legend_label"], colors=theme_cfg["legend_text_color"])
    cbar.outline.set_edgecolor(theme_cfg["axis_color"])
    cbar.outline.set_linewidth(0.5)
    if not transparent_bg:
        cbar_ax.set_facecolor(theme_cfg["bg_color"])

    # --- Gene Count Bubble Legend ---
    legend_counts = _choose_legend_sizes(gene_counts)
    legend_bubbles = []
    legend_labels_text = []
    for count in legend_counts:
        if gc_max > gc_min:
            sn = (count - gc_min) / (gc_max - gc_min)
        else:
            sn = 0.5
        size = (
            PLOT_DEFAULTS["bubble_min_size"]
            + sn * (PLOT_DEFAULTS["bubble_max_size"] - PLOT_DEFAULTS["bubble_min_size"])
        ) * (bubble_scale / 8)
        legend_bubbles.append(
            Line2D(
                [0], [0],
                marker="o",
                color="w",
                markerfacecolor="#888888",
                markersize=np.sqrt(size),
                markeredgecolor="white" if theme == "light" else "#444444",
                markeredgewidth=0.5,
                linestyle="None",
            )
        )
        legend_labels_text.append(str(count))

    gene_legend = ax.legend(
        legend_bubbles,
        legend_labels_text,
        title="Gene count",
        loc="lower right",
        bbox_to_anchor=(1.22, 0.0),
        frameon=False,
        fontsize=fs["legend_label"],
        title_fontsize=fs["legend_title"],
        labelspacing=1.2,
        handletextpad=0.5,
    )
    gene_legend.get_title().set_color(theme_cfg["legend_text_color"])
    for text in gene_legend.get_texts():
        text.set_color(theme_cfg["legend_text_color"])
    ax.add_artist(gene_legend)

    # --- Adjust layout ---
    plt.subplots_adjust(
        left=PLOT_DEFAULTS["margin_left"],
        right=0.82,
        top=0.92,
        bottom=0.1,
    )

    # --- Save ---
    os.makedirs(output_dir, exist_ok=True)
    for fmt in output_formats:
        outpath = os.path.join(output_dir, f"{output_prefix}.{fmt}")
        fig.savefig(
            outpath,
            format=fmt,
            dpi=dpi,
            transparent=transparent_bg,
            bbox_inches="tight",
        )
        print(f"Saved: {outpath}")

    return fig


def _choose_legend_sizes(gene_counts: np.ndarray, n_legend: int = 4) -> list:
    """Choose representative gene count values for the bubble size legend."""
    gc_min = int(gene_counts.min())
    gc_max = int(gene_counts.max())

    if gc_max == gc_min:
        return [gc_min]

    # Use "nice" round numbers
    raw_steps = np.linspace(gc_min, gc_max, n_legend)
    nice = []
    for v in raw_steps:
        v = int(v)
        if v <= 5:
            nice.append(v)
        elif v <= 20:
            nice.append(round(v / 5) * 5)
        elif v <= 100:
            nice.append(round(v / 10) * 10)
        else:
            nice.append(round(v / 50) * 50)

    # Remove duplicates, keep order
    seen = set()
    result = []
    for v in nice:
        if v not in seen:
            seen.add(v)
            result.append(v)
    return result


# =============================================================================
# Convenience Functions
# =============================================================================

def list_categories(data: Union[str, pd.DataFrame]) -> list:
    """
    List available enrichment categories in the data.

    Parameters
    ----------
    data : str or pd.DataFrame

    Returns
    -------
    list of category names.
    """
    if isinstance(data, str):
        df = load_enrichment_data(data)
    else:
        df = data
    if "category" in df.columns:
        return df["category"].unique().tolist()
    return []


def plot_all_categories(
    data: Union[str, pd.DataFrame],
    output_dir: str = ".",
    **kwargs,
) -> dict:
    """
    Generate plots for all enrichment categories in the data.

    Parameters
    ----------
    data : str or pd.DataFrame
    output_dir : str
    **kwargs : passed to plot_enrichment()

    Returns
    -------
    dict mapping category name to Figure.
    """
    categories = list_categories(data)
    figures = {}
    for cat in categories:
        safe_name = cat.replace(" ", "_").replace("/", "-")
        prefix = kwargs.pop("output_prefix", f"enrichment_{safe_name}")
        fig = plot_enrichment(
            data,
            category=cat,
            output_dir=output_dir,
            output_prefix=prefix,
            **kwargs,
        )
        figures[cat] = fig
        kwargs["output_prefix"] = kwargs.get("output_prefix", prefix)  # Restore

    return figures


# =============================================================================
# CLI Entry Point
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate STRING-style GO enrichment lollipop plots."
    )
    parser.add_argument("input", help="Path to enrichment TSV/CSV file")
    parser.add_argument("-o", "--output-dir", default=".", help="Output directory")
    parser.add_argument("-p", "--prefix", default="enrichment_plot", help="Output filename prefix")
    parser.add_argument("-c", "--category", default=None, help="Enrichment category to plot")
    parser.add_argument("--all-categories", action="store_true", help="Plot all categories")
    parser.add_argument("-n", "--max-terms", type=int, default=10, help="Max terms to show")
    parser.add_argument("-s", "--sort-by", default="signal",
                        choices=["signal", "fdr", "strength", "gene_count"])
    parser.add_argument("--similarity", type=float, default=0.5,
                        help="Jaccard similarity threshold for grouping")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--theme", default="light", choices=["light", "dark"])
    parser.add_argument("--palette", default="mint - blue", help="Color palette name")
    parser.add_argument("--transparent", action="store_true", help="Transparent background")
    parser.add_argument("--dpi", type=int, default=150, help="Output DPI")
    parser.add_argument("--formats", nargs="+", default=["svg", "png"],
                        help="Output formats (e.g., svg png pdf)")

    args = parser.parse_args()

    if args.all_categories:
        plot_all_categories(
            args.input,
            output_dir=args.output_dir,
            max_terms=args.max_terms,
            sort_by=args.sort_by,
            similarity_threshold=args.similarity,
            fdr_threshold=args.fdr,
            theme=args.theme,
            color_palette=args.palette,
            transparent_bg=args.transparent,
            dpi=args.dpi,
            output_formats=args.formats,
        )
    else:
        plot_enrichment(
            args.input,
            category=args.category,
            output_dir=args.output_dir,
            output_prefix=args.prefix,
            max_terms=args.max_terms,
            sort_by=args.sort_by,
            similarity_threshold=args.similarity,
            fdr_threshold=args.fdr,
            theme=args.theme,
            color_palette=args.palette,
            transparent_bg=args.transparent,
            dpi=args.dpi,
            output_formats=args.formats,
        )
