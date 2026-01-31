"""Dot plot visualizations with statistical annotations.

This module provides dot plots for visualizing reaction scores
across groups with size indicating expression fraction and
color indicating mean activity.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes


def plot_reaction_dotplot(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 20,
    size_scale: Literal["fraction_active", "std", "n_cells"] = "fraction_active",
    color_scale: Literal["mean_score", "median_score"] = "mean_score",
    threshold: float = 0.5,
    cmap: str = "RdBu_r",
    size_range: tuple[float, float] = (10, 200),
    figsize: tuple[float, float] | None = None,
    title: str = "Reaction Activity Dotplot",
    xlabel: str = "Group",
    ylabel: str = "Reaction",
    show_legend: bool = True,
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Create dot plot of reaction activity.

    Dot plots show two dimensions of information:
    - Size: fraction of cells with activity below threshold (more active)
    - Color: mean/median activity score

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells). Lower scores = more active.
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Reactions to display. If None, uses top variable reactions.
    n_top : int
        Number of top variable reactions if reactions not specified.
    size_scale : str
        What dot size represents: 'fraction_active', 'std', or 'n_cells'.
    color_scale : str
        What color represents: 'mean_score' or 'median_score'.
    threshold : float
        Threshold for considering a reaction "active" (score below this).
    cmap : str
        Colormap for color scale.
    size_range : tuple
        Min and max dot sizes.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    show_legend : bool
        Whether to show size legend.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to scatter.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Select reactions
    if reactions is not None:
        plot_reactions = [r for r in reactions if r in reaction_scores.index]
    else:
        variances = reaction_scores.var(axis=1)
        plot_reactions = list(variances.nlargest(n_top).index)

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]

    # Compute statistics per group
    unique_groups = sorted(group_labels.unique())
    stats = []

    for reaction in plot_reactions:
        for group in unique_groups:
            group_cells = group_labels[group_labels == group].index
            reaction_values = scores.loc[reaction, group_cells]

            mean_score = reaction_values.mean()
            median_score = reaction_values.median()
            std_score = reaction_values.std()
            fraction_active = (reaction_values < threshold).mean()
            n_cells = len(group_cells)

            stats.append(
                {
                    "reaction": reaction,
                    "group": group,
                    "mean_score": mean_score,
                    "median_score": median_score,
                    "std_score": std_score,
                    "fraction_active": fraction_active,
                    "n_cells": n_cells,
                }
            )

    stats_df = pd.DataFrame(stats)

    # Map to numeric positions
    reaction_to_idx = {r: i for i, r in enumerate(plot_reactions)}
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}

    stats_df["x"] = stats_df["group"].map(group_to_idx)
    stats_df["y"] = stats_df["reaction"].map(reaction_to_idx)

    # Determine sizes
    if size_scale == "fraction_active":
        sizes = stats_df["fraction_active"]
    elif size_scale == "std":
        sizes = 1 - (stats_df["std_score"] / stats_df["std_score"].max())
    else:  # n_cells
        sizes = stats_df["n_cells"] / stats_df["n_cells"].max()

    # Normalize sizes to range
    sizes_norm = (sizes - sizes.min()) / (sizes.max() - sizes.min() + 1e-10)
    sizes_plot = size_range[0] + sizes_norm * (size_range[1] - size_range[0])

    # Determine colors
    if color_scale == "mean_score":
        colors = stats_df["mean_score"]
    else:
        colors = stats_df["median_score"]

    # Auto figure size
    if figsize is None:
        width = max(6, len(unique_groups) * 1.2 + 3)
        height = max(4, len(plot_reactions) * 0.4 + 2)
        figsize = (width, height)

    # Create plot
    fig, ax = plt.subplots(figsize=figsize)

    scatter = ax.scatter(
        stats_df["x"],
        stats_df["y"],
        s=sizes_plot,
        c=colors,
        cmap=cmap,
        alpha=0.8,
        edgecolors="white",
        linewidths=0.5,
        **kwargs,
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(color_scale.replace("_", " ").title())

    # Axis labels
    ax.set_xticks(range(len(unique_groups)))
    ax.set_xticklabels(unique_groups, rotation=45, ha="right")
    ax.set_yticks(range(len(plot_reactions)))
    ax.set_yticklabels(plot_reactions)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    # Size legend
    if show_legend:
        legend_sizes = [0.25, 0.5, 0.75, 1.0]
        legend_points = []
        legend_labels = []

        for ls in legend_sizes:
            size = size_range[0] + ls * (size_range[1] - size_range[0])
            legend_points.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor="gray",
                    markersize=np.sqrt(size),
                    label=f"{ls:.0%}",
                )
            )
            legend_labels.append(f"{ls:.0%}")

        ax.legend(
            handles=legend_points,
            title=size_scale.replace("_", " ").title(),
            loc="upper left",
            bbox_to_anchor=(1.15, 1),
            frameon=False,
        )

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_enrichment_dotplot(
    enrichment_results: pd.DataFrame,
    term_col: str = "go_term",
    name_col: str = "go_name",
    pvalue_col: str = "padj",
    fold_col: str = "fold_enrichment",
    n_top: int = 15,
    pvalue_threshold: float = 0.05,
    size_range: tuple[float, float] = (50, 300),
    cmap: str = "Reds",
    figsize: tuple[float, float] | None = None,
    title: str = "GO Enrichment",
    save: str | None = None,
) -> Axes:
    """Create dot plot for enrichment results.

    Parameters
    ----------
    enrichment_results : pd.DataFrame
        Enrichment results from GOEnrichmentAnalyzer.
    term_col : str
        Column with term IDs.
    name_col : str
        Column with term names.
    pvalue_col : str
        Column with p-values.
    fold_col : str
        Column with fold enrichment.
    n_top : int
        Number of top terms to show.
    pvalue_threshold : float
        Only show terms below this threshold.
    size_range : tuple
        Min and max dot sizes.
    cmap : str
        Colormap.
    figsize : tuple, optional
        Figure size.
    title : str
        Plot title.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib axes.
    """
    import matplotlib.pyplot as plt

    # Filter and sort
    df = enrichment_results[enrichment_results[pvalue_col] < pvalue_threshold].copy()
    df = df.nsmallest(n_top, pvalue_col)

    if len(df) == 0:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No significant terms", ha="center", va="center")
        return ax

    # Prepare labels
    if name_col in df.columns and df[name_col].notna().any():
        labels = df[name_col].fillna(df[term_col])
    else:
        labels = df[term_col]

    # Compute -log10 p-value for size
    neg_log_p = -np.log10(df[pvalue_col].clip(lower=1e-300))
    p_range = neg_log_p.max() - neg_log_p.min() + 1e-10
    sizes_norm = (neg_log_p - neg_log_p.min()) / p_range
    sizes = size_range[0] + sizes_norm * (size_range[1] - size_range[0])

    # Auto figure size
    if figsize is None:
        height = max(4, len(df) * 0.4 + 2)
        figsize = (8, height)

    fig, ax = plt.subplots(figsize=figsize)

    y_positions = range(len(df))

    scatter = ax.scatter(
        df[fold_col],
        y_positions,
        s=sizes,
        c=neg_log_p,
        cmap=cmap,
        alpha=0.8,
        edgecolors="black",
        linewidths=0.5,
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("-log10(adjusted p-value)")

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Fold Enrichment")
    ax.set_title(title)

    # Invert y-axis so most significant is at top
    ax.invert_yaxis()

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax
