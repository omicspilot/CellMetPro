"""Advanced visualization types for metabolic analysis.

This module provides specialized visualizations including:
- Stacked bar charts for subsystem composition
- Ridge plots for distributions across groups
- Radar/spider plots for profile comparisons
- Waterfall plots for ranked differential results
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_stacked_bar(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    subsystem_mapping: dict[str, list[str]],
    n_top_subsystems: int = 10,
    normalize: bool = True,
    palette: str | list | None = None,
    figsize: tuple[float, float] = (12, 6),
    title: str = "Subsystem Composition by Group",
    xlabel: str = "Group",
    ylabel: str | None = None,
    legend_loc: str = "right",
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot stacked bar chart of subsystem composition per group.

    Shows the relative contribution of metabolic subsystems to
    overall activity in each cell group.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    subsystem_mapping : dict
        Mapping of subsystem names to lists of reaction IDs.
    n_top_subsystems : int
        Number of top subsystems to show (rest grouped as "Other").
    normalize : bool
        If True, normalize bars to 100%. If False, show absolute values.
    palette : str or list, optional
        Color palette for subsystems.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str, optional
        Y-axis label. Auto-set based on normalize.
    legend_loc : str
        Legend location: 'right', 'bottom', or 'none'.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to pandas.DataFrame.plot.bar.

    Returns
    -------
    Axes
        Matplotlib axes.

    Example
    -------
    >>> from cellmetpro.models import get_subsystem_reactions
    >>> subsystems = get_subsystem_reactions(model)
    >>> ax = plot_stacked_bar(scores, groups, subsystems)
    """
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores[common_cells]
    group_labels = groups[common_cells]

    # Compute subsystem-level scores per group
    subsystem_scores = {}
    for subsystem, reactions in subsystem_mapping.items():
        available = [r for r in reactions if r in scores.index]
        if len(available) > 0:
            subsystem_data = scores.loc[available]
            # Aggregate across reactions (mean) then by group
            cell_means = subsystem_data.mean(axis=0)
            group_means = cell_means.groupby(group_labels).mean()
            subsystem_scores[subsystem] = group_means

    if len(subsystem_scores) == 0:
        raise ValueError("No subsystems found with valid reactions")

    # Create DataFrame
    df = pd.DataFrame(subsystem_scores).T
    df = df.fillna(0)

    # Select top subsystems by overall activity
    total_activity = df.sum(axis=1)
    top_subsystems = total_activity.nlargest(n_top_subsystems).index.tolist()

    # Group remaining as "Other"
    other_subsystems = [s for s in df.index if s not in top_subsystems]
    if len(other_subsystems) > 0:
        other_row = df.loc[other_subsystems].sum()
        df = df.loc[top_subsystems]
        df.loc["Other"] = other_row
    else:
        df = df.loc[top_subsystems]

    # Transpose for plotting (groups as x-axis)
    plot_df = df.T

    # Normalize if requested
    if normalize:
        row_sums = plot_df.sum(axis=1)
        # Avoid division by zero - replace 0 sums with 1 to get 0% instead of NaN
        row_sums = row_sums.replace(0, 1)
        plot_df = plot_df.div(row_sums, axis=0) * 100

    # Create figure
    if legend_loc == "right":
        fig, ax = plt.subplots(figsize=figsize)
    elif legend_loc == "bottom":
        fig, ax = plt.subplots(figsize=(figsize[0], figsize[1] + 2))
    else:
        fig, ax = plt.subplots(figsize=figsize)

    # Set colors
    if palette is None:
        colors = sns.color_palette("tab20", len(plot_df.columns))
    elif isinstance(palette, str):
        colors = sns.color_palette(palette, len(plot_df.columns))
    else:
        colors = palette

    # Create stacked bar plot
    plot_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=colors,
        width=0.8,
        **kwargs,
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    else:
        ax.set_ylabel("Proportion (%)" if normalize else "Mean Activity")

    ax.tick_params(axis="x", rotation=45)

    # Handle legend
    if legend_loc == "right":
        ax.legend(
            title="Subsystem",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=8,
        )
    elif legend_loc == "bottom":
        ax.legend(
            title="Subsystem",
            bbox_to_anchor=(0.5, -0.15),
            loc="upper center",
            ncol=3,
            fontsize=8,
        )
    else:
        ax.legend().remove()

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_ridge(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 10,
    overlap: float = 0.5,
    palette: str | list | None = None,
    figsize: tuple[float, float] | None = None,
    title: str = "Reaction Score Distributions",
    xlabel: str = "Score",
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot ridge plot (joy plot) of reaction distributions across groups.

    Ridge plots show the distribution of scores for each reaction
    stacked vertically with slight overlap, making it easy to compare
    many distributions.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Specific reactions to plot. If None, uses top variable.
    n_top : int
        Number of top variable reactions if reactions not specified.
    overlap : float
        Amount of overlap between ridges (0-1).
    palette : str or list, optional
        Color palette for groups.
    figsize : tuple, optional
        Figure size.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to kdeplot.

    Returns
    -------
    Figure
        Matplotlib figure.

    Example
    -------
    >>> fig = plot_ridge(scores, groups, n_top=15)
    """

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select reactions
    if reactions is not None:
        plot_reactions = [r for r in reactions if r in reaction_scores.index]
    else:
        variances = reaction_scores.var(axis=1)
        plot_reactions = list(variances.nlargest(n_top).index)

    if len(plot_reactions) == 0:
        raise ValueError("No valid reactions to plot")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]

    # Prepare long format data
    data_long = scores.T.melt(
        var_name="reaction",
        value_name="score",
        ignore_index=False,
    )
    data_long["group"] = group_labels[data_long.index].values
    data_long = data_long.reset_index(drop=True)

    # Figure size
    n_reactions = len(plot_reactions)
    if figsize is None:
        figsize = (10, max(6, n_reactions * 0.8))

    # Set up figure with FacetGrid for ridge plot
    if palette is None:
        palette = "Set2"

    g = sns.FacetGrid(
        data_long,
        row="reaction",
        hue="group",
        aspect=3,
        height=0.7,
        palette=palette,
    )

    # Draw the densities
    g.map(
        sns.kdeplot,
        "score",
        bw_adjust=0.5,
        clip_on=False,
        fill=True,
        alpha=0.6,
        linewidth=1.5,
        **kwargs,
    )
    g.map(sns.kdeplot, "score", clip_on=False, color="w", lw=2, bw_adjust=0.5)

    # Overlap adjustment
    g.figure.subplots_adjust(hspace=-overlap)

    # Remove axes details for clean look
    g.set_titles("")
    g.set(yticks=[], ylabel="")
    g.despine(bottom=True, left=True)

    # Add reaction labels on the left
    for i, ax in enumerate(g.axes.flat):
        ax.text(
            -0.1,
            0.5,
            plot_reactions[i],
            transform=ax.transAxes,
            fontsize=9,
            ha="right",
            va="center",
        )

    # Add legend
    g.add_legend(title="Group")

    g.figure.suptitle(title, y=1.02)
    g.set_xlabels(xlabel)

    if save:
        g.savefig(save, dpi=150, bbox_inches="tight")

    return g.figure


def plot_radar(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 8,
    agg_method: Literal["mean", "median"] = "mean",
    normalize: bool = True,
    palette: str | list | None = None,
    figsize: tuple[float, float] = (8, 8),
    title: str = "Metabolic Profile Comparison",
    fill_alpha: float = 0.25,
    save: str | None = None,
) -> Axes:
    """Plot radar/spider chart comparing metabolic profiles across groups.

    Radar plots are useful for comparing multiple groups across
    several metabolic features simultaneously.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Specific reactions to plot. If None, uses top variable.
    n_top : int
        Number of top variable reactions if reactions not specified.
        Recommended: 5-12 for readability.
    agg_method : str
        Aggregation method: 'mean' or 'median'.
    normalize : bool
        Whether to normalize values to 0-1 range per reaction.
    palette : str or list, optional
        Color palette for groups.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    fill_alpha : float
        Transparency of filled area.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib polar axes.

    Example
    -------
    >>> ax = plot_radar(scores, groups, reactions=["R1", "R2", "R3", "R4", "R5"])
    """
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select reactions
    if reactions is not None:
        plot_reactions = [r for r in reactions if r in reaction_scores.index]
    else:
        variances = reaction_scores.var(axis=1)
        plot_reactions = list(variances.nlargest(n_top).index)

    if len(plot_reactions) < 3:
        raise ValueError("Need at least 3 reactions for radar plot")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]
    unique_groups = sorted(group_labels.unique())

    # Aggregate by group
    agg_func = np.mean if agg_method == "mean" else np.median
    group_data = {}
    for g in unique_groups:
        group_cells = group_labels[group_labels == g].index
        group_scores = scores[group_cells]
        group_data[g] = agg_func(group_scores, axis=1)

    df = pd.DataFrame(group_data)

    # Normalize if requested (per-reaction, so each axis has 0-1 range)
    if normalize:
        df_min = df.min(axis=1)
        df_max = df.max(axis=1)
        df_range = df_max - df_min
        # Avoid division by zero for constant reactions
        df_range = df_range.replace(0, 1)
        df = df.sub(df_min, axis=0).div(df_range, axis=0)

    # Set up radar chart
    categories = plot_reactions
    n_cats = len(categories)

    # Compute angle for each category
    angles = [n / float(n_cats) * 2 * np.pi for n in range(n_cats)]
    angles += angles[:1]  # Close the polygon

    # Create figure
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(polar=True))

    # Set colors
    if palette is None:
        colors = sns.color_palette("Set2", len(unique_groups))
    elif isinstance(palette, str):
        colors = sns.color_palette(palette, len(unique_groups))
    else:
        colors = palette

    # Plot each group
    for i, group in enumerate(unique_groups):
        values = df[group].values.tolist()
        values += values[:1]  # Close the polygon

        ax.plot(
            angles,
            values,
            linewidth=2,
            linestyle="solid",
            label=group,
            color=colors[i],
        )
        ax.fill(angles, values, alpha=fill_alpha, color=colors[i])

    # Set category labels
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, size=9)

    # Add legend and title
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.0))
    ax.set_title(title, size=12, y=1.08)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_waterfall(
    differential_results: pd.DataFrame,
    n_top: int = 30,
    log2fc_col: str = "log2fc",
    reaction_col: str = "reaction",
    pvalue_col: str = "padj_bh",
    pvalue_threshold: float = 0.05,
    up_color: str = "#E74C3C",
    down_color: str = "#3498DB",
    ns_color: str = "#95A5A6",
    figsize: tuple[float, float] | None = None,
    title: str = "Differential Reactions (Waterfall)",
    xlabel: str = "log2 Fold Change",
    show_labels: bool = True,
    label_fontsize: int = 8,
    save: str | None = None,
) -> Axes:
    """Plot waterfall chart of differential analysis results.

    Shows reactions ranked by log2 fold change with colors indicating
    direction and significance.

    Parameters
    ----------
    differential_results : pd.DataFrame
        Results from differential analysis with log2fc column.
    n_top : int
        Number of top reactions to show (by absolute log2fc).
    log2fc_col : str
        Column name for log2 fold change.
    reaction_col : str
        Column name for reaction IDs.
    pvalue_col : str
        Column name for adjusted p-values.
    pvalue_threshold : float
        Threshold for significance.
    up_color : str
        Color for upregulated reactions.
    down_color : str
        Color for downregulated reactions.
    ns_color : str
        Color for non-significant reactions.
    figsize : tuple, optional
        Figure size.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    show_labels : bool
        Whether to show reaction labels.
    label_fontsize : int
        Font size for labels.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib axes.

    Example
    -------
    >>> ax = plot_waterfall(diff_results, n_top=40)
    """
    import matplotlib.pyplot as plt

    # Validate input is not empty
    if len(differential_results) == 0:
        raise ValueError("differential_results DataFrame is empty")

    # Validate columns
    required_cols = [log2fc_col, reaction_col]
    for col in required_cols:
        if col not in differential_results.columns:
            raise ValueError(f"Column '{col}' not found in results")

    # Sort by absolute log2fc and take top n
    df = differential_results.copy()
    df["abs_log2fc"] = df[log2fc_col].abs()
    df = df.nlargest(n_top, "abs_log2fc").sort_values(log2fc_col, ascending=True)

    # Determine colors based on significance and direction
    colors = []
    for _, row in df.iterrows():
        is_sig = True
        if pvalue_col in df.columns:
            is_sig = row[pvalue_col] < pvalue_threshold

        if not is_sig:
            colors.append(ns_color)
        elif row[log2fc_col] > 0:
            colors.append(up_color)
        else:
            colors.append(down_color)

    # Figure size
    if figsize is None:
        figsize = (10, max(6, n_top * 0.25))

    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal bar plot
    y_pos = range(len(df))
    ax.barh(y_pos, df[log2fc_col], color=colors, edgecolor="white", height=0.8)

    # Add reaction labels
    if show_labels:
        ax.set_yticks(y_pos)
        ax.set_yticklabels(df[reaction_col], fontsize=label_fontsize)
    else:
        ax.set_yticks([])

    # Add vertical line at 0
    ax.axvline(x=0, color="black", linewidth=0.8)

    # Labels and title
    ax.set_xlabel(xlabel)
    ax.set_title(title)

    # Add legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=up_color, label="Up (significant)"),
        Patch(facecolor=down_color, label="Down (significant)"),
    ]
    if ns_color in colors:
        legend_elements.append(Patch(facecolor=ns_color, label="Not significant"))

    ax.legend(handles=legend_elements, loc="lower right")

    # Adjust layout
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_subsystem_waterfall(
    enrichment_results: pd.DataFrame,
    n_top: int = 20,
    fold_col: str = "fold_enrichment",
    term_col: str = "pathway",
    pvalue_col: str = "padj",
    pvalue_threshold: float = 0.05,
    sig_color: str = "#2ECC71",
    ns_color: str = "#95A5A6",
    figsize: tuple[float, float] | None = None,
    title: str = "Subsystem Enrichment",
    xlabel: str = "Fold Enrichment",
    save: str | None = None,
) -> Axes:
    """Plot waterfall chart for enrichment results.

    Similar to plot_waterfall but designed for pathway/subsystem
    enrichment results where all values are positive.

    Parameters
    ----------
    enrichment_results : pd.DataFrame
        Results from enrichment analysis.
    n_top : int
        Number of top terms to show.
    fold_col : str
        Column name for fold enrichment.
    term_col : str
        Column name for term/pathway names.
    pvalue_col : str
        Column name for adjusted p-values.
    pvalue_threshold : float
        Threshold for significance.
    sig_color : str
        Color for significant terms.
    ns_color : str
        Color for non-significant terms.
    figsize : tuple, optional
        Figure size.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib axes.
    """
    import matplotlib.pyplot as plt

    # Validate input is not empty
    if len(enrichment_results) == 0:
        raise ValueError("enrichment_results DataFrame is empty")

    # Validate columns
    if fold_col not in enrichment_results.columns:
        raise ValueError(f"Column '{fold_col}' not found in results")
    if term_col not in enrichment_results.columns:
        raise ValueError(f"Column '{term_col}' not found in results")

    # Sort by fold enrichment and take top n
    df = enrichment_results.copy()
    df = df.nlargest(n_top, fold_col).sort_values(fold_col, ascending=True)

    # Determine colors
    colors = []
    for _, row in df.iterrows():
        if pvalue_col in df.columns and row[pvalue_col] < pvalue_threshold:
            colors.append(sig_color)
        else:
            colors.append(ns_color)

    # Figure size
    if figsize is None:
        figsize = (10, max(6, n_top * 0.3))

    fig, ax = plt.subplots(figsize=figsize)

    y_pos = range(len(df))
    ax.barh(y_pos, df[fold_col], color=colors, edgecolor="white", height=0.8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(df[term_col], fontsize=9)

    ax.axvline(x=1, color="black", linewidth=0.8, linestyle="--")

    ax.set_xlabel(xlabel)
    ax.set_title(title)

    # Legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor=sig_color, label=f"Significant (p < {pvalue_threshold})"),
        Patch(facecolor=ns_color, label="Not significant"),
    ]
    ax.legend(handles=legend_elements, loc="lower right")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return ax
