"""Heatmap visualizations for pathway and reaction activity.

This module provides heatmap plots comparing metabolic activity
across cell clusters, conditions, or samples.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_pathway_heatmap(
    pathway_scores: pd.DataFrame,
    groups: pd.Series,
    pathways: list[str] | None = None,
    n_top: int = 30,
    agg_method: Literal["mean", "median"] = "mean",
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    scale: bool = True,
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (10, 8),
    annot: bool = False,
    title: str = "Pathway Activity by Group",
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot heatmap of pathway activity across groups.

    Parameters
    ----------
    pathway_scores : pd.DataFrame
        Pathway scores (pathways x cells).
    groups : pd.Series
        Group labels for cells.
    pathways : list[str], optional
        Specific pathways to show. If None, shows top variable.
    n_top : int
        Number of top variable pathways if pathways not specified.
    agg_method : str
        Aggregation method for groups.
    cluster_rows : bool
        Whether to cluster pathways.
    cluster_cols : bool
        Whether to cluster groups.
    scale : bool
        Whether to z-score normalize rows.
    cmap : str
        Colormap.
    figsize : tuple
        Figure size.
    annot : bool
        Whether to show values in cells.
    title : str
        Plot title.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to clustermap.

    Returns
    -------
    Figure
        Matplotlib figure with the plot.
    """

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select pathways
    if pathways is not None:
        plot_data = pathway_scores.loc[pathway_scores.index.isin(pathways)]
    else:
        variances = pathway_scores.var(axis=1)
        top_pathways = variances.nlargest(n_top).index
        plot_data = pathway_scores.loc[top_pathways]

    # Get common cells
    common_cells = plot_data.columns.intersection(groups.index)
    plot_data = plot_data[common_cells]
    group_labels = groups[common_cells]

    # Aggregate by group
    if agg_method == "mean":
        grouped = plot_data.T.groupby(group_labels).mean().T
    else:
        grouped = plot_data.T.groupby(group_labels).median().T

    # Scale if requested
    if scale:
        grouped = grouped.subtract(grouped.mean(axis=1), axis=0)
        grouped = grouped.divide(grouped.std(axis=1).clip(lower=1e-10), axis=0)

    # Create clustermap
    g = sns.clustermap(
        grouped,
        cmap=cmap,
        center=0 if scale else None,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        figsize=figsize,
        annot=annot,
        fmt=".2f" if annot else None,
        **kwargs,
    )

    g.fig.suptitle(title, y=1.02)
    g.ax_cbar.set_ylabel("Z-score" if scale else "Score", rotation=270, labelpad=15)

    if save:
        g.savefig(save, dpi=150, bbox_inches="tight")

    return g.fig  # type: ignore[no-any-return]


def plot_reaction_heatmap(
    reaction_scores: pd.DataFrame,
    groups: pd.Series | None = None,
    reactions: list[str] | None = None,
    n_top: int = 50,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    scale: Literal["row", "column", None] = "row",
    cmap: str = "RdBu_r",
    vmin: float | None = None,
    vmax: float | None = None,
    figsize: tuple[float, float] | None = None,
    show_row_labels: bool = True,
    show_col_labels: bool = False,
    row_label_size: int = 8,
    title: str = "Reaction Activity Heatmap",
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot heatmap of reaction activity across cells.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series, optional
        Group labels for cells. Used for column coloring.
    reactions : list[str], optional
        Specific reactions to show.
    n_top : int
        Number of top variable reactions if reactions not specified.
    cluster_rows : bool
        Whether to cluster reactions.
    cluster_cols : bool
        Whether to cluster cells.
    scale : str or None
        Scale data by 'row', 'column', or None.
    cmap : str
        Colormap.
    vmin : float, optional
        Minimum value for color scale.
    vmax : float, optional
        Maximum value for color scale.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    show_row_labels : bool
        Whether to show reaction names.
    show_col_labels : bool
        Whether to show cell names.
    row_label_size : int
        Font size for row labels.
    title : str
        Plot title.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to clustermap.

    Returns
    -------
    Figure
        Matplotlib figure with the plot.
    """

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select reactions
    if reactions is not None:
        plot_data = reaction_scores.loc[reaction_scores.index.isin(reactions)]
    else:
        variances = reaction_scores.var(axis=1)
        top_reactions = variances.nlargest(n_top).index
        plot_data = reaction_scores.loc[top_reactions]

    # Handle NaN
    plot_data = plot_data.fillna(0)

    # Scale data
    if scale == "row":
        plot_data = plot_data.subtract(plot_data.mean(axis=1), axis=0)
        plot_data = plot_data.divide(plot_data.std(axis=1).clip(lower=1e-10), axis=0)
    elif scale == "column":
        plot_data = plot_data.subtract(plot_data.mean(axis=0), axis=1)
        plot_data = plot_data.divide(plot_data.std(axis=0).clip(lower=1e-10), axis=1)

    # Auto figure size
    if figsize is None:
        n_rows, n_cols = plot_data.shape
        width = max(8, min(20, n_cols * 0.1 + 3))
        height = max(6, min(15, n_rows * 0.15 + 2))
        figsize = (width, height)

    # Prepare column colors
    col_colors = None
    if groups is not None:
        common_cells = plot_data.columns.intersection(groups.index)
        if len(common_cells) > 0:
            plot_data = plot_data[common_cells]
            group_labels = groups[common_cells]

            unique_groups = group_labels.unique()
            palette = sns.color_palette("husl", len(unique_groups))
            color_map = dict(zip(unique_groups, palette))
            col_colors = group_labels.map(color_map)

    # Create clustermap
    g = sns.clustermap(
        plot_data,
        cmap=cmap,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        col_colors=col_colors,
        figsize=figsize,
        vmin=vmin,
        vmax=vmax,
        xticklabels=show_col_labels,
        yticklabels=show_row_labels,
        **kwargs,
    )

    g.ax_heatmap.set_ylabel("")
    g.ax_heatmap.set_xlabel("")

    if show_row_labels:
        g.ax_heatmap.tick_params(axis="y", labelsize=row_label_size)

    g.fig.suptitle(title, y=1.02)
    g.ax_cbar.set_ylabel("Z-score" if scale else "Score", rotation=270, labelpad=15)

    if save:
        g.savefig(save, dpi=150, bbox_inches="tight")

    return g.fig  # type: ignore[no-any-return]


def plot_grouped_heatmap(
    data: pd.DataFrame,
    groups: pd.Series,
    features: list[str] | None = None,
    n_top: int = 30,
    agg_method: Literal["mean", "median"] = "mean",
    scale: bool = True,
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (10, 8),
    annot: bool = False,
    title: str = "Grouped Heatmap",
    xlabel: str = "Group",
    ylabel: str = "Feature",
    save: str | None = None,
) -> Axes:
    """Create heatmap of scores aggregated by group.

    Parameters
    ----------
    data : pd.DataFrame
        Scores (features x cells).
    groups : pd.Series
        Group labels for cells.
    features : list[str], optional
        Specific features to include.
    n_top : int
        Number of top variable features.
    agg_method : str
        Aggregation method: 'mean' or 'median'.
    scale : bool
        Whether to z-score normalize rows.
    cmap : str
        Colormap.
    figsize : tuple
        Figure size.
    annot : bool
        Whether to annotate cells with values.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib axes.
    """
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select features
    if features is not None:
        plot_data = data.loc[data.index.isin(features)]
    else:
        variances = data.var(axis=1)
        top_features = variances.nlargest(n_top).index
        plot_data = data.loc[top_features]

    # Get common cells
    common_cells = plot_data.columns.intersection(groups.index)
    plot_data = plot_data[common_cells]
    group_labels = groups[common_cells]

    # Aggregate by group
    if agg_method == "mean":
        grouped = plot_data.T.groupby(group_labels).mean().T
    else:
        grouped = plot_data.T.groupby(group_labels).median().T

    # Scale if requested
    if scale:
        grouped = grouped.subtract(grouped.mean(axis=1), axis=0)
        grouped = grouped.divide(grouped.std(axis=1).clip(lower=1e-10), axis=0)

    # Create heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        grouped,
        cmap=cmap,
        center=0 if scale else None,
        annot=annot,
        fmt=".2f" if annot else None,
        ax=ax,
        cbar_kws={"label": "Z-score" if scale else "Score"},
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_correlation_heatmap(
    data: pd.DataFrame,
    method: Literal["pearson", "spearman", "kendall"] = "pearson",
    features: list[str] | None = None,
    n_top: int = 50,
    cluster: bool = True,
    cmap: str = "coolwarm",
    figsize: tuple[float, float] = (10, 10),
    title: str = "Feature Correlation",
    save: str | None = None,
) -> Figure:
    """Create correlation heatmap between features.

    Parameters
    ----------
    data : pd.DataFrame
        Scores (features x samples).
    method : str
        Correlation method.
    features : list[str], optional
        Specific features to include.
    n_top : int
        Number of most variable features to include.
    cluster : bool
        Whether to cluster rows/columns.
    cmap : str
        Colormap.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.
    """

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Select features
    if features is not None:
        plot_data = data.loc[data.index.isin(features)]
    else:
        variances = data.var(axis=1)
        top_features = variances.nlargest(n_top).index
        plot_data = data.loc[top_features]

    # Compute correlation
    corr = plot_data.T.corr(method=method)

    # Create clustermap
    g = sns.clustermap(
        corr,
        cmap=cmap,
        center=0,
        row_cluster=cluster,
        col_cluster=cluster,
        figsize=figsize,
        vmin=-1,
        vmax=1,
    )

    g.fig.suptitle(title, y=1.02)

    if save:
        g.savefig(save, dpi=150, bbox_inches="tight")

    return g.fig  # type: ignore[no-any-return]
