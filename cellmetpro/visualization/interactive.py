"""Interactive visualizations using Plotly.

This module provides interactive versions of CellMetPro visualizations
with hover tooltips, zoom/pan, and HTML export capabilities.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

if TYPE_CHECKING:
    from plotly.graph_objects import Figure


def plot_volcano_interactive(
    results: pd.DataFrame,
    log2fc_col: str = "log2fc",
    pvalue_col: str = "pvalue",
    padj_col: str = "padj_bh",
    reaction_col: str = "reaction",
    log2fc_threshold: float = 0.5,
    pvalue_threshold: float = 0.05,
    use_adjusted: bool = True,
    title: str = "Volcano Plot",
    height: int = 600,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create an interactive volcano plot with Plotly.

    Parameters
    ----------
    results : pd.DataFrame
        Differential analysis results with log2fc and p-value columns.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for raw p-values.
    padj_col : str
        Column name for adjusted p-values.
    reaction_col : str
        Column name for reaction identifiers.
    log2fc_threshold : float
        Threshold for log2 fold change significance.
    pvalue_threshold : float
        Threshold for p-value significance.
    use_adjusted : bool
        Whether to use adjusted p-values for significance.
    title : str
        Plot title.
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    df = results.copy()

    # Determine which p-value column to use
    pval_col = padj_col if use_adjusted and padj_col in df.columns else pvalue_col

    # Calculate -log10(pvalue)
    df["neg_log10_pval"] = -np.log10(df[pval_col].clip(lower=1e-300))

    # Classify points
    conditions = [
        (df[log2fc_col] >= log2fc_threshold) & (df[pval_col] < pvalue_threshold),
        (df[log2fc_col] <= -log2fc_threshold) & (df[pval_col] < pvalue_threshold),
    ]
    choices = ["Up-regulated", "Down-regulated"]
    df["regulation"] = np.select(conditions, choices, default="Not significant")

    # Color mapping
    color_map = {
        "Up-regulated": "#e74c3c",
        "Down-regulated": "#3498db",
        "Not significant": "#95a5a6",
    }

    # Create figure
    fig = px.scatter(
        df,
        x=log2fc_col,
        y="neg_log10_pval",
        color="regulation",
        color_discrete_map=color_map,
        hover_name=reaction_col if reaction_col in df.columns else None,
        hover_data={
            log2fc_col: ":.3f",
            pval_col: ":.2e",
            "neg_log10_pval": False,
            "regulation": False,
        },
        title=title,
        labels={
            log2fc_col: "Log2 Fold Change",
            "neg_log10_pval": f"-log10({'Adjusted ' if use_adjusted else ''}P-value)",
        },
    )

    # Add threshold lines
    fig.add_hline(
        y=-np.log10(pvalue_threshold),
        line_dash="dash",
        line_color="gray",
        annotation_text=f"p = {pvalue_threshold}",
    )
    fig.add_vline(
        x=log2fc_threshold,
        line_dash="dash",
        line_color="gray",
    )
    fig.add_vline(
        x=-log2fc_threshold,
        line_dash="dash",
        line_color="gray",
    )

    # Update layout
    fig.update_layout(
        height=height,
        width=width,
        legend_title_text="Regulation",
        hovermode="closest",
    )

    fig.update_traces(marker={"size": 8, "opacity": 0.7})

    if save:
        fig.write_html(save)

    return fig


def plot_embedding_interactive(
    embedding: np.ndarray,
    color: pd.Series | np.ndarray | None = None,
    labels: list[str] | np.ndarray | None = None,
    title: str = "Embedding",
    xlabel: str = "Dimension 1",
    ylabel: str = "Dimension 2",
    colorbar_title: str = "Value",
    height: int = 600,
    width: int = 700,
    save: str | None = None,
) -> Figure:
    """Create an interactive 2D embedding plot with Plotly.

    Parameters
    ----------
    embedding : np.ndarray
        2D embedding coordinates (n_samples x 2).
    color : pd.Series or np.ndarray, optional
        Values for coloring points (categorical or continuous).
    labels : list or np.ndarray, optional
        Labels for hover text.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    colorbar_title : str
        Title for colorbar (continuous) or legend (categorical).
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    df = pd.DataFrame({"x": embedding[:, 0], "y": embedding[:, 1]})

    if labels is not None:
        df["label"] = labels
    else:
        df["label"] = [f"Point {i}" for i in range(len(embedding))]

    # Determine if color is categorical or continuous
    if color is not None:
        df["color"] = color
        is_categorical = (
            isinstance(color, pd.Series)
            and color.dtype == "object"
            or (hasattr(color, "dtype") and np.issubdtype(color.dtype, np.str_))
            or len(np.unique(color)) < 20
        )
    else:
        df["color"] = "Points"
        is_categorical = True

    if is_categorical:
        fig = px.scatter(
            df,
            x="x",
            y="y",
            color="color",
            hover_name="label",
            title=title,
            labels={"x": xlabel, "y": ylabel, "color": colorbar_title},
            color_discrete_sequence=px.colors.qualitative.Set2,
        )
    else:
        fig = px.scatter(
            df,
            x="x",
            y="y",
            color="color",
            hover_name="label",
            title=title,
            labels={"x": xlabel, "y": ylabel, "color": colorbar_title},
            color_continuous_scale="Viridis",
        )

    fig.update_layout(
        height=height,
        width=width,
        hovermode="closest",
    )

    fig.update_traces(marker={"size": 6, "opacity": 0.8})

    if save:
        fig.write_html(save)

    return fig


def plot_umap_interactive(
    embedding: np.ndarray,
    color: pd.Series | np.ndarray | None = None,
    labels: list[str] | np.ndarray | None = None,
    title: str = "UMAP",
    **kwargs,
) -> Figure:
    """Create an interactive UMAP plot.

    Wrapper around plot_embedding_interactive with UMAP-specific defaults.

    Parameters
    ----------
    embedding : np.ndarray
        UMAP coordinates (n_samples x 2).
    color : pd.Series or np.ndarray, optional
        Values for coloring points.
    labels : list or np.ndarray, optional
        Labels for hover text.
    title : str
        Plot title.
    **kwargs
        Additional arguments passed to plot_embedding_interactive.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    return plot_embedding_interactive(
        embedding,
        color=color,
        labels=labels,
        title=title,
        xlabel="UMAP 1",
        ylabel="UMAP 2",
        **kwargs,
    )


def plot_tsne_interactive(
    embedding: np.ndarray,
    color: pd.Series | np.ndarray | None = None,
    labels: list[str] | np.ndarray | None = None,
    title: str = "t-SNE",
    **kwargs,
) -> Figure:
    """Create an interactive t-SNE plot.

    Wrapper around plot_embedding_interactive with t-SNE-specific defaults.

    Parameters
    ----------
    embedding : np.ndarray
        t-SNE coordinates (n_samples x 2).
    color : pd.Series or np.ndarray, optional
        Values for coloring points.
    labels : list or np.ndarray, optional
        Labels for hover text.
    title : str
        Plot title.
    **kwargs
        Additional arguments passed to plot_embedding_interactive.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    return plot_embedding_interactive(
        embedding,
        color=color,
        labels=labels,
        title=title,
        xlabel="t-SNE 1",
        ylabel="t-SNE 2",
        **kwargs,
    )


def plot_heatmap_interactive(
    data: pd.DataFrame,
    groups: pd.Series | None = None,
    title: str = "Heatmap",
    xlabel: str = "Samples",
    ylabel: str = "Features",
    colorscale: str = "RdBu_r",
    cluster_rows: bool = False,
    cluster_cols: bool = False,
    height: int = 600,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create an interactive heatmap with Plotly.

    Parameters
    ----------
    data : pd.DataFrame
        Data matrix (features x samples).
    groups : pd.Series, optional
        Group labels for columns (for annotation).
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    colorscale : str
        Plotly colorscale name.
    cluster_rows : bool
        Whether to cluster rows (requires scipy).
    cluster_cols : bool
        Whether to cluster columns (requires scipy).
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    plot_data = data.copy()

    # Cluster if requested
    if cluster_rows or cluster_cols:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist

        if cluster_rows and len(plot_data) > 1:
            row_linkage = linkage(pdist(plot_data.values), method="average")
            row_order = leaves_list(row_linkage)
            plot_data = plot_data.iloc[row_order]

        if cluster_cols and len(plot_data.columns) > 1:
            col_linkage = linkage(pdist(plot_data.values.T), method="average")
            col_order = leaves_list(col_linkage)
            plot_data = plot_data.iloc[:, col_order]

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=plot_data.values,
            x=plot_data.columns.tolist(),
            y=plot_data.index.tolist(),
            colorscale=colorscale,
            hovertemplate="Sample: %{x}<br>Feature: %{y}<br>Value: %{z:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        height=height,
        width=width,
    )

    # Reverse y-axis to show first row at top
    fig.update_yaxes(autorange="reversed")

    if save:
        fig.write_html(save)

    return fig


def plot_dotplot_interactive(
    scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    title: str = "Reaction Activity",
    size_metric: str = "fraction",
    color_metric: str = "mean",
    threshold: float = 0.0,
    height: int = 600,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create an interactive dotplot with Plotly.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list, optional
        Subset of reactions to plot.
    title : str
        Plot title.
    size_metric : str
        Metric for dot size: "fraction" (cells > threshold) or "std".
    color_metric : str
        Metric for dot color: "mean" or "median".
    threshold : float
        Threshold for fraction calculation.
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    # Subset reactions if specified
    if reactions is not None:
        scores = scores.loc[scores.index.intersection(reactions)]

    # Align groups with scores
    common_cells = scores.columns.intersection(groups.index)
    scores = scores[common_cells]
    groups = groups[common_cells]

    # Compute statistics per group
    unique_groups = sorted(groups.unique())
    data_records = []

    for rxn in scores.index:
        for grp in unique_groups:
            mask = groups == grp
            values = scores.loc[rxn, mask]

            if color_metric == "mean":
                color_val = values.mean()
            else:
                color_val = values.median()

            if size_metric == "fraction":
                size_val = (values > threshold).mean()
            else:
                size_val = values.std()

            data_records.append(
                {
                    "reaction": rxn,
                    "group": grp,
                    "color": color_val,
                    "size": size_val,
                    "n_cells": mask.sum(),
                }
            )

    df = pd.DataFrame(data_records)

    # Normalize size for plotting
    df["size_scaled"] = (df["size"] - df["size"].min()) / (
        df["size"].max() - df["size"].min() + 1e-10
    )
    df["size_scaled"] = df["size_scaled"] * 30 + 5  # Scale to reasonable marker sizes

    fig = px.scatter(
        df,
        x="group",
        y="reaction",
        size="size_scaled",
        color="color",
        hover_data={
            "size": ":.3f",
            "color": ":.3f",
            "n_cells": True,
            "size_scaled": False,
        },
        title=title,
        labels={
            "group": "Group",
            "reaction": "Reaction",
            "color": f"{color_metric.capitalize()} Score",
            "size": f"{size_metric.capitalize()}",
        },
        color_continuous_scale="Viridis",
    )

    fig.update_layout(
        height=height,
        width=width,
        yaxis={"categoryorder": "category ascending"},
    )

    if save:
        fig.write_html(save)

    return fig


def plot_violin_interactive(
    scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    title: str = "Reaction Score Distribution",
    show_points: bool = False,
    height: int = 500,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create interactive violin plots with Plotly.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list, optional
        Subset of reactions to plot. If None, plots top 10 by variance.
    title : str
        Plot title.
    show_points : bool
        Whether to show individual data points.
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    # Subset reactions
    if reactions is None:
        # Select top 10 by variance
        variances = scores.var(axis=1).nlargest(10)
        reactions = variances.index.tolist()
    else:
        reactions = [r for r in reactions if r in scores.index]

    # Align data
    common_cells = scores.columns.intersection(groups.index)
    scores = scores[common_cells]
    groups = groups[common_cells]

    # Reshape for plotly
    data_records = []
    for rxn in reactions:
        for cell in common_cells:
            data_records.append(
                {"reaction": rxn, "group": groups[cell], "score": scores.loc[rxn, cell]}
            )

    df = pd.DataFrame(data_records)

    fig = px.violin(
        df,
        x="reaction",
        y="score",
        color="group",
        title=title,
        labels={"reaction": "Reaction", "score": "Score", "group": "Group"},
        points="all" if show_points else False,
        box=True,
    )

    fig.update_layout(
        height=height,
        width=width,
        xaxis_tickangle=-45,
        legend_title_text="Group",
    )

    if save:
        fig.write_html(save)

    return fig


def plot_box_interactive(
    scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    title: str = "Reaction Score Distribution",
    show_points: bool = False,
    height: int = 500,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create interactive box plots with Plotly.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list, optional
        Subset of reactions to plot. If None, plots top 10 by variance.
    title : str
        Plot title.
    show_points : bool
        Whether to show individual data points.
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    # Subset reactions
    if reactions is None:
        variances = scores.var(axis=1).nlargest(10)
        reactions = variances.index.tolist()
    else:
        reactions = [r for r in reactions if r in scores.index]

    # Align data
    common_cells = scores.columns.intersection(groups.index)
    scores = scores[common_cells]
    groups = groups[common_cells]

    # Reshape for plotly
    data_records = []
    for rxn in reactions:
        for cell in common_cells:
            data_records.append(
                {"reaction": rxn, "group": groups[cell], "score": scores.loc[rxn, cell]}
            )

    df = pd.DataFrame(data_records)

    fig = px.box(
        df,
        x="reaction",
        y="score",
        color="group",
        title=title,
        labels={"reaction": "Reaction", "score": "Score", "group": "Group"},
        points="all" if show_points else "outliers",
    )

    fig.update_layout(
        height=height,
        width=width,
        xaxis_tickangle=-45,
        legend_title_text="Group",
    )

    if save:
        fig.write_html(save)

    return fig


def plot_enrichment_interactive(
    results: pd.DataFrame,
    pvalue_col: str = "padj",
    term_col: str = "pathway",
    fold_col: str = "fold_enrichment",
    pvalue_threshold: float = 0.05,
    n_top: int = 20,
    title: str = "Pathway Enrichment",
    height: int = 600,
    width: int = 800,
    save: str | None = None,
) -> Figure:
    """Create an interactive enrichment bar plot with Plotly.

    Parameters
    ----------
    results : pd.DataFrame
        Enrichment analysis results.
    pvalue_col : str
        Column name for p-values.
    term_col : str
        Column name for pathway/term names.
    fold_col : str
        Column name for fold enrichment.
    pvalue_threshold : float
        P-value threshold for significance.
    n_top : int
        Number of top terms to display.
    title : str
        Plot title.
    height : int
        Plot height in pixels.
    width : int
        Plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    df = results.copy()

    # Filter significant and get top N
    if pvalue_col in df.columns:
        df = df[df[pvalue_col] < pvalue_threshold]
        df = df.nsmallest(n_top, pvalue_col)
        df["neg_log10_p"] = -np.log10(df[pvalue_col].clip(lower=1e-300))
    else:
        df = df.head(n_top)
        df["neg_log10_p"] = 1

    if len(df) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="No significant terms found", xref="paper", yref="paper", x=0.5, y=0.5
        )
        fig.update_layout(title=title, height=height, width=width)
        return fig

    # Sort by fold enrichment
    df = df.sort_values(fold_col, ascending=True)

    fig = px.bar(
        df,
        x=fold_col,
        y=term_col,
        color="neg_log10_p",
        orientation="h",
        title=title,
        labels={
            fold_col: "Fold Enrichment",
            term_col: "Pathway",
            "neg_log10_p": "-log10(p-value)",
        },
        color_continuous_scale="Viridis",
        hover_data={pvalue_col: ":.2e"} if pvalue_col in df.columns else None,
    )

    fig.update_layout(
        height=height,
        width=width,
        yaxis={"categoryorder": "total ascending"},
    )

    if save:
        fig.write_html(save)

    return fig


def plot_feature_expression_interactive(
    embedding: np.ndarray,
    features: pd.DataFrame,
    feature_names: list[str] | None = None,
    labels: list[str] | np.ndarray | None = None,
    ncols: int = 3,
    title: str = "Feature Expression",
    height_per_row: int = 350,
    width: int = 1000,
    save: str | None = None,
) -> Figure:
    """Create interactive multi-panel feature expression plots.

    Parameters
    ----------
    embedding : np.ndarray
        2D embedding coordinates (n_samples x 2).
    features : pd.DataFrame
        Feature values (features x samples).
    feature_names : list, optional
        Subset of features to plot. If None, plots first 9.
    labels : list or np.ndarray, optional
        Labels for hover text.
    ncols : int
        Number of columns in subplot grid.
    title : str
        Overall plot title.
    height_per_row : int
        Height per row in pixels.
    width : int
        Total plot width in pixels.
    save : str, optional
        Path to save HTML file.

    Returns
    -------
    Figure
        Plotly figure object.
    """
    if feature_names is None:
        feature_names = features.index[:9].tolist()
    else:
        feature_names = [f for f in feature_names if f in features.index]

    n_features = len(feature_names)
    nrows = (n_features + ncols - 1) // ncols

    fig = make_subplots(
        rows=nrows,
        cols=ncols,
        subplot_titles=feature_names,
        horizontal_spacing=0.05,
        vertical_spacing=0.1,
    )

    for i, feat in enumerate(feature_names):
        row = i // ncols + 1
        col = i % ncols + 1

        values = features.loc[feat].values

        hover_text = (
            labels
            if labels is not None
            else [f"Cell {j}" for j in range(len(embedding))]
        )

        fig.add_trace(
            go.Scatter(
                x=embedding[:, 0],
                y=embedding[:, 1],
                mode="markers",
                marker={
                    "color": values,
                    "colorscale": "Viridis",
                    "size": 5,
                    "colorbar": {"title": feat} if i == 0 else None,
                },
                text=hover_text,
                hovertemplate="%{text}<br>Value: %{marker.color:.3f}<extra></extra>",
                showlegend=False,
            ),
            row=row,
            col=col,
        )

    fig.update_layout(
        title=title,
        height=height_per_row * nrows,
        width=width,
    )

    if save:
        fig.write_html(save)

    return fig
