"""Volcano plot for differential metabolic analysis.

This module provides volcano plots for visualizing differential
reaction activity between conditions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes


def plot_volcano(
    results: pd.DataFrame,
    log2fc_col: str = "log2fc",
    pvalue_col: str = "padj_bh",
    reaction_col: str = "reaction",
    log2fc_threshold: float = 1.0,
    pvalue_threshold: float = 0.05,
    highlight: list[str] | None = None,
    label_top: int = 10,
    colors: tuple[str, str, str] = ("#3182bd", "#969696", "#e6550d"),
    point_size: float = 20,
    alpha: float = 0.7,
    figsize: tuple[float, float] = (8, 6),
    ax: Axes | None = None,
    title: str = "Volcano Plot",
    xlabel: str = "log2 Fold Change",
    ylabel: str = "-log10(adjusted p-value)",
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Create volcano plot of differential analysis results.

    Parameters
    ----------
    results : pd.DataFrame
        Differential analysis results with log2fc and p-values.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for (adjusted) p-value.
    reaction_col : str
        Column name for reaction identifiers.
    log2fc_threshold : float
        Threshold for significant fold change.
    pvalue_threshold : float
        Threshold for significant p-value.
    highlight : list[str], optional
        Specific reactions to highlight with labels.
    label_top : int
        Number of top significant reactions to label.
    colors : tuple
        Colors for (down, not significant, up) points.
    point_size : float
        Size of scatter points.
    alpha : float
        Transparency of points.
    figsize : tuple
        Figure size.
    ax : Axes, optional
        Matplotlib axes to plot on.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
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

    # Create figure if needed
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Prepare data
    df = results.copy()
    df["neg_log10_pval"] = -np.log10(df[pvalue_col].clip(lower=1e-300))

    # Classify points
    is_significant = df[pvalue_col] < pvalue_threshold
    is_down = df[log2fc_col] <= -log2fc_threshold
    is_up = df[log2fc_col] >= log2fc_threshold
    conditions = [
        is_down & is_significant,  # Down
        is_up & is_significant,    # Up
    ]
    choices = ["down", "up"]
    df["category"] = np.select(conditions, choices, default="ns")

    # Color mapping
    color_map = {"down": colors[0], "ns": colors[1], "up": colors[2]}
    df["color"] = df["category"].map(color_map)

    # Plot non-significant first (background)
    ns_mask = df["category"] == "ns"
    ax.scatter(
        df.loc[ns_mask, log2fc_col],
        df.loc[ns_mask, "neg_log10_pval"],
        c=colors[1],
        s=point_size,
        alpha=alpha * 0.5,
        label="Not significant",
        **kwargs,
    )

    # Plot significant points
    for cat, color, label in [
        ("down", colors[0], "Down-regulated"),
        ("up", colors[2], "Up-regulated"),
    ]:
        mask = df["category"] == cat
        if mask.any():
            ax.scatter(
                df.loc[mask, log2fc_col],
                df.loc[mask, "neg_log10_pval"],
                c=color,
                s=point_size,
                alpha=alpha,
                label=label,
                **kwargs,
            )

    # Add threshold lines
    ax.axhline(-np.log10(pvalue_threshold), color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(-log2fc_threshold, color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(log2fc_threshold, color="gray", linestyle="--", linewidth=0.8)

    # Label top significant reactions
    sig_df = df[df["category"] != "ns"].nsmallest(label_top, pvalue_col)

    # Add highlighted reactions
    if highlight:
        highlight_df = df[df[reaction_col].isin(highlight)]
        sig_df = pd.concat([sig_df, highlight_df]).drop_duplicates()

    # Add labels
    texts = []
    for _, row in sig_df.iterrows():
        texts.append(ax.annotate(
            row[reaction_col],
            (row[log2fc_col], row["neg_log10_pval"]),
            fontsize=8,
            alpha=0.8,
        ))

    # Try to adjust text positions to avoid overlap
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax)
    except ImportError:
        pass  # adjustText not available

    # Labels and styling
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Save if requested
    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_ma(
    results: pd.DataFrame,
    mean_col: str = "group1_mean",
    log2fc_col: str = "log2fc",
    pvalue_col: str = "padj_bh",
    reaction_col: str = "reaction",
    pvalue_threshold: float = 0.05,
    figsize: tuple[float, float] = (8, 6),
    ax: Axes | None = None,
    save: str | None = None,
) -> Axes:
    """Create MA plot (mean vs log2fc).

    Parameters
    ----------
    results : pd.DataFrame
        Differential analysis results.
    mean_col : str
        Column name for mean expression/score.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for p-value.
    reaction_col : str
        Column name for reaction identifiers.
    pvalue_threshold : float
        Threshold for significance.
    figsize : tuple
        Figure size.
    ax : Axes, optional
        Matplotlib axes.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df = results.copy()
    significant = df[pvalue_col] < pvalue_threshold

    # Plot non-significant
    ax.scatter(
        df.loc[~significant, mean_col],
        df.loc[~significant, log2fc_col],
        c="gray",
        alpha=0.5,
        s=20,
        label="Not significant",
    )

    # Plot significant
    ax.scatter(
        df.loc[significant, mean_col],
        df.loc[significant, log2fc_col],
        c="red",
        alpha=0.7,
        s=20,
        label="Significant",
    )

    ax.axhline(0, color="black", linestyle="-", linewidth=0.5)
    ax.set_xlabel("Mean Score")
    ax.set_ylabel("log2 Fold Change")
    ax.set_title("MA Plot")
    ax.legend()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax
