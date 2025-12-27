"""Heatmap visualizations for pathway and reaction activity.

This module provides heatmap plots comparing metabolic activity
across cell clusters, conditions, or samples.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes


def plot_pathway_heatmap(
    pathway_scores: pd.DataFrame,
    groups: pd.Series,
    pathways: list[str] | None = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (10, 8),
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot heatmap of pathway activity across groups.

    Parameters
    ----------
    pathway_scores : pd.DataFrame
        Pathway scores (pathways x cells).
    groups : pd.Series
        Group labels for cells.
    pathways : list[str], optional
        Specific pathways to show. If None, shows top variable.
    cluster_rows : bool
        Whether to cluster pathways.
    cluster_cols : bool
        Whether to cluster groups.
    cmap : str
        Colormap.
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to clustermap.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    raise NotImplementedError


def plot_reaction_heatmap(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 50,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    cmap: str = "RdBu_r",
    figsize: tuple[float, float] = (12, 10),
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot heatmap of reaction activity across groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Specific reactions to show.
    n_top : int
        Number of top variable reactions if reactions not specified.
    cluster_rows : bool
        Whether to cluster reactions.
    cluster_cols : bool
        Whether to cluster groups.
    cmap : str
        Colormap.
    figsize : tuple
        Figure size.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to clustermap.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    raise NotImplementedError
