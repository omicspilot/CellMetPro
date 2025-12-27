"""Dot plot visualizations with statistical annotations.

This module provides dot plots for visualizing reaction scores
across groups with size indicating expression fraction and
color indicating mean activity.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes


def plot_reaction_dotplot(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str],
    size_scale: str = "fraction_active",
    color_scale: str = "mean_score",
    threshold: float = 0.0,
    cmap: str = "Reds",
    figsize: tuple[float, float] | None = None,
    ax: Axes | None = None,
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Create dot plot of reaction activity.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str]
        Reactions to display.
    size_scale : str
        What dot size represents: 'fraction_active' or 'std'.
    color_scale : str
        What color represents: 'mean_score' or 'median_score'.
    threshold : float
        Threshold for considering a reaction "active".
    cmap : str
        Colormap for color scale.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    ax : Axes, optional
        Matplotlib axes to plot on.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    raise NotImplementedError
