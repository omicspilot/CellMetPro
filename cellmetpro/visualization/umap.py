"""Dimensionality reduction plots for metabolic profiles.

This module provides UMAP and t-SNE visualizations of cells
based on their metabolic activity rather than gene expression.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import matplotlib.pyplot as plt
    from matplotlib.axes import Axes


def plot_metabolic_umap(
    embedding: np.ndarray,
    color: np.ndarray | pd.Series | None = None,
    labels: np.ndarray | None = None,
    title: str = "Metabolic UMAP",
    cmap: str = "viridis",
    figsize: tuple[float, float] = (8, 6),
    ax: Axes | None = None,
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot UMAP embedding colored by metabolic features.

    Parameters
    ----------
    embedding : np.ndarray
        UMAP coordinates (n_cells x 2).
    color : np.ndarray or pd.Series, optional
        Values to color points by (continuous or categorical).
    labels : np.ndarray, optional
        Cluster labels for legend.
    title : str
        Plot title.
    cmap : str
        Colormap for continuous values.
    figsize : tuple
        Figure size.
    ax : Axes, optional
        Matplotlib axes to plot on.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to scatter.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    raise NotImplementedError


def plot_metabolic_tsne(
    embedding: np.ndarray,
    color: np.ndarray | pd.Series | None = None,
    labels: np.ndarray | None = None,
    title: str = "Metabolic t-SNE",
    cmap: str = "viridis",
    figsize: tuple[float, float] = (8, 6),
    ax: Axes | None = None,
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot t-SNE embedding colored by metabolic features.

    Parameters
    ----------
    embedding : np.ndarray
        t-SNE coordinates (n_cells x 2).
    color : np.ndarray or pd.Series, optional
        Values to color points by.
    labels : np.ndarray, optional
        Cluster labels for legend.
    title : str
        Plot title.
    cmap : str
        Colormap for continuous values.
    figsize : tuple
        Figure size.
    ax : Axes, optional
        Matplotlib axes to plot on.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to scatter.

    Returns
    -------
    Axes
        Matplotlib axes with the plot.
    """
    raise NotImplementedError
