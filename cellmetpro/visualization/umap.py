"""Dimensionality reduction plots for metabolic profiles.

This module provides UMAP, t-SNE, and PCA visualizations of cells
based on their metabolic activity rather than gene expression.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_embedding(
    embedding: np.ndarray,
    color: np.ndarray | pd.Series | None = None,
    labels: list[str] | None = None,
    title: str = "Embedding",
    xlabel: str = "Dimension 1",
    ylabel: str = "Dimension 2",
    cmap: str = "viridis",
    categorical_palette: str = "husl",
    point_size: float = 20,
    alpha: float = 0.7,
    figsize: tuple[float, float] = (10, 7),
    ax: Axes | None = None,
    legend: bool = True,
    legend_loc: str = "best",
    colorbar: bool = True,
    colorbar_label: str = "",
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot 2D embedding colored by values.

    Parameters
    ----------
    embedding : np.ndarray
        2D coordinates (n_samples x 2).
    color : np.ndarray or pd.Series, optional
        Values to color points by (continuous or categorical).
    labels : list[str], optional
        Sample labels for hover/legend.
    title : str
        Plot title.
    xlabel : str
        X-axis label.
    ylabel : str
        Y-axis label.
    cmap : str
        Colormap for continuous values.
    categorical_palette : str
        Seaborn palette for categorical values.
    point_size : float
        Size of scatter points.
    alpha : float
        Transparency of points.
    figsize : tuple
        Figure size.
    ax : Axes, optional
        Matplotlib axes to plot on.
    legend : bool
        Whether to show legend for categorical data.
    legend_loc : str
        Legend location.
    colorbar : bool
        Whether to show colorbar for continuous data.
    colorbar_label : str
        Label for colorbar.
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

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Determine if color is categorical or continuous
    if color is not None:
        if isinstance(color, pd.Series):
            color_values = color.values
        else:
            color_values = color

        # Check if categorical
        is_categorical = (
            isinstance(color_values[0], str) or len(np.unique(color_values)) < 20
        )

        if is_categorical:
            _plot_categorical(
                ax,
                embedding,
                color_values,
                categorical_palette,
                point_size,
                alpha,
                legend,
                legend_loc,
                **kwargs,
            )
        else:
            scatter = ax.scatter(
                embedding[:, 0],
                embedding[:, 1],
                c=color_values,
                cmap=cmap,
                s=point_size,
                alpha=alpha,
                **kwargs,
            )
            if colorbar:
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label(colorbar_label)
    else:
        ax.scatter(
            embedding[:, 0],
            embedding[:, 1],
            c="steelblue",
            s=point_size,
            alpha=alpha,
            **kwargs,
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def _plot_categorical(
    ax: Axes,
    embedding: np.ndarray,
    categories: np.ndarray,
    palette: str,
    point_size: float,
    alpha: float,
    legend: bool,
    legend_loc: str,
    **kwargs,
) -> None:
    """Plot embedding with categorical coloring."""
    try:
        import seaborn as sns

        colors = sns.color_palette(palette, len(np.unique(categories)))
    except ImportError:
        import matplotlib

        n_categories = len(np.unique(categories))
        colors = matplotlib.colormaps["tab10"](np.linspace(0, 1, n_categories))

    unique_categories = np.unique(categories)
    color_map = dict(zip(unique_categories, colors))

    for cat in unique_categories:
        mask = categories == cat
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            c=[color_map[cat]],
            s=point_size,
            alpha=alpha,
            label=str(cat),
            **kwargs,
        )

    if legend:
        n_cats = len(unique_categories)
        if legend_loc == "best" and n_cats > 10:
            ax.legend(
                loc="upper right",
                bbox_to_anchor=(1.02, 1),
                borderaxespad=0,
                frameon=False,
                fontsize=8,
            )
        else:
            ax.legend(loc=legend_loc, frameon=False)  # type: ignore[call-overload]


def plot_metabolic_umap(
    embedding: np.ndarray,
    color: np.ndarray | pd.Series | None = None,
    labels: list[str] | None = None,
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
    labels : list[str], optional
        Cell labels.
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
    return plot_embedding(
        embedding=embedding,
        color=color,
        labels=labels,
        title=title,
        xlabel="UMAP 1",
        ylabel="UMAP 2",
        cmap=cmap,
        figsize=figsize,
        ax=ax,
        save=save,
        **kwargs,
    )


def plot_metabolic_tsne(
    embedding: np.ndarray,
    color: np.ndarray | pd.Series | None = None,
    labels: list[str] | None = None,
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
    labels : list[str], optional
        Cell labels.
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
    return plot_embedding(
        embedding=embedding,
        color=color,
        labels=labels,
        title=title,
        xlabel="t-SNE 1",
        ylabel="t-SNE 2",
        cmap=cmap,
        figsize=figsize,
        ax=ax,
        save=save,
        **kwargs,
    )


def plot_pca_variance(
    pca_model,
    n_components: int = 20,
    cumulative: bool = True,
    figsize: tuple[float, float] = (8, 4),
    title: str = "PCA Variance Explained",
    save: str | None = None,
) -> Axes:
    """Plot PCA variance explained.

    Parameters
    ----------
    pca_model : PCA
        Fitted sklearn PCA model.
    n_components : int
        Number of components to show.
    cumulative : bool
        Whether to show cumulative variance.
    figsize : tuple
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

    fig, ax = plt.subplots(figsize=figsize)

    n_components = min(n_components, len(pca_model.explained_variance_ratio_))
    variance = pca_model.explained_variance_ratio_[:n_components]
    components = np.arange(1, n_components + 1)

    ax.bar(components, variance, alpha=0.7, label="Individual")

    if cumulative:
        cumvar = np.cumsum(variance)
        ax.plot(components, cumvar, "r-o", markersize=4, label="Cumulative")

    ax.set_xlabel("Principal Component")
    ax.set_ylabel("Variance Explained")
    ax.set_title(title)
    ax.legend()
    ax.set_xticks(components)

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_feature_on_embedding(
    embedding: np.ndarray,
    features: pd.DataFrame,
    feature_names: list[str],
    ncols: int = 3,
    figsize_per_plot: tuple[float, float] = (4, 3),
    cmap: str = "viridis",
    point_size: float = 10,
    title: str | None = None,
    save: str | None = None,
) -> Figure:
    """Plot multiple features on embedding in a grid.

    Parameters
    ----------
    embedding : np.ndarray
        2D coordinates (n_samples x 2).
    features : pd.DataFrame
        Feature values (features x samples).
    feature_names : list[str]
        Names of features to plot.
    ncols : int
        Number of columns in grid.
    figsize_per_plot : tuple
        Size of each subplot.
    cmap : str
        Colormap.
    point_size : float
        Size of points.
    title : str, optional
        Overall title.
    save : str, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    import matplotlib.pyplot as plt

    n_features = len(feature_names)
    nrows = int(np.ceil(n_features / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows)
    )
    axes = np.atleast_2d(axes)

    for idx, feature_name in enumerate(feature_names):
        row = idx // ncols
        col = idx % ncols
        ax = axes[row, col]

        if feature_name in features.index:
            values = features.loc[feature_name].values
            scatter = ax.scatter(
                embedding[:, 0],
                embedding[:, 1],
                c=values,
                cmap=cmap,
                s=point_size,
                alpha=0.7,
            )
            plt.colorbar(scatter, ax=ax)
        else:
            ax.text(
                0.5,
                0.5,
                f"'{feature_name}' not found",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )

        ax.set_title(feature_name)
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide empty subplots
    for idx in range(n_features, nrows * ncols):
        row = idx // ncols
        col = idx % ncols
        axes[row, col].set_visible(False)

    if title:
        fig.suptitle(title, y=1.02)

    plt.tight_layout()

    if save:
        plt.savefig(save, dpi=150, bbox_inches="tight")

    return fig
