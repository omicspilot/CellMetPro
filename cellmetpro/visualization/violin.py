"""Violin plot visualizations for metabolic data.

This module provides violin and box plots for comparing reaction
activity distributions across cell groups.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import pandas as pd

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure


def plot_reaction_violin(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 10,
    orient: Literal["v", "h"] = "v",
    inner: Literal["box", "quartile", "point", "stick", None] = "box",
    density_norm: Literal["area", "count", "width"] = "width",
    palette: str | list | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot violin plots comparing reaction scores across groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Specific reactions to plot. If None, uses top variable reactions.
    n_top : int
        Number of top variable reactions if reactions not specified.
    orient : str
        Orientation: 'v' for vertical, 'h' for horizontal.
    inner : str or None
        Inner annotation: 'box', 'quartile', 'point', 'stick', or None.
    density_norm : str
        Method to scale violin width: 'area', 'count', or 'width'.
    palette : str or list, optional
        Color palette for groups.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to seaborn.violinplot.

    Returns
    -------
    Figure
        Matplotlib figure with violin plots.

    Example
    -------
    >>> fig = plot_reaction_violin(scores, groups, reactions=["R1", "R2", "R3"])
    >>> fig = plot_reaction_violin(scores, groups, n_top=5, inner="quartile")
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

    if len(plot_reactions) == 0:
        raise ValueError("No valid reactions to plot")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]

    # Prepare data in long format for seaborn
    data_long = scores.T.melt(
        var_name="reaction",
        value_name="score",
        ignore_index=False,
    )
    data_long["group"] = group_labels[data_long.index].values
    data_long = data_long.reset_index(drop=True)

    # Auto figure size
    n_reactions = len(plot_reactions)
    n_groups = len(group_labels.unique())
    if figsize is None:
        if orient == "v":
            width: float = max(8, min(20, n_reactions * n_groups * 0.5 + 2))
            height: float = 6
        else:
            width = 10
            height = max(6, min(15, n_reactions * 0.8 + 2))
        figsize = (width, height)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Set default palette if not specified
    if palette is None:
        palette = "Set2"

    # Create violin plot
    if orient == "v":
        sns.violinplot(
            data=data_long,
            x="reaction",
            y="score",
            hue="group",
            inner=inner,
            density_norm=density_norm,
            palette=palette,
            ax=ax,
            **kwargs,
        )
    else:
        sns.violinplot(
            data=data_long,
            x="score",
            y="reaction",
            hue="group",
            inner=inner,
            density_norm=density_norm,
            palette=palette,
            ax=ax,
            **kwargs,
        )

    # Set labels
    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    elif orient == "v":
        ax.set_xlabel("Reaction")
    else:
        ax.set_xlabel("Score")

    if ylabel:
        ax.set_ylabel(ylabel)
    elif orient == "v":
        ax.set_ylabel("Score")
    else:
        ax.set_ylabel("Reaction")

    # Rotate x labels if many reactions
    if orient == "v" and n_reactions > 5:
        ax.tick_params(axis="x", rotation=45)

    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return fig


def plot_reaction_boxplot(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str] | None = None,
    n_top: int = 10,
    orient: Literal["v", "h"] = "v",
    showfliers: bool = True,
    palette: str | list | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot box plots comparing reaction scores across groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str], optional
        Specific reactions to plot. If None, uses top variable reactions.
    n_top : int
        Number of top variable reactions if reactions not specified.
    orient : str
        Orientation: 'v' for vertical, 'h' for horizontal.
    showfliers : bool
        Whether to show outliers.
    palette : str or list, optional
        Color palette for groups.
    figsize : tuple, optional
        Figure size.
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to seaborn.boxplot.

    Returns
    -------
    Figure
        Matplotlib figure with box plots.
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

    if len(plot_reactions) == 0:
        raise ValueError("No valid reactions to plot")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]

    # Prepare data in long format
    data_long = scores.T.melt(
        var_name="reaction",
        value_name="score",
        ignore_index=False,
    )
    data_long["group"] = group_labels[data_long.index].values
    data_long = data_long.reset_index(drop=True)

    # Auto figure size
    n_reactions = len(plot_reactions)
    n_groups = len(group_labels.unique())
    if figsize is None:
        if orient == "v":
            width: float = max(8, min(20, n_reactions * n_groups * 0.4 + 2))
            height: float = 6
        else:
            width = 10
            height = max(6, min(15, n_reactions * 0.6 + 2))
        figsize = (width, height)

    fig, ax = plt.subplots(figsize=figsize)

    if palette is None:
        palette = "Set2"

    if orient == "v":
        sns.boxplot(
            data=data_long,
            x="reaction",
            y="score",
            hue="group",
            palette=palette,
            showfliers=showfliers,
            ax=ax,
            **kwargs,
        )
    else:
        sns.boxplot(
            data=data_long,
            x="score",
            y="reaction",
            hue="group",
            palette=palette,
            showfliers=showfliers,
            orient="h",
            ax=ax,
            **kwargs,
        )

    if title:
        ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    elif orient == "v":
        ax.set_xlabel("Reaction")
    else:
        ax.set_xlabel("Score")

    if ylabel:
        ax.set_ylabel(ylabel)
    elif orient == "v":
        ax.set_ylabel("Score")
    else:
        ax.set_ylabel("Reaction")

    if orient == "v" and n_reactions > 5:
        ax.tick_params(axis="x", rotation=45)

    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return fig


def plot_single_reaction_violin(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reaction: str,
    split: bool = False,
    inner: Literal["box", "quartile", "point", "stick", None] = "box",
    palette: str | list | None = None,
    figsize: tuple[float, float] = (8, 6),
    title: str | None = None,
    show_points: bool = False,
    point_size: float = 3,
    save: str | None = None,
    **kwargs,
) -> Axes:
    """Plot violin plot for a single reaction across groups.

    This is useful for detailed visualization of one reaction's
    distribution across cell groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reaction : str
        Reaction ID to plot.
    split : bool
        If True and exactly 2 groups, split violins.
    inner : str or None
        Inner annotation style.
    palette : str or list, optional
        Color palette.
    figsize : tuple
        Figure size.
    title : str, optional
        Plot title. Defaults to reaction name.
    show_points : bool
        Whether to overlay individual data points.
    point_size : float
        Size of data points if shown.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to seaborn.violinplot.

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

    if reaction not in reaction_scores.index:
        raise ValueError(f"Reaction {reaction} not found in data")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[reaction, common_cells]
    group_labels = groups[common_cells]

    # Prepare data
    data = pd.DataFrame(
        {
            "score": scores.values,
            "group": group_labels.values,
        }
    )

    fig, ax = plt.subplots(figsize=figsize)

    if palette is None:
        palette = "Set2"

    # Check if split is appropriate
    n_groups = len(group_labels.unique())
    use_split = split and n_groups == 2

    sns.violinplot(
        data=data,
        x="group",
        y="score",
        hue="group",
        inner=inner,
        palette=palette,
        split=use_split,
        legend=False,
        ax=ax,
        **kwargs,
    )

    if show_points:
        sns.stripplot(
            data=data,
            x="group",
            y="score",
            color="black",
            alpha=0.3,
            size=point_size,
            ax=ax,
        )

    ax.set_xlabel("Group")
    ax.set_ylabel("Score")
    ax.set_title(title if title else reaction)

    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return ax


def plot_multi_reaction_stripplot(
    reaction_scores: pd.DataFrame,
    groups: pd.Series,
    reactions: list[str],
    dodge: bool = True,
    alpha: float = 0.6,
    size: float = 4,
    palette: str | list | None = None,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    save: str | None = None,
    **kwargs,
) -> Figure:
    """Plot strip plot (jittered points) for multiple reactions.

    Useful for visualizing individual cell values rather than
    summarized distributions.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction scores (reactions x cells).
    groups : pd.Series
        Group labels for cells.
    reactions : list[str]
        Reactions to plot.
    dodge : bool
        Whether to separate points by group.
    alpha : float
        Point transparency.
    size : float
        Point size.
    palette : str or list, optional
        Color palette.
    figsize : tuple, optional
        Figure size.
    title : str, optional
        Plot title.
    save : str, optional
        Path to save figure.
    **kwargs
        Additional arguments passed to seaborn.stripplot.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("Seaborn is required. Install with: pip install seaborn")

    # Filter valid reactions
    plot_reactions = [r for r in reactions if r in reaction_scores.index]
    if len(plot_reactions) == 0:
        raise ValueError("No valid reactions to plot")

    # Get common cells
    common_cells = reaction_scores.columns.intersection(groups.index)
    if len(common_cells) == 0:
        raise ValueError("No common cells between scores and groups")

    scores = reaction_scores.loc[plot_reactions, common_cells]
    group_labels = groups[common_cells]

    # Prepare long format
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
        width = max(8, min(20, n_reactions * 2 + 2))
        figsize = (width, 6)

    fig, ax = plt.subplots(figsize=figsize)

    if palette is None:
        palette = "Set2"

    sns.stripplot(
        data=data_long,
        x="reaction",
        y="score",
        hue="group",
        dodge=dodge,
        alpha=alpha,
        size=size,
        palette=palette,
        ax=ax,
        **kwargs,
    )

    if title:
        ax.set_title(title)
    ax.set_xlabel("Reaction")
    ax.set_ylabel("Score")

    if n_reactions > 5:
        ax.tick_params(axis="x", rotation=45)

    ax.legend(title="Group", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()

    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")

    return fig
