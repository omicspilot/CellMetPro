"""Tests for visualization module."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt

from cellmetpro.visualization import (
    plot_embedding,
    plot_enrichment_dotplot,
    plot_grouped_heatmap,
    plot_ma,
    plot_metabolic_tsne,
    plot_metabolic_umap,
    plot_reaction_dotplot,
    plot_reaction_heatmap,
    plot_volcano,
)


@pytest.fixture
def differential_results():
    """Create mock differential analysis results."""
    np.random.seed(42)
    n_reactions = 50

    return pd.DataFrame({
        "reaction": [f"R{i}" for i in range(n_reactions)],
        "log2fc": np.random.randn(n_reactions) * 2,
        "pvalue": np.random.uniform(0.0001, 0.5, n_reactions),
        "padj_bh": np.random.uniform(0.001, 0.6, n_reactions),
        "group1_mean": np.random.rand(n_reactions),
        "group2_mean": np.random.rand(n_reactions),
    })


@pytest.fixture
def embedding():
    """Create mock 2D embedding."""
    np.random.seed(42)
    return np.random.randn(100, 2)


@pytest.fixture
def reaction_scores():
    """Create mock reaction scores."""
    np.random.seed(42)
    reactions = [f"R{i}" for i in range(20)]
    cells = [f"cell_{i}" for i in range(30)]
    data = np.random.rand(20, 30)
    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def groups():
    """Create mock group labels."""
    cells = [f"cell_{i}" for i in range(30)]
    labels = ["A"] * 10 + ["B"] * 10 + ["C"] * 10
    return pd.Series(labels, index=cells)


@pytest.fixture
def enrichment_results():
    """Create mock enrichment results."""
    return pd.DataFrame({
        "go_term": [f"GO:{i:07d}" for i in range(10)],
        "go_name": [f"Process {i}" for i in range(10)],
        "fold_enrichment": np.random.uniform(1.5, 5, 10),
        "pvalue": np.random.uniform(0.0001, 0.1, 10),
        "padj": np.random.uniform(0.001, 0.15, 10),
    })


# =============================================================================
# TESTS FOR VOLCANO PLOT
# =============================================================================


def test_plot_volcano_returns_axes(differential_results):
    """Test that plot_volcano returns matplotlib axes."""
    ax = plot_volcano(differential_results)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_volcano_custom_columns(differential_results):
    """Test volcano plot with custom column names."""
    ax = plot_volcano(
        differential_results,
        log2fc_col="log2fc",
        pvalue_col="pvalue",
        reaction_col="reaction"
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_volcano_custom_thresholds(differential_results):
    """Test volcano plot with custom thresholds."""
    ax = plot_volcano(
        differential_results,
        log2fc_threshold=0.5,
        pvalue_threshold=0.1
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_volcano_with_highlight(differential_results):
    """Test volcano plot with highlighted reactions."""
    ax = plot_volcano(
        differential_results,
        highlight=["R0", "R1", "R2"]
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_ma_returns_axes(differential_results):
    """Test that plot_ma returns matplotlib axes."""
    ax = plot_ma(differential_results)

    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR EMBEDDING PLOTS
# =============================================================================


def test_plot_embedding_returns_axes(embedding):
    """Test that plot_embedding returns matplotlib axes."""
    ax = plot_embedding(embedding)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_embedding_with_continuous_color(embedding):
    """Test embedding plot with continuous color."""
    colors = np.random.rand(100)
    ax = plot_embedding(embedding, color=colors)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_embedding_with_categorical_color(embedding):
    """Test embedding plot with categorical color."""
    categories = np.array(["A"] * 50 + ["B"] * 50)
    ax = plot_embedding(embedding, color=categories)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_metabolic_umap(embedding):
    """Test UMAP-specific plot function."""
    ax = plot_metabolic_umap(embedding)

    assert isinstance(ax, plt.Axes)
    # Check axis labels
    assert "UMAP" in ax.get_xlabel()
    plt.close()


def test_plot_metabolic_tsne(embedding):
    """Test t-SNE-specific plot function."""
    ax = plot_metabolic_tsne(embedding)

    assert isinstance(ax, plt.Axes)
    assert "t-SNE" in ax.get_xlabel()
    plt.close()


# =============================================================================
# TESTS FOR HEATMAPS
# =============================================================================


def test_plot_reaction_heatmap_returns_figure(reaction_scores):
    """Test that plot_reaction_heatmap returns a figure."""
    fig = plot_reaction_heatmap(reaction_scores, n_top=10)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_heatmap_with_groups(reaction_scores, groups):
    """Test heatmap with group coloring."""
    fig = plot_reaction_heatmap(reaction_scores, groups=groups, n_top=10)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_heatmap_specific_reactions(reaction_scores):
    """Test heatmap with specific reactions."""
    fig = plot_reaction_heatmap(
        reaction_scores,
        reactions=["R0", "R1", "R2"]
    )

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_grouped_heatmap(reaction_scores, groups):
    """Test grouped heatmap."""
    ax = plot_grouped_heatmap(reaction_scores, groups, n_top=10)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_grouped_heatmap_specific_features(reaction_scores, groups):
    """Test grouped heatmap with specific features."""
    ax = plot_grouped_heatmap(
        reaction_scores,
        groups,
        features=["R0", "R1", "R2"]
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR DOT PLOTS
# =============================================================================


def test_plot_reaction_dotplot_returns_axes(reaction_scores, groups):
    """Test that plot_reaction_dotplot returns matplotlib axes."""
    ax = plot_reaction_dotplot(reaction_scores, groups, n_top=10)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_reaction_dotplot_specific_reactions(reaction_scores, groups):
    """Test dotplot with specific reactions."""
    ax = plot_reaction_dotplot(
        reaction_scores,
        groups,
        reactions=["R0", "R1", "R2"]
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_enrichment_dotplot(enrichment_results):
    """Test enrichment dotplot."""
    ax = plot_enrichment_dotplot(enrichment_results)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_enrichment_dotplot_custom_columns(enrichment_results):
    """Test enrichment dotplot with custom columns."""
    ax = plot_enrichment_dotplot(
        enrichment_results,
        term_col="go_term",
        name_col="go_name",
        pvalue_col="padj",
        fold_col="fold_enrichment"
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR EDGE CASES
# =============================================================================


def test_volcano_empty_significant():
    """Test volcano plot with no significant points."""
    df = pd.DataFrame({
        "reaction": ["R1", "R2"],
        "log2fc": [0.1, 0.2],
        "padj_bh": [0.9, 0.8],  # All non-significant
    })

    ax = plot_volcano(df)
    assert isinstance(ax, plt.Axes)
    plt.close()


def test_heatmap_single_group(reaction_scores):
    """Test heatmap with single group."""
    single_group = pd.Series(["A"] * 30, index=reaction_scores.columns)
    fig = plot_reaction_heatmap(reaction_scores, groups=single_group, n_top=5)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_enrichment_dotplot_no_significant():
    """Test enrichment dotplot with no significant terms."""
    df = pd.DataFrame({
        "go_term": ["GO:0000001"],
        "go_name": ["Test"],
        "fold_enrichment": [1.5],
        "pvalue": [0.5],
        "padj": [0.8],  # Not significant
    })

    ax = plot_enrichment_dotplot(df, pvalue_threshold=0.05)
    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR SAVE FUNCTIONALITY
# =============================================================================


def test_volcano_save(differential_results, tmp_path):
    """Test volcano plot save functionality."""
    save_path = tmp_path / "volcano.png"
    plot_volcano(differential_results, save=str(save_path))

    assert save_path.exists()
    plt.close()


def test_embedding_save(embedding, tmp_path):
    """Test embedding plot save functionality."""
    save_path = tmp_path / "embedding.png"
    plot_embedding(embedding, save=str(save_path))

    assert save_path.exists()
    plt.close()


# =============================================================================
# TESTS FOR VIOLIN PLOTS
# =============================================================================


def test_plot_reaction_violin_returns_figure(reaction_scores, groups):
    """Test that plot_reaction_violin returns a figure."""
    from cellmetpro.visualization import plot_reaction_violin

    fig = plot_reaction_violin(reaction_scores, groups, n_top=5)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_violin_specific_reactions(reaction_scores, groups):
    """Test violin plot with specific reactions."""
    from cellmetpro.visualization import plot_reaction_violin

    fig = plot_reaction_violin(
        reaction_scores, groups, reactions=["R0", "R1", "R2"]
    )

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_violin_horizontal(reaction_scores, groups):
    """Test horizontal violin plot."""
    from cellmetpro.visualization import plot_reaction_violin

    fig = plot_reaction_violin(
        reaction_scores, groups, n_top=3, orient="h"
    )

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_boxplot_returns_figure(reaction_scores, groups):
    """Test that plot_reaction_boxplot returns a figure."""
    from cellmetpro.visualization import plot_reaction_boxplot

    fig = plot_reaction_boxplot(reaction_scores, groups, n_top=5)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_single_reaction_violin(reaction_scores, groups):
    """Test single reaction violin plot."""
    from cellmetpro.visualization import plot_single_reaction_violin

    ax = plot_single_reaction_violin(reaction_scores, groups, "R0")

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_single_reaction_violin_with_points(reaction_scores, groups):
    """Test single reaction violin with overlaid points."""
    from cellmetpro.visualization import plot_single_reaction_violin

    ax = plot_single_reaction_violin(
        reaction_scores, groups, "R0", show_points=True
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_multi_reaction_stripplot(reaction_scores, groups):
    """Test multi-reaction strip plot."""
    from cellmetpro.visualization import plot_multi_reaction_stripplot

    fig = plot_multi_reaction_stripplot(
        reaction_scores, groups, reactions=["R0", "R1", "R2"]
    )

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_violin_save(reaction_scores, groups, tmp_path):
    """Test violin plot save functionality."""
    from cellmetpro.visualization import plot_reaction_violin

    save_path = tmp_path / "violin.png"
    plot_reaction_violin(reaction_scores, groups, n_top=3, save=str(save_path))

    assert save_path.exists()
    plt.close()


# =============================================================================
# TESTS FOR CORRELATION HEATMAP
# =============================================================================


def test_plot_correlation_heatmap_returns_figure(reaction_scores):
    """Test that plot_correlation_heatmap returns a figure."""
    from cellmetpro.visualization import plot_correlation_heatmap

    fig = plot_correlation_heatmap(reaction_scores, n_top=10)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_correlation_heatmap_specific_features(reaction_scores):
    """Test correlation heatmap with specific features."""
    from cellmetpro.visualization import plot_correlation_heatmap

    fig = plot_correlation_heatmap(
        reaction_scores, features=["R0", "R1", "R2", "R3", "R4"]
    )

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_correlation_heatmap_spearman(reaction_scores):
    """Test correlation heatmap with Spearman method."""
    from cellmetpro.visualization import plot_correlation_heatmap

    fig = plot_correlation_heatmap(reaction_scores, method="spearman", n_top=10)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_correlation_heatmap_no_cluster(reaction_scores):
    """Test correlation heatmap without clustering."""
    from cellmetpro.visualization import plot_correlation_heatmap

    fig = plot_correlation_heatmap(reaction_scores, cluster=False, n_top=10)

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_correlation_heatmap_save(reaction_scores, tmp_path):
    """Test correlation heatmap save functionality."""
    from cellmetpro.visualization import plot_correlation_heatmap

    save_path = tmp_path / "correlation.png"
    plot_correlation_heatmap(reaction_scores, n_top=10, save=str(save_path))

    assert save_path.exists()
    plt.close()
