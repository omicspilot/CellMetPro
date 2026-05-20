"""Tests for visualization module."""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")  # Use non-interactive backend for testing
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

    return pd.DataFrame(
        {
            "reaction": [f"R{i}" for i in range(n_reactions)],
            "log2fc": np.random.randn(n_reactions) * 2,
            "pvalue": np.random.uniform(0.0001, 0.5, n_reactions),
            "padj_bh": np.random.uniform(0.001, 0.6, n_reactions),
            "group1_mean": np.random.rand(n_reactions),
            "group2_mean": np.random.rand(n_reactions),
        }
    )


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
    return pd.DataFrame(
        {
            "go_term": [f"GO:{i:07d}" for i in range(10)],
            "go_name": [f"Process {i}" for i in range(10)],
            "fold_enrichment": np.random.uniform(1.5, 5, 10),
            "pvalue": np.random.uniform(0.0001, 0.1, 10),
            "padj": np.random.uniform(0.001, 0.15, 10),
        }
    )


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
        reaction_col="reaction",
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_volcano_custom_thresholds(differential_results):
    """Test volcano plot with custom thresholds."""
    ax = plot_volcano(differential_results, log2fc_threshold=0.5, pvalue_threshold=0.1)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_volcano_with_highlight(differential_results):
    """Test volcano plot with highlighted reactions."""
    ax = plot_volcano(differential_results, highlight=["R0", "R1", "R2"])

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
    fig = plot_reaction_heatmap(reaction_scores, reactions=["R0", "R1", "R2"])

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_grouped_heatmap(reaction_scores, groups):
    """Test grouped heatmap."""
    ax = plot_grouped_heatmap(reaction_scores, groups, n_top=10)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_grouped_heatmap_specific_features(reaction_scores, groups):
    """Test grouped heatmap with specific features."""
    ax = plot_grouped_heatmap(reaction_scores, groups, features=["R0", "R1", "R2"])

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
    ax = plot_reaction_dotplot(reaction_scores, groups, reactions=["R0", "R1", "R2"])

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
        fold_col="fold_enrichment",
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR EDGE CASES
# =============================================================================


def test_volcano_empty_significant():
    """Test volcano plot with no significant points."""
    df = pd.DataFrame(
        {
            "reaction": ["R1", "R2"],
            "log2fc": [0.1, 0.2],
            "padj_bh": [0.9, 0.8],  # All non-significant
        }
    )

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
    df = pd.DataFrame(
        {
            "go_term": ["GO:0000001"],
            "go_name": ["Test"],
            "fold_enrichment": [1.5],
            "pvalue": [0.5],
            "padj": [0.8],  # Not significant
        }
    )

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

    fig = plot_reaction_violin(reaction_scores, groups, reactions=["R0", "R1", "R2"])

    assert isinstance(fig, plt.Figure)
    plt.close()


def test_plot_reaction_violin_horizontal(reaction_scores, groups):
    """Test horizontal violin plot."""
    from cellmetpro.visualization import plot_reaction_violin

    fig = plot_reaction_violin(reaction_scores, groups, n_top=3, orient="h")

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

    ax = plot_single_reaction_violin(reaction_scores, groups, "R0", show_points=True)

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


# =============================================================================
# TESTS FOR ADVANCED VISUALIZATIONS
# =============================================================================


@pytest.fixture
def subsystem_mapping():
    """Create mock subsystem mapping."""
    return {
        "Glycolysis": ["R0", "R1", "R2", "R3"],
        "TCA Cycle": ["R4", "R5", "R6"],
        "Oxidative Phosphorylation": ["R7", "R8", "R9", "R10"],
        "Pentose Phosphate": ["R11", "R12", "R13"],
        "Fatty Acid Synthesis": ["R14", "R15", "R16", "R17", "R18", "R19"],
    }


def test_plot_stacked_bar_returns_axes(reaction_scores, groups, subsystem_mapping):
    """Test that plot_stacked_bar returns matplotlib axes."""
    from cellmetpro.visualization import plot_stacked_bar

    ax = plot_stacked_bar(reaction_scores, groups, subsystem_mapping)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_stacked_bar_not_normalized(reaction_scores, groups, subsystem_mapping):
    """Test stacked bar chart without normalization."""
    from cellmetpro.visualization import plot_stacked_bar

    ax = plot_stacked_bar(reaction_scores, groups, subsystem_mapping, normalize=False)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_stacked_bar_custom_n_top(reaction_scores, groups, subsystem_mapping):
    """Test stacked bar chart with custom n_top subsystems."""
    from cellmetpro.visualization import plot_stacked_bar

    ax = plot_stacked_bar(
        reaction_scores, groups, subsystem_mapping, n_top_subsystems=3
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_ridge_returns_figure(reaction_scores, groups):
    """Test that plot_ridge returns a figure."""
    from cellmetpro.visualization import plot_ridge

    fig = plot_ridge(reaction_scores, groups, n_top=5)

    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_plot_ridge_specific_reactions(reaction_scores, groups):
    """Test ridge plot with specific reactions."""
    from cellmetpro.visualization import plot_ridge

    fig = plot_ridge(reaction_scores, groups, reactions=["R0", "R1", "R2"])

    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_plot_ridge_custom_overlap(reaction_scores, groups):
    """Test ridge plot with custom overlap."""
    from cellmetpro.visualization import plot_ridge

    fig = plot_ridge(reaction_scores, groups, n_top=3, overlap=0.3)

    assert isinstance(fig, plt.Figure)
    plt.close("all")


def test_plot_radar_returns_axes(reaction_scores, groups):
    """Test that plot_radar returns polar axes."""
    from cellmetpro.visualization import plot_radar

    ax = plot_radar(reaction_scores, groups, n_top=5)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_radar_specific_reactions(reaction_scores, groups):
    """Test radar plot with specific reactions."""
    from cellmetpro.visualization import plot_radar

    ax = plot_radar(reaction_scores, groups, reactions=["R0", "R1", "R2", "R3", "R4"])

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_radar_median_aggregation(reaction_scores, groups):
    """Test radar plot with median aggregation."""
    from cellmetpro.visualization import plot_radar

    ax = plot_radar(reaction_scores, groups, n_top=5, agg_method="median")

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_radar_not_normalized(reaction_scores, groups):
    """Test radar plot without normalization."""
    from cellmetpro.visualization import plot_radar

    ax = plot_radar(reaction_scores, groups, n_top=5, normalize=False)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_radar_minimum_reactions():
    """Test radar plot requires at least 3 reactions."""
    from cellmetpro.visualization import plot_radar

    # Create minimal data
    reactions = ["R0", "R1"]  # Only 2 reactions
    cells = [f"cell_{i}" for i in range(10)]
    data = np.random.rand(2, 10)
    scores = pd.DataFrame(data, index=reactions, columns=cells)
    groups = pd.Series(["A"] * 5 + ["B"] * 5, index=cells)

    with pytest.raises(ValueError, match="at least 3"):
        plot_radar(scores, groups, reactions=reactions)
    plt.close()


def test_plot_waterfall_returns_axes(differential_results):
    """Test that plot_waterfall returns matplotlib axes."""
    from cellmetpro.visualization import plot_waterfall

    ax = plot_waterfall(differential_results)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_waterfall_custom_n_top(differential_results):
    """Test waterfall plot with custom n_top."""
    from cellmetpro.visualization import plot_waterfall

    ax = plot_waterfall(differential_results, n_top=15)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_waterfall_no_labels(differential_results):
    """Test waterfall plot without labels."""
    from cellmetpro.visualization import plot_waterfall

    ax = plot_waterfall(differential_results, show_labels=False)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_waterfall_custom_colors(differential_results):
    """Test waterfall plot with custom colors."""
    from cellmetpro.visualization import plot_waterfall

    ax = plot_waterfall(
        differential_results,
        up_color="#FF0000",
        down_color="#0000FF",
        ns_color="#CCCCCC",
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_subsystem_waterfall_returns_axes(enrichment_results):
    """Test that plot_subsystem_waterfall returns matplotlib axes."""
    from cellmetpro.visualization import plot_subsystem_waterfall

    # Add pathway column expected by the function
    enrichment_results["pathway"] = enrichment_results["go_name"]

    ax = plot_subsystem_waterfall(enrichment_results)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_subsystem_waterfall_custom_columns(enrichment_results):
    """Test subsystem waterfall with custom columns."""
    from cellmetpro.visualization import plot_subsystem_waterfall

    ax = plot_subsystem_waterfall(
        enrichment_results,
        fold_col="fold_enrichment",
        term_col="go_name",
        pvalue_col="padj",
    )

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_subsystem_waterfall_custom_n_top(enrichment_results):
    """Test subsystem waterfall with custom n_top."""
    from cellmetpro.visualization import plot_subsystem_waterfall

    enrichment_results["pathway"] = enrichment_results["go_name"]

    ax = plot_subsystem_waterfall(enrichment_results, n_top=5)

    assert isinstance(ax, plt.Axes)
    plt.close()


def test_stacked_bar_save(reaction_scores, groups, subsystem_mapping, tmp_path):
    """Test stacked bar save functionality."""
    from cellmetpro.visualization import plot_stacked_bar

    save_path = tmp_path / "stacked_bar.png"
    plot_stacked_bar(reaction_scores, groups, subsystem_mapping, save=str(save_path))

    assert save_path.exists()
    plt.close()


def test_radar_save(reaction_scores, groups, tmp_path):
    """Test radar plot save functionality."""
    from cellmetpro.visualization import plot_radar

    save_path = tmp_path / "radar.png"
    plot_radar(reaction_scores, groups, n_top=5, save=str(save_path))

    assert save_path.exists()
    plt.close()


def test_waterfall_save(differential_results, tmp_path):
    """Test waterfall plot save functionality."""
    from cellmetpro.visualization import plot_waterfall

    save_path = tmp_path / "waterfall.png"
    plot_waterfall(differential_results, save=str(save_path))

    assert save_path.exists()
    plt.close()


# =============================================================================
# TESTS FOR EDGE CASES IN ADVANCED VISUALIZATIONS
# =============================================================================


def test_plot_waterfall_empty_input():
    """Test waterfall plot with empty input raises error."""
    from cellmetpro.visualization import plot_waterfall

    empty_df = pd.DataFrame(columns=["reaction", "log2fc", "padj_bh"])

    with pytest.raises(ValueError, match="empty"):
        plot_waterfall(empty_df)
    plt.close()


def test_plot_subsystem_waterfall_empty_input():
    """Test subsystem waterfall with empty input raises error."""
    from cellmetpro.visualization import plot_subsystem_waterfall

    empty_df = pd.DataFrame(columns=["pathway", "fold_enrichment", "padj"])

    with pytest.raises(ValueError, match="empty"):
        plot_subsystem_waterfall(empty_df)
    plt.close()


def test_plot_stacked_bar_zero_values(groups, subsystem_mapping):
    """Test stacked bar handles all-zero values gracefully."""
    from cellmetpro.visualization import plot_stacked_bar

    # Create reaction scores with all zeros for some groups
    cells = [f"cell_{i}" for i in range(30)]
    reactions = [f"R{i}" for i in range(20)]
    data = np.zeros((20, 30))  # All zeros
    scores = pd.DataFrame(data, index=reactions, columns=cells)

    # Should not raise division by zero
    ax = plot_stacked_bar(scores, groups, subsystem_mapping)
    assert isinstance(ax, plt.Axes)
    plt.close()


def test_plot_radar_constant_reactions(reaction_scores, groups):
    """Test radar plot handles constant (no variance) reactions."""
    from cellmetpro.visualization import plot_radar

    # Create reactions with constant values
    cells = [f"cell_{i}" for i in range(30)]
    reactions = ["R0", "R1", "R2", "R3", "R4"]
    data = np.ones((5, 30))  # All ones - no variance
    scores = pd.DataFrame(data, index=reactions, columns=cells)

    # Should not raise division by zero
    ax = plot_radar(scores, groups, reactions=reactions)
    assert isinstance(ax, plt.Axes)
    plt.close()


# =============================================================================
# TESTS FOR INTERACTIVE PLOTS (PLOTLY)
# =============================================================================


class TestInteractiveVolcano:
    """Tests for interactive volcano plot."""

    def test_plot_volcano_interactive_returns_figure(self, differential_results):
        """Test that plot_volcano_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_volcano_interactive

        fig = plot_volcano_interactive(differential_results)
        assert isinstance(fig, Figure)

    def test_plot_volcano_interactive_custom_thresholds(self, differential_results):
        """Test interactive volcano with custom thresholds."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_volcano_interactive

        fig = plot_volcano_interactive(
            differential_results,
            log2fc_threshold=1.0,
            pvalue_threshold=0.01,
        )
        assert isinstance(fig, Figure)


class TestInteractiveEmbedding:
    """Tests for interactive embedding plots."""

    def test_plot_embedding_interactive_returns_figure(self, embedding):
        """Test that plot_embedding_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_embedding_interactive

        fig = plot_embedding_interactive(embedding)
        assert isinstance(fig, Figure)

    def test_plot_embedding_interactive_with_color(self, embedding):
        """Test interactive embedding with color values."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_embedding_interactive

        color = np.random.rand(len(embedding))
        fig = plot_embedding_interactive(embedding, color=color)
        assert isinstance(fig, Figure)

    def test_plot_embedding_interactive_with_categorical(self, embedding):
        """Test interactive embedding with categorical colors."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_embedding_interactive

        color = pd.Series(["A"] * 50 + ["B"] * 50)
        fig = plot_embedding_interactive(embedding, color=color)
        assert isinstance(fig, Figure)

    def test_plot_umap_interactive(self, embedding):
        """Test plot_umap_interactive wrapper."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_umap_interactive

        fig = plot_umap_interactive(embedding)
        assert isinstance(fig, Figure)

    def test_plot_tsne_interactive(self, embedding):
        """Test plot_tsne_interactive wrapper."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_tsne_interactive

        fig = plot_tsne_interactive(embedding)
        assert isinstance(fig, Figure)


class TestInteractiveHeatmap:
    """Tests for interactive heatmap."""

    def test_plot_heatmap_interactive_returns_figure(self, reaction_scores):
        """Test that plot_heatmap_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_heatmap_interactive

        fig = plot_heatmap_interactive(reaction_scores)
        assert isinstance(fig, Figure)

    def test_plot_heatmap_interactive_with_clustering(self, reaction_scores):
        """Test interactive heatmap with clustering."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_heatmap_interactive

        fig = plot_heatmap_interactive(
            reaction_scores, cluster_rows=True, cluster_cols=True
        )
        assert isinstance(fig, Figure)


class TestInteractiveDotplot:
    """Tests for interactive dotplot."""

    def test_plot_dotplot_interactive_returns_figure(self, reaction_scores, groups):
        """Test that plot_dotplot_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_dotplot_interactive

        fig = plot_dotplot_interactive(reaction_scores, groups)
        assert isinstance(fig, Figure)

    def test_plot_dotplot_interactive_subset_reactions(self, reaction_scores, groups):
        """Test interactive dotplot with subset of reactions."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_dotplot_interactive

        reactions = ["R0", "R1", "R2", "R3", "R4"]
        fig = plot_dotplot_interactive(reaction_scores, groups, reactions=reactions)
        assert isinstance(fig, Figure)


class TestInteractiveViolin:
    """Tests for interactive violin and box plots."""

    def test_plot_violin_interactive_returns_figure(self, reaction_scores, groups):
        """Test that plot_violin_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_violin_interactive

        fig = plot_violin_interactive(reaction_scores, groups)
        assert isinstance(fig, Figure)

    def test_plot_violin_interactive_with_points(self, reaction_scores, groups):
        """Test interactive violin with data points."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_violin_interactive

        fig = plot_violin_interactive(reaction_scores, groups, show_points=True)
        assert isinstance(fig, Figure)

    def test_plot_box_interactive_returns_figure(self, reaction_scores, groups):
        """Test that plot_box_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_box_interactive

        fig = plot_box_interactive(reaction_scores, groups)
        assert isinstance(fig, Figure)


class TestInteractiveEnrichment:
    """Tests for interactive enrichment plot."""

    def test_plot_enrichment_interactive_returns_figure(self, enrichment_results):
        """Test that plot_enrichment_interactive returns a Plotly figure."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_enrichment_interactive

        fig = plot_enrichment_interactive(
            enrichment_results,
            term_col="go_name",
            fold_col="fold_enrichment",
            pvalue_col="padj",
        )
        assert isinstance(fig, Figure)

    def test_plot_enrichment_interactive_empty_results(self):
        """Test interactive enrichment with no significant results."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_enrichment_interactive

        # Create results with no significant terms
        results = pd.DataFrame(
            {
                "pathway": ["P1", "P2"],
                "fold_enrichment": [1.5, 2.0],
                "padj": [0.5, 0.6],  # All above threshold
            }
        )
        fig = plot_enrichment_interactive(results, pvalue_threshold=0.05)
        assert isinstance(fig, Figure)


class TestInteractiveFeatureExpression:
    """Tests for interactive feature expression plots."""

    def test_plot_feature_expression_interactive(self, embedding, reaction_scores):
        """Test multi-panel feature expression plot."""
        from plotly.graph_objects import Figure

        from cellmetpro.visualization import plot_feature_expression_interactive

        # Ensure columns match
        reaction_scores.columns = [f"cell_{i}" for i in range(30)]

        fig = plot_feature_expression_interactive(
            embedding[:30],  # Match size
            reaction_scores,
            feature_names=["R0", "R1", "R2"],
        )
        assert isinstance(fig, Figure)


# =============================================================================
# TESTS FOR PCA VARIANCE AND FEATURE-ON-EMBEDDING PLOTS
# =============================================================================


class TestPcaVariancePlot:
    """Tests for plot_pca_variance."""

    @pytest.fixture
    def pca_model(self):
        """Create a fitted PCA model."""
        from sklearn.decomposition import PCA

        np.random.seed(42)
        X = np.random.rand(50, 20)
        pca = PCA(n_components=10)
        pca.fit(X)
        return pca

    def test_returns_axes(self, pca_model):
        """Test that plot_pca_variance returns Axes."""
        from cellmetpro.visualization.umap import plot_pca_variance

        ax = plot_pca_variance(pca_model)
        assert isinstance(ax, plt.Axes)
        plt.close()

    def test_cumulative_line_shown(self, pca_model):
        """Test with cumulative=True (default)."""
        from cellmetpro.visualization.umap import plot_pca_variance

        ax = plot_pca_variance(pca_model, cumulative=True)
        # Should have at least 2 artists (bars + line)
        assert len(ax.lines) > 0
        plt.close()

    def test_no_cumulative(self, pca_model):
        """Test with cumulative=False."""
        from cellmetpro.visualization.umap import plot_pca_variance

        ax = plot_pca_variance(pca_model, cumulative=False)
        assert isinstance(ax, plt.Axes)
        assert len(ax.lines) == 0
        plt.close()

    def test_n_components_capped(self, pca_model):
        """Test that n_components is capped by the model's actual components."""
        from cellmetpro.visualization.umap import plot_pca_variance

        # Request more components than the model has (10)
        ax = plot_pca_variance(pca_model, n_components=50)
        assert isinstance(ax, plt.Axes)
        plt.close()

    def test_save(self, pca_model, tmp_path):
        """Test saving the figure."""
        from cellmetpro.visualization.umap import plot_pca_variance

        save_path = tmp_path / "pca_variance.png"
        plot_pca_variance(pca_model, save=str(save_path))
        assert save_path.exists()
        plt.close()


class TestPlotFeatureOnEmbedding:
    """Tests for plot_feature_on_embedding."""

    @pytest.fixture
    def small_embedding(self):
        np.random.seed(42)
        return np.random.randn(30, 2)

    @pytest.fixture
    def small_reaction_scores(self):
        np.random.seed(42)
        reactions = [f"R{i}" for i in range(6)]
        cells = [f"cell_{i}" for i in range(30)]
        return pd.DataFrame(np.random.rand(6, 30), index=reactions, columns=cells)

    def test_returns_figure(self, small_embedding, small_reaction_scores):
        """Test that plot_feature_on_embedding returns a Figure."""
        from cellmetpro.visualization.umap import plot_feature_on_embedding

        fig = plot_feature_on_embedding(
            small_embedding,
            small_reaction_scores,
            feature_names=["R0", "R1", "R2"],
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_grid_layout(self, small_embedding, small_reaction_scores):
        """Test that the grid accommodates all requested features."""
        from cellmetpro.visualization.umap import plot_feature_on_embedding

        fig = plot_feature_on_embedding(
            small_embedding,
            small_reaction_scores,
            feature_names=["R0", "R1", "R2", "R3"],
            ncols=2,
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_missing_feature_handled(self, small_embedding, small_reaction_scores):
        """Test graceful handling of feature names not in the DataFrame."""
        from cellmetpro.visualization.umap import plot_feature_on_embedding

        fig = plot_feature_on_embedding(
            small_embedding,
            small_reaction_scores,
            feature_names=["R0", "NONEXISTENT"],
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_single_feature(self, small_embedding, small_reaction_scores):
        """Test with a single feature."""
        from cellmetpro.visualization.umap import plot_feature_on_embedding

        fig = plot_feature_on_embedding(
            small_embedding,
            small_reaction_scores,
            feature_names=["R0"],
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_save(self, small_embedding, small_reaction_scores, tmp_path):
        """Test saving the figure."""
        from cellmetpro.visualization.umap import plot_feature_on_embedding

        save_path = tmp_path / "features.png"
        plot_feature_on_embedding(
            small_embedding,
            small_reaction_scores,
            feature_names=["R0", "R1"],
            save=str(save_path),
        )
        assert save_path.exists()
        plt.close()


class TestCategoricalLegendPlacement:
    """Tests for _plot_categorical legend logic (>10 categories → outside)."""

    def test_many_categories_legend_outside(self):
        """With >10 categories the legend should be placed outside the axes."""
        from cellmetpro.visualization.umap import plot_embedding

        np.random.seed(42)
        n = 150
        embedding = np.random.randn(n, 2)
        # 15 distinct string categories
        categories = np.array([f"cluster_{i}" for i in range(15)] * 10)

        ax = plot_embedding(embedding, color=categories, legend=True)
        assert isinstance(ax, plt.Axes)
        legend = ax.get_legend()
        assert legend is not None
        plt.close()

    def test_few_categories_legend_inside(self):
        """With ≤10 categories the legend stays inside (default loc)."""
        from cellmetpro.visualization.umap import plot_embedding

        np.random.seed(42)
        embedding = np.random.randn(50, 2)
        categories = np.array(["A"] * 25 + ["B"] * 25)

        ax = plot_embedding(embedding, color=categories, legend=True)
        legend = ax.get_legend()
        assert legend is not None
        plt.close()


# =============================================================================
# TESTS FOR PATHWAY HEATMAP
# =============================================================================


class TestPlotPathwayHeatmap:
    """Tests for plot_pathway_heatmap."""

    @pytest.fixture
    def pathway_scores(self):
        np.random.seed(42)
        pathways = [f"PATH_{i}" for i in range(15)]
        cells = [f"cell_{i}" for i in range(20)]
        return pd.DataFrame(np.random.rand(15, 20), index=pathways, columns=cells)

    @pytest.fixture
    def groups_20(self):
        cells = [f"cell_{i}" for i in range(20)]
        return pd.Series(["A"] * 10 + ["B"] * 10, index=cells)

    def test_returns_figure(self, pathway_scores, groups_20):
        """Test that plot_pathway_heatmap returns a Figure."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        fig = plot_pathway_heatmap(pathway_scores, groups_20, n_top=5)
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_specific_pathways(self, pathway_scores, groups_20):
        """Test with explicit pathway list."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        fig = plot_pathway_heatmap(
            pathway_scores, groups_20, pathways=["PATH_0", "PATH_1"]
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_no_scale(self, pathway_scores, groups_20):
        """Test without z-score scaling."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        fig = plot_pathway_heatmap(pathway_scores, groups_20, scale=False, n_top=5)
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_median_aggregation(self, pathway_scores, groups_20):
        """Test median aggregation method."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        fig = plot_pathway_heatmap(
            pathway_scores, groups_20, agg_method="median", n_top=5
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_no_clustering(self, pathway_scores, groups_20):
        """Test with row/col clustering disabled."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        fig = plot_pathway_heatmap(
            pathway_scores,
            groups_20,
            cluster_rows=False,
            cluster_cols=False,
            n_top=5,
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_save(self, pathway_scores, groups_20, tmp_path):
        """Test saving the figure."""
        from cellmetpro.visualization.heatmap import plot_pathway_heatmap

        save_path = tmp_path / "pathway_heatmap.png"
        plot_pathway_heatmap(pathway_scores, groups_20, n_top=5, save=str(save_path))
        assert save_path.exists()
        plt.close()


# =============================================================================
# EDGE CASE TESTS FOR VIOLIN / BOXPLOT / STRIP PLOT
# =============================================================================


class TestViolinEdgeCases:
    """Targeted tests to cover uncovered branches in violin.py."""

    def test_violin_with_title_xlabel_ylabel(self, reaction_scores, groups):
        """title, xlabel, ylabel args exercise the label-setting paths."""
        from cellmetpro.visualization import plot_reaction_violin

        fig = plot_reaction_violin(
            reaction_scores, groups, n_top=3,
            title="My Violin", xlabel="Rxn", ylabel="Activity",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_violin_explicit_figsize(self, reaction_scores, groups):
        """Explicit figsize skips the auto-size branch."""
        from cellmetpro.visualization import plot_reaction_violin

        fig = plot_reaction_violin(
            reaction_scores, groups, n_top=3, figsize=(12, 5)
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_violin_explicit_palette(self, reaction_scores, groups):
        """Explicit palette skips the default-palette branch."""
        from cellmetpro.visualization import plot_reaction_violin

        fig = plot_reaction_violin(
            reaction_scores, groups, n_top=3, palette="Blues"
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_violin_many_reactions_rotates_labels(self, reaction_scores, groups):
        """n_top > 5 in vertical orientation triggers label rotation."""
        from cellmetpro.visualization import plot_reaction_violin

        fig = plot_reaction_violin(
            reaction_scores, groups, n_top=10, orient="v"
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_violin_no_valid_reactions_raises(self, reaction_scores, groups):
        """Empty reactions list raises ValueError."""
        from cellmetpro.visualization import plot_reaction_violin

        import pytest
        with pytest.raises(ValueError, match="No valid reactions"):
            plot_reaction_violin(
                reaction_scores, groups, reactions=["nonexistent_rxn"]
            )
        plt.close("all")

    def test_violin_no_common_cells_raises(self, reaction_scores):
        """No common cells between scores and groups raises ValueError."""
        from cellmetpro.visualization import plot_reaction_violin

        import pytest
        disjoint_groups = pd.Series(["A", "B"], index=["fake_cell_1", "fake_cell_2"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_reaction_violin(reaction_scores, disjoint_groups, n_top=3)
        plt.close("all")

    def test_boxplot_explicit_reactions(self, reaction_scores, groups):
        """Explicit reactions list covers the reactions-is-not-None branch."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, reactions=["R0", "R1", "R2"]
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_horizontal(self, reaction_scores, groups):
        """orient='h' covers the horizontal boxplot code path."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=3, orient="h"
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_with_title_labels(self, reaction_scores, groups):
        """title, xlabel, ylabel cover label-setting lines."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=3,
            title="My Box", xlabel="Rxn", ylabel="Score",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_horizontal_labels(self, reaction_scores, groups):
        """Horizontal orient exercises the horizontal-label code paths."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=3, orient="h",
            xlabel="Score", ylabel="Reaction",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_no_valid_reactions_raises(self, reaction_scores, groups):
        """No valid reactions raises ValueError."""
        from cellmetpro.visualization import plot_reaction_boxplot

        import pytest
        with pytest.raises(ValueError, match="No valid reactions"):
            plot_reaction_boxplot(
                reaction_scores, groups, reactions=["nonexistent"]
            )
        plt.close("all")

    def test_boxplot_no_common_cells_raises(self, reaction_scores):
        """No common cells raises ValueError."""
        from cellmetpro.visualization import plot_reaction_boxplot

        import pytest
        disjoint = pd.Series(["A"], index=["fake_cell"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_reaction_boxplot(reaction_scores, disjoint, n_top=3)
        plt.close("all")

    def test_boxplot_many_reactions_rotates_labels(self, reaction_scores, groups):
        """n_top=10 with vertical orient triggers x-label rotation."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=10, orient="v"
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_save(self, reaction_scores, groups, tmp_path):
        """save argument writes the figure to disk."""
        from cellmetpro.visualization import plot_reaction_boxplot

        save_path = tmp_path / "boxplot.png"
        plot_reaction_boxplot(
            reaction_scores, groups, n_top=3, save=str(save_path)
        )
        assert save_path.exists()
        plt.close()

    def test_boxplot_explicit_figsize(self, reaction_scores, groups):
        """Explicit figsize skips the auto-size branch in boxplot."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=3, figsize=(10, 8)
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_boxplot_explicit_palette(self, reaction_scores, groups):
        """Explicit palette skips the default-palette branch in boxplot."""
        from cellmetpro.visualization import plot_reaction_boxplot

        fig = plot_reaction_boxplot(
            reaction_scores, groups, n_top=3, palette="Pastel1"
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_single_violin_invalid_reaction_raises(self, reaction_scores, groups):
        """Reaction not in data raises ValueError."""
        from cellmetpro.visualization import plot_single_reaction_violin

        import pytest
        with pytest.raises(ValueError, match="not found"):
            plot_single_reaction_violin(
                reaction_scores, groups, "no_such_reaction"
            )
        plt.close("all")

    def test_single_violin_save(self, reaction_scores, groups, tmp_path):
        """save argument writes single violin to disk."""
        from cellmetpro.visualization import plot_single_reaction_violin

        save_path = tmp_path / "single_violin.png"
        plot_single_reaction_violin(
            reaction_scores, groups, "R0", save=str(save_path)
        )
        assert save_path.exists()
        plt.close()

    def test_strip_with_title(self, reaction_scores, groups):
        """title arg covers the title-setting line in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        fig = plot_multi_reaction_stripplot(
            reaction_scores, groups,
            reactions=["R0", "R1", "R2"],
            title="Strip Title",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_strip_many_reactions_rotation(self, reaction_scores, groups):
        """More than 5 reactions triggers label rotation in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        fig = plot_multi_reaction_stripplot(
            reaction_scores, groups,
            reactions=[f"R{i}" for i in range(8)],
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_strip_save(self, reaction_scores, groups, tmp_path):
        """save argument writes strip plot to disk."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        save_path = tmp_path / "strip.png"
        plot_multi_reaction_stripplot(
            reaction_scores, groups,
            reactions=["R0", "R1"],
            save=str(save_path),
        )
        assert save_path.exists()
        plt.close()

    def test_strip_no_valid_reactions_raises(self, reaction_scores, groups):
        """No valid reactions raises ValueError in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        import pytest
        with pytest.raises(ValueError, match="No valid reactions"):
            plot_multi_reaction_stripplot(
                reaction_scores, groups, reactions=["ghost"]
            )
        plt.close("all")

    def test_strip_no_common_cells_raises(self, reaction_scores):
        """No common cells raises ValueError in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        import pytest
        disjoint = pd.Series(["A"], index=["fake_cell"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_multi_reaction_stripplot(
                reaction_scores, disjoint, reactions=["R0"]
            )
        plt.close("all")

    def test_strip_explicit_figsize(self, reaction_scores, groups):
        """Explicit figsize skips auto-size in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        fig = plot_multi_reaction_stripplot(
            reaction_scores, groups,
            reactions=["R0", "R1"],
            figsize=(12, 5),
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_strip_explicit_palette(self, reaction_scores, groups):
        """Explicit palette skips default-palette in strip plot."""
        from cellmetpro.visualization import plot_multi_reaction_stripplot

        fig = plot_multi_reaction_stripplot(
            reaction_scores, groups,
            reactions=["R0", "R1"],
            palette="Set1",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()


# =============================================================================
# EDGE CASE TESTS FOR ADVANCED VISUALIZATION
# =============================================================================


class TestAdvancedVisualizationEdgeCases:
    """Edge cases for advanced.py functions."""

    @pytest.fixture
    def subsystem_mapping(self, reaction_scores):
        return {
            "Glycolysis": list(reaction_scores.index[:5]),
            "TCA": list(reaction_scores.index[5:10]),
        }

    def test_stacked_bar_legend_bottom(self, reaction_scores, groups, subsystem_mapping):
        """legend_loc='bottom' exercises the bottom-legend code path."""
        from cellmetpro.visualization import plot_stacked_bar

        ax = plot_stacked_bar(
            reaction_scores, groups, subsystem_mapping, legend_loc="bottom"
        )
        assert ax is not None
        plt.close()

    def test_stacked_bar_legend_none(self, reaction_scores, groups, subsystem_mapping):
        """legend_loc='none' removes the legend."""
        from cellmetpro.visualization import plot_stacked_bar

        ax = plot_stacked_bar(
            reaction_scores, groups, subsystem_mapping, legend_loc="none"
        )
        assert ax is not None
        plt.close()

    def test_stacked_bar_string_palette(self, reaction_scores, groups, subsystem_mapping):
        """String palette exercises the isinstance(palette, str) branch."""
        from cellmetpro.visualization import plot_stacked_bar

        ax = plot_stacked_bar(
            reaction_scores, groups, subsystem_mapping, palette="Set2"
        )
        assert ax is not None
        plt.close()

    def test_stacked_bar_list_palette(self, reaction_scores, groups, subsystem_mapping):
        """List palette exercises the else branch of palette handling."""
        from cellmetpro.visualization import plot_stacked_bar

        ax = plot_stacked_bar(
            reaction_scores, groups, subsystem_mapping,
            palette=["#ff0000", "#00ff00"],
        )
        assert ax is not None
        plt.close()

    def test_stacked_bar_ylabel(self, reaction_scores, groups, subsystem_mapping):
        """Explicit ylabel covers the ylabel-is-truthy branch."""
        from cellmetpro.visualization import plot_stacked_bar

        ax = plot_stacked_bar(
            reaction_scores, groups, subsystem_mapping, ylabel="Mean Activity"
        )
        assert ax is not None
        plt.close()

    def test_stacked_bar_no_common_cells_raises(self, reaction_scores, subsystem_mapping):
        """No common cells raises ValueError in stacked bar."""
        from cellmetpro.visualization import plot_stacked_bar

        disjoint = pd.Series(["A"], index=["fake_cell"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_stacked_bar(reaction_scores, disjoint, subsystem_mapping)
        plt.close("all")

    def test_stacked_bar_no_valid_subsystems_raises(self, reaction_scores, groups):
        """All reactions missing from subsystem_mapping raises ValueError."""
        from cellmetpro.visualization import plot_stacked_bar

        bad_mapping = {"Glycolysis": ["ghost_rxn_1", "ghost_rxn_2"]}
        with pytest.raises(ValueError, match="No subsystems"):
            plot_stacked_bar(reaction_scores, groups, bad_mapping)
        plt.close("all")

    def test_ridge_no_valid_reactions_raises(self, reaction_scores, groups):
        """No valid reactions raises ValueError in ridge plot."""
        from cellmetpro.visualization import plot_ridge

        with pytest.raises(ValueError, match="No valid reactions"):
            plot_ridge(reaction_scores, groups, reactions=["ghost"])
        plt.close("all")

    def test_ridge_no_common_cells_raises(self, reaction_scores):
        """No common cells raises ValueError in ridge plot."""
        from cellmetpro.visualization import plot_ridge

        disjoint = pd.Series(["A"], index=["fake_cell"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_ridge(reaction_scores, disjoint, n_top=3)
        plt.close("all")

    def test_radar_no_common_cells_raises(self, reaction_scores):
        """No common cells raises ValueError in radar chart."""
        from cellmetpro.visualization import plot_radar

        disjoint = pd.Series(["A", "B", "C"], index=["f1", "f2", "f3"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_radar(reaction_scores, disjoint, n_top=5)
        plt.close("all")

    def test_radar_string_palette(self, reaction_scores, groups):
        """String palette exercises the isinstance(palette, str) branch in radar."""
        from cellmetpro.visualization import plot_radar

        ax = plot_radar(reaction_scores, groups, n_top=5, palette="husl")
        assert ax is not None
        plt.close()

    def test_radar_list_palette(self, reaction_scores, groups):
        """List palette exercises the else branch of palette in radar."""
        from cellmetpro.visualization import plot_radar

        ax = plot_radar(
            reaction_scores, groups, n_top=5,
            palette=["#ff0000", "#00ff00", "#0000ff"],
        )
        assert ax is not None
        plt.close()

    def test_waterfall_no_pvalue_col(self, differential_results):
        """Waterfall without pvalue_col in data covers the no-pvalue branch."""
        from cellmetpro.visualization import plot_waterfall

        df_no_pval = differential_results.drop(columns=["pvalue"])
        ax = plot_waterfall(df_no_pval, pvalue_col="pvalue")
        assert ax is not None
        plt.close()

    def test_waterfall_nonsignificant_results(self, differential_results):
        """All p-values > threshold makes every bar ns_color (covers 623->626)."""
        from cellmetpro.visualization import plot_waterfall

        df = differential_results.copy()
        df["pvalue"] = 0.99  # all above default threshold
        ax = plot_waterfall(df, pvalue_col="pvalue", pvalue_threshold=0.05)
        assert ax is not None
        plt.close()

    def test_waterfall_explicit_figsize(self, differential_results):
        """Explicit figsize skips the auto-size branch in waterfall."""
        from cellmetpro.visualization import plot_waterfall

        ax = plot_waterfall(differential_results, figsize=(12, 8))
        assert ax is not None
        plt.close()

    def test_subsystem_waterfall_missing_fold_col_raises(self):
        """Missing fold_col raises ValueError in subsystem_waterfall."""
        from cellmetpro.visualization.advanced import plot_subsystem_waterfall

        df = pd.DataFrame({"pathway": ["P1"], "padj": [0.01]})
        with pytest.raises(ValueError, match="fold_enrichment"):
            plot_subsystem_waterfall(df)
        plt.close("all")

    def test_subsystem_waterfall_missing_term_col_raises(self):
        """Missing term_col raises ValueError in subsystem_waterfall."""
        from cellmetpro.visualization.advanced import plot_subsystem_waterfall

        df = pd.DataFrame({"fold_enrichment": [2.5], "padj": [0.01]})
        with pytest.raises(ValueError, match="pathway"):
            plot_subsystem_waterfall(df)
        plt.close("all")

    def test_subsystem_waterfall_explicit_figsize(self):
        """Explicit figsize skips auto-size in subsystem_waterfall."""
        from cellmetpro.visualization.advanced import plot_subsystem_waterfall

        df = pd.DataFrame({
            "pathway": ["P1", "P2", "P3"],
            "fold_enrichment": [2.0, 3.0, 1.5],
            "padj": [0.01, 0.5, 0.02],
        })
        ax = plot_subsystem_waterfall(df, figsize=(10, 6))
        assert ax is not None
        plt.close()

    def test_subsystem_waterfall_save(self, tmp_path):
        """save arg writes subsystem waterfall to disk."""
        from cellmetpro.visualization.advanced import plot_subsystem_waterfall

        df = pd.DataFrame({
            "pathway": ["P1", "P2", "P3"],
            "fold_enrichment": [2.0, 3.0, 1.5],
            "padj": [0.01, 0.5, 0.02],
        })
        save_path = tmp_path / "sw.png"
        plot_subsystem_waterfall(df, save=str(save_path))
        assert save_path.exists()
        plt.close()

    def test_waterfall_missing_required_column_raises(self, differential_results):
        """Waterfall raises ValueError when log2fc_col is absent."""
        from cellmetpro.visualization import plot_waterfall

        df = differential_results.drop(columns=["log2fc"])
        with pytest.raises(ValueError, match="log2fc"):
            plot_waterfall(df, log2fc_col="log2fc")
        plt.close("all")

    def test_ridge_save(self, reaction_scores, groups, tmp_path):
        """save arg writes ridge plot to disk."""
        from cellmetpro.visualization import plot_ridge

        save_path = tmp_path / "ridge.png"
        plot_ridge(reaction_scores, groups, n_top=3, save=str(save_path))
        assert save_path.exists()
        plt.close()

    def test_ridge_explicit_figsize_and_palette(self, reaction_scores, groups):
        """Explicit figsize and palette cover the False branches in ridge."""
        from cellmetpro.visualization import plot_ridge

        fig = plot_ridge(
            reaction_scores, groups, n_top=3,
            figsize=(12, 6), palette="husl",
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_single_violin_no_common_cells_raises(self, reaction_scores):
        """No common cells raises ValueError in single reaction violin."""
        from cellmetpro.visualization import plot_single_reaction_violin

        disjoint = pd.Series(["A", "B"], index=["fake1", "fake2"])
        with pytest.raises(ValueError, match="No common cells"):
            plot_single_reaction_violin(reaction_scores, disjoint, "R0")
        plt.close("all")

    def test_single_violin_explicit_palette(self, reaction_scores, groups):
        """Explicit palette covers the palette-is-not-None branch in single violin."""
        from cellmetpro.visualization import plot_single_reaction_violin

        ax = plot_single_reaction_violin(
            reaction_scores, groups, "R0", palette="Blues"
        )
        assert ax is not None
        plt.close()


# =============================================================================
# EDGE CASE TESTS FOR DOTPLOT
# =============================================================================


class TestDotplotEdgeCases:
    """Targeted tests to cover uncovered branches in dotplot.py."""

    def test_reaction_dotplot_size_scale_std(self, reaction_scores, groups):
        """size_scale='std' covers the std-scaling code path."""
        from cellmetpro.visualization import plot_reaction_dotplot

        ax = plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, size_scale="std"
        )
        assert ax is not None
        plt.close()

    def test_reaction_dotplot_size_scale_n_cells(self, reaction_scores, groups):
        """size_scale='n_cells' covers the n_cells-scaling code path."""
        from cellmetpro.visualization import plot_reaction_dotplot

        ax = plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, size_scale="n_cells"
        )
        assert ax is not None
        plt.close()

    def test_reaction_dotplot_color_median(self, reaction_scores, groups):
        """color_scale='median_score' covers the median-color code path."""
        from cellmetpro.visualization import plot_reaction_dotplot

        ax = plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, color_scale="median_score"
        )
        assert ax is not None
        plt.close()

    def test_reaction_dotplot_no_legend(self, reaction_scores, groups):
        """show_legend=False skips the legend-building block."""
        from cellmetpro.visualization import plot_reaction_dotplot

        ax = plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, show_legend=False
        )
        assert ax is not None
        plt.close()

    def test_reaction_dotplot_save(self, reaction_scores, groups, tmp_path):
        """save writes reaction dotplot to disk."""
        from cellmetpro.visualization import plot_reaction_dotplot

        save_path = tmp_path / "dotplot.png"
        plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, save=str(save_path)
        )
        assert save_path.exists()
        plt.close()

    def test_enrichment_dotplot_no_name_col(self, enrichment_results):
        """Enrichment dotplot without name_col falls back to term_col labels."""
        from cellmetpro.visualization import plot_enrichment_dotplot

        df = enrichment_results.rename(columns={"go_name": "ignore_col"})
        ax = plot_enrichment_dotplot(df, name_col="go_name")
        assert ax is not None
        plt.close()

    def test_enrichment_dotplot_explicit_figsize(self, enrichment_results):
        """Explicit figsize in enrichment dotplot skips auto-size branch."""
        from cellmetpro.visualization import plot_enrichment_dotplot

        ax = plot_enrichment_dotplot(enrichment_results, figsize=(10, 6))
        assert ax is not None
        plt.close()

    def test_enrichment_dotplot_save(self, enrichment_results, tmp_path):
        """save writes enrichment dotplot to disk."""
        from cellmetpro.visualization import plot_enrichment_dotplot

        save_path = tmp_path / "enrich_dotplot.png"
        plot_enrichment_dotplot(enrichment_results, save=str(save_path))
        assert save_path.exists()
        plt.close()


# =============================================================================
# EDGE CASE TESTS FOR HEATMAP
# =============================================================================


class TestHeatmapEdgeCases:
    """Targeted tests to cover uncovered branches in heatmap.py."""

    def test_reaction_heatmap_scale_column(self, reaction_scores, groups):
        """scale='column' covers the column-scaling code path."""
        from cellmetpro.visualization import plot_reaction_heatmap

        fig = plot_reaction_heatmap(reaction_scores, groups, n_top=5, scale="column")
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_reaction_heatmap_explicit_figsize(self, reaction_scores, groups):
        """Explicit figsize skips the auto-size branch in reaction heatmap."""
        from cellmetpro.visualization import plot_reaction_heatmap

        fig = plot_reaction_heatmap(
            reaction_scores, groups, n_top=5, figsize=(10, 8)
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_reaction_heatmap_save(self, reaction_scores, groups, tmp_path):
        """save writes reaction heatmap to disk."""
        from cellmetpro.visualization import plot_reaction_heatmap

        save_path = tmp_path / "heatmap.png"
        plot_reaction_heatmap(
            reaction_scores, groups, n_top=5, save=str(save_path)
        )
        assert save_path.exists()
        plt.close()

    def test_grouped_heatmap_median_agg(self, reaction_scores, groups):
        """agg_method='median' covers the median-aggregation path."""
        from cellmetpro.visualization import plot_grouped_heatmap

        ax = plot_grouped_heatmap(
            reaction_scores, groups, n_top=5, agg_method="median"
        )
        assert ax is not None
        plt.close()

    def test_grouped_heatmap_no_scale(self, reaction_scores, groups):
        """scale=False skips the scaling block in grouped heatmap."""
        from cellmetpro.visualization import plot_grouped_heatmap

        ax = plot_grouped_heatmap(
            reaction_scores, groups, n_top=5, scale=False
        )
        assert ax is not None
        plt.close()

    def test_grouped_heatmap_save(self, reaction_scores, groups, tmp_path):
        """save writes grouped heatmap to disk."""
        from cellmetpro.visualization import plot_grouped_heatmap

        save_path = tmp_path / "grouped_heatmap.png"
        plot_grouped_heatmap(
            reaction_scores, groups, n_top=5, save=str(save_path)
        )
        assert save_path.exists()
        plt.close()

    def test_reaction_heatmap_no_scale(self, reaction_scores, groups):
        """scale=None covers the branch when neither row nor column scaling."""
        from cellmetpro.visualization import plot_reaction_heatmap

        fig = plot_reaction_heatmap(
            reaction_scores, groups, n_top=5, scale=None
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_reaction_heatmap_hide_row_labels(self, reaction_scores, groups):
        """show_row_labels=False covers the False branch of label-tick path."""
        from cellmetpro.visualization import plot_reaction_heatmap

        fig = plot_reaction_heatmap(
            reaction_scores, groups, n_top=5, show_row_labels=False
        )
        assert isinstance(fig, plt.Figure)
        plt.close()

    def test_reaction_heatmap_disjoint_groups(self, reaction_scores):
        """Groups with no common cells leaves col_colors=None."""
        from cellmetpro.visualization import plot_reaction_heatmap

        disjoint = pd.Series(["A", "B"], index=["fake1", "fake2"])
        fig = plot_reaction_heatmap(reaction_scores, disjoint, n_top=5)
        assert isinstance(fig, plt.Figure)
        plt.close()


# =============================================================================
# EDGE CASE TESTS FOR UMAP AND INTERACTIVE VISUALIZATIONS
# =============================================================================


class TestUmapEdgeCases:
    """Targeted tests to cover umap.py and interactive.py uncovered branches."""

    def test_plot_embedding_explicit_ax(self, embedding):
        """Providing an explicit ax skips the fig,ax creation branch."""
        from cellmetpro.visualization import plot_embedding

        fig, ax = plt.subplots()
        result_ax = plot_embedding(embedding, ax=ax)
        assert result_ax is ax
        plt.close()

    def test_plot_embedding_series_color(self, embedding):
        """pd.Series color covers the isinstance(color, pd.Series) branch."""
        from cellmetpro.visualization import plot_embedding

        color = pd.Series(np.linspace(0, 1, len(embedding)))
        ax = plot_embedding(embedding, color=color)
        assert ax is not None
        plt.close()

    def test_plot_embedding_no_colorbar(self, embedding):
        """colorbar=False skips colorbar creation in continuous-color path."""
        from cellmetpro.visualization import plot_embedding

        color = np.random.rand(len(embedding))
        ax = plot_embedding(embedding, color=color, colorbar=False)
        assert ax is not None
        plt.close()

    def test_plot_categorical_no_legend(self, embedding):
        """legend=False skips legend in categorical embed (186->exit branch)."""
        from cellmetpro.visualization import plot_embedding

        categories = pd.Series(["A"] * 50 + ["B"] * 50)
        ax = plot_embedding(embedding, color=categories, legend=False)
        assert ax is not None
        plt.close()

    def test_reaction_dotplot_explicit_figsize(self, reaction_scores, groups):
        """Explicit figsize in reaction dotplot skips auto-size branch (153->159)."""
        from cellmetpro.visualization import plot_reaction_dotplot

        ax = plot_reaction_dotplot(
            reaction_scores, groups, n_top=5, figsize=(10, 8)
        )
        assert ax is not None
        plt.close()

    def test_interactive_volcano_save(self, differential_results, tmp_path):
        """save path in plot_volcano_interactive covers line 142."""
        try:
            from cellmetpro.visualization import plot_volcano_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        save_path = str(tmp_path / "volcano.html")
        fig = plot_volcano_interactive(differential_results, save=save_path)
        assert fig is not None
        import os
        assert os.path.exists(save_path)

    def test_interactive_embedding_save(self, embedding, tmp_path):
        """save path in plot_embedding_interactive covers line 192."""
        try:
            from cellmetpro.visualization import plot_embedding_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        save_path = str(tmp_path / "embedding.html")
        fig = plot_embedding_interactive(embedding, save=save_path)
        assert fig is not None
        import os
        assert os.path.exists(save_path)

    def test_interactive_heatmap_no_clustering(self, reaction_scores):
        """Disabling clustering in interactive heatmap covers the False branches."""
        try:
            from cellmetpro.visualization.interactive import plot_heatmap_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        fig = plot_heatmap_interactive(
            reaction_scores.iloc[:5, :10],
            cluster_rows=False,
            cluster_cols=False,
        )
        assert fig is not None

    def test_interactive_heatmap_save(self, reaction_scores, tmp_path):
        """save path in plot_heatmap_interactive covers its save line."""
        try:
            from cellmetpro.visualization.interactive import plot_heatmap_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        save_path = str(tmp_path / "heatmap.html")
        plot_heatmap_interactive(
            reaction_scores.iloc[:5, :10], save=save_path
        )
        import os
        assert os.path.exists(save_path)

    def test_interactive_embedding_with_labels(self, embedding):
        """Providing labels exercises the 'labels is not None' branch (line 192)."""
        try:
            from cellmetpro.visualization import plot_embedding_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        labels = [f"pt_{i}" for i in range(len(embedding))]
        fig = plot_embedding_interactive(embedding, labels=labels)
        assert fig is not None

    def test_interactive_heatmap_col_cluster_only(self, reaction_scores):
        """cluster_rows=False, cluster_cols=True covers the row-skip branch."""
        try:
            from cellmetpro.visualization.interactive import plot_heatmap_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        fig = plot_heatmap_interactive(
            reaction_scores.iloc[:5, :10],
            cluster_rows=False,
            cluster_cols=True,
        )
        assert fig is not None

    def test_interactive_heatmap_row_cluster_only(self, reaction_scores):
        """cluster_rows=True, cluster_cols=False covers the col-skip branch."""
        try:
            from cellmetpro.visualization.interactive import plot_heatmap_interactive
        except ImportError:
            pytest.skip("plotly not installed")

        fig = plot_heatmap_interactive(
            reaction_scores.iloc[:5, :10],
            cluster_rows=True,
            cluster_cols=False,
        )
        assert fig is not None
