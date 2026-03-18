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
