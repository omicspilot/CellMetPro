"""Visualization module for metabolic profiles.

This module provides plotting functions for:
- Dimensionality reduction plots (UMAP, t-SNE, PCA)
- Heatmaps for pathway/reaction activity
- Dot plots with statistical annotations
- Volcano plots for differential analysis
- Interactive dashboards
"""

from .advanced import (
    plot_radar,
    plot_ridge,
    plot_stacked_bar,
    plot_subsystem_waterfall,
    plot_waterfall,
)
from .dotplot import (
    plot_enrichment_dotplot,
    plot_reaction_dotplot,
)
from .heatmap import (
    plot_correlation_heatmap,
    plot_grouped_heatmap,
    plot_pathway_heatmap,
    plot_reaction_heatmap,
)

# Interactive Plotly visualizations
from .interactive import (
    plot_box_interactive,
    plot_dotplot_interactive,
    plot_embedding_interactive,
    plot_enrichment_interactive,
    plot_feature_expression_interactive,
    plot_heatmap_interactive,
    plot_tsne_interactive,
    plot_umap_interactive,
    plot_violin_interactive,
    plot_volcano_interactive,
)
from .umap import (
    plot_embedding,
    plot_feature_on_embedding,
    plot_metabolic_tsne,
    plot_metabolic_umap,
    plot_pca_variance,
)
from .violin import (
    plot_multi_reaction_stripplot,
    plot_reaction_boxplot,
    plot_reaction_violin,
    plot_single_reaction_violin,
)
from .volcano import (
    plot_ma,
    plot_volcano,
)

__all__ = [
    # Embedding plots
    "plot_embedding",
    "plot_metabolic_umap",
    "plot_metabolic_tsne",
    "plot_pca_variance",
    "plot_feature_on_embedding",
    # Heatmaps
    "plot_pathway_heatmap",
    "plot_reaction_heatmap",
    "plot_grouped_heatmap",
    "plot_correlation_heatmap",
    # Dot plots
    "plot_reaction_dotplot",
    "plot_enrichment_dotplot",
    # Violin/Box plots
    "plot_reaction_violin",
    "plot_reaction_boxplot",
    "plot_single_reaction_violin",
    "plot_multi_reaction_stripplot",
    # Volcano plots
    "plot_volcano",
    "plot_ma",
    # Advanced plots
    "plot_stacked_bar",
    "plot_ridge",
    "plot_radar",
    "plot_waterfall",
    "plot_subsystem_waterfall",
    # Interactive Plotly plots
    "plot_volcano_interactive",
    "plot_embedding_interactive",
    "plot_umap_interactive",
    "plot_tsne_interactive",
    "plot_heatmap_interactive",
    "plot_dotplot_interactive",
    "plot_violin_interactive",
    "plot_box_interactive",
    "plot_enrichment_interactive",
    "plot_feature_expression_interactive",
]
