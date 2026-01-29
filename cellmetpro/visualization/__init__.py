"""Visualization module for metabolic profiles.

This module provides plotting functions for:
- Dimensionality reduction plots (UMAP, t-SNE, PCA)
- Heatmaps for pathway/reaction activity
- Dot plots with statistical annotations
- Volcano plots for differential analysis
- Interactive dashboards
"""

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
from .umap import (
    plot_embedding,
    plot_feature_on_embedding,
    plot_metabolic_tsne,
    plot_metabolic_umap,
    plot_pca_variance,
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
    # Volcano plots
    "plot_volcano",
    "plot_ma",
]
