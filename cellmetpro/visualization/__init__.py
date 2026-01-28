"""Visualization module for metabolic profiles.

This module provides plotting functions for:
- Dimensionality reduction plots (UMAP, t-SNE, PCA)
- Heatmaps for pathway/reaction activity
- Dot plots with statistical annotations
- Volcano plots for differential analysis
- Interactive dashboards
"""

from .umap import (
    plot_embedding,
    plot_metabolic_umap,
    plot_metabolic_tsne,
    plot_pca_variance,
    plot_feature_on_embedding,
)
from .heatmap import (
    plot_pathway_heatmap,
    plot_reaction_heatmap,
    plot_grouped_heatmap,
    plot_correlation_heatmap,
)
from .dotplot import (
    plot_reaction_dotplot,
    plot_enrichment_dotplot,
)
from .volcano import (
    plot_volcano,
    plot_ma,
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
