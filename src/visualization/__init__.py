"""Visualization module for metabolic profiles.

This module provides plotting functions for:
- Dimensionality reduction plots (UMAP, t-SNE)
- Heatmaps for pathway/reaction activity
- Dot plots with statistical annotations
- Volcano plots for differential analysis
- Interactive dashboards
"""

from .umap import plot_metabolic_umap, plot_metabolic_tsne
from .heatmap import plot_pathway_heatmap, plot_reaction_heatmap
from .dotplot import plot_reaction_dotplot
from .volcano import plot_volcano

__all__ = [
    "plot_metabolic_umap",
    "plot_metabolic_tsne",
    "plot_pathway_heatmap",
    "plot_reaction_heatmap",
    "plot_reaction_dotplot",
    "plot_volcano",
]
