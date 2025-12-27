"""Metabolic-based cell clustering.

This module provides clustering methods based on metabolic profiles
rather than gene expression, enabling identification of metabolically
distinct cell populations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


class MetabolicClustering:
    """Cluster cells based on metabolic reaction scores.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    n_clusters : int, optional
        Number of clusters for k-means. If None, determined automatically.

    Attributes
    ----------
    labels : np.ndarray
        Cluster assignments for each cell.
    embedding : np.ndarray
        Dimensionality-reduced representation of cells.
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        n_clusters: int | None = None,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.n_clusters = n_clusters
        self.labels: np.ndarray | None = None
        self.embedding: np.ndarray | None = None

    def compute_pca(self, n_components: int = 50) -> np.ndarray:
        """Compute PCA on reaction scores.

        Parameters
        ----------
        n_components : int
            Number of principal components.

        Returns
        -------
        np.ndarray
            PCA-transformed data.
        """
        raise NotImplementedError

    def compute_umap(
        self,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
    ) -> np.ndarray:
        """Compute UMAP embedding.

        Parameters
        ----------
        n_neighbors : int
            Number of neighbors for UMAP.
        min_dist : float
            Minimum distance parameter.

        Returns
        -------
        np.ndarray
            UMAP embedding coordinates.
        """
        raise NotImplementedError

    def cluster(
        self,
        method: Literal["leiden", "kmeans", "louvain"] = "leiden",
        resolution: float = 1.0,
    ) -> np.ndarray:
        """Perform clustering on metabolic profiles.

        Parameters
        ----------
        method : str
            Clustering method to use.
        resolution : float
            Resolution parameter for community detection methods.

        Returns
        -------
        np.ndarray
            Cluster labels.
        """
        raise NotImplementedError
