"""Metabolic-based cell clustering.

This module provides clustering methods based on metabolic profiles
rather than gene expression, enabling identification of metabolically
distinct cell populations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

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
    pca_components : np.ndarray
        PCA-transformed data.
    cell_ids : list
        Cell identifiers.

    Example
    -------
    >>> mc = MetabolicClustering(reaction_scores)
    >>> mc.compute_pca(n_components=50)
    >>> mc.compute_umap()
    >>> labels = mc.cluster(method="leiden")
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
        self.pca_components: np.ndarray | None = None
        self.cell_ids = list(reaction_scores.columns)

        # Transpose to cells x reactions for sklearn
        self._data = reaction_scores.T.values

    def compute_pca(
        self,
        n_components: int = 50,
        scale: bool = True,
    ) -> np.ndarray:
        """Compute PCA on reaction scores.

        Parameters
        ----------
        n_components : int
            Number of principal components. Will be capped at
            min(n_cells, n_reactions).
        scale : bool
            Whether to standardize features before PCA.

        Returns
        -------
        np.ndarray
            PCA-transformed data (cells x components).
        """
        data = self._data.copy()

        # Handle NaN values
        data = np.nan_to_num(data, nan=0.0)

        # Scale if requested
        if scale:
            scaler = StandardScaler()
            data = scaler.fit_transform(data)

        # Cap n_components
        max_components = min(data.shape[0], data.shape[1])
        n_components = min(n_components, max_components)

        # Compute PCA
        pca = PCA(n_components=n_components)
        self.pca_components = pca.fit_transform(data)
        self._pca_model = pca

        return self.pca_components

    def compute_umap(
        self,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        n_components: int = 2,
        metric: str = "euclidean",
        use_pca: bool = True,
        random_state: int = 42,
    ) -> np.ndarray:
        """Compute UMAP embedding.

        Parameters
        ----------
        n_neighbors : int
            Number of neighbors for UMAP.
        min_dist : float
            Minimum distance parameter.
        n_components : int
            Number of UMAP dimensions (usually 2).
        metric : str
            Distance metric.
        use_pca : bool
            Whether to use PCA components as input (recommended).
        random_state : int
            Random seed for reproducibility.

        Returns
        -------
        np.ndarray
            UMAP embedding coordinates (cells x n_components).
        """
        try:
            import umap
        except ImportError:
            raise ImportError(
                "UMAP is required for compute_umap(). "
                "Install with: pip install umap-learn"
            )

        # Use PCA components if available and requested
        if use_pca and self.pca_components is not None:
            data = self.pca_components
        else:
            data = np.nan_to_num(self._data, nan=0.0)

        # Compute UMAP
        reducer = umap.UMAP(
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            n_components=n_components,
            metric=metric,
            random_state=random_state,
        )
        self.embedding = reducer.fit_transform(data)
        self._umap_model = reducer

        return self.embedding

    def compute_tsne(
        self,
        n_components: int = 2,
        perplexity: float = 30.0,
        learning_rate: float | str = "auto",
        max_iter: int = 1000,
        use_pca: bool = True,
        random_state: int = 42,
    ) -> np.ndarray:
        """Compute t-SNE embedding.

        Parameters
        ----------
        n_components : int
            Number of t-SNE dimensions (usually 2).
        perplexity : float
            Perplexity parameter (typically 5-50).
        learning_rate : float or "auto"
            Learning rate for optimization.
        max_iter : int
            Maximum number of iterations.
        use_pca : bool
            Whether to use PCA components as input (recommended for speed).
        random_state : int
            Random seed for reproducibility.

        Returns
        -------
        np.ndarray
            t-SNE embedding coordinates (cells x n_components).
        """
        # Use PCA components if available and requested
        if use_pca and self.pca_components is not None:
            data = self.pca_components
        else:
            data = np.nan_to_num(self._data, nan=0.0)

        # Adjust perplexity if needed
        n_samples = data.shape[0]
        if perplexity >= n_samples:
            perplexity = max(5.0, n_samples / 4)

        # Compute t-SNE
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            max_iter=max_iter,
            random_state=random_state,
            init="pca" if use_pca else "random",
        )
        self.embedding = tsne.fit_transform(data)
        self._tsne_model = tsne

        return self.embedding

    def cluster(
        self,
        method: Literal["leiden", "kmeans", "louvain"] = "leiden",
        resolution: float = 1.0,
        n_neighbors: int = 15,
        use_pca: bool = True,
    ) -> np.ndarray:
        """Perform clustering on metabolic profiles.

        Parameters
        ----------
        method : str
            Clustering method: 'leiden', 'kmeans', or 'louvain'.
        resolution : float
            Resolution parameter for community detection methods.
            Higher values = more clusters.
        n_neighbors : int
            Number of neighbors for graph-based methods.
        use_pca : bool
            Whether to use PCA components as input.

        Returns
        -------
        np.ndarray
            Cluster labels.
        """
        # Prepare data
        if use_pca and self.pca_components is not None:
            data = self.pca_components
        else:
            data = np.nan_to_num(self._data, nan=0.0)

        if method == "kmeans":
            self.labels = self._cluster_kmeans(data)
        elif method in ("leiden", "louvain"):
            self.labels = self._cluster_graph(data, method, resolution, n_neighbors)
        else:
            raise ValueError(f"Unknown clustering method: {method}")

        return self.labels

    def _cluster_kmeans(self, data: np.ndarray) -> np.ndarray:
        """K-means clustering."""
        n_clusters = self.n_clusters
        if n_clusters is None:
            # Estimate number of clusters using elbow method heuristic
            n_clusters = min(10, max(2, data.shape[0] // 50))

        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        return kmeans.fit_predict(data)

    def _cluster_graph(
        self,
        data: np.ndarray,
        method: str,
        resolution: float,
        n_neighbors: int,
    ) -> np.ndarray:
        """Graph-based clustering (Leiden or Louvain)."""
        try:
            import igraph as ig
            import leidenalg
        except ImportError:
            raise ImportError(
                f"{method.capitalize()} clustering requires leidenalg and igraph. "
                "Install with: pip install leidenalg igraph"
            )

        # Build k-NN graph
        n_neighbors = min(n_neighbors, data.shape[0] - 1)
        nn = NearestNeighbors(n_neighbors=n_neighbors, metric="euclidean")
        nn.fit(data)
        distances, indices = nn.kneighbors(data)

        # Build adjacency list for igraph
        edges = []
        weights = []
        for i in range(data.shape[0]):
            for j, dist in zip(indices[i], distances[i]):
                if i != j:
                    edges.append((i, j))
                    # Convert distance to similarity
                    weights.append(1.0 / (1.0 + dist))

        # Create igraph graph
        g = ig.Graph(n=data.shape[0], edges=edges, directed=False)
        g.es["weight"] = weights

        # Simplify (remove duplicates and self-loops)
        g = g.simplify(combine_edges="max")

        # Run clustering
        if method == "leiden":
            partition = leidenalg.find_partition(
                g,
                leidenalg.RBConfigurationVertexPartition,
                weights="weight",
                resolution_parameter=resolution,
            )
        else:  # louvain
            partition = leidenalg.find_partition(
                g,
                leidenalg.ModularityVertexPartition,
                weights="weight",
            )

        return np.array(partition.membership)

    def get_cluster_markers(
        self,
        n_top: int = 10,
    ) -> pd.DataFrame:
        """Identify top marker reactions for each cluster.

        Parameters
        ----------
        n_top : int
            Number of top markers per cluster.

        Returns
        -------
        pd.DataFrame
            Marker reactions with cluster, reaction, and score columns.
        """
        if self.labels is None:
            raise ValueError("Must run cluster() first")

        results = []
        reaction_ids = list(self.reaction_scores.index)

        for cluster_id in np.unique(self.labels):
            # Get cells in this cluster
            cluster_mask = self.labels == cluster_id
            other_mask = ~cluster_mask

            # Compare mean scores
            cluster_mean = self._data[cluster_mask].mean(axis=0)
            other_mean = self._data[other_mask].mean(axis=0)

            # Score = difference in means (lower penalty = more active)
            # Negative score means more active in cluster
            scores = cluster_mean - other_mean

            # Get top markers (most negative = most active in cluster)
            top_indices = np.argsort(scores)[:n_top]

            for idx in top_indices:
                results.append({
                    "cluster": cluster_id,
                    "reaction": reaction_ids[idx],
                    "cluster_mean": cluster_mean[idx],
                    "other_mean": other_mean[idx],
                    "score_diff": scores[idx],
                })

        return pd.DataFrame(results)

    def to_dataframe(self) -> pd.DataFrame:
        """Export clustering results to DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with cell_id, cluster, and embedding columns.
        """
        data = {"cell_id": self.cell_ids}

        if self.labels is not None:
            data["cluster"] = self.labels

        if self.embedding is not None:
            for i in range(self.embedding.shape[1]):
                data[f"dim_{i+1}"] = self.embedding[:, i]

        if self.pca_components is not None:
            for i in range(min(3, self.pca_components.shape[1])):
                data[f"PC{i+1}"] = self.pca_components[:, i]

        return pd.DataFrame(data)

    def to_anndata(self, adata: ad.AnnData) -> ad.AnnData:
        """Add clustering results to AnnData object.

        Parameters
        ----------
        adata : ad.AnnData
            AnnData object to add results to.

        Returns
        -------
        ad.AnnData
            Updated AnnData object.
        """
        # Add cluster labels
        if self.labels is not None:
            # Map cell IDs to adata
            label_series = pd.Series(self.labels, index=self.cell_ids)
            common_cells = adata.obs_names.intersection(label_series.index)
            adata.obs["metabolic_cluster"] = label_series[common_cells].astype(str)

        # Add embeddings
        if self.embedding is not None:
            embedding_df = pd.DataFrame(
                self.embedding,
                index=self.cell_ids,
                columns=[f"metabolic_{i+1}" for i in range(self.embedding.shape[1])]
            )
            common_cells = adata.obs_names.intersection(embedding_df.index)
            adata.obsm["X_metabolic"] = embedding_df.loc[common_cells].values

        # Add PCA
        if self.pca_components is not None:
            pca_df = pd.DataFrame(
                self.pca_components,
                index=self.cell_ids,
            )
            common_cells = adata.obs_names.intersection(pca_df.index)
            adata.obsm["X_metabolic_pca"] = pca_df.loc[common_cells].values

        return adata


def find_optimal_clusters(
    data: np.ndarray,
    max_clusters: int = 15,
    method: Literal["elbow", "silhouette"] = "silhouette",
) -> int:
    """Find optimal number of clusters using elbow or silhouette method.

    Parameters
    ----------
    data : np.ndarray
        Data matrix (samples x features).
    max_clusters : int
        Maximum number of clusters to try.
    method : str
        Method to use: 'elbow' or 'silhouette'.

    Returns
    -------
    int
        Optimal number of clusters.
    """
    from sklearn.metrics import silhouette_score

    max_clusters = min(max_clusters, data.shape[0] - 1)

    if method == "silhouette":
        scores = []
        for k in range(2, max_clusters + 1):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            labels = kmeans.fit_predict(data)
            score = silhouette_score(data, labels)
            scores.append((k, score))

        # Return k with highest silhouette score
        best_k = max(scores, key=lambda x: x[1])[0]
        return int(best_k)

    else:  # elbow
        inertias = []
        for k in range(1, max_clusters + 1):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            kmeans.fit(data)
            inertias.append(kmeans.inertia_)

        # Find elbow using second derivative
        if len(inertias) < 3:
            return 2

        diffs = np.diff(inertias)
        diffs2 = np.diff(diffs)
        elbow_idx = int(np.argmax(diffs2)) + 2  # +2 because of two diffs

        return max(2, elbow_idx)
