"""Metabolic-based cell clustering.

This module provides clustering methods based on metabolic profiles
rather than gene expression, enabling identification of metabolically
distinct cell populations.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from sklearn.cluster import (
    DBSCAN,
    AgglomerativeClustering,
    KMeans,
    MiniBatchKMeans,
    SpectralClustering,
)
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
        method: Literal[
            "leiden",
            "kmeans",
            "louvain",
            "hierarchical",
            "spectral",
            "dbscan",
            "minibatch_kmeans",
            "hdbscan",
        ] = "leiden",
        resolution: float = 1.0,
        n_neighbors: int = 15,
        use_pca: bool = True,
        linkage: Literal["ward", "complete", "average", "single"] = "ward",
        eps: float | None = None,
        min_samples: int = 5,
    ) -> np.ndarray:
        """Perform clustering on metabolic profiles.

        Parameters
        ----------
        method : str
            Clustering method:
            - 'leiden': Graph-based community detection (default)
            - 'louvain': Modularity-based community detection
            - 'kmeans': K-means clustering
            - 'minibatch_kmeans': Memory-efficient K-means for large datasets
            - 'hierarchical': Agglomerative hierarchical clustering
            - 'spectral': Spectral clustering using graph Laplacian
            - 'dbscan': Density-based clustering (auto-detects n_clusters)
            - 'hdbscan': Hierarchical DBSCAN (requires hdbscan package)
        resolution : float
            Resolution parameter for community detection methods.
            Higher values = more clusters.
        n_neighbors : int
            Number of neighbors for graph-based methods.
        use_pca : bool
            Whether to use PCA components as input.
        linkage : str
            Linkage criterion for hierarchical clustering:
            'ward', 'complete', 'average', or 'single'.
        eps : float, optional
            The maximum distance between samples for DBSCAN.
            If None, estimated from data.
        min_samples : int
            Minimum samples in a neighborhood for DBSCAN core points.

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
        elif method == "minibatch_kmeans":
            self.labels = self._cluster_minibatch_kmeans(data)
        elif method in ("leiden", "louvain"):
            self.labels = self._cluster_graph(data, method, resolution, n_neighbors)
        elif method == "hierarchical":
            self.labels = self._cluster_hierarchical(data, linkage)
        elif method == "spectral":
            self.labels = self._cluster_spectral(data, n_neighbors)
        elif method == "dbscan":
            self.labels = self._cluster_dbscan(data, eps, min_samples)
        elif method == "hdbscan":
            self.labels = self._cluster_hdbscan(data, min_samples)
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
        return np.asarray(kmeans.fit_predict(data))

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
                seed=42,
            )
        else:  # louvain
            partition = leidenalg.find_partition(
                g,
                leidenalg.ModularityVertexPartition,
                weights="weight",
                seed=42,
            )

        return np.array(partition.membership)

    def _cluster_minibatch_kmeans(self, data: np.ndarray) -> np.ndarray:
        """MiniBatch K-means clustering (memory efficient for large datasets)."""
        n_clusters = self.n_clusters
        if n_clusters is None:
            n_clusters = min(10, max(2, data.shape[0] // 50))

        kmeans = MiniBatchKMeans(
            n_clusters=n_clusters,
            random_state=42,
            batch_size=min(1024, data.shape[0]),
            n_init=10,
        )
        return np.asarray(kmeans.fit_predict(data))

    def _cluster_hierarchical(
        self,
        data: np.ndarray,
        linkage: str = "ward",
    ) -> np.ndarray:
        """Agglomerative hierarchical clustering."""
        n_clusters = self.n_clusters
        if n_clusters is None:
            n_clusters = min(10, max(2, data.shape[0] // 50))

        clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage=linkage,
        )
        return np.asarray(clustering.fit_predict(data))

    def _cluster_spectral(
        self,
        data: np.ndarray,
        n_neighbors: int = 15,
    ) -> np.ndarray:
        """Spectral clustering using graph Laplacian."""
        n_clusters = self.n_clusters
        if n_clusters is None:
            n_clusters = min(10, max(2, data.shape[0] // 50))

        # Cap n_neighbors for small datasets
        n_neighbors = min(n_neighbors, data.shape[0] - 1)

        clustering = SpectralClustering(
            n_clusters=n_clusters,
            affinity="nearest_neighbors",
            n_neighbors=n_neighbors,
            random_state=42,
        )
        return np.asarray(clustering.fit_predict(data))

    def _cluster_dbscan(
        self,
        data: np.ndarray,
        eps: float | None = None,
        min_samples: int = 5,
    ) -> np.ndarray:
        """DBSCAN density-based clustering.

        Automatically detects number of clusters based on density.
        Points labeled as -1 are noise/outliers.
        """
        if eps is None:
            # Estimate eps using k-nearest neighbors distances
            k = min(min_samples, data.shape[0] - 1)
            nn = NearestNeighbors(n_neighbors=k)
            nn.fit(data)
            distances, _ = nn.kneighbors(data)
            # Use the median of k-th nearest neighbor distances
            eps = float(np.median(distances[:, -1]) * 1.5)

        clustering = DBSCAN(
            eps=eps,
            min_samples=min_samples,
            metric="euclidean",
        )
        labels: np.ndarray = clustering.fit_predict(data)

        # Assign noise points (-1) to nearest cluster rather than grouping
        # them together, which would create a spurious metabolic state
        if -1 in labels:
            noise_mask = labels == -1
            n_noise = noise_mask.sum()
            if labels.max() >= 0:
                # Find nearest non-noise neighbor for each noise point
                non_noise_mask = ~noise_mask
                nn = NearestNeighbors(n_neighbors=1, metric="euclidean")
                nn.fit(data[non_noise_mask])
                _, indices = nn.kneighbors(data[noise_mask])
                # Map back to original labels
                non_noise_labels = labels[non_noise_mask]
                labels[noise_mask] = non_noise_labels[indices.flatten()]

                import logging

                logging.getLogger(__name__).info(
                    f"DBSCAN: assigned {n_noise} noise points to nearest clusters"
                )
            else:
                # All points are noise - assign to cluster 0
                labels[:] = 0

        return np.asarray(labels)

    def _cluster_hdbscan(
        self,
        data: np.ndarray,
        min_samples: int = 5,
    ) -> np.ndarray:
        """HDBSCAN hierarchical density-based clustering.

        Automatically detects number of clusters and handles varying densities.
        """
        try:
            import hdbscan
        except ImportError:
            raise ImportError(
                "HDBSCAN clustering requires the hdbscan package. "
                "Install with: pip install hdbscan"
            )

        min_cluster_size = max(5, data.shape[0] // 50)

        clustering = hdbscan.HDBSCAN(
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            metric="euclidean",
        )
        labels: np.ndarray = clustering.fit_predict(data)

        # Assign noise points (-1) to nearest cluster rather than grouping
        # them together, which would create a spurious metabolic state
        if -1 in labels:
            noise_mask = labels == -1
            if labels.max() >= 0:
                non_noise_mask = ~noise_mask
                nn = NearestNeighbors(n_neighbors=1, metric="euclidean")
                nn.fit(data[non_noise_mask])
                _, indices = nn.kneighbors(data[noise_mask])
                non_noise_labels = labels[non_noise_mask]
                labels[noise_mask] = non_noise_labels[indices.flatten()]
            else:
                labels[:] = 0

        return np.asarray(labels)

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
                results.append(
                    {
                        "cluster": cluster_id,
                        "reaction": reaction_ids[idx],
                        "cluster_mean": cluster_mean[idx],
                        "other_mean": other_mean[idx],
                        "score_diff": scores[idx],
                    }
                )

        return pd.DataFrame(results)

    def to_dataframe(self) -> pd.DataFrame:
        """Export clustering results to DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with cell_id, cluster, and embedding columns.
        """
        data: dict[str, list | np.ndarray] = {"cell_id": self.cell_ids}

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
                columns=[f"metabolic_{i+1}" for i in range(self.embedding.shape[1])],
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


def evaluate_clustering(
    data: np.ndarray,
    labels: np.ndarray,
) -> dict[str, float]:
    """Evaluate clustering quality using multiple metrics.

    Parameters
    ----------
    data : np.ndarray
        Data matrix (samples x features).
    labels : np.ndarray
        Cluster labels for each sample.

    Returns
    -------
    dict
        Dictionary with metric names and scores:
        - silhouette: Silhouette score (-1 to 1, higher is better)
        - calinski_harabasz: Calinski-Harabasz index (higher is better)
        - davies_bouldin: Davies-Bouldin index (lower is better)
        - n_clusters: Number of clusters found

    Example
    -------
    >>> mc = MetabolicClustering(reaction_scores)
    >>> mc.compute_pca()
    >>> labels = mc.cluster(method="leiden")
    >>> metrics = evaluate_clustering(mc.pca_components, labels)
    >>> print(f"Silhouette: {metrics['silhouette']:.3f}")
    """
    from sklearn.metrics import (
        calinski_harabasz_score,
        davies_bouldin_score,
        silhouette_score,
    )

    n_clusters = len(np.unique(labels))
    n_samples = len(labels)

    # Need at least 2 clusters and more samples than clusters
    if n_clusters < 2 or n_clusters >= n_samples:
        return {
            "silhouette": np.nan,
            "calinski_harabasz": np.nan,
            "davies_bouldin": np.nan,
            "n_clusters": n_clusters,
        }

    return {
        "silhouette": float(silhouette_score(data, labels)),
        "calinski_harabasz": float(calinski_harabasz_score(data, labels)),
        "davies_bouldin": float(davies_bouldin_score(data, labels)),
        "n_clusters": n_clusters,
    }


def compare_clusterings(
    labels1: np.ndarray,
    labels2: np.ndarray,
) -> dict[str, float]:
    """Compare two clustering solutions using multiple metrics.

    Useful for comparing results from different algorithms or parameter settings.

    Parameters
    ----------
    labels1 : np.ndarray
        First set of cluster labels.
    labels2 : np.ndarray
        Second set of cluster labels.

    Returns
    -------
    dict
        Dictionary with comparison metrics:
        - adjusted_rand: Adjusted Rand Index (-1 to 1, 1 = identical)
        - normalized_mutual_info: Normalized MI (0 to 1, 1 = identical)
        - fowlkes_mallows: Fowlkes-Mallows index (0 to 1, 1 = identical)
        - adjusted_mutual_info: Adjusted MI (0 to 1, 1 = identical)

    Example
    -------
    >>> labels_leiden = mc.cluster(method="leiden")
    >>> labels_kmeans = mc.cluster(method="kmeans")
    >>> comparison = compare_clusterings(labels_leiden, labels_kmeans)
    >>> print(f"ARI: {comparison['adjusted_rand']:.3f}")
    """
    from sklearn.metrics import (
        adjusted_mutual_info_score,
        adjusted_rand_score,
        fowlkes_mallows_score,
        normalized_mutual_info_score,
    )

    if len(labels1) != len(labels2):
        raise ValueError(
            f"Label arrays must have same length: {len(labels1)} vs {len(labels2)}"
        )

    return {
        "adjusted_rand": float(adjusted_rand_score(labels1, labels2)),
        "normalized_mutual_info": float(normalized_mutual_info_score(labels1, labels2)),
        "fowlkes_mallows": float(fowlkes_mallows_score(labels1, labels2)),
        "adjusted_mutual_info": float(adjusted_mutual_info_score(labels1, labels2)),
    }


def benchmark_clustering_methods(
    data: np.ndarray,
    methods: list[str] | None = None,
    n_clusters: int | None = None,
) -> pd.DataFrame:
    """Benchmark multiple clustering methods on the same data.

    Parameters
    ----------
    data : np.ndarray
        Data matrix (samples x features).
    methods : list of str, optional
        Clustering methods to benchmark. Defaults to
        ['kmeans', 'hierarchical', 'spectral'].
    n_clusters : int, optional
        Number of clusters. If None, uses auto-estimation.

    Returns
    -------
    pd.DataFrame
        Results with columns: method, silhouette, calinski_harabasz,
        davies_bouldin, n_clusters.

    Example
    -------
    >>> mc = MetabolicClustering(reaction_scores)
    >>> mc.compute_pca()
    >>> results = benchmark_clustering_methods(mc.pca_components, n_clusters=5)
    >>> print(results.sort_values("silhouette", ascending=False))
    """
    if methods is None:
        methods = ["kmeans", "hierarchical", "spectral"]

    results: list[dict[str, str | float | int]] = []

    for method in methods:
        try:
            # Create temporary clustering object with the data transposed back
            # We need a DataFrame, so create a dummy one
            dummy_df = pd.DataFrame(data.T)
            mc = MetabolicClustering(dummy_df, n_clusters=n_clusters)
            mc._data = data  # Set the data directly
            mc.pca_components = data  # Use data as PCA components

            labels = mc.cluster(method=method, use_pca=True)  # type: ignore[arg-type]
            metrics: dict[str, str | float | int] = dict(
                evaluate_clustering(data, labels)
            )
            metrics["method"] = method
            results.append(metrics)
        except (ImportError, ValueError) as e:
            # Skip methods with missing dependencies or errors
            results.append(
                {
                    "method": method,
                    "silhouette": float("nan"),
                    "calinski_harabasz": float("nan"),
                    "davies_bouldin": float("nan"),
                    "n_clusters": 0,
                    "error": str(e),
                }
            )

    return pd.DataFrame(results)
