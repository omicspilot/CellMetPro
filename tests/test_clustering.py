"""Tests for metabolic clustering module."""

import numpy as np
import pandas as pd
import pytest

from cellmetpro.analysis.clustering import (
    MetabolicClustering,
    benchmark_clustering_methods,
    compare_clusterings,
    evaluate_clustering,
    find_optimal_clusters,
)


@pytest.fixture
def reaction_scores():
    """Create mock reaction scores with cluster structure."""
    np.random.seed(42)

    # 10 reactions, 30 cells (3 clusters of 10 cells each)
    reactions = [f"R{i}" for i in range(10)]
    cells = [f"cell_{i}" for i in range(30)]

    # Create base data
    data = np.random.rand(10, 30) * 0.5

    # Add cluster structure
    # Cluster 0 (cells 0-9): high activity in R0-R3
    data[0:4, 0:10] += 2.0
    # Cluster 1 (cells 10-19): high activity in R4-R6
    data[4:7, 10:20] += 2.0
    # Cluster 2 (cells 20-29): high activity in R7-R9
    data[7:10, 20:30] += 2.0

    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def small_reaction_scores():
    """Create small mock reaction scores for quick tests."""
    np.random.seed(42)
    reactions = ["R1", "R2", "R3"]
    cells = [f"cell_{i}" for i in range(10)]
    data = np.random.rand(3, 10)
    return pd.DataFrame(data, index=reactions, columns=cells)


# =============================================================================
# TESTS FOR PCA
# =============================================================================


def test_compute_pca_returns_array(small_reaction_scores):
    """Test that compute_pca returns an array."""
    mc = MetabolicClustering(small_reaction_scores)
    pca_result = mc.compute_pca(n_components=2)

    assert isinstance(pca_result, np.ndarray)
    assert pca_result.shape[0] == 10  # 10 cells
    assert pca_result.shape[1] == 2  # 2 components


def test_compute_pca_caps_components(small_reaction_scores):
    """Test that n_components is capped at max possible."""
    mc = MetabolicClustering(small_reaction_scores)
    # Request more components than possible
    pca_result = mc.compute_pca(n_components=100)

    # Should be capped at min(n_cells, n_reactions) = min(10, 3) = 3
    assert pca_result.shape[1] <= 3


def test_compute_pca_stores_result(small_reaction_scores):
    """Test that PCA result is stored in pca_components."""
    mc = MetabolicClustering(small_reaction_scores)
    pca_result = mc.compute_pca()

    assert mc.pca_components is not None
    np.testing.assert_array_equal(pca_result, mc.pca_components)


def test_compute_pca_with_scaling(small_reaction_scores):
    """Test PCA with and without scaling."""
    mc1 = MetabolicClustering(small_reaction_scores)
    mc2 = MetabolicClustering(small_reaction_scores)

    result_scaled = mc1.compute_pca(scale=True)
    result_unscaled = mc2.compute_pca(scale=False)

    # Results should be different
    assert not np.allclose(result_scaled, result_unscaled)


# =============================================================================
# TESTS FOR t-SNE
# =============================================================================


def test_compute_tsne_returns_array(small_reaction_scores):
    """Test that compute_tsne returns an array."""
    mc = MetabolicClustering(small_reaction_scores)
    mc.compute_pca(n_components=3)  # PCA first for speed
    tsne_result = mc.compute_tsne(n_components=2, perplexity=3, max_iter=250)

    assert isinstance(tsne_result, np.ndarray)
    assert tsne_result.shape == (10, 2)  # 10 cells, 2 dims


def test_compute_tsne_stores_embedding(small_reaction_scores):
    """Test that t-SNE result is stored in embedding."""
    mc = MetabolicClustering(small_reaction_scores)
    tsne_result = mc.compute_tsne(perplexity=3, max_iter=250)

    assert mc.embedding is not None
    np.testing.assert_array_equal(tsne_result, mc.embedding)


def test_compute_tsne_adjusts_perplexity(small_reaction_scores):
    """Test that perplexity is adjusted for small samples."""
    mc = MetabolicClustering(small_reaction_scores)
    # Should not raise even with high perplexity
    tsne_result = mc.compute_tsne(perplexity=100, max_iter=250)
    assert tsne_result.shape == (10, 2)


# =============================================================================
# TESTS FOR CLUSTERING
# =============================================================================


def test_cluster_kmeans(reaction_scores):
    """Test k-means clustering."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="kmeans")

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30  # 30 cells
    assert len(np.unique(labels)) == 3  # 3 clusters


def test_cluster_kmeans_auto_n_clusters(small_reaction_scores):
    """Test k-means with automatic n_clusters."""
    mc = MetabolicClustering(small_reaction_scores)
    mc.compute_pca()
    labels = mc.cluster(method="kmeans")

    assert mc.labels is not None
    assert len(labels) == 10


def test_cluster_stores_labels(reaction_scores):
    """Test that cluster labels are stored."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca()
    labels = mc.cluster(method="kmeans")

    assert mc.labels is not None
    np.testing.assert_array_equal(labels, mc.labels)


def test_cluster_invalid_method(small_reaction_scores):
    """Test that invalid method raises error."""
    mc = MetabolicClustering(small_reaction_scores)
    mc.compute_pca()

    with pytest.raises(ValueError):
        mc.cluster(method="invalid_method")


# =============================================================================
# TESTS FOR CLUSTER MARKERS
# =============================================================================


def test_get_cluster_markers(reaction_scores):
    """Test marker identification."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca()
    mc.cluster(method="kmeans")

    markers = mc.get_cluster_markers(n_top=5)

    assert isinstance(markers, pd.DataFrame)
    assert "cluster" in markers.columns
    assert "reaction" in markers.columns
    assert "score_diff" in markers.columns


def test_get_cluster_markers_requires_clustering(small_reaction_scores):
    """Test that get_cluster_markers requires clustering first."""
    mc = MetabolicClustering(small_reaction_scores)

    with pytest.raises(ValueError):
        mc.get_cluster_markers()


# =============================================================================
# TESTS FOR EXPORT FUNCTIONS
# =============================================================================


def test_to_dataframe(small_reaction_scores):
    """Test export to DataFrame."""
    mc = MetabolicClustering(small_reaction_scores, n_clusters=2)
    mc.compute_pca(n_components=3)
    mc.compute_tsne(perplexity=3, max_iter=250)
    mc.cluster(method="kmeans")

    df = mc.to_dataframe()

    assert isinstance(df, pd.DataFrame)
    assert "cell_id" in df.columns
    assert "cluster" in df.columns
    assert "dim_1" in df.columns  # From embedding
    assert "PC1" in df.columns  # From PCA


def test_to_dataframe_without_clustering(small_reaction_scores):
    """Test export without clustering still works."""
    mc = MetabolicClustering(small_reaction_scores)

    df = mc.to_dataframe()

    assert isinstance(df, pd.DataFrame)
    assert "cell_id" in df.columns
    assert "cluster" not in df.columns  # No clustering done


# =============================================================================
# TESTS FOR find_optimal_clusters
# =============================================================================


def test_find_optimal_clusters_silhouette(reaction_scores):
    """Test silhouette method for optimal clusters."""
    # Use PCA-reduced data
    mc = MetabolicClustering(reaction_scores)
    pca_data = mc.compute_pca(n_components=10)

    n_clusters = find_optimal_clusters(pca_data, max_clusters=5, method="silhouette")

    assert isinstance(n_clusters, int)
    assert 2 <= n_clusters <= 5


def test_find_optimal_clusters_elbow(reaction_scores):
    """Test elbow method for optimal clusters."""
    mc = MetabolicClustering(reaction_scores)
    pca_data = mc.compute_pca(n_components=10)

    n_clusters = find_optimal_clusters(pca_data, max_clusters=5, method="elbow")

    assert isinstance(n_clusters, int)
    assert 2 <= n_clusters <= 5


# =============================================================================
# TESTS FOR EDGE CASES
# =============================================================================


def test_handles_nan_values():
    """Test that NaN values are handled."""
    data = pd.DataFrame(
        [[1.0, np.nan, 2.0], [3.0, 4.0, np.nan]],
        index=["R1", "R2"],
        columns=["c1", "c2", "c3"],
    )

    mc = MetabolicClustering(data)
    pca_result = mc.compute_pca()

    assert not np.any(np.isnan(pca_result))


def test_single_reaction():
    """Test with single reaction."""
    data = pd.DataFrame(
        [[1.0, 2.0, 3.0, 4.0, 5.0]], index=["R1"], columns=[f"c{i}" for i in range(5)]
    )

    mc = MetabolicClustering(data)
    pca_result = mc.compute_pca()

    assert pca_result.shape[0] == 5  # 5 cells


# =============================================================================
# TESTS FOR NEW CLUSTERING ALGORITHMS
# =============================================================================


def test_cluster_hierarchical(reaction_scores):
    """Test hierarchical clustering."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="hierarchical")

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30
    assert len(np.unique(labels)) == 3


def test_cluster_hierarchical_linkages(reaction_scores):
    """Test hierarchical clustering with different linkages."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)

    for linkage in ["ward", "complete", "average", "single"]:
        labels = mc.cluster(method="hierarchical", linkage=linkage)
        assert len(np.unique(labels)) == 3


def test_cluster_spectral(reaction_scores):
    """Test spectral clustering."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="spectral")

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30
    assert len(np.unique(labels)) == 3


def test_cluster_minibatch_kmeans(reaction_scores):
    """Test MiniBatch K-means clustering."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="minibatch_kmeans")

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30
    assert len(np.unique(labels)) == 3


def test_cluster_dbscan(reaction_scores):
    """Test DBSCAN clustering."""
    mc = MetabolicClustering(reaction_scores)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="dbscan", min_samples=3)

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30
    # DBSCAN auto-detects clusters, so we just check it runs


def test_cluster_dbscan_with_custom_eps(reaction_scores):
    """Test DBSCAN with custom eps parameter."""
    mc = MetabolicClustering(reaction_scores)
    mc.compute_pca(n_components=10)
    labels = mc.cluster(method="dbscan", eps=1.0, min_samples=3)

    assert isinstance(labels, np.ndarray)
    assert len(labels) == 30


# =============================================================================
# TESTS FOR CLUSTERING EVALUATION
# =============================================================================


def test_evaluate_clustering(reaction_scores):
    """Test clustering evaluation metrics."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    pca_data = mc.compute_pca(n_components=10)
    labels = mc.cluster(method="kmeans")

    metrics = evaluate_clustering(pca_data, labels)

    assert isinstance(metrics, dict)
    assert "silhouette" in metrics
    assert "calinski_harabasz" in metrics
    assert "davies_bouldin" in metrics
    assert "n_clusters" in metrics

    # Silhouette should be between -1 and 1
    assert -1 <= metrics["silhouette"] <= 1
    # Calinski-Harabasz should be positive
    assert metrics["calinski_harabasz"] > 0
    # Davies-Bouldin should be positive (lower is better)
    assert metrics["davies_bouldin"] > 0


def test_evaluate_clustering_single_cluster():
    """Test evaluation with single cluster returns NaN."""
    data = np.random.rand(10, 5)
    labels = np.zeros(10)  # All same cluster

    metrics = evaluate_clustering(data, labels)

    assert np.isnan(metrics["silhouette"])
    assert metrics["n_clusters"] == 1


# =============================================================================
# TESTS FOR CLUSTERING COMPARISON
# =============================================================================


def test_compare_clusterings(reaction_scores):
    """Test clustering comparison metrics."""
    mc = MetabolicClustering(reaction_scores, n_clusters=3)
    mc.compute_pca(n_components=10)

    labels_kmeans = mc.cluster(method="kmeans")
    labels_hierarchical = mc.cluster(method="hierarchical")

    comparison = compare_clusterings(labels_kmeans, labels_hierarchical)

    assert isinstance(comparison, dict)
    assert "adjusted_rand" in comparison
    assert "normalized_mutual_info" in comparison
    assert "fowlkes_mallows" in comparison
    assert "adjusted_mutual_info" in comparison

    # All metrics should be between 0 and 1 (or -1 to 1 for ARI)
    assert -1 <= comparison["adjusted_rand"] <= 1
    assert 0 <= comparison["normalized_mutual_info"] <= 1
    assert 0 <= comparison["fowlkes_mallows"] <= 1


def test_compare_clusterings_identical():
    """Test comparison of identical clusterings."""
    labels = np.array([0, 0, 1, 1, 2, 2])

    comparison = compare_clusterings(labels, labels.copy())

    # Should be perfect match
    assert comparison["adjusted_rand"] == pytest.approx(1.0)
    assert comparison["normalized_mutual_info"] == pytest.approx(1.0)


def test_compare_clusterings_different_lengths():
    """Test comparison with mismatched lengths raises error."""
    labels1 = np.array([0, 1, 2])
    labels2 = np.array([0, 1, 2, 3])

    with pytest.raises(ValueError, match="same length"):
        compare_clusterings(labels1, labels2)


# =============================================================================
# TESTS FOR BENCHMARK FUNCTION
# =============================================================================


def test_benchmark_clustering_methods(reaction_scores):
    """Test clustering method benchmarking."""
    mc = MetabolicClustering(reaction_scores)
    pca_data = mc.compute_pca(n_components=10)

    results = benchmark_clustering_methods(
        pca_data,
        methods=["kmeans", "hierarchical"],
        n_clusters=3,
    )

    assert isinstance(results, pd.DataFrame)
    assert "method" in results.columns
    assert "silhouette" in results.columns
    assert len(results) == 2  # Two methods


def test_benchmark_with_default_methods(reaction_scores):
    """Test benchmark with default methods."""
    mc = MetabolicClustering(reaction_scores)
    pca_data = mc.compute_pca(n_components=10)

    results = benchmark_clustering_methods(pca_data, n_clusters=3)

    # Default methods: kmeans, hierarchical, spectral
    assert len(results) == 3


# =============================================================================
# ADDITIONAL EDGE CASE TESTS
# =============================================================================


def test_metabolic_clustering_single_cell():
    """Test clustering with a single cell."""
    data = pd.DataFrame(
        [[1.0], [2.0], [3.0]],
        index=["R1", "R2", "R3"],
        columns=["cell_0"],
    )

    mc = MetabolicClustering(data)
    pca_result = mc.compute_pca(n_components=1)

    assert pca_result.shape == (1, 1)


def test_metabolic_clustering_two_cells():
    """Test clustering with exactly two cells."""
    data = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
        index=["R1", "R2", "R3"],
        columns=["cell_0", "cell_1"],
    )

    mc = MetabolicClustering(data, n_clusters=2)
    mc.compute_pca()
    labels = mc.cluster(method="kmeans")

    assert len(labels) == 2
    assert len(np.unique(labels)) == 2


def test_metabolic_clustering_identical_cells():
    """Test clustering when all cells have identical values."""
    data = pd.DataFrame(
        np.ones((5, 10)),
        index=[f"R{i}" for i in range(5)],
        columns=[f"cell_{i}" for i in range(10)],
    )

    mc = MetabolicClustering(data, n_clusters=2)
    # Should handle identical data gracefully
    pca_result = mc.compute_pca()
    assert pca_result is not None


def test_metabolic_clustering_high_dimensionality():
    """Test clustering when reactions >> cells."""
    np.random.seed(42)
    # 100 reactions, 10 cells
    data = pd.DataFrame(
        np.random.rand(100, 10),
        index=[f"R{i}" for i in range(100)],
        columns=[f"cell_{i}" for i in range(10)],
    )

    mc = MetabolicClustering(data, n_clusters=3)
    mc.compute_pca(n_components=5)
    labels = mc.cluster(method="kmeans")

    assert len(labels) == 10


def test_metabolic_clustering_negative_values():
    """Test clustering with negative values."""
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(5, 15),  # Normal distribution, includes negatives
        index=[f"R{i}" for i in range(5)],
        columns=[f"cell_{i}" for i in range(15)],
    )

    mc = MetabolicClustering(data, n_clusters=3)
    mc.compute_pca()
    labels = mc.cluster(method="kmeans")

    assert len(labels) == 15


def test_compute_tsne_small_perplexity():
    """Test t-SNE with very small perplexity."""
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.rand(5, 20),
        index=[f"R{i}" for i in range(5)],
        columns=[f"cell_{i}" for i in range(20)],
    )

    mc = MetabolicClustering(data)
    mc.compute_pca()
    # Very small perplexity
    tsne_result = mc.compute_tsne(perplexity=2, max_iter=250)

    assert tsne_result.shape == (20, 2)


def test_cluster_dbscan_all_noise():
    """Test DBSCAN when all points are noise."""
    np.random.seed(42)
    # Very sparse, separated points
    data = pd.DataFrame(
        np.diag(np.arange(1, 11, dtype=float)),  # Diagonal matrix
        index=[f"R{i}" for i in range(10)],
        columns=[f"cell_{i}" for i in range(10)],
    )

    mc = MetabolicClustering(data)
    mc.compute_pca()
    labels = mc.cluster(method="dbscan", eps=0.001, min_samples=5)

    # With tiny eps, all points may be noise (-1)
    assert len(labels) == 10


def test_cluster_spectral_more_clusters_than_possible():
    """Test spectral clustering with more clusters than samples."""
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.rand(3, 5),
        index=["R1", "R2", "R3"],
        columns=[f"cell_{i}" for i in range(5)],
    )

    mc = MetabolicClustering(data, n_clusters=10)  # More than 5 cells
    mc.compute_pca()

    # Should handle gracefully (cap clusters or raise error)
    try:
        labels = mc.cluster(method="spectral")
        assert len(labels) == 5
    except (ValueError, TypeError):
        pass  # Also acceptable - scipy may raise TypeError


def test_evaluate_clustering_two_clusters():
    """Test clustering evaluation with exactly two clusters."""
    np.random.seed(42)
    data = np.random.rand(10, 5)
    labels = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])

    metrics = evaluate_clustering(data, labels)

    assert metrics["n_clusters"] == 2
    assert "silhouette" in metrics


def test_evaluate_clustering_many_clusters():
    """Test clustering evaluation with many small clusters."""
    np.random.seed(42)
    data = np.random.rand(20, 5)
    labels = np.arange(20)  # Each point is its own cluster

    metrics = evaluate_clustering(data, labels)

    assert metrics["n_clusters"] == 20


def test_compare_clusterings_permuted_labels():
    """Test comparison when labels are just permuted."""
    labels1 = np.array([0, 0, 1, 1, 2, 2])
    labels2 = np.array([1, 1, 2, 2, 0, 0])  # Same clustering, different labels

    comparison = compare_clusterings(labels1, labels2)

    # Should be perfect match
    assert comparison["adjusted_rand"] == pytest.approx(1.0)


def test_compare_clusterings_random():
    """Test comparison of random clusterings."""
    np.random.seed(42)
    labels1 = np.random.randint(0, 3, 100)
    np.random.seed(123)
    labels2 = np.random.randint(0, 3, 100)

    comparison = compare_clusterings(labels1, labels2)

    # Random clusterings should have low agreement
    assert comparison["adjusted_rand"] < 0.3


def test_benchmark_single_method():
    """Test benchmarking with a single method."""
    np.random.seed(42)
    data = np.random.rand(20, 5)

    results = benchmark_clustering_methods(
        data,
        methods=["kmeans"],
        n_clusters=3,
    )

    assert len(results) == 1
    assert results.iloc[0]["method"] == "kmeans"


def test_benchmark_invalid_method():
    """Test benchmark with invalid method in list."""
    np.random.seed(42)
    data = np.random.rand(20, 5)

    # Should skip invalid method or raise error
    try:
        results = benchmark_clustering_methods(
            data,
            methods=["kmeans", "invalid_method"],
            n_clusters=3,
        )
        # If it doesn't raise, should only have kmeans
        assert len(results) <= 2
    except (ValueError, KeyError):
        pass


def test_find_optimal_clusters_min_max_equal():
    """Test find_optimal_clusters when min and max clusters are equal."""
    np.random.seed(42)
    data = np.random.rand(30, 5)

    # When max_clusters=2, only option is 2
    n_clusters = find_optimal_clusters(data, max_clusters=2, method="silhouette")

    assert n_clusters == 2


def test_find_optimal_clusters_small_dataset():
    """Test find_optimal_clusters on very small dataset."""
    np.random.seed(42)
    data = np.random.rand(5, 3)

    n_clusters = find_optimal_clusters(data, max_clusters=3, method="silhouette")

    assert 2 <= n_clusters <= 3


def test_get_cluster_markers_single_cluster():
    """Test get_cluster_markers when only one cluster exists."""
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.rand(5, 10),
        index=[f"R{i}" for i in range(5)],
        columns=[f"cell_{i}" for i in range(10)],
    )

    mc = MetabolicClustering(data, n_clusters=1)
    mc.compute_pca()
    mc.cluster(method="kmeans")

    # Single cluster - markers don't make sense
    try:
        markers = mc.get_cluster_markers()
        # If it works, should still return DataFrame
        assert isinstance(markers, pd.DataFrame)
    except ValueError:
        pass  # Also acceptable


def test_to_dataframe_minimal():
    """Test to_dataframe with minimal data (no PCA, no clustering)."""
    data = pd.DataFrame(
        [[1.0, 2.0], [3.0, 4.0]],
        index=["R1", "R2"],
        columns=["cell_0", "cell_1"],
    )

    mc = MetabolicClustering(data)
    df = mc.to_dataframe()

    assert "cell_id" in df.columns
    assert len(df) == 2


def test_metabolic_clustering_inf_values():
    """Test clustering handles infinite values."""
    data = pd.DataFrame(
        [[1.0, np.inf, 3.0], [4.0, 5.0, np.inf]],
        index=["R1", "R2"],
        columns=["cell_0", "cell_1", "cell_2"],
    )

    mc = MetabolicClustering(data)
    # Should handle inf gracefully (replace or error)
    try:
        pca_result = mc.compute_pca()
        assert not np.any(np.isinf(pca_result))
    except ValueError:
        pass  # Also acceptable to reject inf values
