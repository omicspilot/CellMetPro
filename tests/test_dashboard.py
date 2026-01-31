"""Tests for the Streamlit dashboard module.

These tests focus on the data loading and computation functions,
not the Streamlit UI components which require a running Streamlit server.
"""

import importlib.util
import json

import numpy as np
import pandas as pd
import pytest

# Skip all tests in this module if streamlit is not installed
pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("streamlit") is None,
    reason="Streamlit not installed (optional dependency)",
)


class TestLoadResults:
    """Tests for load_results function."""

    def test_load_scores(self, tmp_path):
        """Test loading reaction scores."""
        from cellmetpro.visualization.dashboard import load_results

        # Create test scores file
        scores = pd.DataFrame(
            np.random.rand(10, 20),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(20)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        data = load_results(tmp_path)

        assert "scores" in data
        assert data["scores"].shape == (10, 20)

    def test_load_config(self, tmp_path):
        """Test loading config file."""
        from cellmetpro.visualization.dashboard import load_results

        # Create test config file
        config = {
            "n_cells": 100,
            "n_genes": 500,
            "n_reactions": 50,
            "model": "human",
        }
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        data = load_results(tmp_path)

        assert "config" in data
        assert data["config"]["n_cells"] == 100
        assert data["config"]["model"] == "human"

    def test_load_differential_results(self, tmp_path):
        """Test loading differential analysis results."""
        from cellmetpro.visualization.dashboard import load_results

        # Create test differential results
        diff = pd.DataFrame(
            {
                "reaction": ["RXN1", "RXN2", "RXN3"],
                "log2fc": [1.5, -0.8, 0.2],
                "pvalue": [0.001, 0.05, 0.5],
            }
        )
        diff.to_csv(tmp_path / "differential_A_vs_B.csv", index=False)

        data = load_results(tmp_path)

        assert "differential" in data
        assert "A_vs_B" in data["differential"]
        assert len(data["differential"]["A_vs_B"]) == 3

    def test_load_clustering_results(self, tmp_path):
        """Test loading clustering results."""
        from cellmetpro.visualization.dashboard import load_results

        # Create test clustering results
        clusters = pd.DataFrame(
            {
                "cell_id": [f"cell{i}" for i in range(20)],
                "cluster": np.random.randint(0, 3, 20),
            }
        )
        clusters.to_csv(tmp_path / "clustering_results.csv", index=False)

        data = load_results(tmp_path)

        assert "clusters" in data
        assert "cluster" in data["clusters"].columns

    def test_load_enrichment_results(self, tmp_path):
        """Test loading enrichment results."""
        from cellmetpro.visualization.dashboard import load_results

        # Create test enrichment results
        enrichment = pd.DataFrame(
            {
                "pathway": ["Glycolysis", "TCA cycle", "OXPHOS"],
                "pvalue": [0.001, 0.01, 0.1],
                "padj": [0.005, 0.03, 0.2],
            }
        )
        enrichment.to_csv(tmp_path / "subsystem_enrichment.csv", index=False)

        data = load_results(tmp_path)

        assert "enrichment" in data
        assert "subsystem" in data["enrichment"]

    def test_load_empty_directory(self, tmp_path):
        """Test loading from empty directory."""
        from cellmetpro.visualization.dashboard import load_results

        data = load_results(tmp_path)
        assert data == {}


class TestComputeEmbedding:
    """Tests for compute_embedding function."""

    def test_pca_embedding(self):
        """Test PCA embedding computation."""
        from cellmetpro.visualization.dashboard import compute_embedding

        scores = pd.DataFrame(
            np.random.rand(20, 50),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(50)],
        )

        embedding = compute_embedding(scores, method="pca", n_pcs=10)

        assert embedding.shape == (50, 2)
        assert not np.any(np.isnan(embedding))

    def test_tsne_embedding(self):
        """Test t-SNE embedding computation."""
        from cellmetpro.visualization.dashboard import compute_embedding

        scores = pd.DataFrame(
            np.random.rand(20, 30),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(30)],
        )

        embedding = compute_embedding(scores, method="tsne", n_pcs=10)

        assert embedding.shape == (30, 2)
        assert not np.any(np.isnan(embedding))

    def test_embedding_with_nan_values(self):
        """Test embedding computation with NaN values."""
        from cellmetpro.visualization.dashboard import compute_embedding

        scores = pd.DataFrame(
            np.random.rand(20, 30),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(30)],
        )
        # Add some NaN values
        scores.iloc[0, 0] = np.nan
        scores.iloc[5, 10] = np.nan

        embedding = compute_embedding(scores, method="pca", n_pcs=10)

        assert embedding.shape == (30, 2)
        assert not np.any(np.isnan(embedding))

    def test_embedding_with_small_dataset(self):
        """Test embedding with very small dataset."""
        from cellmetpro.visualization.dashboard import compute_embedding

        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )

        embedding = compute_embedding(scores, method="pca", n_pcs=50)

        assert embedding.shape == (10, 2)


class TestDashboardIntegration:
    """Integration tests for dashboard data loading."""

    def test_full_results_directory(self, tmp_path):
        """Test loading a complete results directory."""
        from cellmetpro.visualization.dashboard import load_results

        # Create all result files
        n_cells, n_reactions = 50, 30

        # Scores
        scores = pd.DataFrame(
            np.random.rand(n_reactions, n_cells),
            index=[f"RXN{i}" for i in range(n_reactions)],
            columns=[f"cell{i}" for i in range(n_cells)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        # Penalties
        penalties = pd.DataFrame(
            np.random.rand(n_reactions, n_cells),
            index=[f"RXN{i}" for i in range(n_reactions)],
            columns=[f"cell{i}" for i in range(n_cells)],
        )
        penalties.to_csv(tmp_path / "reaction_penalties.csv")

        # Config
        config = {
            "n_cells": n_cells,
            "n_genes": 100,
            "n_reactions": n_reactions,
            "model": "human",
            "beta": 0.95,
        }
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        # Differential
        diff = pd.DataFrame(
            {
                "reaction": [f"RXN{i}" for i in range(10)],
                "log2fc": np.random.randn(10),
                "pvalue": np.random.rand(10) * 0.1,
                "padj_bh": np.random.rand(10) * 0.2,
            }
        )
        diff.to_csv(tmp_path / "differential_A_vs_B.csv", index=False)

        # Clustering
        clusters = pd.DataFrame(
            {
                "cell_id": [f"cell{i}" for i in range(n_cells)],
                "cluster": np.random.randint(0, 4, n_cells),
                "UMAP1": np.random.randn(n_cells),
                "UMAP2": np.random.randn(n_cells),
            }
        )
        clusters.to_csv(tmp_path / "clustering_results.csv", index=False)

        # Enrichment
        enrichment = pd.DataFrame(
            {
                "pathway": ["Glycolysis", "TCA", "OXPHOS", "PPP", "Glutamine"],
                "n_overlap": [5, 3, 8, 2, 4],
                "n_pathway": [10, 8, 15, 5, 7],
                "odds_ratio": [2.5, 1.8, 3.2, 1.2, 2.1],
                "pvalue": [0.001, 0.01, 0.0001, 0.1, 0.05],
                "padj": [0.005, 0.03, 0.001, 0.2, 0.1],
            }
        )
        enrichment.to_csv(tmp_path / "subsystem_enrichment.csv", index=False)

        # Load all data
        data = load_results(tmp_path)

        # Verify all components loaded
        assert "scores" in data
        assert "penalties" in data
        assert "config" in data
        assert "differential" in data
        assert "clusters" in data
        assert "enrichment" in data

        # Verify dimensions
        assert data["scores"].shape == (n_reactions, n_cells)
        assert data["penalties"].shape == (n_reactions, n_cells)
        assert data["config"]["n_cells"] == n_cells
        assert len(data["differential"]["A_vs_B"]) == 10
        assert len(data["clusters"]) == n_cells
        assert len(data["enrichment"]["subsystem"]) == 5
