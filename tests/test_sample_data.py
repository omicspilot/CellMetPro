"""Tests for sample data module."""

import numpy as np
import pandas as pd
import pytest


class TestLoadSampleExpression:
    """Tests for load_sample_expression function."""

    def test_default_dimensions(self):
        """Test default expression matrix dimensions."""
        from cellmetpro.data import load_sample_expression

        expression = load_sample_expression()

        assert expression.shape == (50, 100)
        assert isinstance(expression, pd.DataFrame)

    def test_custom_dimensions(self):
        """Test custom dimensions."""
        from cellmetpro.data import load_sample_expression

        expression = load_sample_expression(n_cells=50, n_genes=25)

        assert expression.shape == (25, 50)

    def test_reproducibility(self):
        """Test that same seed produces same data."""
        from cellmetpro.data import load_sample_expression

        expr1 = load_sample_expression(seed=123)
        expr2 = load_sample_expression(seed=123)

        pd.testing.assert_frame_equal(expr1, expr2)

    def test_different_seeds(self):
        """Test that different seeds produce different data."""
        from cellmetpro.data import load_sample_expression

        expr1 = load_sample_expression(seed=1)
        expr2 = load_sample_expression(seed=2)

        assert not expr1.equals(expr2)

    def test_non_negative_values(self):
        """Test that all values are non-negative."""
        from cellmetpro.data import load_sample_expression

        expression = load_sample_expression()

        assert (expression >= 0).all().all()

    def test_has_metabolic_genes(self):
        """Test that expression includes metabolic pathway genes."""
        from cellmetpro.data import load_sample_expression

        expression = load_sample_expression()

        # Should include common metabolic genes
        metabolic_genes = ["HK1", "HK2", "GAPDH", "PKM", "CS", "NDUFA1"]
        found_genes = [g for g in metabolic_genes if g in expression.index]

        assert len(found_genes) > 0

    def test_cell_names_have_cluster_info(self):
        """Test that cell names include cluster information."""
        from cellmetpro.data import load_sample_expression

        expression = load_sample_expression()

        # Check for expected cluster prefixes
        expected_clusters = ["Proliferating", "Quiescent", "Hypoxic", "Oxidative"]
        cell_prefixes = [c.split("_")[0] for c in expression.columns]

        for cluster in expected_clusters:
            assert cluster in cell_prefixes


class TestLoadSampleGroups:
    """Tests for load_sample_groups function."""

    def test_default_dimensions(self):
        """Test default group annotations dimensions."""
        from cellmetpro.data import load_sample_groups

        groups = load_sample_groups()

        assert len(groups) == 100
        assert isinstance(groups, pd.DataFrame)

    def test_required_columns(self):
        """Test that required columns are present."""
        from cellmetpro.data import load_sample_groups

        groups = load_sample_groups()

        required_cols = ["cell", "group", "cell_type", "treatment"]
        for col in required_cols:
            assert col in groups.columns

    def test_matching_cells(self):
        """Test that group cell names match expression cell names."""
        from cellmetpro.data import load_sample_expression, load_sample_groups

        expression = load_sample_expression(n_cells=100, seed=42)
        groups = load_sample_groups(n_cells=100, seed=42)

        assert set(expression.columns) == set(groups["cell"])

    def test_balanced_treatments(self):
        """Test that treatments are approximately balanced."""
        from cellmetpro.data import load_sample_groups

        groups = load_sample_groups()
        treatment_counts = groups["treatment"].value_counts()

        # Should be roughly 50/50
        assert abs(treatment_counts["control"] - treatment_counts["treatment"]) <= 10

    def test_cell_types_present(self):
        """Test that expected cell types are present."""
        from cellmetpro.data import load_sample_groups

        groups = load_sample_groups()
        cell_types = groups["cell_type"].unique()

        expected_types = ["Proliferating", "Quiescent", "Hypoxic", "Oxidative"]
        for ct in expected_types:
            assert ct in cell_types


class TestLoadSampleReactionScores:
    """Tests for load_sample_reaction_scores function."""

    def test_default_dimensions(self):
        """Test default reaction scores dimensions."""
        from cellmetpro.data import load_sample_reaction_scores

        scores = load_sample_reaction_scores()

        assert scores.shape == (30, 100)
        assert isinstance(scores, pd.DataFrame)

    def test_custom_dimensions(self):
        """Test custom dimensions."""
        from cellmetpro.data import load_sample_reaction_scores

        scores = load_sample_reaction_scores(n_cells=50, n_reactions=20)

        assert scores.shape == (20, 50)

    def test_score_range(self):
        """Test that scores are in valid range [0, 1]."""
        from cellmetpro.data import load_sample_reaction_scores

        scores = load_sample_reaction_scores()

        assert (scores >= 0).all().all()
        assert (scores <= 1).all().all()

    def test_has_metabolic_reactions(self):
        """Test that scores include metabolic reactions."""
        from cellmetpro.data import load_sample_reaction_scores

        scores = load_sample_reaction_scores()

        # Should include common reactions
        expected_rxns = ["HEX1", "PGI", "PFK", "CS", "GLCt1"]
        found_rxns = [r for r in expected_rxns if r in scores.index]

        assert len(found_rxns) > 0

    def test_matching_cells(self):
        """Test that cells match between scores and groups."""
        from cellmetpro.data import load_sample_groups, load_sample_reaction_scores

        scores = load_sample_reaction_scores(n_cells=100, seed=42)
        groups = load_sample_groups(n_cells=100, seed=42)

        assert set(scores.columns) == set(groups["cell"])


class TestCreateSampleModel:
    """Tests for create_sample_model function."""

    def test_model_creation(self):
        """Test that model is created successfully."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()

        assert model is not None
        assert model.id == "sample_metabolic_model"

    def test_model_has_reactions(self):
        """Test that model has reactions."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()

        assert len(model.reactions) > 5

    def test_model_has_genes(self):
        """Test that model has genes with GPR rules."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()

        assert len(model.genes) > 0

        # Check that some reactions have GPR rules
        rxns_with_gpr = [r for r in model.reactions if r.gene_reaction_rule]
        assert len(rxns_with_gpr) > 0

    def test_model_has_metabolites(self):
        """Test that model has metabolites."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()

        assert len(model.metabolites) > 0

    def test_model_optimizes(self):
        """Test that model can be optimized."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()
        solution = model.optimize()

        assert solution.status == "optimal"
        assert solution.objective_value >= 0

    def test_model_has_exchange_reactions(self):
        """Test that model has exchange reactions."""
        from cellmetpro.data import create_sample_model

        model = create_sample_model()

        exchanges = model.exchanges
        assert len(exchanges) > 0


class TestGetSampleDataPath:
    """Tests for get_sample_data_path function."""

    def test_returns_path(self):
        """Test that function returns a Path object."""
        from pathlib import Path

        from cellmetpro.data import get_sample_data_path

        path = get_sample_data_path()

        assert isinstance(path, Path)
        assert path.exists()


class TestSampleDataIntegration:
    """Integration tests for sample data with CellMetPro modules."""

    def test_expression_with_compass(self):
        """Test that sample expression works with CompassScorer."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer
        from cellmetpro.data import create_sample_model, load_sample_expression

        expression = load_sample_expression(n_cells=10, n_genes=30)
        model = create_sample_model()

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(model, expression, config)
        penalties = scorer.compute_reaction_penalties()

        assert penalties is not None
        assert not penalties.empty

    def test_scores_with_differential_analysis(self):
        """Test that sample scores work with DifferentialAnalysis."""
        from cellmetpro.analysis.differential import DifferentialAnalysis
        from cellmetpro.data import load_sample_groups, load_sample_reaction_scores

        scores = load_sample_reaction_scores(n_cells=100)
        groups = load_sample_groups(n_cells=100)

        # Create groups Series
        group_series = groups.set_index("cell")["cell_type"]

        da = DifferentialAnalysis(scores, group_series)
        results = da.compare_groups("Proliferating", "Quiescent")

        assert results is not None
        assert len(results) > 0

    def test_scores_with_clustering(self):
        """Test that sample scores work with MetabolicClustering."""
        from cellmetpro.analysis.clustering import MetabolicClustering
        from cellmetpro.data import load_sample_reaction_scores

        scores = load_sample_reaction_scores(n_cells=50, n_reactions=20)

        mc = MetabolicClustering(scores, n_clusters=3)
        mc.compute_pca(n_components=10)
        labels = mc.cluster(method="kmeans")

        assert labels is not None
        assert len(labels) == 50
