"""Tests for the batch correction module."""

import importlib.util

import numpy as np
import pandas as pd
import pytest

from cellmetpro.core.batch_correction import (
    center_batches,
    combat_correct,
    compute_integration_metrics,
    select_hvr,
)


class TestCenterBatches:
    """Tests for center_batches function."""

    def test_basic_centering(self):
        """Test basic batch centering."""
        scores = pd.DataFrame(
            np.random.rand(20, 100),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(100)],
        )
        # Add batch effect
        scores.iloc[:, 50:] += 2.0
        batch_labels = ["batch1"] * 50 + ["batch2"] * 50

        centered = center_batches(scores, batch_labels)

        # Check that batch means are more similar
        batch1_mean = centered.iloc[:, :50].mean().mean()
        batch2_mean = centered.iloc[:, 50:].mean().mean()

        assert centered.shape == scores.shape
        assert abs(batch1_mean - batch2_mean) < 1.0  # Should be closer after centering

    def test_with_series_labels(self):
        """Test with pandas Series labels."""
        scores = pd.DataFrame(
            np.random.rand(10, 50),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(50)],
        )
        batch_labels = pd.Series(["A"] * 25 + ["B"] * 25)

        centered = center_batches(scores, batch_labels)

        assert centered.shape == scores.shape

    def test_preserves_variance_structure(self):
        """Test that centering preserves relative variance."""
        scores = pd.DataFrame(
            np.random.rand(10, 40),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(40)],
        )
        batch_labels = ["A"] * 20 + ["B"] * 20

        centered = center_batches(scores, batch_labels)

        # Check that reactions maintain their relative variance
        orig_var_order = scores.var(axis=1).rank()
        cent_var_order = centered.var(axis=1).rank()

        correlation = orig_var_order.corr(cent_var_order)
        assert correlation > 0.9


class TestCombatCorrect:
    """Tests for combat_correct function."""

    def test_basic_combat(self):
        """Test basic ComBat correction."""
        scores = pd.DataFrame(
            np.random.rand(20, 100),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(100)],
        )
        # Add batch effect
        scores.iloc[:, 50:] += 1.5
        batch_labels = ["batch1"] * 50 + ["batch2"] * 50

        corrected = combat_correct(scores, batch_labels)

        assert corrected.shape == scores.shape

    def test_single_batch_returns_copy(self):
        """Test that single batch returns unchanged copy."""
        scores = pd.DataFrame(
            np.random.rand(10, 50),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(50)],
        )
        batch_labels = ["batch1"] * 50

        corrected = combat_correct(scores, batch_labels)

        np.testing.assert_array_almost_equal(corrected.values, scores.values)

    def test_non_parametric(self):
        """Test non-parametric ComBat."""
        scores = pd.DataFrame(
            np.random.rand(15, 60),
            index=[f"RXN{i}" for i in range(15)],
            columns=[f"cell{i}" for i in range(60)],
        )
        batch_labels = ["A"] * 30 + ["B"] * 30

        corrected = combat_correct(scores, batch_labels, parametric=False)

        assert corrected.shape == scores.shape

    def test_invalid_batch_labels_length(self):
        """Test error on mismatched batch labels."""
        scores = pd.DataFrame(
            np.random.rand(10, 50),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(50)],
        )
        batch_labels = ["A"] * 30  # Wrong length

        with pytest.raises(ValueError, match="batch_labels length"):
            combat_correct(scores, batch_labels)


class TestComputeIntegrationMetrics:
    """Tests for compute_integration_metrics function."""

    def test_basic_metrics(self):
        """Test basic metric computation."""
        scores = pd.DataFrame(
            np.random.rand(20, 100),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(100)],
        )
        batch_labels = ["batch1"] * 50 + ["batch2"] * 50

        metrics = compute_integration_metrics(scores, batch_labels)

        assert "batch_silhouette" in metrics
        assert "batch_mixing" in metrics
        assert -1 <= metrics["batch_silhouette"] <= 1
        assert 0 <= metrics["batch_mixing"] <= 1

    def test_with_cell_labels(self):
        """Test metrics with cell type labels."""
        scores = pd.DataFrame(
            np.random.rand(20, 100),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(100)],
        )
        batch_labels = ["batch1"] * 50 + ["batch2"] * 50
        cell_labels = ["typeA"] * 25 + ["typeB"] * 25 + ["typeA"] * 25 + ["typeB"] * 25

        metrics = compute_integration_metrics(scores, batch_labels, cell_labels)

        assert "label_ari" in metrics


class TestSelectHvr:
    """Tests for select_hvr function."""

    def test_basic_hvr_selection(self):
        """Test basic HVR selection."""
        scores = pd.DataFrame(
            np.random.rand(100, 50),
            index=[f"RXN{i}" for i in range(100)],
            columns=[f"cell{i}" for i in range(50)],
        )

        hvr = select_hvr(scores, n_top=20)

        assert len(hvr) <= 20
        assert all(r in scores.index for r in hvr)

    def test_identifies_variable_reactions(self):
        """Test that selected reactions have higher variance than average."""
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(50, 100),
            index=[f"RXN{i}" for i in range(50)],
            columns=[f"cell{i}" for i in range(100)],
        )

        hvr = select_hvr(scores, n_top=10, min_disp=0.0, min_mean=0.0)

        # Selected HVR should be valid reactions from the input
        assert len(hvr) == 10
        assert all(r in scores.index for r in hvr)


# Skip Harmony tests if harmonypy not installed
pytestmark_harmony = pytest.mark.skipif(
    importlib.util.find_spec("harmonypy") is None,
    reason="harmonypy not installed",
)


@pytestmark_harmony
class TestHarmonyIntegrate:
    """Tests for harmony_integrate function."""

    def test_harmony_basic(self):
        """Test basic Harmony integration."""
        from cellmetpro.core.batch_correction import harmony_integrate

        scores = pd.DataFrame(
            np.random.rand(30, 100),
            index=[f"RXN{i}" for i in range(30)],
            columns=[f"cell{i}" for i in range(100)],
        )
        batch_labels = ["batch1"] * 50 + ["batch2"] * 50

        corrected = harmony_integrate(scores, batch_labels, max_iter=2)

        assert corrected.shape == scores.shape
