"""Tests for the trajectory analysis module."""

import numpy as np
import pandas as pd
import pytest

from cellmetpro.analysis.trajectory import (
    cluster_trajectory_patterns,
    compute_metabolic_velocity,
    compute_pseudotime,
    fit_trajectory_genes,
    identify_branch_points,
    trajectory_differential,
)


@pytest.fixture
def trajectory_scores():
    """Create test scores with trajectory-like structure."""
    np.random.seed(42)
    n_reactions = 30
    n_cells = 100

    # Create pseudotime-dependent data
    true_pt = np.linspace(0, 1, n_cells)

    scores = pd.DataFrame(
        np.random.rand(n_reactions, n_cells) * 0.1,
        index=[f"RXN{i}" for i in range(n_reactions)],
        columns=[f"cell{i}" for i in range(n_cells)],
    )

    # Make some reactions increase with pseudotime
    for i in range(5):
        scores.iloc[i] = true_pt + np.random.rand(n_cells) * 0.1

    # Make some reactions decrease
    for i in range(5, 10):
        scores.iloc[i] = (1 - true_pt) + np.random.rand(n_cells) * 0.1

    return scores


class TestComputePseudotime:
    """Tests for compute_pseudotime function."""

    def test_basic_pseudotime(self, trajectory_scores):
        """Test basic pseudotime computation."""
        pseudotime = compute_pseudotime(trajectory_scores)

        assert isinstance(pseudotime, pd.Series)
        assert len(pseudotime) == trajectory_scores.shape[1]
        assert pseudotime.min() >= 0
        assert pseudotime.max() <= 1

    def test_with_root_cell(self, trajectory_scores):
        """Test pseudotime with specified root cell."""
        pseudotime = compute_pseudotime(trajectory_scores, root_cell="cell0")

        assert pseudotime["cell0"] == pytest.approx(0, abs=0.1)

    def test_with_root_cell_index(self, trajectory_scores):
        """Test pseudotime with root cell as index."""
        pseudotime = compute_pseudotime(trajectory_scores, root_cell=0)

        assert len(pseudotime) == trajectory_scores.shape[1]

    def test_dpt_method(self, trajectory_scores):
        """Test DPT method explicitly."""
        pseudotime = compute_pseudotime(trajectory_scores, method="dpt")

        assert isinstance(pseudotime, pd.Series)
        assert not np.any(np.isnan(pseudotime))

    def test_principal_curve_method(self, trajectory_scores):
        """Test principal curve method."""
        pseudotime = compute_pseudotime(trajectory_scores, method="principal_curve")

        assert isinstance(pseudotime, pd.Series)

    def test_correlation_method(self, trajectory_scores):
        """Test correlation-based method."""
        pseudotime = compute_pseudotime(trajectory_scores, method="correlation")

        assert isinstance(pseudotime, pd.Series)

    def test_invalid_method(self, trajectory_scores):
        """Test error on invalid method."""
        with pytest.raises(ValueError, match="Unknown method"):
            compute_pseudotime(trajectory_scores, method="invalid")


class TestComputeMetabolicVelocity:
    """Tests for compute_metabolic_velocity function."""

    def test_basic_velocity(self, trajectory_scores):
        """Test basic velocity computation."""
        pseudotime = compute_pseudotime(trajectory_scores)
        velocity = compute_metabolic_velocity(trajectory_scores, pseudotime)

        assert isinstance(velocity, pd.DataFrame)
        assert velocity.shape[0] == trajectory_scores.shape[0]

    def test_velocity_sign(self, trajectory_scores):
        """Test that velocity correctly captures direction."""
        # Use a known pseudotime ordering (by cell index) to test velocity
        n_cells = trajectory_scores.shape[1]
        true_pt = pd.Series(
            np.linspace(0, 1, n_cells),
            index=trajectory_scores.columns,
            name="pseudotime",
        )

        velocity = compute_metabolic_velocity(
            trajectory_scores, true_pt, window_size=30
        )

        # Reactions 0-4 increase with pseudotime, should have positive mean velocity
        # Reactions 5-9 decrease, should have negative mean velocity
        increasing_velocities = [velocity.iloc[i].dropna().mean() for i in range(5)]
        decreasing_velocities = [velocity.iloc[i].dropna().mean() for i in range(5, 10)]

        # At least some increasing reactions should have positive velocity
        assert any(v > 0 for v in increasing_velocities)
        # At least some decreasing reactions should have negative velocity
        assert any(v < 0 for v in decreasing_velocities)

    def test_small_window(self, trajectory_scores):
        """Test with small window size."""
        pseudotime = compute_pseudotime(trajectory_scores)
        velocity = compute_metabolic_velocity(
            trajectory_scores, pseudotime, window_size=15, min_cells=5
        )

        assert isinstance(velocity, pd.DataFrame)


class TestIdentifyBranchPoints:
    """Tests for identify_branch_points function."""

    def test_basic_branch_detection(self, trajectory_scores):
        """Test basic branch point identification."""
        pseudotime = compute_pseudotime(trajectory_scores)
        branches = identify_branch_points(trajectory_scores, pseudotime)

        assert isinstance(branches, pd.DataFrame)
        assert "cell" in branches.columns
        assert "branch_score" in branches.columns
        assert "is_branch" in branches.columns

    def test_branch_score_range(self, trajectory_scores):
        """Test that branch scores are normalized."""
        pseudotime = compute_pseudotime(trajectory_scores)
        branches = identify_branch_points(trajectory_scores, pseudotime)

        assert branches["branch_score"].min() >= 0
        assert branches["branch_score"].max() <= 1

    def test_custom_threshold(self, trajectory_scores):
        """Test with custom threshold."""
        pseudotime = compute_pseudotime(trajectory_scores)
        branches = identify_branch_points(trajectory_scores, pseudotime, threshold=0.8)

        # Higher threshold should result in fewer branch points
        n_branches_high = branches["is_branch"].sum()

        branches_low = identify_branch_points(
            trajectory_scores, pseudotime, threshold=0.2
        )
        n_branches_low = branches_low["is_branch"].sum()

        assert n_branches_high <= n_branches_low


class TestTrajectoryDifferential:
    """Tests for trajectory_differential function."""

    def test_basic_differential(self, trajectory_scores):
        """Test basic trajectory differential analysis."""
        pseudotime = compute_pseudotime(trajectory_scores)
        diff = trajectory_differential(trajectory_scores, pseudotime)

        assert isinstance(diff, pd.DataFrame)
        assert "reaction" in diff.columns
        assert "correlation" in diff.columns
        assert "pvalue" in diff.columns
        assert "trend" in diff.columns

    def test_identifies_increasing_reactions(self, trajectory_scores):
        """Test that reactions with trajectory-correlated patterns are identified."""
        # Use true pseudotime that matches the fixture data structure
        # (reactions 0-4 increase with cell index, 5-9 decrease)
        n_cells = trajectory_scores.shape[1]
        true_pt = pd.Series(
            np.linspace(0, 1, n_cells),
            index=trajectory_scores.columns,
            name="pseudotime",
        )
        diff = trajectory_differential(trajectory_scores, true_pt)

        # Check that reactions with significant trends are identified
        significant_trends = diff[diff["trend"].isin(["increasing", "decreasing"])][
            "reaction"
        ].tolist()

        # With true pseudotime matching the data, we expect strong correlations
        assert len(significant_trends) > 0

    def test_log2fc_computation(self, trajectory_scores):
        """Test log2 fold change computation."""
        pseudotime = compute_pseudotime(trajectory_scores)
        diff = trajectory_differential(trajectory_scores, pseudotime)

        assert "log2fc" in diff.columns
        assert "early_mean" in diff.columns
        assert "late_mean" in diff.columns

    def test_adjusted_pvalues(self, trajectory_scores):
        """Test adjusted p-value computation."""
        pseudotime = compute_pseudotime(trajectory_scores)
        diff = trajectory_differential(trajectory_scores, pseudotime)

        assert "padj" in diff.columns
        assert (diff["padj"] >= diff["pvalue"]).all()


class TestFitTrajectoryGenes:
    """Tests for fit_trajectory_genes function."""

    def test_basic_fitting(self, trajectory_scores):
        """Test basic trajectory fitting."""
        pseudotime = compute_pseudotime(trajectory_scores)
        fit = fit_trajectory_genes(trajectory_scores, pseudotime, n_knots=3)

        assert "fitted_values" in fit
        assert "pseudotime_grid" in fit
        assert "r_squared" in fit

    def test_specific_reactions(self, trajectory_scores):
        """Test fitting specific reactions."""
        pseudotime = compute_pseudotime(trajectory_scores)
        reactions = ["RXN0", "RXN5", "RXN10"]
        fit = fit_trajectory_genes(trajectory_scores, pseudotime, reactions=reactions)

        assert len(fit["fitted_values"].columns) <= len(reactions)

    def test_r_squared_values(self, trajectory_scores):
        """Test R-squared value computation."""
        pseudotime = compute_pseudotime(trajectory_scores)
        fit = fit_trajectory_genes(trajectory_scores, pseudotime)

        r_squared = fit["r_squared"]
        assert isinstance(r_squared, pd.Series)
        # R-squared should be between 0 and 1 for well-fitted data
        assert r_squared.max() <= 1.0


class TestClusterTrajectoryPatterns:
    """Tests for cluster_trajectory_patterns function."""

    def test_basic_clustering(self, trajectory_scores):
        """Test basic pattern clustering."""
        pseudotime = compute_pseudotime(trajectory_scores)
        fit = fit_trajectory_genes(trajectory_scores, pseudotime)
        clusters = cluster_trajectory_patterns(fit, n_clusters=3)

        assert isinstance(clusters, pd.Series)
        assert clusters.nunique() <= 3

    def test_cluster_labels(self, trajectory_scores):
        """Test that cluster labels are numeric."""
        pseudotime = compute_pseudotime(trajectory_scores)
        fit = fit_trajectory_genes(trajectory_scores, pseudotime)
        clusters = cluster_trajectory_patterns(fit, n_clusters=4)

        assert all(isinstance(c, (int, np.integer)) for c in clusters.values)
