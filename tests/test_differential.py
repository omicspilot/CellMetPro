"""Tests for differential analysis module."""

import numpy as np
import pandas as pd
import pytest

from cellmetpro.analysis.differential import DifferentialAnalysis, PseudoBulkAnalysis


@pytest.fixture
def reaction_scores():
    """Create mock reaction scores (reactions x cells)."""
    np.random.seed(42)

    # 5 reactions, 10 cells
    reactions = ["R1", "R2", "R3", "R4", "R5"]
    cells = [f"cell_{i}" for i in range(10)]

    # Create data where some reactions differ between groups
    data = np.random.rand(5, 10)

    # Make R1 clearly higher in group B (cells 5-9)
    data[0, 5:10] += 2.0

    # Make R2 clearly higher in group A (cells 0-4)
    data[1, 0:5] += 2.0

    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def group_labels():
    """Create group labels for cells."""
    cells = [f"cell_{i}" for i in range(10)]
    groups = ["A"] * 5 + ["B"] * 5
    return pd.Series(groups, index=cells)


## TEST A: Basic output structure


def test_compare_groups_returns_dataframe(reaction_scores, group_labels):
    """Test that compare_groups returns a DataFrame with expected columns."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Check it's a DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check expected columns exist
    expected_cols = [
        "reaction",
        "group1_mean",
        "group2_mean",
        "log2fc",
        "statistic",
        "pvalue",
        "padj_bh",
        "padj_bonf",
    ]
    for col in expected_cols:
        assert col in result.columns

    # Check number of rows equals number of reactions
    assert len(result) == len(reaction_scores)


## TEST B: Results sorted by p-value


def test_compare_groups_sorted_by_pvalue(reaction_scores, group_labels):
    """Test that results are sorted by p-value ascending."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Check sorting: each p-value should be <= the next one
    pvalues = result["pvalue"].values
    assert all(pvalues[i] <= pvalues[i + 1] for i in range(len(pvalues) - 1))


## TEST C: All three statistical methods work


@pytest.mark.parametrize("method", ["wilcoxon", "ttest", "mannwhitneyu"])
def test_compare_groups_all_methods(reaction_scores, group_labels, method):
    """Test that all statistical methods work without errors."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B", method=method)

    # Should return results for all reactions
    assert len(result) == len(reaction_scores)
    # No NaN p-values
    assert not result["pvalue"].isna().any()


## TEST D: Invalid method raises error


def test_compare_groups_invalid_method(reaction_scores, group_labels):
    """Test that invalid method raises AssertionError."""
    da = DifferentialAnalysis(reaction_scores, group_labels)

    with pytest.raises(ValueError):
        da.compare_groups("A", "B", method="invalid_method")


## TEST E: Log2FC direction is correct


def test_compare_groups_log2fc_direction(reaction_scores, group_labels):
    """Test that log2fc sign correctly reflects direction of change.

    log2fc = log2(group2_mean / group1_mean)
    - Positive log2fc means group2 > group1
    - Negative log2fc means group2 < group1
    """
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")  # group1=A, group2=B

    # R1 was made higher in group B (cells 5-9), so log2fc should be positive
    r1_row = result[result["reaction"] == "R1"].iloc[0]
    assert r1_row["log2fc"] > 0, "R1 should have positive log2fc (higher in B)"

    # R2 was made higher in group A (cells 0-4), so log2fc should be negative
    r2_row = result[result["reaction"] == "R2"].iloc[0]
    assert r2_row["log2fc"] < 0, "R2 should have negative log2fc (higher in A)"


## TEST F: P-value adjustments are valid


def test_compare_groups_padj_bounds(reaction_scores, group_labels):
    """Test that all p-values are in valid range [0, 1]."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Raw p-values in [0, 1]
    assert (result["pvalue"] >= 0).all()
    assert (result["pvalue"] <= 1).all()

    # BH-adjusted p-values in [0, 1]
    assert (result["padj_bh"] >= 0).all()
    assert (result["padj_bh"] <= 1).all()

    # Bonferroni-adjusted p-values in [0, 1]
    assert (result["padj_bonf"] >= 0).all()
    assert (result["padj_bonf"] <= 1).all()


def test_compare_groups_bonferroni_more_conservative(reaction_scores, group_labels):
    """Test that Bonferroni correction is more conservative than BH.

    Bonferroni multiplies p-values by n (number of tests).
    BH controls FDR and is less strict.
    Therefore: padj_bonf >= padj_bh (with floating point tolerance).
    """
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Bonferroni should always be >= BH (allow tiny floating point error)
    assert (result["padj_bonf"] >= result["padj_bh"] - 1e-10).all()


## TEST G: Means are calculated correctly


def test_compare_groups_means_correct(reaction_scores, group_labels):
    """Test that group means match manual calculation."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Manually calculate expected means for R1
    cells_a = [f"cell_{i}" for i in range(5)]  # group A
    cells_b = [f"cell_{i}" for i in range(5, 10)]  # group B

    expected_mean_a = reaction_scores.loc["R1", cells_a].mean()
    expected_mean_b = reaction_scores.loc["R1", cells_b].mean()

    r1_row = result[result["reaction"] == "R1"].iloc[0]

    np.testing.assert_almost_equal(r1_row["group1_mean"], expected_mean_a)
    np.testing.assert_almost_equal(r1_row["group2_mean"], expected_mean_b)


## TEST H: Partial cell overlap between scores and labels


def test_compare_groups_partial_cell_overlap():
    """Test behavior when cells don't fully overlap between scores and labels.

    Real data may have cells in expression data that aren't in metadata,
    or vice versa. The function should use only the intersection.
    """
    np.random.seed(42)

    # Reaction scores have cells 0-7
    scores = pd.DataFrame(
        np.random.rand(3, 8),
        index=["R1", "R2", "R3"],
        columns=[f"cell_{i}" for i in range(8)],
    )

    # Labels have cells 2-9 (overlap is cells 2-7, i.e., 6 cells)
    labels = pd.Series(["A"] * 4 + ["B"] * 4, index=[f"cell_{i}" for i in range(2, 10)])

    da = DifferentialAnalysis(scores, labels)
    result = da.compare_groups("A", "B")

    # Should still produce results for all 3 reactions
    assert len(result) == 3
    # No NaN values
    assert not result["pvalue"].isna().any()


## TEST I: Significant reactions are detected


def test_compare_groups_detects_significant_difference(reaction_scores, group_labels):
    """Test that reactions with clear differences have low p-values.

    R1 and R2 were artificially made different between groups.
    They should have much lower p-values than the random R3, R4, R5.
    """
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compare_groups("A", "B")

    # Get p-values for manipulated reactions
    r1_pval = result[result["reaction"] == "R1"]["pvalue"].iloc[0]
    r2_pval = result[result["reaction"] == "R2"]["pvalue"].iloc[0]

    # These should be significant (p < 0.05)
    assert r1_pval < 0.05, f"R1 should be significant, got p={r1_pval}"
    assert r2_pval < 0.05, f"R2 should be significant, got p={r2_pval}"


# =============================================================================
# TESTS FOR rank_reactions()
# =============================================================================


## TEST J: rank_reactions returns correct structure


def test_rank_reactions_returns_dataframe(reaction_scores, group_labels):
    """Test that rank_reactions returns a DataFrame with expected columns."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.rank_reactions("A", n_top=5)

    assert isinstance(result, pd.DataFrame)

    expected_cols = [
        "reaction",
        "mean_score",
        "std_score",
        "min_score",
        "max_score",
        "n_cells",
        "rank",
    ]
    for col in expected_cols:
        assert col in result.columns


## TEST K: rank_reactions respects n_top parameter


def test_rank_reactions_n_top(reaction_scores, group_labels):
    """Test that rank_reactions returns at most n_top reactions."""
    da = DifferentialAnalysis(reaction_scores, group_labels)

    result_3 = da.rank_reactions("A", n_top=3)
    result_10 = da.rank_reactions("A", n_top=10)

    assert len(result_3) == 3
    # Only 5 reactions exist, so n_top=10 should return 5
    assert len(result_10) == 5


## TEST L: rank_reactions sorted by mean_score ascending


def test_rank_reactions_sorted(reaction_scores, group_labels):
    """Test that results are sorted by mean_score ascending (lowest = most active)."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.rank_reactions("A")

    mean_scores = result["mean_score"].values
    is_sorted = all(
        mean_scores[i] <= mean_scores[i + 1] for i in range(len(mean_scores) - 1)
    )
    assert is_sorted


## TEST M: rank_reactions rank column is correct


def test_rank_reactions_rank_column(reaction_scores, group_labels):
    """Test that rank column starts at 1 and increments."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.rank_reactions("A", n_top=5)

    expected_ranks = list(range(1, len(result) + 1))
    assert list(result["rank"]) == expected_ranks


## TEST N: rank_reactions n_cells is correct


def test_rank_reactions_n_cells(reaction_scores, group_labels):
    """Test that n_cells correctly counts cells in the group."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.rank_reactions("A")

    # Group A has 5 cells (cell_0 to cell_4)
    assert (result["n_cells"] == 5).all()


# =============================================================================
# TESTS FOR compute_effect_size()
# =============================================================================


## TEST O: compute_effect_size returns correct structure


def test_compute_effect_size_returns_dataframe(reaction_scores, group_labels):
    """Test that compute_effect_size returns a DataFrame with expected columns."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    assert isinstance(result, pd.DataFrame)

    expected_cols = [
        "reaction",
        "cohens_d",
        "abs_cohens_d",
        "interpretation",
        "group1_mean",
        "group2_mean",
        "pooled_std",
    ]
    for col in expected_cols:
        assert col in result.columns


## TEST P: compute_effect_size sorted by absolute effect size


def test_compute_effect_size_sorted(reaction_scores, group_labels):
    """Test results are sorted by abs_cohens_d descending (largest first)."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    abs_d = result["abs_cohens_d"].values
    assert all(abs_d[i] >= abs_d[i + 1] for i in range(len(abs_d) - 1))


## TEST Q: compute_effect_size interpretation is valid


def test_compute_effect_size_interpretation(reaction_scores, group_labels):
    """Test that interpretation values are valid categories."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    valid_interpretations = {"negligible", "small", "medium", "large"}
    assert set(result["interpretation"].unique()).issubset(valid_interpretations)


## TEST R: compute_effect_size detects large effects


def test_compute_effect_size_detects_large_effects(reaction_scores, group_labels):
    """Test that reactions with large differences have large effect sizes.

    R1 and R2 were artificially made different between groups.
    They should have large effect sizes.
    """
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    r1_row = result[result["reaction"] == "R1"].iloc[0]
    r2_row = result[result["reaction"] == "R2"].iloc[0]

    # These should have large effect sizes
    r1_interp = r1_row["interpretation"]
    r2_interp = r2_row["interpretation"]
    assert r1_interp == "large", f"R1 should have large effect, got {r1_interp}"
    assert r2_interp == "large", f"R2 should have large effect, got {r2_interp}"


## TEST S: compute_effect_size direction matches Cohen's d sign


def test_compute_effect_size_direction(reaction_scores, group_labels):
    """Test that Cohen's d sign reflects which group has higher values.

    Cohen's d = (mean1 - mean2) / pooled_std
    Positive d means group1 > group2
    Negative d means group1 < group2
    """
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    # R1 is higher in B, so group1(A) < group2(B), d should be negative
    r1_row = result[result["reaction"] == "R1"].iloc[0]
    assert r1_row["cohens_d"] < 0, "R1 should have negative Cohen's d (higher in B)"

    # R2 is higher in A, so group1(A) > group2(B), d should be positive
    r2_row = result[result["reaction"] == "R2"].iloc[0]
    assert r2_row["cohens_d"] > 0, "R2 should have positive Cohen's d (higher in A)"


## TEST T: compute_effect_size abs_cohens_d matches cohens_d


def test_compute_effect_size_abs_matches(reaction_scores, group_labels):
    """Test that abs_cohens_d is the absolute value of cohens_d."""
    da = DifferentialAnalysis(reaction_scores, group_labels)
    result = da.compute_effect_size("A", "B")

    for _, row in result.iterrows():
        np.testing.assert_almost_equal(row["abs_cohens_d"], abs(row["cohens_d"]))


# =============================================================================
# FIXTURES FOR MULTI-GROUP TESTS
# =============================================================================


@pytest.fixture
def multi_group_reaction_scores():
    """Create mock reaction scores with 3 groups."""
    np.random.seed(42)

    # 5 reactions, 15 cells (5 per group)
    reactions = ["R1", "R2", "R3", "R4", "R5"]
    cells = [f"cell_{i}" for i in range(15)]

    # Create base data
    data = np.random.rand(5, 15)

    # Make R1 clearly different across groups
    # Group A (cells 0-4): low values
    data[0, 0:5] = np.random.rand(5) * 0.5
    # Group B (cells 5-9): medium values
    data[0, 5:10] = np.random.rand(5) * 0.5 + 1.5
    # Group C (cells 10-14): high values
    data[0, 10:15] = np.random.rand(5) * 0.5 + 3.0

    # Make R2 different between A and C but not B
    data[1, 0:5] = np.random.rand(5) * 0.5
    data[1, 5:10] = np.random.rand(5) * 0.5 + 0.5
    data[1, 10:15] = np.random.rand(5) * 0.5 + 2.0

    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def multi_group_labels():
    """Create group labels for 3 groups."""
    cells = [f"cell_{i}" for i in range(15)]
    groups = ["A"] * 5 + ["B"] * 5 + ["C"] * 5
    return pd.Series(groups, index=cells)


# =============================================================================
# TESTS FOR compare_multiple_groups()
# =============================================================================


def test_compare_multiple_groups_returns_dataframe(
    multi_group_reaction_scores, multi_group_labels
):
    """Test compare_multiple_groups returns DataFrame with expected columns."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups()

    assert isinstance(result, pd.DataFrame)

    # Check expected columns
    expected_cols = [
        "reaction",
        "statistic",
        "pvalue",
        "n_groups",
        "padj_bh",
        "padj_bonf",
    ]
    for col in expected_cols:
        assert col in result.columns

    # Check group-specific columns exist
    for group in ["A", "B", "C"]:
        assert f"{group}_mean" in result.columns
        assert f"{group}_std" in result.columns
        assert f"{group}_n" in result.columns


def test_compare_multiple_groups_kruskal(
    multi_group_reaction_scores, multi_group_labels
):
    """Test Kruskal-Wallis method."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups(method="kruskal")

    assert len(result) == 5  # 5 reactions
    assert not result["pvalue"].isna().any()


def test_compare_multiple_groups_anova(multi_group_reaction_scores, multi_group_labels):
    """Test ANOVA method."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups(method="anova")

    assert len(result) == 5
    assert not result["pvalue"].isna().any()


def test_compare_multiple_groups_detects_difference(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that reactions with clear differences have low p-values."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups(method="kruskal")

    # R1 was made to differ clearly across groups
    r1_pval = result[result["reaction"] == "R1"]["pvalue"].iloc[0]
    assert r1_pval < 0.05, f"R1 should be significant, got p={r1_pval}"


def test_compare_multiple_groups_sorted_by_pvalue(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that results are sorted by p-value."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups()

    pvalues = result["pvalue"].values
    assert all(pvalues[i] <= pvalues[i + 1] for i in range(len(pvalues) - 1))


def test_compare_multiple_groups_subset_of_groups(
    multi_group_reaction_scores, multi_group_labels
):
    """Test comparison with subset of groups."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.compare_multiple_groups(groups=["A", "C"])

    assert result["n_groups"].iloc[0] == 2
    assert "A_mean" in result.columns
    assert "C_mean" in result.columns
    # B should not be included
    assert "B_mean" not in result.columns


def test_compare_multiple_groups_invalid_method(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that invalid method raises error."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    with pytest.raises(ValueError):
        da.compare_multiple_groups(method="invalid")


def test_compare_multiple_groups_too_few_groups(
    multi_group_reaction_scores, multi_group_labels
):
    """Test error when only one group specified."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    with pytest.raises(ValueError):
        da.compare_multiple_groups(groups=["A"])


# =============================================================================
# TESTS FOR posthoc_tests()
# =============================================================================


def test_posthoc_dunn_returns_dataframe(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that Dunn's post-hoc test returns DataFrame."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.posthoc_tests("R1", method="dunn")

    assert isinstance(result, pd.DataFrame)

    # Should have n*(n-1)/2 = 3 pairwise comparisons for 3 groups
    assert len(result) == 3

    # Check expected columns
    expected_cols = ["group1", "group2", "pvalue", "padj", "significant"]
    for col in expected_cols:
        assert col in result.columns


def test_posthoc_tukey_returns_dataframe(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that Tukey's HSD post-hoc test returns DataFrame."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.posthoc_tests("R1", method="tukey")

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3  # 3 pairwise comparisons


def test_posthoc_conover_returns_dataframe(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that Conover's post-hoc test returns DataFrame."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.posthoc_tests("R1", method="conover")

    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3


def test_posthoc_detects_significant_pairs(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that post-hoc test identifies significant pairwise differences."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    # R1 was made to differ between A and C
    result = da.posthoc_tests("R1", method="dunn")

    # Find A vs C comparison
    a_vs_c = result[(result["group1"] == "A") & (result["group2"] == "C")]
    if len(a_vs_c) == 0:
        a_vs_c = result[(result["group1"] == "C") & (result["group2"] == "A")]

    assert len(a_vs_c) == 1
    assert a_vs_c["padj"].iloc[0] < 0.05, "A vs C should be significant"


def test_posthoc_invalid_reaction(multi_group_reaction_scores, multi_group_labels):
    """Test error when reaction not found."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    with pytest.raises(ValueError):
        da.posthoc_tests("nonexistent_reaction")


def test_posthoc_invalid_method(multi_group_reaction_scores, multi_group_labels):
    """Test error when invalid method specified."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    with pytest.raises(ValueError):
        da.posthoc_tests("R1", method="invalid")


def test_posthoc_pvalues_in_range(multi_group_reaction_scores, multi_group_labels):
    """Test that p-values are in valid range [0, 1]."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    for method in ["dunn", "tukey", "conover"]:
        result = da.posthoc_tests("R1", method=method)
        assert (result["pvalue"] >= 0).all()
        assert (result["pvalue"] <= 1).all()


# =============================================================================
# TESTS FOR all_pairwise_comparisons()
# =============================================================================


def test_all_pairwise_comparisons_returns_dataframe(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that all_pairwise_comparisons returns DataFrame."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.all_pairwise_comparisons()

    assert isinstance(result, pd.DataFrame)

    # Check expected columns
    expected_cols = ["comparison", "reaction", "group1", "group2", "log2fc", "pvalue"]
    for col in expected_cols:
        assert col in result.columns


def test_all_pairwise_comparisons_correct_count(
    multi_group_reaction_scores, multi_group_labels
):
    """Test correct number of comparisons."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.all_pairwise_comparisons()

    # 5 reactions x 3 pairwise comparisons = 15 rows
    assert len(result) == 15


def test_all_pairwise_comparisons_subset_groups(
    multi_group_reaction_scores, multi_group_labels
):
    """Test with subset of groups."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)
    result = da.all_pairwise_comparisons(groups=["A", "B"])

    # 5 reactions x 1 comparison = 5 rows
    assert len(result) == 5
    assert result["comparison"].unique()[0] == "A_vs_B"


def test_all_pairwise_comparisons_methods(
    multi_group_reaction_scores, multi_group_labels
):
    """Test that different methods work."""
    da = DifferentialAnalysis(multi_group_reaction_scores, multi_group_labels)

    for method in ["wilcoxon", "ttest", "mannwhitneyu"]:
        result = da.all_pairwise_comparisons(method=method)
        assert len(result) == 15
        assert not result["pvalue"].isna().any()


# =============================================================================
# EDGE CASE TESTS
# =============================================================================


class TestDifferentialEdgeCases:
    """Edge case tests for differential analysis."""

    def test_compare_groups_single_reaction(self):
        """Test comparison with a single reaction."""
        np.random.seed(42)

        scores = pd.DataFrame(
            np.random.rand(1, 10),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(10)],
        )
        labels = pd.Series(
            ["A"] * 5 + ["B"] * 5, index=[f"cell_{i}" for i in range(10)]
        )

        da = DifferentialAnalysis(scores, labels)
        result = da.compare_groups("A", "B")

        assert len(result) == 1
        assert result.iloc[0]["reaction"] == "R1"

    def test_compare_groups_single_cell_per_group(self):
        """Test comparison with minimal cells per group."""
        np.random.seed(42)

        scores = pd.DataFrame(
            np.array([[1.0, 5.0], [2.0, 6.0]]),
            index=["R1", "R2"],
            columns=["cell_0", "cell_1"],
        )
        labels = pd.Series(["A", "B"], index=["cell_0", "cell_1"])

        da = DifferentialAnalysis(scores, labels)
        # Should handle but may produce NaN p-values for some methods
        result = da.compare_groups("A", "B", method="ttest")

        assert len(result) == 2

    def test_compare_groups_constant_values_one_group(self):
        """Test when one group has constant values (zero variance)."""
        scores = pd.DataFrame(
            np.array([[1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(10)],
        )
        labels = pd.Series(
            ["A"] * 5 + ["B"] * 5, index=[f"cell_{i}" for i in range(10)]
        )

        da = DifferentialAnalysis(scores, labels)
        # Should handle zero variance gracefully
        result = da.compare_groups("A", "B")

        assert len(result) == 1

    def test_compare_groups_identical_groups(self):
        """Test when both groups have identical values."""
        np.random.seed(42)
        values = np.random.rand(3, 10)

        scores = pd.DataFrame(
            values,
            index=["R1", "R2", "R3"],
            columns=[f"cell_{i}" for i in range(10)],
        )
        # Same data in both groups (no actual difference)
        labels = pd.Series(
            ["A"] * 5 + ["B"] * 5, index=[f"cell_{i}" for i in range(10)]
        )

        # Make values identical between groups for R1
        scores.loc["R1", "cell_0":"cell_4"] = 1.0
        scores.loc["R1", "cell_5":"cell_9"] = 1.0

        da = DifferentialAnalysis(scores, labels)
        result = da.compare_groups("A", "B")

        # R1 should have log2fc close to 0
        r1_log2fc = result[result["reaction"] == "R1"]["log2fc"].iloc[0]
        assert abs(r1_log2fc) < 1e-10

    def test_compare_groups_with_nan_values(self):
        """Test handling of NaN values in reaction scores."""
        scores = pd.DataFrame(
            np.array([[1.0, 2.0, np.nan, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]]),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(10)],
        )
        labels = pd.Series(
            ["A"] * 5 + ["B"] * 5, index=[f"cell_{i}" for i in range(10)]
        )

        da = DifferentialAnalysis(scores, labels)
        result = da.compare_groups("A", "B")

        # Should produce valid result despite NaN
        assert len(result) == 1

    def test_compare_groups_nonexistent_group(self):
        """Test behavior when comparing with non-existent group.

        The code may return NaN p-values or empty results when group doesn't exist.
        """
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(3, 6),
            index=["R1", "R2", "R3"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)

        # Code handles missing group gracefully (returns NaN p-values)
        result = da.compare_groups("A", "NONEXISTENT")
        # Result will have NaN p-values due to empty group
        assert len(result) == 3
        assert result["pvalue"].isna().all()

    def test_compute_effect_size_zero_pooled_std(self):
        """Test effect size when pooled std is zero (constant values)."""
        scores = pd.DataFrame(
            np.array([[1.0, 1.0, 1.0, 2.0, 2.0, 2.0]]),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)
        result = da.compute_effect_size("A", "B")

        # Should handle division by zero gracefully
        assert len(result) == 1
        # Effect size may be inf or handled differently

    def test_rank_reactions_nonexistent_group(self):
        """Test rank_reactions with non-existent group.

        Code handles missing group gracefully with NaN statistics.
        """
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(3, 6),
            index=["R1", "R2", "R3"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)

        result = da.rank_reactions("NONEXISTENT")
        # Result will have NaN means due to empty group
        assert len(result) == 3
        assert result["mean_score"].isna().all()

    def test_compare_multiple_groups_two_groups_only(self):
        """Test multi-group comparison with exactly two groups."""
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(3, 6),
            index=["R1", "R2", "R3"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)
        result = da.compare_multiple_groups()

        assert len(result) == 3
        assert result["n_groups"].iloc[0] == 2

    def test_posthoc_two_groups_only(self):
        """Test post-hoc test with only two groups (degenerates to pairwise)."""
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(3, 6),
            index=["R1", "R2", "R3"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)
        result = da.posthoc_tests("R1", method="dunn")

        # Should have exactly 1 pairwise comparison
        assert len(result) == 1

    def test_all_pairwise_comparisons_single_group_cells(self):
        """Test pairwise comparisons when one group has very few cells."""
        np.random.seed(42)
        scores = pd.DataFrame(
            np.random.rand(2, 7),
            index=["R1", "R2"],
            columns=[f"cell_{i}" for i in range(7)],
        )
        # Group A has 1 cell, B has 3, C has 3
        labels = pd.Series(
            ["A"] + ["B"] * 3 + ["C"] * 3,
            index=[f"cell_{i}" for i in range(7)],
        )

        da = DifferentialAnalysis(scores, labels)
        # Should still work, though statistical power is limited
        result = da.all_pairwise_comparisons()

        # 2 reactions x 3 comparisons = 6 rows
        assert len(result) == 6

    def test_differential_empty_intersection(self):
        """Test when there's no cell overlap between scores and labels.

        Code handles this by returning NaN p-values due to empty groups.
        """
        scores = pd.DataFrame(
            np.random.rand(3, 5),
            index=["R1", "R2", "R3"],
            columns=["cell_0", "cell_1", "cell_2", "cell_3", "cell_4"],
        )
        labels = pd.Series(
            ["A", "A", "B", "B", "B"],
            index=["cell_100", "cell_101", "cell_102", "cell_103", "cell_104"],
        )

        da = DifferentialAnalysis(scores, labels)

        # Code handles this gracefully - produces NaN p-values
        result = da.compare_groups("A", "B")
        assert len(result) == 3
        assert result["pvalue"].isna().all()

    def test_compare_groups_large_log2fc_handling(self):
        """Test handling of very large fold changes (near-zero values)."""
        scores = pd.DataFrame(
            np.array([[1e-10, 1e-10, 1e-10, 1.0, 1.0, 1.0]]),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(6)],
        )
        labels = pd.Series(["A"] * 3 + ["B"] * 3, index=[f"cell_{i}" for i in range(6)])

        da = DifferentialAnalysis(scores, labels)
        result = da.compare_groups("A", "B")

        # log2fc should be finite (not inf)
        assert np.isfinite(result["log2fc"].iloc[0])


# =============================================================================
# FIXTURES FOR PSEUDO-BULK TESTS
# =============================================================================
# Design: 3 reactions, 4 samples (2 per group), 5 cells each = 20 cells total.
# R1 is spiked higher in group B to give a detectable signal.


@pytest.fixture
def pb_scores():
    """Reaction scores (reactions x cells): 3 reactions, 20 cells."""
    np.random.seed(0)
    reactions = ["R1", "R2", "R3"]
    cells = [f"cell_{i}" for i in range(20)]
    data = np.random.rand(3, 20)
    # R1 clearly higher in group B (cells 10-19, samples s3/s4)
    data[0, 10:] += 3.0
    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def pb_groups():
    """Group labels: cells 0-9 → group A, cells 10-19 → group B."""
    cells = [f"cell_{i}" for i in range(20)]
    return pd.Series(["A"] * 10 + ["B"] * 10, index=cells)


@pytest.fixture
def pb_samples():
    """Sample labels: s1/s2 in group A (5 cells each), s3/s4 in group B."""
    cells = [f"cell_{i}" for i in range(20)]
    return pd.Series(
        ["s1"] * 5 + ["s2"] * 5 + ["s3"] * 5 + ["s4"] * 5, index=cells
    )


# =============================================================================
# TESTS FOR PseudoBulkAnalysis.aggregate()
# =============================================================================


def test_aggregate_returns_dataframe(pb_scores, pb_groups, pb_samples):
    """aggregate() should return a DataFrame of shape (reactions x samples)."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.aggregate()

    assert isinstance(result, pd.DataFrame)
    assert result.shape == (3, 4)  # 3 reactions, 4 samples
    assert set(result.columns) == {"s1", "s2", "s3", "s4"}
    assert list(result.index) == ["R1", "R2", "R3"]


def test_aggregate_mean_correct(pb_scores, pb_groups, pb_samples):
    """Mean aggregation should match manual per-sample average."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.aggregate(method="mean")

    # s1 contains cells 0-4
    s1_cells = [f"cell_{i}" for i in range(5)]
    expected_s1_r1 = pb_scores.loc["R1", s1_cells].mean()
    np.testing.assert_almost_equal(result.loc["R1", "s1"], expected_s1_r1)


def test_aggregate_sum_correct(pb_scores, pb_groups, pb_samples):
    """Sum aggregation should match manual per-sample sum."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.aggregate(method="sum")

    s3_cells = [f"cell_{i}" for i in range(10, 15)]
    expected_s3_r2 = pb_scores.loc["R2", s3_cells].sum()
    np.testing.assert_almost_equal(result.loc["R2", "s3"], expected_s3_r2)


def test_aggregate_sample_in_multiple_groups_raises(pb_scores, pb_groups, pb_samples):
    """aggregate() must raise if one sample has cells in two different groups."""
    # Corrupt one cell's group so s1 spans both A and B
    bad_groups = pb_groups.copy()
    bad_groups["cell_0"] = "B"

    pb = PseudoBulkAnalysis(pb_scores, bad_groups, pb_samples)
    with pytest.raises(ValueError, match="multiple groups"):
        pb.aggregate()


def test_aggregate_partial_cell_overlap(pb_scores, pb_groups, pb_samples):
    """aggregate() uses only cells present in all three inputs."""
    # Drop a few cells from groups so they're missing
    trimmed_groups = pb_groups.drop(["cell_0", "cell_1"])
    pb = PseudoBulkAnalysis(pb_scores, trimmed_groups, pb_samples)
    result = pb.aggregate()

    # Should still return 4 samples (s1 just has 3 cells instead of 5)
    assert result.shape[1] == 4


# =============================================================================
# TESTS FOR PseudoBulkAnalysis.compare_groups()
# =============================================================================


def test_compare_groups_pb_returns_dataframe(pb_scores, pb_groups, pb_samples):
    """compare_groups() should return a DataFrame with expected columns."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    assert isinstance(result, pd.DataFrame)
    expected_cols = [
        "reaction",
        "group1_mean",
        "group2_mean",
        "log2fc",
        "statistic",
        "pvalue",
        "padj_bh",
        "padj_bonf",
        "group1_n_samples",
        "group2_n_samples",
    ]
    for col in expected_cols:
        assert col in result.columns


def test_compare_groups_pb_n_samples_correct(pb_scores, pb_groups, pb_samples):
    """n_samples columns should reflect the number of replicates per group."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    assert (result["group1_n_samples"] == 2).all()
    assert (result["group2_n_samples"] == 2).all()


def test_compare_groups_pb_log2fc_direction(pb_scores, pb_groups, pb_samples):
    """R1 was spiked in group B → log2fc should be positive (B > A)."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    r1_row = result[result["reaction"] == "R1"].iloc[0]
    assert r1_row["log2fc"] > 0, f"Expected positive log2fc for R1, got {r1_row['log2fc']}"


def test_compare_groups_pb_detects_signal(pb_scores, pb_groups, pb_samples):
    """R1 has a large spike so its p-value should be the lowest."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    r1_pval = result[result["reaction"] == "R1"]["pvalue"].iloc[0]
    other_pvals = result[result["reaction"] != "R1"]["pvalue"]
    assert r1_pval < other_pvals.min(), "R1 should have the smallest p-value"


def test_compare_groups_pb_sorted_by_pvalue(pb_scores, pb_groups, pb_samples):
    """Results should be sorted by p-value ascending."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    pvalues = result["pvalue"].values
    assert all(pvalues[i] <= pvalues[i + 1] for i in range(len(pvalues) - 1))


def test_compare_groups_pb_padj_bounds(pb_scores, pb_groups, pb_samples):
    """All p-values and adjusted p-values must be in [0, 1]."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B")

    for col in ("pvalue", "padj_bh", "padj_bonf"):
        assert (result[col] >= 0).all()
        assert (result[col] <= 1).all()


@pytest.mark.parametrize("method", ["ttest", "wilcoxon", "mannwhitneyu"])
def test_compare_groups_pb_all_methods(pb_scores, pb_groups, pb_samples, method):
    """All three statistical methods should run without errors."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    result = pb.compare_groups("A", "B", method=method)

    assert len(result) == 3
    assert not result["pvalue"].isna().any()


def test_compare_groups_pb_invalid_method(pb_scores, pb_groups, pb_samples):
    """An unknown method name should raise ValueError."""
    pb = PseudoBulkAnalysis(pb_scores, pb_groups, pb_samples)
    with pytest.raises(ValueError, match="Invalid method"):
        pb.compare_groups("A", "B", method="bad_method")


def test_compare_groups_pb_too_few_replicates(pb_scores, pb_groups, pb_samples):
    """Groups with only one replicate must raise ValueError."""
    # Remap all group-B cells to a single sample
    single_rep_samples = pb_samples.copy()
    single_rep_samples[single_rep_samples == "s4"] = "s3"

    pb = PseudoBulkAnalysis(pb_scores, pb_groups, single_rep_samples)
    with pytest.raises(ValueError, match="at least 2 biological replicates"):
        pb.compare_groups("A", "B")


# =============================================================================
# TESTS FOR PseudoBulkAnalysis.compare_multiple_groups()
# =============================================================================


@pytest.fixture
def pb_three_group_scores():
    """Reaction scores for a 3-group experiment (6 samples, 5 cells each)."""
    np.random.seed(1)
    reactions = ["R1", "R2", "R3"]
    cells = [f"cell_{i}" for i in range(30)]
    data = np.random.rand(3, 30)
    # R1 increases monotonically across groups A < B < C
    data[0, 0:10] = np.random.rand(10) * 0.5          # group A: low
    data[0, 10:20] = np.random.rand(10) * 0.5 + 1.5   # group B: mid
    data[0, 20:30] = np.random.rand(10) * 0.5 + 3.0   # group C: high
    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def pb_three_group_groups():
    cells = [f"cell_{i}" for i in range(30)]
    return pd.Series(["A"] * 10 + ["B"] * 10 + ["C"] * 10, index=cells)


@pytest.fixture
def pb_three_group_samples():
    cells = [f"cell_{i}" for i in range(30)]
    labels = (
        ["s1"] * 5 + ["s2"] * 5     # group A
        + ["s3"] * 5 + ["s4"] * 5   # group B
        + ["s5"] * 5 + ["s6"] * 5   # group C
    )
    return pd.Series(labels, index=cells)


def test_compare_multiple_groups_pb_returns_dataframe(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
):
    """compare_multiple_groups() should return a DataFrame with expected columns."""
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    result = pb.compare_multiple_groups()

    assert isinstance(result, pd.DataFrame)
    for col in ("reaction", "statistic", "pvalue", "n_groups", "padj_bh", "padj_bonf"):
        assert col in result.columns
    for g in ("A", "B", "C"):
        assert f"{g}_mean" in result.columns
        assert f"{g}_n_samples" in result.columns


def test_compare_multiple_groups_pb_detects_signal(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
):
    """R1 has a clear gradient across groups and should be the most significant."""
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    result = pb.compare_multiple_groups(method="anova")

    r1_pval = result[result["reaction"] == "R1"]["pvalue"].iloc[0]
    other_pvals = result[result["reaction"] != "R1"]["pvalue"]
    assert r1_pval < other_pvals.min()


@pytest.mark.parametrize("method", ["anova", "kruskal"])
def test_compare_multiple_groups_pb_methods(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples, method
):
    """Both anova and kruskal methods should run without errors."""
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    result = pb.compare_multiple_groups(method=method)

    assert len(result) == 3
    assert not result["pvalue"].isna().any()


def test_compare_multiple_groups_pb_invalid_method(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
):
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    with pytest.raises(ValueError, match="Invalid method"):
        pb.compare_multiple_groups(method="bad")


def test_compare_multiple_groups_pb_subset(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
):
    """Passing a subset of groups should restrict the comparison."""
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    result = pb.compare_multiple_groups(groups=["A", "C"])

    assert result["n_groups"].iloc[0] == 2
    assert "A_mean" in result.columns
    assert "C_mean" in result.columns
    assert "B_mean" not in result.columns


def test_compare_multiple_groups_pb_too_few_groups(
    pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
):
    pb = PseudoBulkAnalysis(
        pb_three_group_scores, pb_three_group_groups, pb_three_group_samples
    )
    with pytest.raises(ValueError, match="at least 2 groups"):
        pb.compare_multiple_groups(groups=["A"])
