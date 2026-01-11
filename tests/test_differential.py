"""Tests for differential analysis module."""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from cellmetpro.analysis.differential import DifferentialAnalysis


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
    expected_cols = ["reaction", "group1_mean", "group2_mean",
                     "log2fc", "statistic", "pvalue", "padj_bh", "padj_bonf"]
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

    with pytest.raises(AssertionError):
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
    cells_a = [f"cell_{i}" for i in range(5)]      # group A
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
        columns=[f"cell_{i}" for i in range(8)]
    )

    # Labels have cells 2-9 (overlap is cells 2-7, i.e., 6 cells)
    labels = pd.Series(
        ["A"] * 4 + ["B"] * 4,
        index=[f"cell_{i}" for i in range(2, 10)]
    )

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
