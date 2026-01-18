"""Statistical comparison of metabolic states.

This module provides methods for identifying differentially active
reactions and pathways between cell populations or conditions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


if TYPE_CHECKING:
    import anndata as ad


class DifferentialAnalysis:
    """Compare metabolic activity between groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    groups : pd.Series
        Group labels for each cell.

    Attributes
    ----------
    results : pd.DataFrame
        Differential analysis results with statistics.
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        groups: pd.Series,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.groups = groups
        self.results: pd.DataFrame | None = None

    def compare_groups(
        self,
        group1: str,
        group2: str,
        method: Literal["wilcoxon", "ttest", "mannwhitneyu"] = "wilcoxon",
    ) -> pd.DataFrame:
        """Compare two groups statistically.

        Parameters
        ----------
        group1 : str
            Name of first group.
        group2 : str
            Name of second group.
        method : str
            Statistical test to use.

        Returns
        -------
        pd.DataFrame
            Results with columns: reaction, group1_mean, group2_mean, log2fc, pvalue, padj_bh, padj_bonf.
        """
        
        assert method in {"wilcoxon", "ttest", "mannwhitneyu"}, "Invalid method"
    
        # Ensure columns in reaction_scores match index in groups
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        # Filter cells by group
        cells1 = group_labels[group_labels == group1].index
        cells2 = group_labels[group_labels == group2].index

        results = []
        for reaction_id in scores.index:
            # reaction scores
            scores1 = scores.loc[reaction_id, cells1]
            scores2 = scores.loc[reaction_id, cells2]
            
            # means
            mean1 = np.mean(scores1)
            mean2 = np.mean(scores2)
            
            # log2 for change
            epsilon = 1e-9
            log2fc = np.log2((mean2 + epsilon) / (mean1 + epsilon))
            
            # run statistical test
            if method == "wilcoxon":
                stat, pval = stats.ranksums(scores1, scores2)
            elif method == "ttest":
                stat, pval = stats.ttest_ind(scores1, scores2)
            elif method == "mannwhitneyu":
                stat, pval = stats.mannwhitneyu(scores1, scores2, alternative='two-sided')
        
            results.append({
                "reaction": reaction_id,
                "group1_mean": mean1,
                "group2_mean": mean2,
                "log2fc": log2fc,
                "statistic": stat,
                "pvalue": pval,
            })
            
        # FDR correction
        results_df = pd.DataFrame(results)
        _, padj_bh, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")      # Less strict
        _, padj_bonf, _, _ = multipletests(results_df["pvalue"], method="bonferroni") # More strict

        results_df["padj_bh"] = padj_bh
        results_df["padj_bonf"] = padj_bonf
        
        # sort results based on p-val
        return results_df.sort_values("pvalue")

    def rank_reactions(
        self,
        group: str,
        n_top: int = 50,
    ) -> pd.DataFrame:
        """Rank reactions by activity in a group.

        Parameters
        ----------
        group : str
            Group to analyze.
        n_top : int
            Number of top reactions to return.

        Returns
        -------
        pd.DataFrame
            Top reactions ranked by mean activity (lowest penalty = highest activity).
        """
        # Get cells belonging to this group
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        group_cells = group_labels[group_labels == group].index

        # Compute statistics for each reaction
        results = []
        for reaction_id in scores.index:
            reaction_scores = scores.loc[reaction_id, group_cells]
            results.append({
                "reaction": reaction_id,
                "mean_score": np.mean(reaction_scores),
                "std_score": np.std(reaction_scores),
                "min_score": np.min(reaction_scores),
                "max_score": np.max(reaction_scores),
                "n_cells": len(group_cells),
            })

        # Sort by mean score (ascending = lowest penalty = highest activity)
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values("mean_score", ascending=True)

        # Add rank column
        results_df["rank"] = range(1, len(results_df) + 1)

        # Return top N
        return results_df.head(n_top).reset_index(drop=True)

    def compute_effect_size(
        self,
        group1: str,
        group2: str,
    ) -> pd.DataFrame:
        """Compute effect size (Cohen's d) between groups.

        Cohen's d measures the standardized difference between two means:
            d = (mean1 - mean2) / pooled_std

        Interpretation:
            |d| < 0.2  : negligible
            |d| 0.2-0.5: small
            |d| 0.5-0.8: medium
            |d| > 0.8  : large

        Parameters
        ----------
        group1 : str
            Name of first group.
        group2 : str
            Name of second group.

        Returns
        -------
        pd.DataFrame
            Effect sizes for each reaction with interpretation.
        """
        # Get common cells
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        # Filter cells by group
        cells1 = group_labels[group_labels == group1].index
        cells2 = group_labels[group_labels == group2].index

        n1 = len(cells1)
        n2 = len(cells2)

        results = []
        for reaction_id in scores.index:
            scores1 = scores.loc[reaction_id, cells1].values
            scores2 = scores.loc[reaction_id, cells2].values

            mean1 = np.mean(scores1)
            mean2 = np.mean(scores2)
            std1 = np.std(scores1, ddof=1)
            std2 = np.std(scores2, ddof=1)

            # Pooled standard deviation
            pooled_std = np.sqrt(
                ((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2)
            )

            # Cohen's d (avoid division by zero)
            if pooled_std > 0:
                cohens_d = (mean1 - mean2) / pooled_std
            else:
                cohens_d = 0.0

            # Interpretation
            abs_d = abs(cohens_d)
            if abs_d < 0.2:
                interpretation = "negligible"
            elif abs_d < 0.5:
                interpretation = "small"
            elif abs_d < 0.8:
                interpretation = "medium"
            else:
                interpretation = "large"

            results.append({
                "reaction": reaction_id,
                "cohens_d": cohens_d,
                "abs_cohens_d": abs_d,
                "interpretation": interpretation,
                "group1_mean": mean1,
                "group2_mean": mean2,
                "pooled_std": pooled_std,
            })

        results_df = pd.DataFrame(results)
        return results_df.sort_values("abs_cohens_d", ascending=False)
