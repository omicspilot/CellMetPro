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
            Top reactions ranked by mean activity.
        """
        raise NotImplementedError

    def compute_effect_size(
        self,
        group1: str,
        group2: str,
    ) -> pd.Series:
        """Compute effect size (Cohen's d) between groups.

        Parameters
        ----------
        group1 : str
            Name of first group.
        group2 : str
            Name of second group.

        Returns
        -------
        pd.Series
            Effect sizes for each reaction.
        """
        raise NotImplementedError
