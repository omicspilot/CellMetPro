"""Statistical comparison of metabolic states.

This module provides methods for identifying differentially active
reactions and pathways between cell populations or conditions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd

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
            Results with columns: reaction, log2fc, pvalue, padj.
        """
        raise NotImplementedError

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
