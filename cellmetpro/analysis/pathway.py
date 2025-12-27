"""Pathway-level aggregation and analysis.

This module aggregates reaction-level scores into pathway scores
and provides pathway enrichment analysis.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import cobra


class PathwayAnalyzer:
    """Aggregate reactions into pathway-level scores.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    model : cobra.Model
        Metabolic model for pathway annotations.

    Attributes
    ----------
    pathway_scores : pd.DataFrame
        Pathway-level aggregated scores.
    pathway_mapping : dict
        Mapping of reactions to pathways.
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        model: cobra.Model,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.model = model
        self.pathway_scores: pd.DataFrame | None = None
        self.pathway_mapping: dict | None = None

    def get_pathway_mapping(self) -> dict[str, list[str]]:
        """Extract pathway-to-reaction mapping from model.

        Returns
        -------
        dict
            Mapping of pathway names to reaction IDs.
        """
        raise NotImplementedError

    def aggregate_to_pathways(
        self,
        method: str = "mean",
    ) -> pd.DataFrame:
        """Aggregate reaction scores to pathway level.

        Parameters
        ----------
        method : str
            Aggregation method: 'mean', 'median', 'sum', or 'max'.

        Returns
        -------
        pd.DataFrame
            Pathway scores (pathways x cells).
        """
        raise NotImplementedError

    def pathway_enrichment(
        self,
        differential_results: pd.DataFrame,
        method: str = "gsea",
    ) -> pd.DataFrame:
        """Perform pathway enrichment analysis.

        Parameters
        ----------
        differential_results : pd.DataFrame
            Results from differential analysis.
        method : str
            Enrichment method: 'gsea', 'ora', or 'fgsea'.

        Returns
        -------
        pd.DataFrame
            Enrichment results with pathway names and statistics.
        """
        raise NotImplementedError
