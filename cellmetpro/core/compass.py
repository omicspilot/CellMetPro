"""COMPASS algorithm implementation for metabolic scoring.

COMPASS (Characterizing Cell states through metabolic Profiling of the
Transcriptome) integrates scRNA-seq data with genome-scale metabolic
models to infer metabolic activity at single-cell resolution.

References:
    Wagner et al. (2021) "Metabolic modeling of single Th17 cells reveals
    regulators of autoimmunity" Cell.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad
    import cobra


class CompassScorer:
    """Compute metabolic reaction scores using the COMPASS algorithm.

    Parameters
    ----------
    model : cobra.Model
        A genome-scale metabolic model (GEM).
    gene_expression : pd.DataFrame
        Gene expression matrix (genes x cells).

    Attributes
    ----------
    reaction_scores : pd.DataFrame
        Computed reaction activity scores (reactions x cells).
    """

    def __init__(
        self,
        model: cobra.Model,
        gene_expression: pd.DataFrame,
    ) -> None:
        self.model = model
        self.gene_expression = gene_expression
        self.reaction_scores: pd.DataFrame | None = None

    def compute_reaction_penalties(self) -> pd.DataFrame:
        """Compute reaction penalties from gene expression.

        Returns
        -------
        pd.DataFrame
            Reaction penalties matrix (reactions x cells).
        """
        raise NotImplementedError

    def optimize_reactions(self) -> pd.DataFrame:
        """Run COMPASS optimization for each reaction.

        Returns
        -------
        pd.DataFrame
            Reaction consistency scores (reactions x cells).
        """
        raise NotImplementedError

    def score(self) -> pd.DataFrame:
        """Compute full COMPASS scores.

        Returns
        -------
        pd.DataFrame
            Final reaction activity scores.
        """
        raise NotImplementedError
