"""Interactive Streamlit dashboard for metabolic exploration.

This module provides a Streamlit-based interactive dashboard
for exploring metabolic profiles with filtering and drill-down.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad
    import pandas as pd


def create_dashboard(
    reaction_scores: pd.DataFrame,
    pathway_scores: pd.DataFrame | None = None,
    metadata: pd.DataFrame | None = None,
    embedding: dict[str, pd.DataFrame] | None = None,
) -> None:
    """Launch interactive Streamlit dashboard.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores.
    pathway_scores : pd.DataFrame, optional
        Pathway-level scores.
    metadata : pd.DataFrame, optional
        Cell metadata (clusters, conditions, etc.).
    embedding : dict, optional
        Dict of embeddings like {'umap': df, 'tsne': df}.

    Notes
    -----
    This function launches a Streamlit app. Run with:
        streamlit run -m cellmetpro.visualization.dashboard
    """
    raise NotImplementedError


def run_app() -> None:
    """Entry point for running the dashboard as a module."""
    raise NotImplementedError


if __name__ == "__main__":
    run_app()
