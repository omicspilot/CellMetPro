"""AnnData round-trip utilities.

Functions for storing CellMetPro analysis results back into an AnnData
object, following scanpy naming conventions so outputs are immediately
usable in downstream Scanpy / Squidpy workflows.

Conventions used
----------------
- Cell-level score matrices  → ``adata.obsm["X_<key>*"]``
- Cell-level scalar values   → ``adata.obs["<key>*"]``
- Analysis result tables     → ``adata.uns["<key>*"]`` (stored as dicts)

All functions operate **in-place** and return the same ``adata`` for
convenient chaining.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import anndata as ad
import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from cellmetpro.core.compass import CompassResult

logger = logging.getLogger(__name__)


def store_compass_result(
    adata: ad.AnnData,
    result: CompassResult,
    *,
    key: str = "compass",
) -> ad.AnnData:
    """Store COMPASS reaction scores and penalties in an AnnData object.

    Writes the following entries in-place:

    - ``adata.obsm[f"X_{key}"]`` — reaction scores (cells × reactions, float32)
    - ``adata.obsm[f"X_{key}_penalties"]`` — reaction penalties (same shape)
    - ``adata.uns[key]`` — reaction names and config summary
    - ``adata.obsm[f"X_{key}_uptake"]`` — uptake scores, if available
    - ``adata.obsm[f"X_{key}_secretion"]`` — secretion scores, if available

    Parameters
    ----------
    adata : anndata.AnnData
        Target object. Cell order is preserved; ``adata.obs_names`` are used
        as the reference — any cells missing from the COMPASS result receive
        NaN values.
    result : CompassResult
        Output from :class:`~cellmetpro.core.compass.CompassScorer`.
    key : str, default ``"compass"``
        Prefix for all keys written into ``adata``.

    Returns
    -------
    anndata.AnnData
        The same ``adata``, modified in-place.

    Raises
    ------
    ValueError
        If no cells overlap between ``adata`` and the COMPASS result.
    """
    result_cells = result.reaction_scores.columns
    common = adata.obs_names.intersection(result_cells)

    if len(common) == 0:
        raise ValueError(
            "No cells in common between adata and the COMPASS result. "
            "Make sure adata.obs_names match the cell IDs used during scoring."
        )

    if len(common) < adata.n_obs:
        logger.warning(
            "%d cells in adata have no COMPASS scores and will contain NaN.",
            adata.n_obs - len(common),
        )

    # Align to adata cell order (reactions × cells → cells × reactions)
    scores = result.reaction_scores.reindex(columns=adata.obs_names).T
    penalties = result.reaction_penalties.reindex(columns=adata.obs_names).T

    adata.obsm[f"X_{key}"] = scores.values.astype(np.float32)
    adata.obsm[f"X_{key}_penalties"] = penalties.values.astype(np.float32)

    cfg = result.config
    adata.uns[key] = {
        "reaction_names": list(result.reaction_scores.index),
        "n_reactions": len(result.reaction_scores),
        "config": {
            "beta": cfg.beta,
            "and_function": cfg.and_function,
            "or_function": cfg.or_function,
            "lambda_penalty": cfg.lambda_penalty,
            "n_neighbors": cfg.n_neighbors,
            "solver": cfg.solver,
        },
    }

    if result.uptake_scores is not None:
        uptake = result.uptake_scores.reindex(columns=adata.obs_names).T
        adata.obsm[f"X_{key}_uptake"] = uptake.values.astype(np.float32)

    if result.secretion_scores is not None:
        secretion = result.secretion_scores.reindex(columns=adata.obs_names).T
        adata.obsm[f"X_{key}_secretion"] = secretion.values.astype(np.float32)

    logger.info(
        "Stored %d reaction scores for %d cells in adata (key='%s').",
        len(result.reaction_scores),
        adata.n_obs,
        key,
    )
    return adata


def store_differential_result(
    adata: ad.AnnData,
    result: pd.DataFrame,
    comparison: str,
    *,
    key: str = "compass_differential",
) -> ad.AnnData:
    """Store differential analysis results in an AnnData object.

    Results are stored under ``adata.uns[key][comparison]``.  The
    ``comparison`` string becomes the sub-key so multiple comparisons
    can co-exist (e.g. ``"control_vs_treated"``, ``"A_vs_B"``).

    Parameters
    ----------
    adata : anndata.AnnData
        Target object.
    result : pd.DataFrame
        Output of :meth:`~cellmetpro.analysis.DifferentialAnalysis.compare_groups`
        or :meth:`~cellmetpro.analysis.PseudoBulkAnalysis.compare_groups`.
        Must contain a ``"reaction"`` column.
    comparison : str
        Label for this comparison, used as the dict key inside
        ``adata.uns[key]``.
    key : str, default ``"compass_differential"``
        Top-level key in ``adata.uns``.

    Returns
    -------
    anndata.AnnData
        The same ``adata``, modified in-place.

    Raises
    ------
    ValueError
        If ``result`` is missing the required ``"reaction"`` column.
    """
    if "reaction" not in result.columns:
        raise ValueError(
            "result DataFrame must have a 'reaction' column. "
            f"Got columns: {list(result.columns)}"
        )

    if key not in adata.uns:
        adata.uns[key] = {}

    adata.uns[key][comparison] = result.to_dict(orient="list")
    logger.info(
        "Stored differential result '%s' (%d reactions) in adata.uns['%s'].",
        comparison,
        len(result),
        key,
    )
    return adata


def store_pathway_result(
    adata: ad.AnnData,
    result: pd.DataFrame,
    *,
    key: str = "compass_pathway",
    analysis_name: str = "enrichment",
) -> ad.AnnData:
    """Store pathway enrichment results in an AnnData object.

    Results are stored under ``adata.uns[key][analysis_name]``.

    Parameters
    ----------
    adata : anndata.AnnData
        Target object.
    result : pd.DataFrame
        Output of :class:`~cellmetpro.analysis.PathwayAnalyzer` or
        :class:`~cellmetpro.analysis.GOEnrichmentAnalyzer`.
    key : str, default ``"compass_pathway"``
        Top-level key in ``adata.uns``.
    analysis_name : str, default ``"enrichment"``
        Sub-key for this result (use distinct names for multiple runs).

    Returns
    -------
    anndata.AnnData
        The same ``adata``, modified in-place.
    """
    if key not in adata.uns:
        adata.uns[key] = {}

    adata.uns[key][analysis_name] = result.to_dict(orient="list")
    logger.info(
        "Stored pathway result '%s' (%d entries) in adata.uns['%s'].",
        analysis_name,
        len(result),
        key,
    )
    return adata


def store_pseudotime(
    adata: ad.AnnData,
    pseudotime: pd.Series,
    *,
    key: str = "compass_pseudotime",
) -> ad.AnnData:
    """Store pseudotime values in an AnnData object.

    Values are aligned to ``adata.obs_names``; cells not present in
    ``pseudotime`` receive NaN.

    Parameters
    ----------
    adata : anndata.AnnData
        Target object.
    pseudotime : pd.Series
        Pseudotime values indexed by cell names.
    key : str, default ``"compass_pseudotime"``
        Column name created in ``adata.obs``.

    Returns
    -------
    anndata.AnnData
        The same ``adata``, modified in-place.
    """
    adata.obs[key] = pseudotime.reindex(adata.obs_names).values
    logger.info("Stored pseudotime in adata.obs['%s'].", key)
    return adata
