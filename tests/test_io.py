"""Tests for AnnData round-trip utilities (cellmetpro/io.py)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from cellmetpro.io import (
    store_compass_result,
    store_differential_result,
    store_pathway_result,
    store_pseudotime,
)


# ---------------------------------------------------------------------------
# Minimal stubs so tests don't need to import the full COMPASS engine
# ---------------------------------------------------------------------------


@dataclass
class _StubConfig:
    beta: float = 0.95
    and_function: str = "min"
    or_function: str = "sum"
    lambda_penalty: float = 0.0
    n_neighbors: int = 30
    solver: str = "glpk"


@dataclass
class _StubCompassResult:
    reaction_scores: pd.DataFrame
    reaction_penalties: pd.DataFrame
    config: _StubConfig = field(default_factory=_StubConfig)
    uptake_scores: Optional[pd.DataFrame] = None
    secretion_scores: Optional[pd.DataFrame] = None


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def cells():
    return [f"cell_{i}" for i in range(10)]


@pytest.fixture
def reactions():
    return ["R1", "R2", "R3", "R4", "R5"]


@pytest.fixture
def adata(cells):
    """Minimal AnnData: 10 cells × 4 genes."""
    rng = np.random.default_rng(42)
    X = rng.random((len(cells), 4))
    return ad.AnnData(X, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=["G1", "G2", "G3", "G4"]))


@pytest.fixture
def compass_result(cells, reactions):
    """Stub CompassResult with 5 reactions × 10 cells."""
    rng = np.random.default_rng(0)
    scores = pd.DataFrame(rng.random((len(reactions), len(cells))), index=reactions, columns=cells)
    penalties = pd.DataFrame(rng.random((len(reactions), len(cells))), index=reactions, columns=cells)
    return _StubCompassResult(reaction_scores=scores, reaction_penalties=penalties)


# ---------------------------------------------------------------------------
# Tests for store_compass_result
# ---------------------------------------------------------------------------


def test_store_compass_result_adds_obsm_keys(adata, compass_result):
    """store_compass_result must write X_compass and X_compass_penalties into obsm."""
    store_compass_result(adata, compass_result)

    assert "X_compass" in adata.obsm
    assert "X_compass_penalties" in adata.obsm


def test_store_compass_result_shape(adata, compass_result, cells, reactions):
    """obsm arrays must be (cells × reactions)."""
    store_compass_result(adata, compass_result)

    assert adata.obsm["X_compass"].shape == (len(cells), len(reactions))
    assert adata.obsm["X_compass_penalties"].shape == (len(cells), len(reactions))


def test_store_compass_result_dtype(adata, compass_result):
    """Scores should be stored as float32 to save memory."""
    store_compass_result(adata, compass_result)

    assert adata.obsm["X_compass"].dtype == np.float32
    assert adata.obsm["X_compass_penalties"].dtype == np.float32


def test_store_compass_result_values_correct(adata, compass_result, cells, reactions):
    """Values must match the transposed input scores (cell order preserved)."""
    store_compass_result(adata, compass_result)

    for i, cell in enumerate(cells):
        for j, rxn in enumerate(reactions):
            expected = compass_result.reaction_scores.loc[rxn, cell]
            np.testing.assert_almost_equal(
                adata.obsm["X_compass"][i, j], expected, decimal=5
            )


def test_store_compass_result_uns_keys(adata, compass_result, reactions):
    """uns['compass'] must contain reaction names, n_reactions, and config."""
    store_compass_result(adata, compass_result)

    uns = adata.uns["compass"]
    assert uns["reaction_names"] == reactions
    assert uns["n_reactions"] == len(reactions)
    assert "config" in uns
    assert uns["config"]["beta"] == 0.95
    assert uns["config"]["solver"] == "glpk"


def test_store_compass_result_custom_key(adata, compass_result):
    """Custom key prefix must be applied to all written entries."""
    store_compass_result(adata, compass_result, key="mymodel")

    assert "X_mymodel" in adata.obsm
    assert "X_mymodel_penalties" in adata.obsm
    assert "mymodel" in adata.uns


def test_store_compass_result_returns_adata(adata, compass_result):
    """Function must return the same adata object (in-place, chainable)."""
    returned = store_compass_result(adata, compass_result)
    assert returned is adata


def test_store_compass_result_no_overlap_raises(adata, reactions):
    """Raise ValueError when no cells are shared."""
    bad_cells = ["cell_x", "cell_y"]
    rng = np.random.default_rng(0)
    scores = pd.DataFrame(rng.random((len(reactions), 2)), index=reactions, columns=bad_cells)
    penalties = scores.copy()
    result = _StubCompassResult(reaction_scores=scores, reaction_penalties=penalties)

    with pytest.raises(ValueError, match="No cells in common"):
        store_compass_result(adata, result)


def test_store_compass_result_partial_overlap(cells, reactions):
    """Cells missing from the result get NaN in obsm, others are filled correctly."""
    # adata has cells 0-9; result only covers cells 0-4
    rng = np.random.default_rng(1)
    subset = cells[:5]
    scores = pd.DataFrame(rng.random((len(reactions), 5)), index=reactions, columns=subset)
    penalties = scores.copy()
    result = _StubCompassResult(reaction_scores=scores, reaction_penalties=penalties)

    X = rng.random((len(cells), 2))
    adata = ad.AnnData(X, obs=pd.DataFrame(index=cells))

    store_compass_result(adata, result)

    stored = adata.obsm["X_compass"]
    # First 5 cells should have real values, last 5 should be NaN
    assert not np.isnan(stored[:5]).any()
    assert np.isnan(stored[5:]).all()


def test_store_compass_result_with_exchange_scores(adata, cells, reactions):
    """uptake and secretion scores are stored when present."""
    rng = np.random.default_rng(2)
    metabolites = ["M1", "M2"]
    scores = pd.DataFrame(rng.random((len(reactions), len(cells))), index=reactions, columns=cells)
    penalties = scores.copy()
    uptake = pd.DataFrame(rng.random((2, len(cells))), index=metabolites, columns=cells)
    secretion = pd.DataFrame(rng.random((2, len(cells))), index=metabolites, columns=cells)

    result = _StubCompassResult(
        reaction_scores=scores,
        reaction_penalties=penalties,
        uptake_scores=uptake,
        secretion_scores=secretion,
    )

    store_compass_result(adata, result)

    assert "X_compass_uptake" in adata.obsm
    assert "X_compass_secretion" in adata.obsm
    assert adata.obsm["X_compass_uptake"].shape == (len(cells), 2)


def test_store_compass_result_no_exchange_scores_absent(adata, compass_result):
    """uptake/secretion keys must NOT be written when scores are None."""
    store_compass_result(adata, compass_result)

    assert "X_compass_uptake" not in adata.obsm
    assert "X_compass_secretion" not in adata.obsm


# ---------------------------------------------------------------------------
# Tests for store_differential_result
# ---------------------------------------------------------------------------


@pytest.fixture
def diff_result(reactions):
    """Minimal differential analysis DataFrame."""
    rng = np.random.default_rng(3)
    return pd.DataFrame(
        {
            "reaction": reactions,
            "log2fc": rng.standard_normal(len(reactions)),
            "pvalue": rng.random(len(reactions)),
            "padj_bh": rng.random(len(reactions)),
        }
    )


def test_store_differential_result_creates_uns_key(adata, diff_result):
    store_differential_result(adata, diff_result, "control_vs_treated")

    assert "compass_differential" in adata.uns
    assert "control_vs_treated" in adata.uns["compass_differential"]


def test_store_differential_result_values_preserved(adata, diff_result, reactions):
    store_differential_result(adata, diff_result, "A_vs_B")

    stored = adata.uns["compass_differential"]["A_vs_B"]
    assert stored["reaction"] == reactions


def test_store_differential_result_multiple_comparisons(adata, diff_result):
    """Multiple comparisons should co-exist under the same top-level key."""
    store_differential_result(adata, diff_result, "A_vs_B")
    store_differential_result(adata, diff_result, "A_vs_C")

    assert "A_vs_B" in adata.uns["compass_differential"]
    assert "A_vs_C" in adata.uns["compass_differential"]


def test_store_differential_result_custom_key(adata, diff_result):
    store_differential_result(adata, diff_result, "A_vs_B", key="my_diff")

    assert "my_diff" in adata.uns
    assert "compass_differential" not in adata.uns


def test_store_differential_result_missing_reaction_column_raises(adata):
    bad_df = pd.DataFrame({"log2fc": [0.1, 0.2]})
    with pytest.raises(ValueError, match="'reaction' column"):
        store_differential_result(adata, bad_df, "A_vs_B")


def test_store_differential_result_returns_adata(adata, diff_result):
    returned = store_differential_result(adata, diff_result, "A_vs_B")
    assert returned is adata


# ---------------------------------------------------------------------------
# Tests for store_pathway_result
# ---------------------------------------------------------------------------


@pytest.fixture
def pathway_result():
    return pd.DataFrame(
        {
            "pathway": ["Glycolysis", "TCA cycle", "Fatty acid oxidation"],
            "pvalue": [0.01, 0.04, 0.2],
            "padj": [0.03, 0.08, 0.4],
            "n_reactions": [10, 8, 12],
        }
    )


def test_store_pathway_result_creates_uns_key(adata, pathway_result):
    store_pathway_result(adata, pathway_result)

    assert "compass_pathway" in adata.uns
    assert "enrichment" in adata.uns["compass_pathway"]


def test_store_pathway_result_values_preserved(adata, pathway_result):
    store_pathway_result(adata, pathway_result)

    stored = adata.uns["compass_pathway"]["enrichment"]
    assert stored["pathway"] == ["Glycolysis", "TCA cycle", "Fatty acid oxidation"]


def test_store_pathway_result_multiple_analyses(adata, pathway_result):
    """Multiple named analyses co-exist under the same key."""
    store_pathway_result(adata, pathway_result, analysis_name="subsystem")
    store_pathway_result(adata, pathway_result, analysis_name="go_bp")

    assert "subsystem" in adata.uns["compass_pathway"]
    assert "go_bp" in adata.uns["compass_pathway"]


def test_store_pathway_result_custom_key(adata, pathway_result):
    store_pathway_result(adata, pathway_result, key="my_pathway")

    assert "my_pathway" in adata.uns


def test_store_pathway_result_returns_adata(adata, pathway_result):
    returned = store_pathway_result(adata, pathway_result)
    assert returned is adata


# ---------------------------------------------------------------------------
# Tests for store_pseudotime
# ---------------------------------------------------------------------------


@pytest.fixture
def pseudotime(cells):
    rng = np.random.default_rng(4)
    return pd.Series(rng.random(len(cells)), index=cells)


def test_store_pseudotime_creates_obs_column(adata, pseudotime):
    store_pseudotime(adata, pseudotime)

    assert "compass_pseudotime" in adata.obs.columns


def test_store_pseudotime_values_correct(adata, pseudotime, cells):
    store_pseudotime(adata, pseudotime)

    for cell in cells:
        np.testing.assert_almost_equal(
            adata.obs.loc[cell, "compass_pseudotime"], pseudotime[cell]
        )


def test_store_pseudotime_custom_key(adata, pseudotime):
    store_pseudotime(adata, pseudotime, key="my_pseudotime")

    assert "my_pseudotime" in adata.obs.columns
    assert "compass_pseudotime" not in adata.obs.columns


def test_store_pseudotime_partial_overlap(cells):
    """Cells not in pseudotime Series receive NaN."""
    rng = np.random.default_rng(5)
    X = rng.random((len(cells), 2))
    adata = ad.AnnData(X, obs=pd.DataFrame(index=cells))

    # Only provide pseudotime for first 5 cells
    pt = pd.Series(rng.random(5), index=cells[:5])
    store_pseudotime(adata, pt)

    stored = adata.obs["compass_pseudotime"].values
    assert not np.isnan(stored[:5]).any()
    assert np.isnan(stored[5:]).all()


def test_store_pseudotime_returns_adata(adata, pseudotime):
    returned = store_pseudotime(adata, pseudotime)
    assert returned is adata


# ---------------------------------------------------------------------------
# Integration: chaining all store functions together
# ---------------------------------------------------------------------------


def test_chaining_all_store_functions(adata, compass_result, diff_result, pathway_result, pseudotime):
    """All four store functions can be called on the same adata."""
    store_compass_result(adata, compass_result)
    store_differential_result(adata, diff_result, "A_vs_B")
    store_pathway_result(adata, pathway_result)
    store_pseudotime(adata, pseudotime)

    assert "X_compass" in adata.obsm
    assert "compass_differential" in adata.uns
    assert "compass_pathway" in adata.uns
    assert "compass_pseudotime" in adata.obs.columns
