"""Tests for core module."""

import importlib.util

import numpy as np
import pandas as pd
import pytest

_has_parquet = (
    importlib.util.find_spec("pyarrow") is not None
    or importlib.util.find_spec("fastparquet") is not None
)

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def simple_model():
    """Create a simple test metabolic model."""
    import cobra

    model = cobra.Model("test_model")

    # Create metabolites
    A = cobra.Metabolite("A", compartment="c")
    B = cobra.Metabolite("B", compartment="c")
    C = cobra.Metabolite("C", compartment="c")

    # Create reactions with GPR rules
    r1 = cobra.Reaction("R1")
    r1.add_metabolites({A: -1, B: 1})
    r1.bounds = (0, 1000)
    r1.gene_reaction_rule = "gene1"

    r2 = cobra.Reaction("R2")
    r2.add_metabolites({B: -1, C: 1})
    r2.bounds = (0, 1000)
    r2.gene_reaction_rule = "gene2 and gene3"

    r3 = cobra.Reaction("R3")
    r3.add_metabolites({A: -1, C: 1})
    r3.bounds = (0, 1000)
    r3.gene_reaction_rule = "gene1 or gene4"

    # Exchange reactions
    ex_A = cobra.Reaction("EX_A")
    ex_A.add_metabolites({A: -1})
    ex_A.bounds = (-10, 0)

    ex_C = cobra.Reaction("EX_C")
    ex_C.add_metabolites({C: -1})
    ex_C.bounds = (0, 1000)

    model.add_reactions([r1, r2, r3, ex_A, ex_C])
    model.objective = "EX_C"

    return model


@pytest.fixture
def expression_df():
    """Create a simple expression DataFrame (genes x cells)."""
    genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
    cells = ["cell1", "cell2", "cell3", "cell4", "cell5"]

    np.random.seed(42)
    data = np.random.rand(5, 5) * 100

    return pd.DataFrame(data, index=genes, columns=cells)


@pytest.fixture
def adata(expression_df):
    """Create a simple AnnData object."""
    import anndata as ad

    # AnnData is cells x genes
    return ad.AnnData(expression_df.T)


@pytest.fixture
def tmp_csv(tmp_path, expression_df):
    """Create a temporary CSV file."""
    csv_path = tmp_path / "expression.csv"
    expression_df.to_csv(csv_path)
    return csv_path


@pytest.fixture
def tmp_h5ad(tmp_path, adata):
    """Create a temporary h5ad file."""
    h5ad_path = tmp_path / "expression.h5ad"
    adata.write(h5ad_path)
    return h5ad_path


# ============================================================================
# DataLoader Tests
# ============================================================================


class TestDataLoader:
    """Tests for DataLoader class."""

    def test_load_csv(self, tmp_csv):
        """Test loading CSV file."""
        from cellmetpro.core.preprocessing import DataLoader

        loader = DataLoader(tmp_csv)
        adata = loader.load()

        assert adata is not None
        assert adata.n_obs == 5  # cells
        assert adata.n_vars == 5  # genes

    def test_load_h5ad(self, tmp_h5ad):
        """Test loading h5ad file."""
        from cellmetpro.core.preprocessing import DataLoader

        loader = DataLoader(tmp_h5ad)
        adata = loader.load()

        assert adata is not None
        assert adata.n_obs == 5
        assert adata.n_vars == 5

    def test_load_nonexistent_file(self):
        """Test loading non-existent file raises error."""
        from cellmetpro.core.preprocessing import DataLoader

        loader = DataLoader("/nonexistent/path/file.csv")

        with pytest.raises(FileNotFoundError):
            loader.load()

    def test_load_unsupported_format(self, tmp_path):
        """Test loading unsupported format raises error."""
        from cellmetpro.core.preprocessing import DataLoader

        bad_file = tmp_path / "data.xyz"
        bad_file.write_text("some data")

        loader = DataLoader(bad_file)

        with pytest.raises(ValueError, match="Unsupported file format"):
            loader.load()


# ============================================================================
# Seurat Loader Tests
# ============================================================================


# Check if rpy2 is available
HAS_RPY2 = importlib.util.find_spec("rpy2") is not None


class TestSeuratLoader:
    """Tests for Seurat object loading."""

    def test_load_seurat_import_error(self):
        """Test that import error is raised when rpy2 not installed."""
        from cellmetpro.core.preprocessing import load_seurat_rds

        # This test verifies the error message is correct
        # If rpy2 IS installed, it will try to load a non-existent file
        if not HAS_RPY2:
            with pytest.raises(ImportError, match="rpy2 is required"):
                load_seurat_rds("fake_file.rds")

    def test_load_seurat_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        if not HAS_RPY2:
            pytest.skip("rpy2 not installed")

        from cellmetpro.core.preprocessing import load_seurat_rds

        with pytest.raises(FileNotFoundError):
            load_seurat_rds("/nonexistent/seurat_object.rds")

    def test_dataloader_recognizes_rds_extension(self, tmp_path):
        """Test that DataLoader recognizes .rds extension."""
        from cellmetpro.core.preprocessing import DataLoader

        # Create a dummy .rds file
        rds_file = tmp_path / "test.rds"
        rds_file.write_bytes(b"fake rds content")

        loader = DataLoader(rds_file)

        # Should try to load as Seurat (and fail with appropriate error)
        if not HAS_RPY2:
            with pytest.raises(ImportError, match="rpy2 is required"):
                loader.load()
        else:
            # If rpy2 is installed, it will fail because the file is not a real RDS
            with pytest.raises(Exception):  # RuntimeError or R error
                loader.load()


# ============================================================================
# Preprocessing Tests
# ============================================================================


class TestPreprocessing:
    """Tests for preprocessing functions."""

    def test_normalize_expression(self, adata):
        """Test expression normalization."""
        from cellmetpro.core.preprocessing import normalize_expression

        normalized = normalize_expression(adata, target_sum=1e4, log_transform=True)

        # Check normalization was applied
        assert "normalization" in normalized.uns
        assert normalized.uns["normalization"]["target_sum"] == 1e4
        assert normalized.uns["normalization"]["log_transform"] is True

    def test_normalize_expression_no_log(self, adata):
        """Test normalization without log transform."""
        from cellmetpro.core.preprocessing import normalize_expression

        normalized = normalize_expression(adata, target_sum=1e6, log_transform=False)

        assert normalized.uns["normalization"]["log_transform"] is False

    def test_filter_cells(self, adata):
        """Test cell filtering."""
        from cellmetpro.core.preprocessing import filter_cells

        filtered = filter_cells(adata, min_genes=1)
        assert filtered.n_obs <= adata.n_obs

    def test_filter_genes(self, adata):
        """Test gene filtering."""
        from cellmetpro.core.preprocessing import filter_genes

        filtered = filter_genes(adata, min_cells=1)
        assert filtered.n_vars <= adata.n_vars

    def test_to_dataframe(self, adata):
        """Test AnnData to DataFrame conversion."""
        from cellmetpro.core.preprocessing import to_dataframe

        # Genes as rows
        df = to_dataframe(adata, genes_as_rows=True)
        assert df.shape == (adata.n_vars, adata.n_obs)

        # Cells as rows
        df = to_dataframe(adata, genes_as_rows=False)
        assert df.shape == (adata.n_obs, adata.n_vars)


# ============================================================================
# CompassScorer Tests
# ============================================================================


class TestCompassScorer:
    """Tests for CompassScorer class."""

    def test_init(self, simple_model, expression_df):
        """Test CompassScorer initialization."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)

        assert scorer.model is not None
        assert len(scorer.cell_names) == 5
        assert len(scorer.gene_names) == 5

    def test_init_with_adata(self, simple_model, adata):
        """Test CompassScorer initialization with AnnData."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, adata)

        assert scorer.model is not None
        assert len(scorer.cell_names) == 5

    def test_compute_reaction_penalties(self, simple_model, expression_df):
        """Test reaction penalty computation."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        assert penalties is not None
        assert penalties.shape[1] == 5  # 5 cells
        assert penalties.shape[0] > 0  # Some reactions with GPR

    def test_penalty_values_valid(self, simple_model, expression_df):
        """Test that penalty values are in valid range."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # Penalties should be between 0 and 1 after normalization
        assert penalties.min().min() >= 0
        assert penalties.max().max() <= 1

    def test_gpr_evaluation_and(self, simple_model, expression_df):
        """Test GPR evaluation with AND operator."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # R2 has "gene2 and gene3" - should have penalty
        assert "R2" in penalties.index

    def test_gpr_evaluation_or(self, simple_model, expression_df):
        """Test GPR evaluation with OR operator."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # R3 has "gene1 or gene4"
        assert "R3" in penalties.index


class TestCompassConfig:
    """Tests for CompassConfig dataclass."""

    def test_default_config(self):
        """Test default configuration values."""
        from cellmetpro.core.compass import CompassConfig

        config = CompassConfig()

        assert config.beta == 0.95
        assert config.exchange_limit == 1000.0
        assert config.and_function == "min"
        assert config.or_function == "sum"
        assert config.lambda_penalty == 0.0
        assert config.n_neighbors == 30
        assert config.n_processes == 1
        assert config.show_progress is True

    def test_custom_config(self):
        """Test custom configuration."""
        from cellmetpro.core.compass import CompassConfig

        config = CompassConfig(
            beta=0.9,
            lambda_penalty=0.5,
            n_processes=4,
            show_progress=False,
        )

        assert config.beta == 0.9
        assert config.lambda_penalty == 0.5
        assert config.n_processes == 4
        assert config.show_progress is False


# ============================================================================
# FluxBalanceAnalyzer Tests
# ============================================================================


class TestFluxBalanceAnalyzer:
    """Tests for FluxBalanceAnalyzer class."""

    def test_init(self, simple_model):
        """Test FBA initialization."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        assert fba.model is not None
        assert fba.solution is None

    def test_optimize(self, simple_model):
        """Test FBA optimization."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.optimize()

        assert fluxes is not None
        assert isinstance(fluxes, pd.Series)
        assert fba.solution is not None
        assert fba.solution.status == "optimal"

    def test_set_bounds(self, simple_model):
        """Test setting reaction bounds."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        fba.set_bounds("R1", lower=0, upper=100)

        assert fba.model.reactions.get_by_id("R1").upper_bound == 100
        assert "R1" in fba._original_bounds

    def test_reset_bounds(self, simple_model):
        """Test resetting reaction bounds."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        original_ub = fba.model.reactions.get_by_id("R1").upper_bound

        fba.set_bounds("R1", upper=100)
        fba.reset_bounds("R1")

        assert fba.model.reactions.get_by_id("R1").upper_bound == original_ub

    def test_knockout(self, simple_model):
        """Test reaction knockout."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.knockout("R1")

        assert fluxes is not None
        # After knockout, R1 should have zero flux
        # (model is restored after knockout due to context manager)

    def test_flux_variability(self, simple_model):
        """Test flux variability analysis."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fva = fba.flux_variability(reactions=["R1", "R2"])

        assert fva is not None
        assert "minimum" in fva.columns
        assert "maximum" in fva.columns

    def test_summary(self, simple_model):
        """Test FBA summary."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        # Before optimization
        summary = fba.summary()
        assert "No solution" in summary

        # After optimization
        fba.optimize()
        summary = fba.summary()
        assert "optimal" in summary.lower()


# ============================================================================
# Microclustering Tests
# ============================================================================


class TestMicroclustering:
    """Tests for microclustering functions."""

    def test_microcluster_config(self):
        """Test MicroclusterConfig defaults."""
        from cellmetpro.core.microclustering import MicroclusterConfig

        config = MicroclusterConfig()

        assert config.cells_per_cluster == 100
        assert config.n_neighbors == 30
        assert config.n_pcs == 20
        assert config.method == "leiden"

    def test_microcluster_small_dataset(self, expression_df):
        """Test microclustering on small dataset."""
        from cellmetpro.core.microclustering import MicroclusterConfig, microcluster

        config = MicroclusterConfig(cells_per_cluster=2, method="kmeans")
        result = microcluster(expression_df, config)

        assert result is not None
        assert result.n_clusters > 0
        assert len(result.cluster_labels) == 5
        assert result.pooled_expression.shape[1] == result.n_clusters

    def test_unpool_results(self, expression_df):
        """Test unpooling cluster results."""
        from cellmetpro.core.microclustering import (
            MicroclusterConfig,
            microcluster,
            unpool_results,
        )

        config = MicroclusterConfig(cells_per_cluster=2, method="kmeans")
        mc_result = microcluster(expression_df, config)

        # Create fake cluster-level results
        cluster_results = pd.DataFrame(
            np.random.rand(10, mc_result.n_clusters),
            index=[f"rxn_{i}" for i in range(10)],
            columns=mc_result.pooled_expression.columns,
        )

        unpooled = unpool_results(cluster_results, mc_result)

        assert unpooled.shape[1] == 5  # Back to 5 cells

    def test_filter_genes_fano(self, expression_df):
        """Test Fano factor gene filtering."""
        from cellmetpro.core.microclustering import filter_genes_fano

        filtered = filter_genes_fano(expression_df, n_genes=3)

        assert filtered.shape[0] == 3
        assert filtered.shape[1] == 5


# ============================================================================
# Cache Tests
# ============================================================================


class TestCache:
    """Tests for caching functionality."""

    def test_compass_cache_init(self, tmp_path):
        """Test CompassCache initialization."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test_model")

        assert cache.cache_dir == tmp_path
        assert cache.model_id == "test_model"

    def test_save_load_max_fluxes(self, tmp_path):
        """Test saving and loading max fluxes."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")

        max_fluxes = {"R1": 100.0, "R2": 50.0, "R3": 75.0}
        cache.save_max_fluxes(max_fluxes)

        loaded = cache.load_max_fluxes()

        assert loaded is not None
        assert loaded["R1"] == 100.0
        assert loaded["R2"] == 50.0

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_save_load_reaction_scores(self, tmp_path):
        """Test saving and loading reaction scores."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")

        scores = pd.DataFrame(
            {"cell1": [1.0, 2.0], "cell2": [3.0, 4.0]},
            index=["R1", "R2"],
        )
        cache.save_reaction_scores(scores, "sample1")

        loaded = cache.load_reaction_scores("sample1")

        assert loaded is not None
        assert loaded.shape == (2, 2)

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_has_sample(self, tmp_path):
        """Test checking for cached samples."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")

        assert not cache.has_sample("sample1")

        scores = pd.DataFrame({"cell1": [1.0]}, index=["R1"])
        cache.save_reaction_scores(scores, "sample1")

        assert cache.has_sample("sample1")

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_clear_cache(self, tmp_path):
        """Test clearing cache."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")

        scores = pd.DataFrame({"cell1": [1.0]}, index=["R1"])
        cache.save_reaction_scores(scores, "sample1")

        cache.clear("sample1")

        assert not cache.has_sample("sample1")

    def test_memory_cache(self):
        """Test in-memory cache."""
        from cellmetpro.core.cache import MemoryCache

        cache = MemoryCache(max_size=3)

        cache.set("key1", "value1")
        cache.set("key2", "value2")
        cache.set("key3", "value3")

        assert cache.get("key1") == "value1"
        assert len(cache) == 3

        # Adding 4th item should evict oldest
        cache.set("key4", "value4")
        assert len(cache) == 3
        assert "key2" not in cache  # key1 was accessed, so key2 is oldest


# ============================================================================
# Model Loading Tests
# ============================================================================


class TestModelLoading:
    """Tests for model loading functions."""

    def test_load_model_from_file_sbml(self, simple_model, tmp_path):
        """Test loading model from SBML file."""
        import cobra

        from cellmetpro.models import load_model_from_file

        # Save model as SBML
        sbml_path = tmp_path / "model.xml"
        cobra.io.write_sbml_model(simple_model, str(sbml_path))

        # Load it back
        loaded = load_model_from_file(sbml_path)

        assert loaded is not None
        assert len(loaded.reactions) == len(simple_model.reactions)

    def test_load_model_from_file_json(self, simple_model, tmp_path):
        """Test loading model from JSON file."""
        import cobra

        from cellmetpro.models import load_model_from_file

        # Save model as JSON
        json_path = tmp_path / "model.json"
        cobra.io.save_json_model(simple_model, str(json_path))

        # Load it back
        loaded = load_model_from_file(json_path)

        assert loaded is not None
        assert len(loaded.reactions) == len(simple_model.reactions)

    def test_load_model_nonexistent(self):
        """Test loading non-existent model file."""
        from cellmetpro.models import load_model_from_file

        with pytest.raises(FileNotFoundError):
            load_model_from_file("/nonexistent/model.xml")

    def test_load_gem_invalid_organism(self):
        """Test loading model with invalid organism name."""
        from cellmetpro.models import load_gem

        with pytest.raises((ValueError, FileNotFoundError)):
            load_gem("invalid_organism_xyz")

    def test_get_reaction_gene_mapping(self, simple_model):
        """Test getting reaction-gene mapping."""
        from cellmetpro.models import get_reaction_gene_mapping

        mapping = get_reaction_gene_mapping(simple_model)

        assert "R1" in mapping
        assert "gene1" in mapping["R1"]

    def test_get_subsystem_reactions(self, simple_model):
        """Test getting subsystem-reaction mapping."""
        from cellmetpro.models import get_subsystem_reactions

        # Simple model has no subsystems, should return empty dict
        subsystems = get_subsystem_reactions(simple_model)
        assert isinstance(subsystems, dict)


# ============================================================================
# Integration Tests
# ============================================================================


class TestCompassPenaltySmoothing:
    """Tests for COMPASS penalty smoothing functionality."""

    def test_penalty_smoothing(self, simple_model, expression_df):
        """Test KNN-based penalty smoothing."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(lambda_penalty=0.5, n_neighbors=2, show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)
        penalties = scorer.compute_reaction_penalties()

        assert penalties is not None
        # Smoothed penalties should still be in valid range
        assert penalties.min().min() >= 0
        assert penalties.max().max() <= 1

    def test_penalty_smoothing_zero_lambda(self, simple_model, expression_df):
        """Test that zero lambda means no smoothing."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        # With no smoothing
        config1 = CompassConfig(lambda_penalty=0.0, show_progress=False)
        scorer1 = CompassScorer(simple_model, expression_df, config1)
        penalties1 = scorer1.compute_reaction_penalties()

        # Re-run with same config (should be identical)
        config2 = CompassConfig(lambda_penalty=0.0, show_progress=False)
        scorer2 = CompassScorer(simple_model, expression_df, config2)
        penalties2 = scorer2.compute_reaction_penalties()

        pd.testing.assert_frame_equal(penalties1, penalties2)


class TestRunCompassFunction:
    """Tests for the run_compass convenience function."""

    def test_run_compass_basic(self, simple_model, expression_df):
        """Test run_compass with basic inputs."""
        from cellmetpro.core.compass import CompassConfig, run_compass

        config = CompassConfig(show_progress=False)
        result = run_compass(simple_model, expression_df, config)

        assert result is not None
        assert result.reaction_penalties is not None
        assert result.reaction_scores is not None

    def test_run_compass_with_exchange(self, simple_model, expression_df):
        """Test run_compass with exchange score computation."""
        from cellmetpro.core.compass import CompassConfig, run_compass

        config = CompassConfig(show_progress=False)
        result = run_compass(simple_model, expression_df, config, compute_exchange=True)

        assert result is not None
        assert result.uptake_scores is not None
        assert result.secretion_scores is not None

    def test_run_compass_no_config(self, simple_model, expression_df):
        """Test run_compass without explicit config."""
        from cellmetpro.core.compass import run_compass

        result = run_compass(simple_model, expression_df, config=None)

        assert result is not None
        # Should use default config
        assert result.config.beta == 0.95


class TestIntegration:
    """Integration tests combining multiple modules."""

    def test_full_compass_pipeline(self, simple_model, expression_df):
        """Test full COMPASS pipeline."""
        from cellmetpro.core.compass import CompassConfig, run_compass

        config = CompassConfig(beta=0.95, n_processes=1, show_progress=False)
        result = run_compass(simple_model, expression_df, config)

        assert result is not None
        assert result.reaction_penalties is not None
        assert result.reaction_scores is not None
        assert result.config.beta == 0.95

    def test_compass_with_adata(self, simple_model, adata):
        """Test COMPASS with AnnData input."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, adata, config)
        penalties = scorer.compute_reaction_penalties()

        assert penalties is not None
        assert penalties.shape[1] == adata.n_obs


# ============================================================================
# Edge Case Tests - Preprocessing
# ============================================================================


class TestPreprocessingEdgeCases:
    """Edge case tests for preprocessing functions."""

    def test_normalize_expression_with_zeros(self):
        """Test normalization handles cells with all zeros."""
        import anndata as ad

        from cellmetpro.core.preprocessing import normalize_expression

        # Create data with one cell having all zeros
        data = np.array([[0, 0, 0], [1, 2, 3], [4, 5, 6]], dtype=float)
        adata = ad.AnnData(data)

        # Should not raise division by zero
        normalized = normalize_expression(adata, target_sum=1e4)
        assert not np.any(np.isnan(normalized.X))

    def test_normalize_expression_single_cell(self):
        """Test normalization with a single cell."""
        import anndata as ad

        from cellmetpro.core.preprocessing import normalize_expression

        data = np.array([[1, 2, 3]], dtype=float)
        adata = ad.AnnData(data)

        normalized = normalize_expression(adata, target_sum=1e4)
        assert normalized.n_obs == 1

    def test_normalize_expression_single_gene(self):
        """Test normalization with a single gene."""
        import anndata as ad

        from cellmetpro.core.preprocessing import normalize_expression

        data = np.array([[1], [2], [3]], dtype=float)
        adata = ad.AnnData(data)

        normalized = normalize_expression(adata, target_sum=1e4)
        assert normalized.n_vars == 1

    def test_filter_cells_all_filtered(self):
        """Test filtering when all cells would be filtered out."""
        import anndata as ad

        from cellmetpro.core.preprocessing import filter_cells

        # All cells have < 10 genes expressed
        data = np.array([[1, 0, 0], [0, 1, 0]], dtype=float)
        adata = ad.AnnData(data)

        # High threshold should result in no cells passing
        filtered = filter_cells(adata, min_genes=100)
        assert filtered.n_obs == 0

    def test_filter_genes_all_filtered(self):
        """Test filtering when all genes would be filtered out."""
        import anndata as ad

        from cellmetpro.core.preprocessing import filter_genes

        # All genes expressed in < 2 cells
        data = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 1]], dtype=float)
        adata = ad.AnnData(data)

        filtered = filter_genes(adata, min_cells=100)
        assert filtered.n_vars == 0

    def test_to_dataframe_with_nan_values(self):
        """Test DataFrame conversion with NaN values."""
        import anndata as ad

        from cellmetpro.core.preprocessing import to_dataframe

        data = np.array([[1, np.nan, 3], [4, 5, np.nan]], dtype=float)
        adata = ad.AnnData(data)

        df = to_dataframe(adata, genes_as_rows=True)
        assert df.isna().sum().sum() == 2  # Two NaN values preserved


# ============================================================================
# Edge Case Tests - CompassScorer
# ============================================================================


class TestCompassScorerEdgeCases:
    """Edge case tests for CompassScorer."""

    def test_compass_no_matching_genes(self, simple_model):
        """Test COMPASS with no matching genes between expression and model."""
        from cellmetpro.core.compass import CompassScorer

        # Expression data with completely different gene names
        expression_df = pd.DataFrame(
            np.random.rand(5, 5) * 100,
            index=["OTHER1", "OTHER2", "OTHER3", "OTHER4", "OTHER5"],
            columns=[f"cell{i}" for i in range(5)],
        )

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # Should still return a result (with default penalties)
        assert penalties is not None

    def test_compass_single_cell(self, simple_model):
        """Test COMPASS with a single cell."""
        from cellmetpro.core.compass import CompassScorer

        expression_df = pd.DataFrame(
            np.random.rand(5, 1) * 100,
            index=["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            columns=["cell0"],
        )

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        assert penalties.shape[1] == 1

    def test_compass_with_inf_values(self, simple_model):
        """Test COMPASS handles infinite values gracefully."""
        from cellmetpro.core.compass import CompassScorer

        expression_df = pd.DataFrame(
            [[np.inf, 1, 2], [3, np.inf, 4], [5, 6, 7], [8, 9, 10], [11, 12, 13]],
            index=["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            columns=["cell0", "cell1", "cell2"],
        )

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # Should not have any inf or nan in result
        assert not np.any(np.isinf(penalties.values))


# ============================================================================
# COMPASS Performance and Optimization Tests
# ============================================================================


class TestCompassScoreAndOptimize:
    """Tests for COMPASS score() and optimize_reactions() methods."""

    def test_score_basic(self, simple_model, expression_df):
        """Test the score() method returns complete results."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)
        result = scorer.score()

        assert result is not None
        assert result.reaction_penalties is not None
        assert result.reaction_scores is not None
        assert result.config.beta == 0.95

    def test_score_with_exchange(self, simple_model, expression_df):
        """Test score() with exchange reaction scoring."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)
        result = scorer.score(compute_exchange=True)

        assert result is not None
        assert result.uptake_scores is not None
        assert result.secretion_scores is not None
        assert isinstance(result.uptake_scores, pd.DataFrame)
        assert isinstance(result.secretion_scores, pd.DataFrame)

    def test_optimize_reactions_basic(self, simple_model, expression_df):
        """Test optimize_reactions() directly."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False, cache_max_fluxes=True)
        scorer = CompassScorer(simple_model, expression_df, config)

        penalties = scorer.compute_reaction_penalties()
        scores = scorer.optimize_reactions(penalties)

        assert scores is not None
        assert isinstance(scores, pd.DataFrame)
        assert scores.shape[1] == len(expression_df.columns)  # Same number of cells

    def test_optimize_reactions_without_penalties(self, simple_model, expression_df):
        """Test optimize_reactions() computes penalties if not provided."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Call without penalties - should compute them internally
        scores = scorer.optimize_reactions(penalties=None)

        assert scores is not None
        assert isinstance(scores, pd.DataFrame)

    def test_optimize_with_batch_processing(self, simple_model):
        """Test optimization with batch processing."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        # Create larger expression data
        np.random.seed(42)
        genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        cells = [f"cell{i}" for i in range(15)]
        expr_df = pd.DataFrame(
            np.random.rand(5, 15) * 100,
            index=genes,
            columns=cells,
        )

        config = CompassConfig(batch_size=5, show_progress=False)
        scorer = CompassScorer(simple_model, expr_df, config)
        result = scorer.score()

        assert result.reaction_scores.shape[1] == 15

    def test_precompute_max_fluxes(self, simple_model, expression_df):
        """Test that max fluxes are precomputed correctly."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(cache_max_fluxes=True, show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Compute penalties first
        penalties = scorer.compute_reaction_penalties()

        # Get internal reactions
        internal_reactions = [
            rxn.id
            for rxn in simple_model.reactions
            if rxn.id in penalties.index and not rxn.boundary
        ]

        # Precompute max fluxes
        model_copy = simple_model.copy()
        scorer._precompute_max_fluxes(model_copy, internal_reactions)

        # Check cache is populated
        assert len(scorer._max_flux_cache) > 0
        for rxn_id in internal_reactions:
            assert rxn_id in scorer._max_flux_cache
            assert scorer._max_flux_cache[rxn_id] >= 0

    def test_reaction_scores_property(self, simple_model, expression_df):
        """Test the reaction_scores property."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Access property - should trigger computation
        scores = scorer.reaction_scores

        assert scores is not None
        assert isinstance(scores, pd.DataFrame)

    def test_reaction_penalties_property(self, simple_model, expression_df):
        """Test the reaction_penalties property."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Access property - should trigger computation
        penalties = scorer.reaction_penalties

        assert penalties is not None
        assert isinstance(penalties, pd.DataFrame)

    def test_score_config_preserved(self, simple_model, expression_df):
        """Test that config is preserved in result."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(beta=0.9, lambda_penalty=0.1, show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)
        result = scorer.score()

        assert result.config.beta == 0.9
        assert result.config.lambda_penalty == 0.1

    def test_compass_result_dataclass(self, simple_model, expression_df):
        """Test CompassResult dataclass structure."""
        from cellmetpro.core.compass import CompassConfig, CompassResult, CompassScorer

        config = CompassConfig(show_progress=False)
        scorer = CompassScorer(simple_model, expression_df, config)
        result = scorer.score()

        assert isinstance(result, CompassResult)
        assert hasattr(result, "reaction_penalties")
        assert hasattr(result, "reaction_scores")
        assert hasattr(result, "uptake_scores")
        assert hasattr(result, "secretion_scores")
        assert hasattr(result, "config")


class TestCompassOptimizations:
    """Tests for COMPASS performance optimizations."""

    def test_parsed_gpr_basic(self):
        """Test ParsedGPR correctly parses simple rules."""
        from cellmetpro.core.compass import ParsedGPR

        # Test single gene
        gpr = ParsedGPR("GENE1")
        assert gpr.genes == ["GENE1"]
        assert gpr.tree == "GENE1"

        # Test AND rule
        gpr = ParsedGPR("GENE1 and GENE2")
        assert set(gpr.genes) == {"GENE1", "GENE2"}
        assert gpr.tree[0] == "AND"

        # Test OR rule
        gpr = ParsedGPR("GENE1 or GENE2")
        assert set(gpr.genes) == {"GENE1", "GENE2"}
        assert gpr.tree[0] == "OR"

    def test_parsed_gpr_nested(self):
        """Test ParsedGPR handles nested rules."""
        from cellmetpro.core.compass import ParsedGPR

        gpr = ParsedGPR("GENE1 and (GENE2 or GENE3)")
        assert set(gpr.genes) == {"GENE1", "GENE2", "GENE3"}
        assert gpr.tree[0] == "AND"

    def test_parsed_gpr_evaluation(self):
        """Test ParsedGPR vectorized evaluation."""
        from cellmetpro.core.compass import ParsedGPR

        # Create expression matrix (3 genes x 4 cells)
        expr = np.array(
            [
                [1.0, 2.0, 3.0, 4.0],  # GENE1
                [5.0, 6.0, 7.0, 8.0],  # GENE2
                [9.0, 10.0, 11.0, 12.0],  # GENE3
            ]
        )
        gene_to_idx = {"GENE1": 0, "GENE2": 1, "GENE3": 2}

        # Test single gene
        gpr = ParsedGPR("GENE1")
        result = gpr.evaluate(expr, gene_to_idx, np.min, np.max)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0, 4.0])

        # Test AND (min)
        gpr = ParsedGPR("GENE1 and GENE2")
        result = gpr.evaluate(expr, gene_to_idx, np.min, np.max)
        np.testing.assert_array_equal(result, [1.0, 2.0, 3.0, 4.0])  # min of each col

        # Test OR (max)
        gpr = ParsedGPR("GENE1 or GENE2")
        result = gpr.evaluate(expr, gene_to_idx, np.min, np.max)
        np.testing.assert_array_equal(result, [5.0, 6.0, 7.0, 8.0])  # max of each col

    def test_compass_config_defaults(self):
        """Test new CompassConfig options have correct defaults."""
        from cellmetpro.core.compass import CompassConfig

        config = CompassConfig()
        assert config.batch_size == 100
        assert config.cache_max_fluxes is True
        assert config.precompute_gpr is True

    def test_compass_precomputed_gpr(self, simple_model, expression_df):
        """Test COMPASS with precomputed GPR rules."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(precompute_gpr=True)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Should have pre-parsed GPRs
        assert len(scorer._parsed_gprs) > 0

        # Compute penalties should work
        penalties = scorer.compute_reaction_penalties()
        assert not penalties.empty

    def test_compass_cached_max_fluxes(self, simple_model, expression_df):
        """Test COMPASS caches max fluxes correctly."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        config = CompassConfig(cache_max_fluxes=True)
        scorer = CompassScorer(simple_model, expression_df, config)

        # Run penalty computation
        penalties = scorer.compute_reaction_penalties()

        # Pre-compute max fluxes
        model = simple_model.copy()
        internal_reactions = [
            rxn.id
            for rxn in model.reactions
            if rxn.id in penalties.index and not rxn.boundary
        ]
        scorer._precompute_max_fluxes(model, internal_reactions)

        # Max flux cache should be populated
        assert len(scorer._max_flux_cache) > 0

    def test_compass_batch_processing(self, simple_model):
        """Test COMPASS batch processing for larger datasets."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        # Create larger expression data (20 cells)
        np.random.seed(42)
        genes = ["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"]
        cells = [f"cell{i}" for i in range(20)]
        expr_df = pd.DataFrame(
            np.random.rand(5, 20) * 100,
            index=genes,
            columns=cells,
        )

        config = CompassConfig(batch_size=5)  # Small batches
        scorer = CompassScorer(simple_model, expr_df, config)

        penalties = scorer.compute_reaction_penalties()
        assert penalties.shape[1] == 20  # All cells processed

    def test_expression_to_penalty_optimized(self, simple_model, expression_df):
        """Test optimized penalty computation produces valid results."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, expression_df)
        penalties = scorer.compute_reaction_penalties()

        # Penalties should be in [0, 1] range
        assert penalties.min().min() >= 0
        assert penalties.max().max() <= 1

        # Should have no NaN or inf
        assert not penalties.isna().any().any()
        assert not np.any(np.isinf(penalties.values))

    def test_compass_optimizations_consistency(self, simple_model, expression_df):
        """Test that optimized and non-optimized produce similar results."""
        from cellmetpro.core.compass import CompassConfig, CompassScorer

        # With optimizations
        config_opt = CompassConfig(precompute_gpr=True, cache_max_fluxes=True)
        scorer_opt = CompassScorer(simple_model, expression_df, config_opt)
        penalties_opt = scorer_opt.compute_reaction_penalties()

        # Without optimizations
        config_basic = CompassConfig(precompute_gpr=False, cache_max_fluxes=False)
        scorer_basic = CompassScorer(simple_model, expression_df, config_basic)
        penalties_basic = scorer_basic.compute_reaction_penalties()

        # Results should be identical
        pd.testing.assert_frame_equal(penalties_opt, penalties_basic)


# ============================================================================
# Edge Case Tests - FluxBalanceAnalyzer
# ============================================================================


class TestFBAAdvanced:
    """Advanced tests for FluxBalanceAnalyzer methods."""

    def test_parsimonious_fba(self, simple_model):
        """Test parsimonious FBA (pFBA)."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.parsimonious_fba(fraction_of_optimum=1.0)

        assert fluxes is not None
        assert isinstance(fluxes, pd.Series)
        assert fba.solution is not None
        assert fba.solution.status == "optimal"

    def test_parsimonious_fba_partial_optimum(self, simple_model):
        """Test pFBA with fraction of optimum."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.parsimonious_fba(fraction_of_optimum=0.9)

        assert fluxes is not None
        assert isinstance(fluxes, pd.Series)

    def test_get_exchange_fluxes(self, simple_model):
        """Test getting exchange reaction fluxes."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        exchange_fluxes = fba.get_exchange_fluxes()

        assert exchange_fluxes is not None
        assert isinstance(exchange_fluxes, pd.Series)
        # Should include EX_A and EX_C
        assert len(exchange_fluxes) > 0

    def test_get_exchange_fluxes_after_optimize(self, simple_model):
        """Test exchange fluxes after explicit optimization."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fba.optimize()
        exchange_fluxes = fba.get_exchange_fluxes()

        assert exchange_fluxes is not None
        # Exchange reactions should have specific IDs
        exchange_ids = [rxn.id for rxn in fba.model.exchanges]
        assert set(exchange_fluxes.index) == set(exchange_ids)

    def test_gene_knockout(self, simple_model):
        """Test gene knockout analysis."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        # gene1 is used by R1 and R3
        fluxes = fba.gene_knockout("gene1")

        assert fluxes is not None
        assert isinstance(fluxes, pd.Series)

    def test_gene_knockout_multiple(self, simple_model):
        """Test knockout of multiple genes."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.gene_knockout(["gene1", "gene2"])

        assert fluxes is not None
        assert isinstance(fluxes, pd.Series)

    def test_optimize_with_custom_objective(self, simple_model):
        """Test optimization with custom objective."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.optimize(objective="R1", direction="max")

        assert fluxes is not None
        assert fba.solution is not None

    def test_optimize_minimize(self, simple_model):
        """Test minimization."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fluxes = fba.optimize(direction="min")

        assert fluxes is not None
        assert fba.solution is not None

    def test_summary_optimal(self, simple_model):
        """Test summary after optimization."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fba.optimize()
        summary = fba.summary()

        assert "FBA Solution Summary" in summary
        assert "optimal" in summary.lower()
        assert "Objective value" in summary
        assert "Top 10 reactions" in summary

    def test_reset_bounds_all(self, simple_model):
        """Test resetting all reaction bounds."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        # Store original bounds
        r1_orig = fba.model.reactions.get_by_id("R1").upper_bound
        r2_orig = fba.model.reactions.get_by_id("R2").upper_bound

        # Modify multiple reactions
        fba.set_bounds("R1", upper=50)
        fba.set_bounds("R2", upper=75)

        # Reset all
        fba.reset_bounds()

        # Verify all restored
        assert fba.model.reactions.get_by_id("R1").upper_bound == r1_orig
        assert fba.model.reactions.get_by_id("R2").upper_bound == r2_orig
        assert len(fba._original_bounds) == 0

    def test_flux_variability_all_reactions(self, simple_model):
        """Test FVA on all reactions (reactions=None)."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fva = fba.flux_variability(reactions=None, fraction_of_optimum=0.9)

        assert fva is not None
        assert "minimum" in fva.columns
        assert "maximum" in fva.columns
        # Should have entries for all reactions
        assert len(fva) == len(simple_model.reactions)


class TestFBAStandaloneFunctions:
    """Tests for standalone FBA utility functions."""

    def test_compute_yield(self, simple_model):
        """Test theoretical yield computation."""
        from cellmetpro.core.fba import compute_yield

        # EX_C is product, EX_A is substrate
        yield_value = compute_yield(
            simple_model,
            product_reaction="EX_C",
            substrate_reaction="EX_A",
            substrate_uptake=10.0,
        )

        assert yield_value is not None
        assert isinstance(yield_value, float)
        assert yield_value >= 0

    def test_compute_yield_zero_uptake_error(self, simple_model):
        """Test yield computation raises error for zero uptake."""
        from cellmetpro.core.fba import compute_yield

        with pytest.raises(ValueError, match="must be positive"):
            compute_yield(
                simple_model,
                product_reaction="EX_C",
                substrate_reaction="EX_A",
                substrate_uptake=0.0,
            )

    def test_compute_yield_negative_uptake_error(self, simple_model):
        """Test yield computation raises error for negative uptake."""
        from cellmetpro.core.fba import compute_yield

        with pytest.raises(ValueError, match="must be positive"):
            compute_yield(
                simple_model,
                product_reaction="EX_C",
                substrate_reaction="EX_A",
                substrate_uptake=-5.0,
            )

    def test_find_blocked_reactions(self, simple_model):
        """Test finding blocked reactions."""
        from cellmetpro.core.fba import find_blocked_reactions

        blocked = find_blocked_reactions(simple_model)

        assert isinstance(blocked, list)
        # In our simple model, all reactions can carry flux

    def test_find_essential_reactions(self, simple_model):
        """Test finding essential reactions."""
        from cellmetpro.core.fba import find_essential_reactions

        essential = find_essential_reactions(simple_model, threshold=0.01)

        assert isinstance(essential, list)
        # Essential reactions are those whose knockout kills growth

    def test_apply_media(self, simple_model):
        """Test applying media constraints."""
        from cellmetpro.core.fba import apply_media

        media = {"EX_A": 5.0}  # Allow uptake of 5 units of A
        model_with_media = apply_media(simple_model, media, default_uptake=0.0)

        assert model_with_media is not None
        # Check that EX_A has the correct bound
        ex_a = model_with_media.reactions.get_by_id("EX_A")
        assert ex_a.lower_bound == -5.0

    def test_apply_media_default_uptake(self, simple_model):
        """Test applying media with default uptake."""
        from cellmetpro.core.fba import apply_media

        media = {"EX_A": 10.0}
        model_with_media = apply_media(simple_model, media, default_uptake=1.0)

        assert model_with_media is not None
        # Original model is not modified
        assert simple_model.reactions.get_by_id("EX_A").lower_bound == -10


class TestFBAEdgeCases:
    """Edge case tests for FluxBalanceAnalyzer."""

    def test_fba_infeasible_model(self, simple_model):
        """Test FBA with infeasible constraints."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        # Set impossible bounds (lower > upper equivalent effect)
        fba.set_bounds("R1", lower=0, upper=0)
        fba.set_bounds("R2", lower=0, upper=0)
        fba.set_bounds("R3", lower=0, upper=0)

        # Should handle infeasible gracefully
        fluxes = fba.optimize()
        # Depending on implementation, may return None or zeros
        if fluxes is not None:
            assert isinstance(fluxes, pd.Series)

    def test_fba_zero_objective(self, simple_model):
        """Test FBA when objective is already zero."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        # Block all paths to objective
        fba.set_bounds("EX_A", lower=0, upper=0)

        fba.optimize()
        assert fba.solution is not None

    def test_fba_knockout_nonexistent_reaction(self, simple_model):
        """Test knockout with non-existent reaction ID."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)

        with pytest.raises((ValueError, KeyError)):
            fba.knockout("NONEXISTENT_REACTION")

    def test_fba_flux_variability_empty_reactions(self, simple_model):
        """Test FVA with empty reaction list."""
        from cellmetpro.core.fba import FluxBalanceAnalyzer

        fba = FluxBalanceAnalyzer(simple_model)
        fva = fba.flux_variability(reactions=[])

        assert len(fva) == 0


# ============================================================================
# Edge Case Tests - Microclustering
# ============================================================================


class TestMicroclusteringEdgeCases:
    """Edge case tests for microclustering."""

    def test_microcluster_few_cells(self):
        """Test microclustering with few cells."""
        from cellmetpro.core.microclustering import MicroclusterConfig, microcluster

        # Need enough cells for clustering - at least a few per cluster
        expression_df = pd.DataFrame(
            np.random.rand(5, 10),
            index=["GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            columns=[f"cell{i}" for i in range(10)],
        )

        config = MicroclusterConfig(cells_per_cluster=3, method="kmeans")
        result = microcluster(expression_df.T, config)

        assert result is not None
        assert result.n_clusters >= 1

    def test_microcluster_all_identical_values(self):
        """Test microclustering when all expression values are identical."""
        from cellmetpro.core.microclustering import MicroclusterConfig, microcluster

        # All values are the same - no variance
        expression_df = pd.DataFrame(
            np.ones((5, 10)),
            index=[f"GENE{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )

        config = MicroclusterConfig(cells_per_cluster=3, method="kmeans")
        result = microcluster(expression_df.T, config)

        assert result is not None

    def test_filter_genes_fano_zero_variance(self):
        """Test Fano filtering with zero-variance genes."""
        from cellmetpro.core.microclustering import filter_genes_fano

        # Create data where some genes have zero variance
        data = np.array(
            [
                [1, 1, 1, 1, 1],  # Zero variance
                [1, 2, 3, 4, 5],  # Has variance
                [5, 5, 5, 5, 5],  # Zero variance
            ]
        )
        expression_df = pd.DataFrame(
            data,
            index=["GENE1", "GENE2", "GENE3"],
            columns=[f"cell{i}" for i in range(5)],
        )

        filtered = filter_genes_fano(expression_df, n_genes=2)
        # Should handle zero variance gracefully
        assert len(filtered) <= 2


# ============================================================================
# Edge Case Tests - Cache
# ============================================================================


class TestCacheEdgeCases:
    """Edge case tests for caching functionality."""

    def test_memory_cache_size_one(self):
        """Test memory cache with size of 1."""
        from cellmetpro.core.cache import MemoryCache

        cache = MemoryCache(max_size=1)

        cache.set("key1", "value1")
        assert cache.get("key1") == "value1"

        # Adding second should evict first
        cache.set("key2", "value2")
        assert cache.get("key1") is None
        assert cache.get("key2") == "value2"

    def test_memory_cache_get_nonexistent(self):
        """Test getting non-existent key returns None."""
        from cellmetpro.core.cache import MemoryCache

        cache = MemoryCache(max_size=10)
        assert cache.get("nonexistent") is None

    def test_compass_cache_load_nonexistent_sample(self, tmp_path):
        """Test loading non-existent sample returns None."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        loaded = cache.load_reaction_scores("nonexistent_sample")

        assert loaded is None

    def test_compass_cache_clear_nonexistent(self, tmp_path):
        """Test clearing non-existent sample doesn't raise error."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        # Should not raise
        cache.clear("nonexistent_sample")

    def test_load_max_fluxes_cache_miss(self, tmp_path):
        """load_max_fluxes returns None when nothing has been saved."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="fresh")
        result = cache.load_max_fluxes()
        assert result is None

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_save_load_penalties(self, tmp_path):
        """Penalties can be saved and reloaded from disk."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        penalties = pd.DataFrame(
            {"cell1": [0.1, 0.9], "cell2": [0.5, 0.3]}, index=["R1", "R2"]
        )
        cache.save_penalties(penalties, "sample_pen")
        loaded = cache.load_penalties("sample_pen")

        assert loaded is not None
        assert loaded.shape == (2, 2)

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_load_penalties_cache_miss(self, tmp_path):
        """load_penalties returns None when nothing has been saved."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        result = cache.load_penalties("nonexistent")
        assert result is None

    def test_save_load_json(self, tmp_path):
        """Arbitrary JSON data can be saved and reloaded."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        data = {"beta": 0.95, "solver": "glpk", "n_reactions": 100}
        cache.save_json(data, "config")

        loaded = cache.load_json("config")
        assert loaded is not None
        assert loaded["beta"] == 0.95
        assert loaded["solver"] == "glpk"

    def test_load_json_cache_miss(self, tmp_path):
        """load_json returns None when the key has not been saved."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        result = cache.load_json("nonexistent_key")
        assert result is None

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_get_cached_samples(self, tmp_path):
        """get_cached_samples returns IDs of all samples with saved scores."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        assert cache.get_cached_samples() == []

        scores = pd.DataFrame({"cell1": [1.0]}, index=["R1"])
        cache.save_reaction_scores(scores, "sample_A")
        cache.save_reaction_scores(scores, "sample_B")

        cached = cache.get_cached_samples()
        assert set(cached) == {"sample_A", "sample_B"}

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_clear_all_cache(self, tmp_path):
        """clear() with no argument wipes all cached data for the model."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        scores = pd.DataFrame({"cell1": [1.0]}, index=["R1"])
        cache.save_reaction_scores(scores, "sample_X")

        assert cache.has_sample("sample_X")
        cache.clear()
        assert not cache.has_sample("sample_X")

    def test_get_cache_size(self, tmp_path):
        """get_cache_size returns a non-negative integer (bytes)."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")
        cache.save_max_fluxes({"R1": 100.0})

        size = cache.get_cache_size()
        assert isinstance(size, int)
        assert size > 0

    def test_get_cache_info(self, tmp_path):
        """get_cache_info returns a dict with the expected keys."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="mymodel")
        cache.save_max_fluxes({"R1": 1.0})

        info = cache.get_cache_info()
        assert info["model_id"] == "mymodel"
        assert "n_cached_samples" in info
        assert "total_size_bytes" in info
        assert "has_max_fluxes" in info
        assert info["has_max_fluxes"] is True

    def test_memory_cache_update_existing_key(self):
        """set() with an existing key updates the value without growing the cache."""
        from cellmetpro.core.cache import MemoryCache

        cache = MemoryCache(max_size=5)
        cache.set("key1", "first")
        cache.set("key2", "second")
        assert len(cache) == 2

        # Update key1 — should not grow the cache
        cache.set("key1", "updated")
        assert len(cache) == 2
        assert cache.get("key1") == "updated"

    def test_memory_cache_clear(self):
        """clear() empties the cache completely."""
        from cellmetpro.core.cache import MemoryCache

        cache = MemoryCache(max_size=5)
        cache.set("a", 1)
        cache.set("b", 2)
        assert len(cache) == 2

        cache.clear()
        assert len(cache) == 0
        assert cache.get("a") is None


# ============================================================================
# Additional Preprocessing Tests
# ============================================================================


class TestPreprocessingCoverage:
    """Additional coverage for preprocessing edge cases."""

    def _make_adata(self, n_cells: int = 10, n_genes: int = 5):
        """Create a simple dense AnnData for testing."""
        import anndata as ad

        rng = np.random.default_rng(42)
        X = rng.integers(0, 10, (n_cells, n_genes)).astype(float)
        return ad.AnnData(X)

    def test_filter_cells_max_genes(self):
        """filter_cells with max_genes removes cells expressing too many genes."""
        from cellmetpro.core.preprocessing import filter_cells

        adata = self._make_adata(n_cells=10, n_genes=5)
        # Give one cell expression in all 5 genes
        adata.X[0] = 1.0

        # Allow at most 3 expressed genes; cell_0 (5 genes) should be removed
        filtered = filter_cells(adata, min_genes=0, max_genes=4)
        assert filtered.n_obs < adata.n_obs

    def test_filter_cells_min_counts(self):
        """filter_cells with min_counts removes low-count cells."""
        from cellmetpro.core.preprocessing import filter_cells

        adata = self._make_adata(n_cells=6, n_genes=4)
        adata.X[:] = 0.0
        adata.X[0] = 100.0  # Only cell_0 has counts

        filtered = filter_cells(adata, min_genes=0, min_counts=50)
        assert filtered.n_obs == 1

    def test_filter_cells_max_counts(self):
        """filter_cells with max_counts removes high-count cells."""
        from cellmetpro.core.preprocessing import filter_cells

        adata = self._make_adata(n_cells=6, n_genes=4)
        adata.X[:] = 1.0
        adata.X[0] = 1000.0  # cell_0 has very high counts

        filtered = filter_cells(adata, min_genes=0, max_counts=100)
        assert filtered.n_obs == adata.n_obs - 1

    def test_filter_genes_min_counts(self):
        """filter_genes with min_counts removes lowly expressed genes."""
        from cellmetpro.core.preprocessing import filter_genes

        adata = self._make_adata(n_cells=5, n_genes=4)
        adata.X[:] = 0.0
        adata.X[:, 0] = 10.0  # Only gene_0 has counts

        filtered = filter_genes(adata, min_cells=0, min_counts=5)
        assert filtered.n_vars == 1

    def test_normalize_expression_sparse_input(self):
        """normalize_expression handles sparse AnnData correctly."""
        import anndata as ad
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import normalize_expression

        X = csr_matrix(np.array([[10.0, 0.0, 5.0], [0.0, 20.0, 0.0]]))
        adata = ad.AnnData(X)

        normalized = normalize_expression(adata, target_sum=1e4, log_transform=False)
        assert not np.any(np.isnan(normalized.X))

    def test_load_mtx_without_gene_barcode_files(self, tmp_path):
        """load_mtx generates generic gene/cell names when support files are absent."""
        from scipy.io import mmwrite
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import DataLoader

        # Write a 3-genes x 4-cells matrix (DataLoader will transpose to 4 x 3)
        matrix = csr_matrix(
            np.array([[1, 2, 0, 1], [0, 0, 3, 1], [2, 1, 0, 0]], dtype=float)
        )
        mtx_path = tmp_path / "matrix.mtx"
        mmwrite(str(mtx_path), matrix)

        loader = DataLoader(mtx_path)
        adata = loader.load()

        assert adata.n_obs == 4  # cells
        assert adata.n_vars == 3  # genes
        assert adata.obs_names[0].startswith("cell_")
        assert adata.var_names[0].startswith("gene_")

    def test_load_mtx_with_gene_and_barcode_files(self, tmp_path):
        """load_mtx reads gene and barcode names from TSV files."""
        import pandas as pd
        from scipy.io import mmwrite
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import DataLoader

        matrix = csr_matrix(
            np.array([[1, 2, 3], [4, 5, 6]], dtype=float)
        )  # 2 genes x 3 cells
        mmwrite(str(tmp_path / "matrix.mtx"), matrix)

        # genes.tsv: 2-column format (ID, name)
        pd.DataFrame([["ENSG001", "GAPDH"], ["ENSG002", "PKM"]]).to_csv(
            tmp_path / "genes.tsv", sep="\t", header=False, index=False
        )
        # barcodes.tsv
        pd.DataFrame(["ACGT-1", "TGCA-1", "GCTA-1"]).to_csv(
            tmp_path / "barcodes.tsv", sep="\t", header=False, index=False
        )

        loader = DataLoader(tmp_path / "matrix.mtx")
        adata = loader.load()

        assert adata.n_obs == 3
        assert adata.n_vars == 2
        assert "GAPDH" in adata.var_names
        assert "ACGT-1" in adata.obs_names

    def test_load_mtx_single_column_genes_file(self, tmp_path):
        """load_mtx reads gene names when genes.tsv has only one column."""
        import pandas as pd
        from scipy.io import mmwrite
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import DataLoader

        matrix = csr_matrix(
            np.array([[1, 2], [3, 4]], dtype=float)
        )  # 2 genes x 2 cells
        mmwrite(str(tmp_path / "matrix.mtx"), matrix)

        pd.DataFrame(["GENE_A", "GENE_B"]).to_csv(
            tmp_path / "genes.tsv", sep="\t", header=False, index=False
        )

        loader = DataLoader(tmp_path / "matrix.mtx")
        adata = loader.load()

        assert "GENE_A" in adata.var_names

    def test_normalize_expression_copy_false(self):
        """normalize_expression with copy=False modifies adata in-place."""
        import anndata as ad

        from cellmetpro.core.preprocessing import normalize_expression

        X = np.array([[10.0, 0.0, 5.0], [0.0, 20.0, 0.0]])
        adata = ad.AnnData(X)
        original_id = id(adata)

        result = normalize_expression(adata, copy=False)
        assert id(result) == original_id

    def test_filter_cells_sparse_input(self):
        """filter_cells handles sparse AnnData correctly."""
        import anndata as ad
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import filter_cells

        X = csr_matrix(
            np.array(
                [
                    [1, 1, 0, 0],  # 2 genes expressed
                    [1, 1, 1, 0],  # 3 genes expressed
                    [0, 0, 0, 0],  # 0 genes expressed
                ],
                dtype=float,
            )
        )
        adata = ad.AnnData(X)

        filtered = filter_cells(adata, min_genes=1)
        assert filtered.n_obs == 2

    def test_filter_cells_copy_false(self):
        """filter_cells with copy=False works without duplicating the object."""
        import anndata as ad

        from cellmetpro.core.preprocessing import filter_cells

        X = np.ones((5, 3))
        adata = ad.AnnData(X)

        result = filter_cells(adata, min_genes=0, copy=False)
        assert result.n_obs == 5

    def test_filter_genes_sparse_input(self):
        """filter_genes handles sparse AnnData correctly."""
        import anndata as ad
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import filter_genes

        X = csr_matrix(
            np.array(
                [
                    [1, 0, 1],  # gene0 and gene2 expressed in cell0
                    [1, 0, 0],  # only gene0 expressed in cell1
                ],
                dtype=float,
            )
        )
        adata = ad.AnnData(X)

        filtered = filter_genes(adata, min_cells=2)
        assert filtered.n_vars == 1  # only gene0 expressed in both cells

    def test_filter_genes_copy_false(self):
        """filter_genes with copy=False works without duplicating the object."""
        import anndata as ad

        from cellmetpro.core.preprocessing import filter_genes

        X = np.ones((4, 3))
        adata = ad.AnnData(X)

        result = filter_genes(adata, min_cells=0, copy=False)
        assert result.n_vars == 3

    def test_to_dataframe_with_layer(self):
        """to_dataframe reads from a named layer when layer is specified."""
        import anndata as ad

        from cellmetpro.core.preprocessing import to_dataframe

        X = np.zeros((3, 4))
        adata = ad.AnnData(X)
        layer_data = np.ones((3, 4)) * 5.0
        adata.layers["normalized"] = layer_data

        df = to_dataframe(adata, layer="normalized", genes_as_rows=False)
        assert df.shape == (3, 4)
        assert (df.values == 5.0).all()

    def test_to_dataframe_genes_as_rows_false(self):
        """to_dataframe with genes_as_rows=False returns cells x genes."""
        import anndata as ad

        from cellmetpro.core.preprocessing import to_dataframe

        rng = np.random.default_rng(0)
        X = rng.random((5, 10))
        adata = ad.AnnData(X)

        df = to_dataframe(adata, genes_as_rows=False)
        assert df.shape == (5, 10)  # cells x genes

    def test_to_dataframe_sparse_x_converted(self):
        """to_dataframe converts sparse X to dense via toarray()."""
        import anndata as ad
        from scipy.sparse import csr_matrix

        from cellmetpro.core.preprocessing import to_dataframe

        X = csr_matrix(np.array([[1, 0, 3], [0, 5, 0]], dtype=float))
        adata = ad.AnnData(X)

        df = to_dataframe(adata, genes_as_rows=False)
        assert df.shape == (2, 3)
        assert df.iloc[0, 0] == 1.0
        assert df.iloc[0, 1] == 0.0

    def test_load_tsv_file(self, tmp_path):
        """DataLoader handles .tsv files via load_csv(sep='\\t')."""
        import anndata as ad
        import pandas as pd

        from cellmetpro.core.preprocessing import DataLoader

        tsv_path = tmp_path / "expr.tsv"
        df = pd.DataFrame(
            np.array([[1, 2], [3, 4]], dtype=float),
            index=["cell0", "cell1"],
            columns=["geneA", "geneB"],
        )
        df.to_csv(tsv_path, sep="\t")

        loader = DataLoader(str(tsv_path))
        adata = loader.load()
        assert isinstance(adata, ad.AnnData)


class TestGetOrComputeMaxFluxes:
    """Tests for get_or_compute_max_fluxes (requires cobra)."""

    @pytest.fixture
    def tiny_model(self):
        """Minimal cobra model with two non-boundary reactions."""
        import cobra

        model = cobra.Model("tiny")
        A = cobra.Metabolite("A_c", compartment="c")
        B = cobra.Metabolite("B_c", compartment="c")

        r1 = cobra.Reaction("R1")
        r1.add_metabolites({A: -1, B: 1})
        r1.bounds = (0, 10)

        ex_A = cobra.Reaction("EX_A")
        ex_A.add_metabolites({A: -1})
        ex_A.bounds = (-10, 1000)

        ex_B = cobra.Reaction("EX_B")
        ex_B.add_metabolites({B: -1})
        ex_B.bounds = (0, 1000)

        model.add_reactions([r1, ex_A, ex_B])
        model.objective = "EX_B"
        return model

    def test_compute_no_cache(self, tiny_model):
        """Computes max fluxes when cache=None."""
        try:
            from cellmetpro.core.cache import get_or_compute_max_fluxes
        except ImportError:
            pytest.skip("cobra not installed")

        result = get_or_compute_max_fluxes(tiny_model, cache=None)

        assert isinstance(result, dict)
        assert len(result) > 0

    def test_compute_and_save_to_cache(self, tiny_model, tmp_path):
        """Computes when cache has no max_fluxes, then saves to cache."""
        try:
            from cellmetpro.core.cache import CompassCache, get_or_compute_max_fluxes
        except ImportError:
            pytest.skip("cobra not installed")

        cache = CompassCache(cache_dir=tmp_path)
        result = get_or_compute_max_fluxes(tiny_model, cache=cache)

        assert isinstance(result, dict)
        cached = cache.load_max_fluxes()
        assert cached is not None
        assert set(cached.keys()) == set(result.keys())

    def test_returns_cached_value(self, tiny_model, tmp_path):
        """Returns pre-existing cached max fluxes without recomputing."""
        try:
            from cellmetpro.core.cache import CompassCache, get_or_compute_max_fluxes
        except ImportError:
            pytest.skip("cobra not installed")

        cache = CompassCache(cache_dir=tmp_path)
        precomputed = {"R1": 99.0, "EX_A": 42.0}
        cache.save_max_fluxes(precomputed)

        result = get_or_compute_max_fluxes(tiny_model, cache=cache)

        assert result == precomputed

    def test_from_model_creates_deterministic_id(self, tiny_model, tmp_path):
        """CompassCache.from_model derives a stable ID from the model."""
        try:
            from cellmetpro.core.cache import CompassCache
        except ImportError:
            pytest.skip("cobra not installed")

        cache1 = CompassCache.from_model(tiny_model, cache_dir=tmp_path)
        cache2 = CompassCache.from_model(tiny_model, cache_dir=tmp_path)

        assert isinstance(cache1, CompassCache)
        assert "tiny" in cache1.model_id
        assert cache1.model_id == cache2.model_id
