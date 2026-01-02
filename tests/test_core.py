"""Tests for core module."""

import numpy as np
import pandas as pd
import pytest


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

    def test_custom_config(self):
        """Test custom configuration."""
        from cellmetpro.core.compass import CompassConfig

        config = CompassConfig(
            beta=0.9,
            lambda_penalty=0.5,
            n_processes=4,
        )

        assert config.beta == 0.9
        assert config.lambda_penalty == 0.5
        assert config.n_processes == 4


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
        original_ub = fba.model.reactions.get_by_id("R1").upper_bound

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
        from cellmetpro.core.microclustering import microcluster, MicroclusterConfig

        config = MicroclusterConfig(cells_per_cluster=2, method="kmeans")
        result = microcluster(expression_df, config)

        assert result is not None
        assert result.n_clusters > 0
        assert len(result.cluster_labels) == 5
        assert result.pooled_expression.shape[1] == result.n_clusters

    def test_unpool_results(self, expression_df):
        """Test unpooling cluster results."""
        from cellmetpro.core.microclustering import (
            microcluster,
            unpool_results,
            MicroclusterConfig,
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

    def test_has_sample(self, tmp_path):
        """Test checking for cached samples."""
        from cellmetpro.core.cache import CompassCache

        cache = CompassCache(cache_dir=tmp_path, model_id="test")

        assert not cache.has_sample("sample1")

        scores = pd.DataFrame({"cell1": [1.0]}, index=["R1"])
        cache.save_reaction_scores(scores, "sample1")

        assert cache.has_sample("sample1")

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


class TestIntegration:
    """Integration tests combining multiple modules."""

    def test_full_compass_pipeline(self, simple_model, expression_df):
        """Test full COMPASS pipeline."""
        from cellmetpro.core.compass import run_compass, CompassConfig

        config = CompassConfig(beta=0.95, n_processes=1)
        result = run_compass(simple_model, expression_df, config)

        assert result is not None
        assert result.reaction_penalties is not None
        assert result.reaction_scores is not None
        assert result.config.beta == 0.95

    def test_compass_with_adata(self, simple_model, adata):
        """Test COMPASS with AnnData input."""
        from cellmetpro.core.compass import CompassScorer

        scorer = CompassScorer(simple_model, adata)
        penalties = scorer.compute_reaction_penalties()

        assert penalties is not None
        assert penalties.shape[1] == adata.n_obs
