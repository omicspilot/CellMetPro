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
        from cellmetpro.core.compass import CompassConfig, run_compass

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
        expr = np.array([
            [1.0, 2.0, 3.0, 4.0],  # GENE1
            [5.0, 6.0, 7.0, 8.0],  # GENE2
            [9.0, 10.0, 11.0, 12.0],  # GENE3
        ])
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
            rxn.id for rxn in model.reactions
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

        fluxes = fba.optimize()
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
        data = np.array([
            [1, 1, 1, 1, 1],  # Zero variance
            [1, 2, 3, 4, 5],  # Has variance
            [5, 5, 5, 5, 5],  # Zero variance
        ])
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
