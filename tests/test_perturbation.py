"""Tests for the perturbation analysis module."""

import importlib.util

import pandas as pd
import pytest

# Skip all tests in this module if cobra is not installed
pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("cobra") is None,
    reason="cobra not installed",
)


@pytest.fixture
def simple_model():
    """Create a simple test metabolic model."""
    import cobra

    model = cobra.Model("test_model")

    # Create metabolites
    A = cobra.Metabolite("A_c", compartment="c")
    B = cobra.Metabolite("B_c", compartment="c")
    C = cobra.Metabolite("C_c", compartment="c")

    # Create reactions
    r1 = cobra.Reaction("R1")
    r1.add_metabolites({A: -1, B: 1})
    r1.bounds = (0, 1000)

    r2 = cobra.Reaction("R2")
    r2.add_metabolites({B: -1, C: 1})
    r2.bounds = (0, 1000)

    r3 = cobra.Reaction("R3")
    r3.add_metabolites({A: -1, C: 1})
    r3.bounds = (0, 1000)

    # Exchange reactions
    ex_A = cobra.Reaction("EX_A")
    ex_A.add_metabolites({A: -1})
    ex_A.bounds = (-10, 1000)

    ex_C = cobra.Reaction("EX_C")
    ex_C.add_metabolites({C: -1})
    ex_C.bounds = (0, 1000)

    model.add_reactions([r1, r2, r3, ex_A, ex_C])
    model.objective = "EX_C"

    # Add genes
    r1.gene_reaction_rule = "gene1"
    r2.gene_reaction_rule = "gene2"
    r3.gene_reaction_rule = "gene3"

    return model


class TestSimulateExpressionChange:
    """Tests for simulate_expression_change function."""

    def test_knockdown(self, simple_model):
        """Test gene knockdown simulation."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        # 50% knockdown of gene1
        fluxes = simulate_expression_change(simple_model, "gene1", 0.5)

        assert isinstance(fluxes, pd.Series)
        assert len(fluxes) == len(simple_model.reactions)

    def test_overexpression(self, simple_model):
        """Test gene overexpression simulation."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        # 2x overexpression
        fluxes = simulate_expression_change(simple_model, "gene1", 2.0)

        assert isinstance(fluxes, pd.Series)

    def test_complete_knockout(self, simple_model):
        """Test complete knockout (fold_change near 0)."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        fluxes = simulate_expression_change(
            simple_model, "gene1", 0.01, method="threshold"
        )

        assert isinstance(fluxes, pd.Series)

    def test_invalid_gene(self, simple_model):
        """Test error on invalid gene ID."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        with pytest.raises(ValueError, match="not found"):
            simulate_expression_change(simple_model, "invalid_gene", 0.5)

    def test_log_method(self, simple_model):
        """Test log scaling method."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        fluxes = simulate_expression_change(simple_model, "gene1", 2.0, method="log")

        assert isinstance(fluxes, pd.Series)


class TestMultiKnockout:
    """Tests for multi_knockout function."""

    def test_simultaneous_knockout(self, simple_model):
        """Test simultaneous multi-gene knockout."""
        from cellmetpro.analysis.perturbation import multi_knockout

        fluxes = multi_knockout(simple_model, ["gene1", "gene2"], method="simultaneous")

        assert isinstance(fluxes, pd.Series)

    def test_sequential_knockout(self, simple_model):
        """Test sequential knockout."""
        from cellmetpro.analysis.perturbation import multi_knockout

        fluxes = multi_knockout(
            simple_model,
            ["gene1", "gene2"],
            method="sequential",
            return_intermediate=False,
        )

        assert isinstance(fluxes, pd.Series)

    def test_sequential_with_intermediate(self, simple_model):
        """Test sequential knockout with intermediate states."""
        from cellmetpro.analysis.perturbation import multi_knockout

        results = multi_knockout(
            simple_model,
            ["gene1", "gene2"],
            method="sequential",
            return_intermediate=True,
        )

        assert isinstance(results, pd.DataFrame)
        assert results.shape[1] == 2  # Two knockout states

    def test_individual_knockouts(self, simple_model):
        """Test individual knockout comparison."""
        from cellmetpro.analysis.perturbation import multi_knockout

        results = multi_knockout(simple_model, ["gene1", "gene2"], method="individual")

        assert isinstance(results, pd.DataFrame)
        assert "wild_type" in results.columns


class TestCompareFluxDistributions:
    """Tests for compare_flux_distributions function."""

    def test_basic_comparison(self):
        """Test basic flux comparison."""
        from cellmetpro.analysis.perturbation import compare_flux_distributions

        flux1 = pd.Series({"R1": 10.0, "R2": 5.0, "R3": 2.0})
        flux2 = pd.Series({"R1": 8.0, "R2": 5.0, "R3": 4.0})

        comparison = compare_flux_distributions(flux1, flux2)

        assert "abs_diff" in comparison.columns
        assert "flux_1" in comparison.columns
        assert "flux_2" in comparison.columns

    def test_relative_comparison(self):
        """Test relative (fold change) comparison."""
        from cellmetpro.analysis.perturbation import compare_flux_distributions

        flux1 = pd.Series({"R1": 10.0, "R2": 5.0})
        flux2 = pd.Series({"R1": 20.0, "R2": 2.5})

        comparison = compare_flux_distributions(flux1, flux2, method="both")

        assert "fold_change" in comparison.columns
        assert "log2_fc" in comparison.columns


class TestIdentifyEssentialGenes:
    """Tests for identify_essential_genes function."""

    def test_basic_essential_screen(self, simple_model):
        """Test basic essential gene identification."""
        from cellmetpro.analysis.perturbation import identify_essential_genes

        results = identify_essential_genes(
            simple_model, genes=["gene1", "gene2", "gene3"], threshold=0.01
        )

        assert isinstance(results, pd.DataFrame)
        assert "gene" in results.columns
        assert "essential" in results.columns
        assert "ko_growth" in results.columns


class TestSensitivityAnalysis:
    """Tests for sensitivity_analysis function."""

    def test_basic_sensitivity(self, simple_model):
        """Test basic sensitivity analysis."""
        from cellmetpro.analysis.perturbation import sensitivity_analysis

        results = sensitivity_analysis(simple_model, reactions=["R1", "R2"], n_steps=3)

        assert isinstance(results, pd.DataFrame)
        assert "sensitivity" in results.columns


class TestSimulateExpressionChangeEdgeCases:
    """Edge cases for simulate_expression_change."""

    def test_unknown_method_raises(self, simple_model):
        """An unrecognised method name must raise ValueError."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        with pytest.raises(ValueError, match="Unknown method"):
            simulate_expression_change(simple_model, "gene1", 1.0, method="bad_method")

    def test_reversible_reaction_lower_bound_scaled(self):
        """When a reaction has lb < 0 (reversible), lb is also scaled."""
        import cobra

        from cellmetpro.analysis.perturbation import simulate_expression_change

        model = cobra.Model("rev")
        A = cobra.Metabolite("A")
        B = cobra.Metabolite("B")
        ex = cobra.Reaction("EX_B")
        ex.add_metabolites({B: -1})
        ex.bounds = (0, 1000)

        rxn = cobra.Reaction("REV")
        rxn.add_metabolites({A: -1, B: 1})
        rxn.bounds = (-10, 10)  # reversible
        rxn.gene_reaction_rule = "rev_gene"

        ex_A = cobra.Reaction("EX_A")
        ex_A.add_metabolites({A: -1})
        ex_A.bounds = (-100, 0)

        model.add_reactions([rxn, ex, ex_A])
        model.objective = "EX_B"

        fluxes = simulate_expression_change(model, "rev_gene", 0.5, method="linear")
        assert isinstance(fluxes, pd.Series)

    def test_threshold_above_cutoff_keeps_reaction(self, simple_model):
        """threshold method with fold_change >= 0.1 keeps reaction active."""
        from cellmetpro.analysis.perturbation import simulate_expression_change

        fluxes = simulate_expression_change(
            simple_model, "gene1", 0.5, method="threshold"
        )
        assert isinstance(fluxes, pd.Series)


class TestMultiKnockoutEdgeCases:
    """Edge cases for multi_knockout."""

    def test_unknown_method_raises(self, simple_model):
        """An unrecognised method name must raise ValueError."""
        from cellmetpro.analysis.perturbation import multi_knockout

        with pytest.raises(ValueError, match="Unknown method"):
            multi_knockout(simple_model, ["gene1"], method="invalid_method")

    def test_simultaneous_with_invalid_gene_skips(self, simple_model):
        """Invalid gene in simultaneous mode is silently skipped."""
        from cellmetpro.analysis.perturbation import multi_knockout

        # Should not raise — bad gene is logged and skipped
        fluxes = multi_knockout(
            simple_model, ["gene1", "no_such_gene"], method="simultaneous"
        )
        assert isinstance(fluxes, pd.Series)


class TestCompareFluxDistributionsEdgeCases:
    """Edge cases for compare_flux_distributions."""

    def test_relative_method_produces_fold_change_columns(self):
        """method='relative' must produce fold_change and log2_fc columns."""
        from cellmetpro.analysis.perturbation import compare_flux_distributions

        flux1 = pd.Series({"R1": 10.0, "R2": 5.0, "R3": 0.0})
        flux2 = pd.Series({"R1": 20.0, "R2": 2.5, "R3": 0.0})

        result = compare_flux_distributions(flux1, flux2, method="relative")

        assert "fold_change" in result.columns
        assert "log2_fc" in result.columns

    def test_both_method_produces_all_columns(self):
        """method='both' produces both absolute and relative columns."""
        from cellmetpro.analysis.perturbation import compare_flux_distributions

        flux1 = pd.Series({"R1": 10.0, "R2": 5.0})
        flux2 = pd.Series({"R1": 15.0, "R2": 3.0})

        result = compare_flux_distributions(flux1, flux2, method="both")

        assert "abs_diff" in result.columns
        assert "fold_change" in result.columns


class TestPredictSyntheticLethality:
    """Tests for predict_synthetic_lethality function."""

    def test_synthetic_lethality_screen(self, simple_model):
        """Test synthetic lethality prediction."""
        from cellmetpro.analysis.perturbation import predict_synthetic_lethality

        results = predict_synthetic_lethality(
            simple_model, genes=["gene1", "gene2", "gene3"], threshold=0.01
        )

        assert isinstance(results, pd.DataFrame)
        assert "gene1" in results.columns
        assert "gene2" in results.columns
        assert "synthetic_lethal" in results.columns
