"""Tests for pathway analysis module."""

import cobra
import numpy as np
import pandas as pd
import pytest

from cellmetpro.analysis.pathway import (
    GOEnrichmentAnalyzer,
    PathwayAnalyzer,
    create_go_annotations_from_dict,
)


@pytest.fixture
def model_with_subsystems():
    """Create a model with subsystem annotations."""
    model = cobra.Model("test_model")

    # Metabolites
    A = cobra.Metabolite("A", compartment="c")
    B = cobra.Metabolite("B", compartment="c")
    C = cobra.Metabolite("C", compartment="c")
    D = cobra.Metabolite("D", compartment="c")

    # Glycolysis reactions
    r1 = cobra.Reaction("R1")
    r1.add_metabolites({A: -1, B: 1})
    r1.gene_reaction_rule = "gene1"
    r1.subsystem = "Glycolysis"

    r2 = cobra.Reaction("R2")
    r2.add_metabolites({B: -1, C: 1})
    r2.gene_reaction_rule = "gene2"
    r2.subsystem = "Glycolysis"

    # TCA cycle reactions
    r3 = cobra.Reaction("R3")
    r3.add_metabolites({C: -1, D: 1})
    r3.gene_reaction_rule = "gene3"
    r3.subsystem = "TCA cycle"

    r4 = cobra.Reaction("R4")
    r4.add_metabolites({D: -1, A: 1})
    r4.gene_reaction_rule = "gene4"
    r4.subsystem = "TCA cycle"

    # No subsystem reaction
    r5 = cobra.Reaction("R5")
    r5.add_metabolites({A: -1, D: 1})
    r5.gene_reaction_rule = "gene5"
    # No subsystem set

    model.add_reactions([r1, r2, r3, r4, r5])
    return model


@pytest.fixture
def reaction_scores():
    """Create mock reaction scores."""
    np.random.seed(42)
    reactions = ["R1", "R2", "R3", "R4", "R5"]
    cells = [f"cell_{i}" for i in range(10)]
    data = np.random.rand(5, 10)
    return pd.DataFrame(data, index=reactions, columns=cells)


@pytest.fixture
def go_annotations():
    """Create mock GO annotations."""
    gene_to_go = {
        "gene1": [
            ("GO:0006096", "glycolytic process", "BP"),
            ("GO:0005737", "cytoplasm", "CC"),
        ],
        "gene2": [
            ("GO:0006096", "glycolytic process", "BP"),
        ],
        "gene3": [
            ("GO:0006099", "tricarboxylic acid cycle", "BP"),
        ],
        "gene4": [
            ("GO:0006099", "tricarboxylic acid cycle", "BP"),
            ("GO:0005739", "mitochondrion", "CC"),
        ],
        "gene5": [
            ("GO:0008150", "biological_process", "BP"),
        ],
    }
    return create_go_annotations_from_dict(gene_to_go)


# =============================================================================
# TESTS FOR PathwayAnalyzer
# =============================================================================


def test_pathway_analyzer_get_pathway_mapping(model_with_subsystems, reaction_scores):
    """Test that pathway mapping is correctly extracted from model."""
    analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)
    mapping = analyzer.get_pathway_mapping()

    assert "Glycolysis" in mapping
    assert "TCA cycle" in mapping
    assert set(mapping["Glycolysis"]) == {"R1", "R2"}
    assert set(mapping["TCA cycle"]) == {"R3", "R4"}


def test_pathway_analyzer_aggregate_mean(model_with_subsystems, reaction_scores):
    """Test pathway aggregation with mean method."""
    analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)
    pathway_scores = analyzer.aggregate(method="mean")

    assert isinstance(pathway_scores, pd.DataFrame)
    assert "Glycolysis" in pathway_scores.index
    assert "TCA cycle" in pathway_scores.index
    assert pathway_scores.shape[1] == 10  # 10 cells


def test_pathway_analyzer_aggregate_methods(model_with_subsystems, reaction_scores):
    """Test all aggregation methods work."""
    analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)

    for method in ["mean", "median", "sum", "max"]:
        result = analyzer.aggregate(method=method)
        assert len(result) > 0


def test_pathway_analyzer_enrich(model_with_subsystems, reaction_scores):
    """Test pathway enrichment analysis."""
    analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)

    # Create mock differential results
    diff_results = pd.DataFrame({
        "reaction": ["R1", "R2", "R3", "R4", "R5"],
        "pvalue": [0.01, 0.02, 0.5, 0.6, 0.7],
        "padj_bh": [0.02, 0.04, 0.5, 0.6, 0.7],
    })

    enrichment = analyzer.enrich(diff_results)

    assert isinstance(enrichment, pd.DataFrame)
    if len(enrichment) > 0:
        assert "pathway" in enrichment.columns
        assert "pvalue" in enrichment.columns
        assert "fold_enrichment" in enrichment.columns


# =============================================================================
# TESTS FOR GOEnrichmentAnalyzer
# =============================================================================


def test_go_enrichment_init(model_with_subsystems, go_annotations):
    """Test GOEnrichmentAnalyzer initialization."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    assert len(analyzer.gene_to_go) > 0
    assert len(analyzer.reaction_to_genes) > 0


def test_go_enrichment_get_reaction_go_terms(model_with_subsystems, go_annotations):
    """Test getting GO terms for a reaction."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    # R1 has gene1, which has GO:0006096 and GO:0005737
    r1_terms = analyzer.get_reaction_go_terms("R1")
    assert "GO:0006096" in r1_terms
    assert "GO:0005737" in r1_terms


def test_go_enrichment_get_go_to_reactions_mapping(
    model_with_subsystems, go_annotations
):
    """Test GO to reactions mapping."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)
    mapping = analyzer.get_go_to_reactions_mapping()

    # GO:0006096 (glycolysis) should map to R1 and R2
    assert "GO:0006096" in mapping
    assert "R1" in mapping["GO:0006096"]
    assert "R2" in mapping["GO:0006096"]


def test_go_enrichment_enrich_reactions(model_with_subsystems, go_annotations):
    """Test GO enrichment analysis."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    significant = {"R1", "R2"}  # Both glycolysis reactions
    background = {"R1", "R2", "R3", "R4", "R5"}

    results = analyzer.enrich_reactions(
        significant,
        background,
        min_genes=1,  # Lower threshold for test
    )

    assert isinstance(results, pd.DataFrame)
    if len(results) > 0:
        assert "go_term" in results.columns
        assert "go_name" in results.columns
        assert "pvalue" in results.columns
        assert "fold_enrichment" in results.columns


def test_go_enrichment_filter_by_namespace(model_with_subsystems, go_annotations):
    """Test filtering GO enrichment by namespace."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    significant = {"R1", "R2", "R3", "R4"}
    background = {"R1", "R2", "R3", "R4", "R5"}

    # Filter to only biological process
    results_bp = analyzer.enrich_reactions(
        significant,
        background,
        min_genes=1,
        namespace="BP",
    )

    if len(results_bp) > 0:
        assert all(results_bp["namespace"] == "BP")


def test_go_enrichment_from_differential(model_with_subsystems, go_annotations):
    """Test GO enrichment from differential analysis results."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    diff_results = pd.DataFrame({
        "reaction": ["R1", "R2", "R3", "R4", "R5"],
        "pvalue": [0.01, 0.02, 0.5, 0.6, 0.7],
        "padj_bh": [0.02, 0.04, 0.5, 0.6, 0.7],
    })

    results = analyzer.enrich_from_differential(
        diff_results,
        pvalue_threshold=0.05,
        min_genes=1,
    )

    assert isinstance(results, pd.DataFrame)


def test_go_enrichment_fold_enrichment_calculation(
    model_with_subsystems, go_annotations
):
    """Test that fold enrichment is calculated correctly."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    # If all glycolysis reactions are significant, fold enrichment should be > 1
    significant = {"R1", "R2"}
    background = {"R1", "R2", "R3", "R4", "R5"}

    results = analyzer.enrich_reactions(significant, background, min_genes=1)

    # Find glycolysis GO term
    glycolysis_row = results[results["go_term"] == "GO:0006096"]
    if len(glycolysis_row) > 0:
        # All glycolysis reactions (R1, R2) are significant
        # So fold enrichment should be high
        assert glycolysis_row["fold_enrichment"].iloc[0] > 1.0


def test_go_enrichment_padj_correction(model_with_subsystems, go_annotations):
    """Test that p-value adjustment is applied."""
    analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

    significant = {"R1", "R2", "R3"}
    background = {"R1", "R2", "R3", "R4", "R5"}

    results = analyzer.enrich_reactions(significant, background, min_genes=1)

    if len(results) > 0:
        assert "padj" in results.columns
        # Adjusted p-values should be >= raw p-values
        assert (results["padj"] >= results["pvalue"] - 1e-10).all()


# =============================================================================
# TESTS FOR helper functions
# =============================================================================


def test_create_go_annotations_from_dict():
    """Test creating GO annotations from dictionary."""
    gene_to_go = {
        "gene1": [
            ("GO:0006096", "glycolytic process", "BP"),
        ],
        "gene2": [
            ("GO:0006099", "TCA cycle", "BP"),
        ],
    }

    df = create_go_annotations_from_dict(gene_to_go)

    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == {"gene_id", "go_term", "go_name", "namespace"}
    assert len(df) == 2
    assert "gene1" in df["gene_id"].values
    assert "GO:0006096" in df["go_term"].values


# =============================================================================
# EDGE CASE TESTS
# =============================================================================


class TestPathwayAnalyzerEdgeCases:
    """Edge case tests for PathwayAnalyzer."""

    def test_pathway_analyzer_no_subsystems(self, reaction_scores):
        """Test analyzer when model has no subsystems."""
        import cobra

        # Create model without subsystems
        model = cobra.Model("empty_subsystems")
        r1 = cobra.Reaction("R1")
        r1.add_metabolites({cobra.Metabolite("A", compartment="c"): -1})
        model.add_reactions([r1])

        analyzer = PathwayAnalyzer(reaction_scores, model)
        mapping = analyzer.get_pathway_mapping()

        # Should return empty mapping or handle gracefully
        assert isinstance(mapping, dict)

    def test_pathway_analyzer_single_reaction_pathway(self, reaction_scores):
        """Test with pathway containing single reaction."""
        import cobra

        model = cobra.Model("single_rxn_pathway")
        A = cobra.Metabolite("A", compartment="c")
        B = cobra.Metabolite("B", compartment="c")

        r1 = cobra.Reaction("R1")
        r1.add_metabolites({A: -1, B: 1})
        r1.subsystem = "SingleReactionPathway"

        model.add_reactions([r1])

        # Scores must include R1
        scores = pd.DataFrame(
            np.random.rand(1, 5),
            index=["R1"],
            columns=[f"cell_{i}" for i in range(5)],
        )

        analyzer = PathwayAnalyzer(scores, model)
        result = analyzer.aggregate(method="mean")

        assert "SingleReactionPathway" in result.index

    def test_pathway_analyzer_no_matching_reactions(self):
        """Test when reaction scores don't match model reactions."""
        import cobra

        model = cobra.Model("test")
        r1 = cobra.Reaction("MODEL_R1")
        r1.add_metabolites({cobra.Metabolite("A", compartment="c"): -1})
        r1.subsystem = "Pathway1"
        model.add_reactions([r1])

        # Scores have completely different reaction names
        scores = pd.DataFrame(
            np.random.rand(3, 5),
            index=["OTHER1", "OTHER2", "OTHER3"],
            columns=[f"cell_{i}" for i in range(5)],
        )

        analyzer = PathwayAnalyzer(scores, model)
        result = analyzer.aggregate(method="mean")

        # Should return empty or handle gracefully
        assert len(result) == 0

    def test_pathway_analyzer_empty_scores(self, model_with_subsystems):
        """Test with empty reaction scores DataFrame."""
        scores = pd.DataFrame(columns=["cell_0", "cell_1", "cell_2"])

        analyzer = PathwayAnalyzer(scores, model_with_subsystems)
        result = analyzer.aggregate(method="mean")

        assert len(result) == 0

    def test_pathway_analyzer_single_cell(self, model_with_subsystems):
        """Test with single cell."""
        scores = pd.DataFrame(
            np.random.rand(5, 1),
            index=["R1", "R2", "R3", "R4", "R5"],
            columns=["cell_0"],
        )

        analyzer = PathwayAnalyzer(scores, model_with_subsystems)
        result = analyzer.aggregate(method="mean")

        assert result.shape[1] == 1

    def test_pathway_analyzer_nan_values(self, model_with_subsystems):
        """Test aggregation with NaN values in scores."""
        scores = pd.DataFrame(
            np.array([
                [1.0, np.nan, 3.0],
                [np.nan, 2.0, 4.0],
                [5.0, 6.0, np.nan],
                [7.0, 8.0, 9.0],
                [10.0, 11.0, 12.0],
            ]),
            index=["R1", "R2", "R3", "R4", "R5"],
            columns=["cell_0", "cell_1", "cell_2"],
        )

        analyzer = PathwayAnalyzer(scores, model_with_subsystems)
        result = analyzer.aggregate(method="mean")

        # Should handle NaN gracefully (nanmean)
        assert isinstance(result, pd.DataFrame)

    def test_pathway_enrich_no_significant_reactions(self, model_with_subsystems, reaction_scores):
        """Test enrichment when no reactions are significant."""
        analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)

        # All p-values > threshold
        diff_results = pd.DataFrame({
            "reaction": ["R1", "R2", "R3", "R4", "R5"],
            "pvalue": [0.9, 0.8, 0.7, 0.6, 0.5],
            "padj_bh": [0.9, 0.8, 0.7, 0.6, 0.5],
        })

        enrichment = analyzer.enrich(diff_results)

        # Should return empty or results with high p-values
        assert isinstance(enrichment, pd.DataFrame)

    def test_pathway_enrich_all_significant(self, model_with_subsystems, reaction_scores):
        """Test enrichment when all reactions are significant."""
        analyzer = PathwayAnalyzer(reaction_scores, model_with_subsystems)

        diff_results = pd.DataFrame({
            "reaction": ["R1", "R2", "R3", "R4", "R5"],
            "pvalue": [0.001, 0.002, 0.003, 0.004, 0.005],
            "padj_bh": [0.005, 0.006, 0.007, 0.008, 0.009],
        })

        enrichment = analyzer.enrich(diff_results)
        assert isinstance(enrichment, pd.DataFrame)


class TestGOEnrichmentEdgeCases:
    """Edge case tests for GO enrichment analysis."""

    def test_go_enrichment_empty_significant_set(self, model_with_subsystems, go_annotations):
        """Test enrichment with empty significant reaction set."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        significant = set()
        background = {"R1", "R2", "R3", "R4", "R5"}

        results = analyzer.enrich_reactions(significant, background, min_genes=1)

        # Should return empty DataFrame
        assert len(results) == 0

    def test_go_enrichment_significant_equals_background(self, model_with_subsystems, go_annotations):
        """Test when all background reactions are significant."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        all_reactions = {"R1", "R2", "R3", "R4", "R5"}
        results = analyzer.enrich_reactions(all_reactions, all_reactions, min_genes=1)

        # Fold enrichment should be 1.0 for all terms
        if len(results) > 0:
            assert all(results["fold_enrichment"] <= 1.0 + 1e-10)

    def test_go_enrichment_single_reaction(self, model_with_subsystems, go_annotations):
        """Test enrichment with single significant reaction."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        significant = {"R1"}
        background = {"R1", "R2", "R3", "R4", "R5"}

        results = analyzer.enrich_reactions(significant, background, min_genes=1)
        assert isinstance(results, pd.DataFrame)

    def test_go_enrichment_no_go_terms(self, model_with_subsystems):
        """Test enrichment when no GO annotations exist."""
        empty_annotations = pd.DataFrame(columns=["gene_id", "go_term", "go_name", "namespace"])

        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, empty_annotations)

        significant = {"R1", "R2"}
        background = {"R1", "R2", "R3", "R4", "R5"}

        results = analyzer.enrich_reactions(significant, background, min_genes=1)
        assert len(results) == 0

    def test_go_enrichment_invalid_namespace(self, model_with_subsystems, go_annotations):
        """Test filtering with invalid namespace."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        significant = {"R1", "R2"}
        background = {"R1", "R2", "R3", "R4", "R5"}

        # Invalid namespace should return empty results
        results = analyzer.enrich_reactions(
            significant, background, min_genes=1, namespace="INVALID"
        )
        assert len(results) == 0

    def test_go_enrichment_high_min_genes_threshold(self, model_with_subsystems, go_annotations):
        """Test enrichment with very high min_genes threshold."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        significant = {"R1", "R2", "R3"}
        background = {"R1", "R2", "R3", "R4", "R5"}

        # Very high threshold - no term should pass
        results = analyzer.enrich_reactions(significant, background, min_genes=100)
        assert len(results) == 0

    def test_go_reaction_with_no_genes(self, model_with_subsystems, go_annotations):
        """Test getting GO terms for reaction without GPR."""
        import cobra

        # Add reaction without genes
        r_no_gpr = cobra.Reaction("R_NO_GPR")
        r_no_gpr.add_metabolites({cobra.Metabolite("X", compartment="c"): -1})
        model_with_subsystems.add_reactions([r_no_gpr])

        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)
        terms = analyzer.get_reaction_go_terms("R_NO_GPR")

        assert len(terms) == 0

    def test_go_enrichment_from_differential_empty_results(self, model_with_subsystems, go_annotations):
        """Test GO enrichment from differential results with no significant."""
        analyzer = GOEnrichmentAnalyzer(model_with_subsystems, go_annotations)

        # All non-significant
        diff_results = pd.DataFrame({
            "reaction": ["R1", "R2", "R3", "R4", "R5"],
            "pvalue": [0.5, 0.6, 0.7, 0.8, 0.9],
            "padj_bh": [0.5, 0.6, 0.7, 0.8, 0.9],
        })

        results = analyzer.enrich_from_differential(
            diff_results, pvalue_threshold=0.05, min_genes=1
        )
        assert len(results) == 0


class TestCreateGOAnnotationsEdgeCases:
    """Edge case tests for GO annotation creation."""

    def test_create_go_annotations_empty_dict(self):
        """Test creating GO annotations from empty dictionary."""
        gene_to_go = {}
        df = create_go_annotations_from_dict(gene_to_go)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0

    def test_create_go_annotations_gene_with_no_terms(self):
        """Test creating annotations when gene has empty GO term list."""
        gene_to_go = {
            "gene1": [],  # No GO terms
            "gene2": [("GO:0006099", "TCA cycle", "BP")],
        }

        df = create_go_annotations_from_dict(gene_to_go)

        # Only gene2 should be in result
        assert len(df) == 1
        assert "gene2" in df["gene_id"].values

    def test_create_go_annotations_duplicate_terms(self):
        """Test handling of duplicate GO terms for same gene."""
        gene_to_go = {
            "gene1": [
                ("GO:0006096", "glycolytic process", "BP"),
                ("GO:0006096", "glycolytic process", "BP"),  # Duplicate
            ],
        }

        df = create_go_annotations_from_dict(gene_to_go)

        # Depending on implementation, may have 2 rows or deduplicated to 1
        assert len(df) >= 1
