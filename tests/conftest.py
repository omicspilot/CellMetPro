"""Pytest configuration and shared fixtures."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture(scope="session")
def simple_model():
    """Create a simple test metabolic model (session-scoped for efficiency)."""
    import cobra

    model = cobra.Model("test_model")

    # Create metabolites
    A = cobra.Metabolite("A", compartment="c", name="Metabolite A")
    B = cobra.Metabolite("B", compartment="c", name="Metabolite B")
    C = cobra.Metabolite("C", compartment="c", name="Metabolite C")

    # Create reactions with GPR rules
    r1 = cobra.Reaction("R1", name="Reaction 1")
    r1.add_metabolites({A: -1, B: 1})
    r1.bounds = (0, 1000)
    r1.gene_reaction_rule = "gene1"

    r2 = cobra.Reaction("R2", name="Reaction 2")
    r2.add_metabolites({B: -1, C: 1})
    r2.bounds = (0, 1000)
    r2.gene_reaction_rule = "gene2 and gene3"

    r3 = cobra.Reaction("R3", name="Reaction 3")
    r3.add_metabolites({A: -1, C: 1})
    r3.bounds = (0, 1000)
    r3.gene_reaction_rule = "gene1 or gene4"

    # Exchange reactions
    ex_A = cobra.Reaction("EX_A", name="A exchange")
    ex_A.add_metabolites({A: -1})
    ex_A.bounds = (-10, 0)

    ex_C = cobra.Reaction("EX_C", name="C exchange")
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
def large_expression_df():
    """Create a larger expression DataFrame for performance tests."""
    n_genes = 100
    n_cells = 50

    genes = [f"GENE{i}" for i in range(n_genes)]
    cells = [f"cell{i}" for i in range(n_cells)]

    np.random.seed(42)
    data = np.random.rand(n_genes, n_cells) * 100

    return pd.DataFrame(data, index=genes, columns=cells)


@pytest.fixture
def adata(expression_df):
    """Create a simple AnnData object."""
    import anndata as ad

    # AnnData is cells x genes
    adata = ad.AnnData(expression_df.T)
    adata.obs["cluster"] = ["A", "A", "B", "B", "C"]
    return adata


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


@pytest.fixture
def tmp_model_sbml(tmp_path, simple_model):
    """Create a temporary SBML model file."""
    import cobra

    sbml_path = tmp_path / "model.xml"
    cobra.io.write_sbml_model(simple_model, str(sbml_path))
    return sbml_path


@pytest.fixture
def tmp_model_json(tmp_path, simple_model):
    """Create a temporary JSON model file."""
    import cobra

    json_path = tmp_path / "model.json"
    cobra.io.save_json_model(simple_model, str(json_path))
    return json_path
