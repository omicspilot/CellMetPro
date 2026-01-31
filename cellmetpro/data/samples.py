"""Sample data generation and loading utilities.

This module provides synthetic datasets for testing and tutorials.
The data is designed to be biologically plausible while remaining
small enough for quick analysis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def get_sample_data_path() -> Path:
    """Get the path to the sample data directory.

    Returns
    -------
    Path
        Path to the sample data directory.
    """
    return Path(__file__).parent


def load_sample_expression(
    n_cells: int = 100,
    n_genes: int = 50,
    seed: int = 42,
) -> pd.DataFrame:
    """Load or generate sample gene expression data.

    Generates a synthetic scRNA-seq expression matrix with genes
    that include common metabolic pathway genes. The data has
    realistic characteristics including:
    - Sparse expression (many zeros)
    - Varying expression levels across cell types
    - Gene-gene correlations within pathways

    Parameters
    ----------
    n_cells : int, default=100
        Number of cells to generate.
    n_genes : int, default=50
        Number of genes to include.
    seed : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Expression matrix with genes as rows and cells as columns.

    Examples
    --------
    >>> expression = load_sample_expression()
    >>> print(expression.shape)
    (50, 100)
    >>> print(expression.index[:5].tolist())
    ['HK1', 'HK2', 'PFKM', 'PFKL', 'ALDOA']
    """
    np.random.seed(seed)

    # Metabolic pathway genes (glycolysis, TCA, oxidative phosphorylation)
    glycolysis_genes = [
        "HK1", "HK2", "PFKM", "PFKL", "ALDOA", "ALDOB",
        "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "PKLR",
    ]

    tca_genes = [
        "CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A",
        "OGDH", "SUCLA2", "SDHA", "SDHB", "FH", "MDH1", "MDH2",
    ]

    oxphos_genes = [
        "NDUFA1", "NDUFB1", "NDUFS1", "SDHC", "UQCRB",
        "COX4I1", "COX5A", "COX6A1", "ATP5A1", "ATP5B", "ATP5C1",
    ]

    amino_acid_genes = [
        "GLUL", "GLUD1", "GOT1", "GOT2", "GPT", "GPT2",
        "BCAT1", "BCAT2", "ASS1", "ASL",
    ]

    lipid_genes = [
        "FASN", "ACACA", "SCD", "ACLY", "ELOVL1",
    ]

    # Combine genes
    all_genes = (
        glycolysis_genes + tca_genes + oxphos_genes
        + amino_acid_genes + lipid_genes
    )

    # Select genes up to n_genes
    if n_genes <= len(all_genes):
        genes = all_genes[:n_genes]
    else:
        # Add generic genes if needed
        extra_genes = [f"GENE{i}" for i in range(n_genes - len(all_genes))]
        genes = all_genes + extra_genes

    # Generate cell names with cluster prefixes
    n_clusters = 4
    cells_per_cluster = n_cells // n_clusters
    cells = []
    cell_types = []
    for i in range(n_clusters):
        cluster_name = ["Proliferating", "Quiescent", "Hypoxic", "Oxidative"][i]
        n = cells_per_cluster if i < n_clusters - 1 else n_cells - len(cells)
        cells.extend([f"{cluster_name}_{j}" for j in range(n)])
        cell_types.extend([cluster_name] * n)

    # Generate expression data with cluster-specific patterns
    expression = np.zeros((len(genes), n_cells))

    # Base expression (log-normal distribution)
    base_expr = np.random.lognormal(mean=1, sigma=1.5, size=(len(genes), n_cells))

    # Add cluster-specific patterns
    for i, cell in enumerate(cells):
        cluster = cell.split("_")[0]

        if cluster == "Proliferating":
            # High glycolysis
            for j, gene in enumerate(genes):
                if gene in glycolysis_genes:
                    expression[j, i] = base_expr[j, i] * 3
                else:
                    expression[j, i] = base_expr[j, i]

        elif cluster == "Quiescent":
            # Low overall metabolism
            expression[:, i] = base_expr[:, i] * 0.5

        elif cluster == "Hypoxic":
            # High glycolysis, low oxphos
            for j, gene in enumerate(genes):
                if gene in glycolysis_genes:
                    expression[j, i] = base_expr[j, i] * 2.5
                elif gene in oxphos_genes:
                    expression[j, i] = base_expr[j, i] * 0.3
                else:
                    expression[j, i] = base_expr[j, i]

        elif cluster == "Oxidative":
            # High TCA and oxphos
            for j, gene in enumerate(genes):
                if gene in tca_genes or gene in oxphos_genes:
                    expression[j, i] = base_expr[j, i] * 2.5
                else:
                    expression[j, i] = base_expr[j, i]

    # Add sparsity (dropout)
    dropout_mask = np.random.random(expression.shape) < 0.3
    expression[dropout_mask] = 0

    # Round to integers and ensure non-negative
    expression = np.maximum(0, np.round(expression)).astype(float)

    return pd.DataFrame(expression, index=genes, columns=cells)


def load_sample_groups(
    n_cells: int = 100,
    seed: int = 42,
) -> pd.DataFrame:
    """Load or generate sample group annotations.

    Generates group annotations matching the sample expression data,
    with cells assigned to treatment conditions and cell types.

    Parameters
    ----------
    n_cells : int, default=100
        Number of cells (must match expression data).
    seed : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns 'cell', 'group', 'cell_type', 'treatment'.

    Examples
    --------
    >>> groups = load_sample_groups()
    >>> print(groups.head())
              cell         group     cell_type treatment
    0  Proliferating_0  Proliferating  Proliferating   control
    1  Proliferating_1  Proliferating  Proliferating   control
    """
    np.random.seed(seed)

    # Generate matching cell names
    n_clusters = 4
    cells_per_cluster = n_cells // n_clusters
    cells = []
    cell_types = []
    cluster_names = ["Proliferating", "Quiescent", "Hypoxic", "Oxidative"]

    for i in range(n_clusters):
        cluster_name = cluster_names[i]
        n = cells_per_cluster if i < n_clusters - 1 else n_cells - len(cells)
        cells.extend([f"{cluster_name}_{j}" for j in range(n)])
        cell_types.extend([cluster_name] * n)

    # Assign treatments (balanced across cell types)
    treatments = []
    for i, ct in enumerate(cell_types):
        if i % 2 == 0:
            treatments.append("control")
        else:
            treatments.append("treatment")

    return pd.DataFrame({
        "cell": cells,
        "group": cell_types,
        "cell_type": cell_types,
        "treatment": treatments,
    })


def load_sample_reaction_scores(
    n_cells: int = 100,
    n_reactions: int = 30,
    seed: int = 42,
) -> pd.DataFrame:
    """Load or generate sample reaction scores.

    Generates synthetic COMPASS-like reaction scores with
    pathway-correlated patterns. Useful for testing visualization
    and analysis functions without running the full COMPASS pipeline.

    Parameters
    ----------
    n_cells : int, default=100
        Number of cells.
    n_reactions : int, default=30
        Number of reactions.
    seed : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Reaction scores matrix (reactions x cells).

    Examples
    --------
    >>> scores = load_sample_reaction_scores()
    >>> print(scores.shape)
    (30, 100)
    >>> print(scores.index[:5].tolist())
    ['HEX1', 'PGI', 'PFK', 'FBA', 'TPI']
    """
    np.random.seed(seed)

    # Metabolic reactions
    glycolysis_rxns = [
        "HEX1", "PGI", "PFK", "FBA", "TPI",
        "GAPD", "PGK", "PGM", "ENO", "PYK",
    ]

    tca_rxns = [
        "CS", "ACONT", "ICDHx", "AKGD", "SUCOAS",
        "SUCD", "FUM", "MDH", "PC", "ME",
    ]

    transport_rxns = [
        "GLCt1", "PYRt2m", "LACt", "ATPt", "O2t",
        "CO2t", "NH4t", "GLNt", "GLUt", "AKGt",
    ]

    all_rxns = glycolysis_rxns + tca_rxns + transport_rxns

    if n_reactions <= len(all_rxns):
        reactions = all_rxns[:n_reactions]
    else:
        extra = [f"RXN{i}" for i in range(n_reactions - len(all_rxns))]
        reactions = all_rxns + extra

    # Get matching cell names and groups
    groups = load_sample_groups(n_cells, seed)
    cells = groups["cell"].tolist()

    # Generate scores with group-specific patterns
    scores = np.random.random((len(reactions), n_cells)) * 0.3 + 0.2

    for i, cell in enumerate(cells):
        cell_type = groups.iloc[i]["cell_type"]

        if cell_type == "Proliferating":
            # Higher glycolysis scores
            for j, rxn in enumerate(reactions):
                if rxn in glycolysis_rxns:
                    scores[j, i] += 0.3
        elif cell_type == "Oxidative":
            # Higher TCA scores
            for j, rxn in enumerate(reactions):
                if rxn in tca_rxns:
                    scores[j, i] += 0.25
        elif cell_type == "Hypoxic":
            # Higher glycolysis, lower TCA
            for j, rxn in enumerate(reactions):
                if rxn in glycolysis_rxns:
                    scores[j, i] += 0.35
                elif rxn in tca_rxns:
                    scores[j, i] -= 0.1

    # Add noise and clip to valid range
    scores += np.random.normal(0, 0.05, scores.shape)
    scores = np.clip(scores, 0, 1)

    return pd.DataFrame(scores, index=reactions, columns=cells)


def create_sample_model():
    """Create a simple metabolic model for testing.

    Returns a minimal COBRApy model with core metabolic reactions
    suitable for testing COMPASS and FBA functionality.

    Returns
    -------
    cobra.Model
        A simple metabolic model with glycolysis and exchange reactions.

    Examples
    --------
    >>> model = create_sample_model()
    >>> print(f"Model: {model.id}")
    >>> print(f"Reactions: {len(model.reactions)}")
    >>> print(f"Genes: {len(model.genes)}")
    """
    import cobra

    model = cobra.Model("sample_metabolic_model")

    # Create metabolites
    metabolites = {
        "glc_e": cobra.Metabolite("glc_e", name="Glucose (extracellular)", compartment="e"),
        "glc_c": cobra.Metabolite("glc_c", name="Glucose (cytosol)", compartment="c"),
        "g6p_c": cobra.Metabolite("g6p_c", name="Glucose-6-phosphate", compartment="c"),
        "f6p_c": cobra.Metabolite("f6p_c", name="Fructose-6-phosphate", compartment="c"),
        "fbp_c": cobra.Metabolite("fbp_c", name="Fructose-1,6-bisphosphate", compartment="c"),
        "g3p_c": cobra.Metabolite("g3p_c", name="Glyceraldehyde-3-phosphate", compartment="c"),
        "pyr_c": cobra.Metabolite("pyr_c", name="Pyruvate (cytosol)", compartment="c"),
        "lac_c": cobra.Metabolite("lac_c", name="Lactate (cytosol)", compartment="c"),
        "lac_e": cobra.Metabolite("lac_e", name="Lactate (extracellular)", compartment="e"),
        "atp_c": cobra.Metabolite("atp_c", name="ATP", compartment="c"),
        "adp_c": cobra.Metabolite("adp_c", name="ADP", compartment="c"),
    }

    # Create reactions with GPR rules
    reactions = []

    # Glucose transport
    glc_transport = cobra.Reaction("GLCt1")
    glc_transport.name = "Glucose transport"
    glc_transport.add_metabolites({
        metabolites["glc_e"]: -1,
        metabolites["glc_c"]: 1,
    })
    glc_transport.bounds = (-10, 10)
    glc_transport.gene_reaction_rule = "SLC2A1 or SLC2A3"
    reactions.append(glc_transport)

    # Hexokinase
    hex1 = cobra.Reaction("HEX1")
    hex1.name = "Hexokinase"
    hex1.add_metabolites({
        metabolites["glc_c"]: -1,
        metabolites["atp_c"]: -1,
        metabolites["g6p_c"]: 1,
        metabolites["adp_c"]: 1,
    })
    hex1.bounds = (0, 1000)
    hex1.gene_reaction_rule = "HK1 or HK2"
    reactions.append(hex1)

    # Phosphoglucose isomerase
    pgi = cobra.Reaction("PGI")
    pgi.name = "Phosphoglucose isomerase"
    pgi.add_metabolites({
        metabolites["g6p_c"]: -1,
        metabolites["f6p_c"]: 1,
    })
    pgi.bounds = (-1000, 1000)
    pgi.gene_reaction_rule = "GPI"
    reactions.append(pgi)

    # Phosphofructokinase
    pfk = cobra.Reaction("PFK")
    pfk.name = "Phosphofructokinase"
    pfk.add_metabolites({
        metabolites["f6p_c"]: -1,
        metabolites["atp_c"]: -1,
        metabolites["fbp_c"]: 1,
        metabolites["adp_c"]: 1,
    })
    pfk.bounds = (0, 1000)
    pfk.gene_reaction_rule = "PFKM or PFKL"
    reactions.append(pfk)

    # Aldolase (simplified: FBP -> 2 G3P)
    fba = cobra.Reaction("FBA")
    fba.name = "Aldolase"
    fba.add_metabolites({
        metabolites["fbp_c"]: -1,
        metabolites["g3p_c"]: 2,
    })
    fba.bounds = (-1000, 1000)
    fba.gene_reaction_rule = "ALDOA or ALDOB"
    reactions.append(fba)

    # Glycolysis (simplified: G3P -> Pyruvate + ATP)
    gapdh_to_pyk = cobra.Reaction("GAPD_PYK")
    gapdh_to_pyk.name = "Glycolysis (G3P to Pyruvate)"
    gapdh_to_pyk.add_metabolites({
        metabolites["g3p_c"]: -1,
        metabolites["adp_c"]: -2,
        metabolites["pyr_c"]: 1,
        metabolites["atp_c"]: 2,
    })
    gapdh_to_pyk.bounds = (0, 1000)
    gapdh_to_pyk.gene_reaction_rule = "GAPDH and PGK1 and ENO1 and PKM"
    reactions.append(gapdh_to_pyk)

    # Lactate dehydrogenase
    ldh = cobra.Reaction("LDH")
    ldh.name = "Lactate dehydrogenase"
    ldh.add_metabolites({
        metabolites["pyr_c"]: -1,
        metabolites["lac_c"]: 1,
    })
    ldh.bounds = (-1000, 1000)
    ldh.gene_reaction_rule = "LDHA or LDHB"
    reactions.append(ldh)

    # Lactate transport
    lac_transport = cobra.Reaction("LACt")
    lac_transport.name = "Lactate transport"
    lac_transport.add_metabolites({
        metabolites["lac_c"]: -1,
        metabolites["lac_e"]: 1,
    })
    lac_transport.bounds = (-1000, 1000)
    lac_transport.gene_reaction_rule = "SLC16A1 or SLC16A3"
    reactions.append(lac_transport)

    # Exchange reactions
    ex_glc = cobra.Reaction("EX_glc_e")
    ex_glc.name = "Glucose exchange"
    ex_glc.add_metabolites({metabolites["glc_e"]: -1})
    ex_glc.bounds = (-10, 0)
    reactions.append(ex_glc)

    ex_lac = cobra.Reaction("EX_lac_e")
    ex_lac.name = "Lactate exchange"
    ex_lac.add_metabolites({metabolites["lac_e"]: -1})
    ex_lac.bounds = (0, 1000)
    reactions.append(ex_lac)

    # ATP maintenance (sink for ATP)
    atp_sink = cobra.Reaction("ATPM")
    atp_sink.name = "ATP maintenance"
    atp_sink.add_metabolites({
        metabolites["atp_c"]: -1,
        metabolites["adp_c"]: 1,
    })
    atp_sink.bounds = (0, 1000)
    reactions.append(atp_sink)

    # Add reactions to model
    model.add_reactions(reactions)

    # Set objective (maximize lactate production)
    model.objective = "EX_lac_e"

    return model
