"""Generate synthetic scRNA-seq test data for CellMetPro.

Produces a biologically meaningful 200-cell × ~2000-gene AnnData object with:
  - CD4+ T cells (80 cells)   — glycolysis / OXPHOS signature
  - M1 Macrophages (80 cells) — aerobic glycolysis (Warburg-like)
  - Naive B cells (40 cells)  — fatty acid oxidation signature

Output: data/synthetic_scrna.h5ad
"""

import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
RNG = np.random.default_rng(42)

# ---------------------------------------------------------------------------
# Metabolic gene signatures (pathway → genes)
# ---------------------------------------------------------------------------
METABOLIC_GENES: dict[str, list[str]] = {
    "Glycolysis": [
        "HK1",
        "HK2",
        "GPI",
        "PFKL",
        "ALDOA",
        "GAPDH",
        "PGK1",
        "PKM",
        "LDHA",
        "ENO1",
    ],
    "TCA_Cycle": ["IDH1", "IDH2", "OGDH", "SDHA", "FH", "MDH2", "CS", "ACO2"],
    "OXPHOS": ["NDUFS1", "SDHB", "UQCRC1", "COX5A", "ATP5A1"],
    "Fatty_Acid_Oxidation": ["CPT1A", "ACADM", "HADHA", "ACADVL"],
    "Glutaminolysis": ["GLS", "GLUD1", "GOT2"],
    "Pentose_Phosphate": ["G6PD", "PGD", "TKT"],
    "Transport": ["SLC2A1", "SLC1A5"],
}

ALL_METABOLIC = [g for genes in METABOLIC_GENES.values() for g in genes]  # 38 genes

# ---------------------------------------------------------------------------
# Cell-type-specific mean expression per metabolic gene
# ---------------------------------------------------------------------------
# Keys that are NOT listed default to a low baseline (mean = 1).
SIGNATURES: dict[str, dict[str, float]] = {
    "CD4_T": {
        # glycolysis up
        "LDHA": 40,
        "PKM": 35,
        "SLC2A1": 30,
        "GPI": 25,
        "GAPDH": 45,
        "ENO1": 30,
        # OXPHOS moderate
        "COX5A": 20,
        "ATP5A1": 18,
        "NDUFS1": 15,
        # FAO low
        "CPT1A": 3,
        "ACADM": 3,
    },
    "M1_Macrophage": {
        # aerobic glycolysis (Warburg)
        "HK2": 50,
        "PFKL": 40,
        "LDHA": 55,
        "PKM": 45,
        "ALDOA": 35,
        "ENO1": 38,
        # OXPHOS suppressed
        "NDUFS1": 5,
        "SDHB": 5,
        "COX5A": 6,
        "ATP5A1": 6,
        # inflammatory / glutaminolysis
        "GLS": 25,
        "GLUD1": 20,
    },
    "Naive_B": {
        # FAO signature
        "CPT1A": 40,
        "ACADM": 38,
        "HADHA": 35,
        "ACADVL": 30,
        # TCA moderate
        "IDH1": 22,
        "MDH2": 20,
        "CS": 18,
        # glycolysis low-moderate
        "GAPDH": 15,
        "PKM": 12,
        # pentose phosphate
        "G6PD": 20,
        "PGD": 18,
    },
}

BASELINE_MEAN = 1.0  # low mean for genes not in signature
DISPERSION = 0.8  # NB dispersion parameter (r)
N_FILLER = 1962  # non-metabolic filler genes to reach ~2000 total


def _nb_counts(mean: float, dispersion: float, n: int) -> np.ndarray:
    """Sample n counts from a Negative Binomial(r=dispersion, p) with given mean."""
    # NB parameterization: mean = r*(1-p)/p  →  p = r / (r + mean)
    if mean <= 0:
        return np.zeros(n, dtype=np.int32)
    p = dispersion / (dispersion + mean)
    return RNG.negative_binomial(n=dispersion, p=p, size=n).astype(np.int32)


def _build_population(
    cell_type: str, n_cells: int, gene_names: list[str]
) -> np.ndarray:
    """Return (n_cells, n_genes) count matrix for one cell population."""
    sig = SIGNATURES[cell_type]
    matrix = np.zeros((n_cells, len(gene_names)), dtype=np.int32)
    for j, gene in enumerate(gene_names):
        mean = sig.get(gene, BASELINE_MEAN)
        matrix[:, j] = _nb_counts(mean, DISPERSION, n_cells)
    return matrix


def generate(output_path: str = "data/synthetic_scrna.h5ad") -> ad.AnnData:
    # -----------------------------------------------------------------------
    # Gene list
    # -----------------------------------------------------------------------
    filler_genes = [f"GENE_{i:04d}" for i in range(N_FILLER)]
    gene_names = ALL_METABOLIC + filler_genes  # 38 + 1962 = 2000

    is_metabolic = [True] * len(ALL_METABOLIC) + [False] * N_FILLER

    # -----------------------------------------------------------------------
    # Cell populations
    # -----------------------------------------------------------------------
    populations = [
        ("CD4_T", 80),
        ("M1_Macrophage", 80),
        ("Naive_B", 40),
    ]

    blocks: list[np.ndarray] = []
    cell_types: list[str] = []

    for cell_type, n_cells in populations:
        blocks.append(_build_population(cell_type, n_cells, gene_names))
        cell_types.extend([cell_type] * n_cells)

    X = np.vstack(blocks)  # (200, 2000)

    # -----------------------------------------------------------------------
    # Build AnnData
    # -----------------------------------------------------------------------
    obs = pd.DataFrame(
        {"cell_type": cell_types, "n_counts": X.sum(axis=1)},
        index=[f"cell_{i:04d}" for i in range(X.shape[0])],
    )
    var = pd.DataFrame(
        {"gene_name": gene_names, "is_metabolic": is_metabolic},
        index=gene_names,
    )

    adata = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
    adata.raw = adata.copy()  # preserve raw counts

    # -----------------------------------------------------------------------
    # Preprocessing (CellMetPro-compatible)
    # -----------------------------------------------------------------------
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    adata.write_h5ad(output_path)
    print(f"Saved {adata.n_obs} cells × {adata.n_vars} genes → {output_path}")
    print(f"Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")
    return adata


if __name__ == "__main__":
    generate()
