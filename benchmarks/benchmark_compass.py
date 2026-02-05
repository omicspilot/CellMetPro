#!/usr/bin/env python
"""Benchmark COMPASS scoring performance.

Run with: python benchmarks/benchmark_compass.py
"""

import time
from functools import wraps

import cobra
import numpy as np
import pandas as pd

from cellmetpro.core.compass import CompassConfig, CompassScorer


def timer(func):
    """Decorator to time function execution."""

    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        return result, elapsed

    return wrapper


def create_test_model(n_reactions: int = 50) -> cobra.Model:
    """Create a test metabolic model with specified number of reactions."""
    model = cobra.Model("benchmark_model")

    # Create metabolites
    metabolites = [
        cobra.Metabolite(f"M{i}", compartment="c") for i in range(n_reactions + 1)
    ]

    # Create reactions with GPR rules
    genes = [f"GENE{i}" for i in range(n_reactions * 2)]

    for i in range(n_reactions):
        rxn = cobra.Reaction(f"R{i}")
        rxn.add_metabolites({metabolites[i]: -1, metabolites[i + 1]: 1})
        rxn.bounds = (0, 1000)

        # Create varied GPR rules
        if i % 3 == 0:
            rxn.gene_reaction_rule = genes[i]
        elif i % 3 == 1:
            rxn.gene_reaction_rule = f"{genes[i]} and {genes[i + 1]}"
        else:
            rxn.gene_reaction_rule = f"{genes[i]} or {genes[i + 1]}"

        model.add_reactions([rxn])

    # Add exchange reactions
    ex_in = cobra.Reaction("EX_IN")
    ex_in.add_metabolites({metabolites[0]: -1})
    ex_in.bounds = (-10, 0)

    ex_out = cobra.Reaction("EX_OUT")
    ex_out.add_metabolites({metabolites[-1]: -1})
    ex_out.bounds = (0, 1000)

    model.add_reactions([ex_in, ex_out])
    model.objective = "EX_OUT"

    return model


def create_expression_data(
    n_genes: int, n_cells: int, seed: int = 42
) -> pd.DataFrame:
    """Create synthetic expression data."""
    np.random.seed(seed)
    genes = [f"GENE{i}" for i in range(n_genes)]
    cells = [f"cell{i}" for i in range(n_cells)]

    data = np.random.rand(n_genes, n_cells) * 100

    return pd.DataFrame(data, index=genes, columns=cells)


@timer
def benchmark_penalty_computation(model, expression, config):
    """Benchmark reaction penalty computation."""
    scorer = CompassScorer(model, expression, config)
    penalties = scorer.compute_reaction_penalties()
    return penalties


def run_benchmark():
    """Run performance benchmark."""
    print("=" * 60)
    print("COMPASS Performance Benchmark")
    print("=" * 60)

    # Test configurations
    configs = [
        {
            "n_reactions": 50,
            "n_genes": 100,
            "n_cells": 50,
            "description": "Small dataset (50 rxns, 100 genes, 50 cells)",
        },
        {
            "n_reactions": 100,
            "n_genes": 200,
            "n_cells": 100,
            "description": "Medium dataset (100 rxns, 200 genes, 100 cells)",
        },
        {
            "n_reactions": 200,
            "n_genes": 400,
            "n_cells": 200,
            "description": "Large dataset (200 rxns, 400 genes, 200 cells)",
        },
    ]

    for cfg in configs:
        print(f"\n{cfg['description']}")
        print("-" * 50)

        model = create_test_model(cfg["n_reactions"])
        expression = create_expression_data(cfg["n_genes"], cfg["n_cells"])

        # Benchmark with optimizations enabled
        config_opt = CompassConfig(
            precompute_gpr=True,
            cache_max_fluxes=True,
            batch_size=50,
        )
        _, time_opt = benchmark_penalty_computation(model, expression, config_opt)
        print(f"  With optimizations:    {time_opt:.3f}s")

        # Benchmark without optimizations
        config_basic = CompassConfig(
            precompute_gpr=False,
            cache_max_fluxes=False,
            batch_size=0,
        )
        _, time_basic = benchmark_penalty_computation(model, expression, config_basic)
        print(f"  Without optimizations: {time_basic:.3f}s")

        if time_basic > 0:
            speedup = time_basic / time_opt
            print(f"  Speedup: {speedup:.2f}x")

    print("\n" + "=" * 60)
    print("Benchmark complete")
    print("=" * 60)


if __name__ == "__main__":
    run_benchmark()
