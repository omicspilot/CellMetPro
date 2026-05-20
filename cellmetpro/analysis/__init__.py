"""Analysis module for metabolic data interpretation.

This module provides tools for:
- Metabolic-based cell clustering
- Differential metabolic analysis between conditions
- Pathway-level aggregation and scoring
"""

from .clustering import (
    MetabolicClustering,
    benchmark_clustering_methods,
    compare_clusterings,
    evaluate_clustering,
    find_optimal_clusters,
)
from .differential import DifferentialAnalysis, PseudoBulkAnalysis
from .pathway import PathwayAnalyzer
from .perturbation import (
    compare_flux_distributions,
    identify_essential_genes,
    multi_knockout,
    predict_synthetic_lethality,
    sensitivity_analysis,
    simulate_expression_change,
)
from .trajectory import (
    cluster_trajectory_patterns,
    compute_metabolic_velocity,
    compute_pseudotime,
    fit_trajectory_genes,
    identify_branch_points,
    trajectory_differential,
)

__all__ = [
    "MetabolicClustering",
    "DifferentialAnalysis",
    "PseudoBulkAnalysis",
    "PathwayAnalyzer",
    # Clustering utilities
    "find_optimal_clusters",
    "evaluate_clustering",
    "compare_clusterings",
    "benchmark_clustering_methods",
    # Perturbation analysis
    "simulate_expression_change",
    "multi_knockout",
    "compare_flux_distributions",
    "identify_essential_genes",
    "sensitivity_analysis",
    "predict_synthetic_lethality",
    # Trajectory analysis
    "compute_pseudotime",
    "compute_metabolic_velocity",
    "identify_branch_points",
    "trajectory_differential",
    "fit_trajectory_genes",
    "cluster_trajectory_patterns",
]
