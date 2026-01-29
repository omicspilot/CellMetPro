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
from .differential import DifferentialAnalysis
from .pathway import PathwayAnalyzer

__all__ = [
    "MetabolicClustering",
    "DifferentialAnalysis",
    "PathwayAnalyzer",
    # Clustering utilities
    "find_optimal_clusters",
    "evaluate_clustering",
    "compare_clusterings",
    "benchmark_clustering_methods",
]
