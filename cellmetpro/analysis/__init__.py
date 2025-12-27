"""Analysis module for metabolic data interpretation.

This module provides tools for:
- Metabolic-based cell clustering
- Differential metabolic analysis between conditions
- Pathway-level aggregation and scoring
"""

from .clustering import MetabolicClustering
from .differential import DifferentialAnalysis
from .pathway import PathwayAnalyzer

__all__ = [
    "MetabolicClustering",
    "DifferentialAnalysis",
    "PathwayAnalyzer",
]
