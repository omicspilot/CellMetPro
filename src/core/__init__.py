"""Core module for metabolic profiling computations.

This module contains the fundamental algorithms and utilities for:
- COMPASS algorithm implementation for metabolic scoring
- Flux Balance Analysis (FBA) utilities
- Data preprocessing and normalization
"""

from .compass import CompassScorer
from .fba import FluxBalanceAnalyzer
from .preprocessing import DataLoader, normalize_expression

__all__ = [
    "CompassScorer",
    "FluxBalanceAnalyzer",
    "DataLoader",
    "normalize_expression",
]
