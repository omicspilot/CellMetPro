"""cellMetPro: Cellular Metabolic Profiler.

A Python-based tool for analyzing and visualizing metabolic profiles
from single-cell RNA-seq data using Flux Balance Analysis (FBA) and
Genome-Scale Metabolic Models (GEMs).
"""

__version__ = "0.1.0"
__author__ = "Oumar Ndiaye"

from . import analysis, core, data, visualization

__all__ = ["core", "analysis", "visualization", "data"]
