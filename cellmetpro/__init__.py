"""cellMetPro: Cellular Metabolic Profiler.

A Python-based tool for analyzing and visualizing metabolic profiles
from single-cell RNA-seq data using Flux Balance Analysis (FBA) and
Genome-Scale Metabolic Models (GEMs).
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cellmetpro")
except PackageNotFoundError:
    __version__ = "unknown"
__author__ = "Oumar Ndiaye"

__all__ = ["core", "analysis", "visualization", "data", "reporting"]


def __getattr__(name: str):
    if name in __all__:
        import importlib

        mod = importlib.import_module(f".{name}", __name__)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
