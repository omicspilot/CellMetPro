"""Data loading and preprocessing utilities.

This module handles loading scRNA-seq data from various formats
and preprocessing steps required for metabolic analysis.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import anndata as ad


class DataLoader:
    """Load scRNA-seq data from various formats.

    Supports loading from CSV, TSV, and h5ad (AnnData) formats.

    Parameters
    ----------
    filepath : str or Path
        Path to the data file.

    Attributes
    ----------
    adata : anndata.AnnData
        The loaded data as an AnnData object.
    """

    def __init__(self, filepath: str | Path) -> None:
        self.filepath = Path(filepath)
        self.adata: ad.AnnData | None = None

    def load(self) -> ad.AnnData:
        """Load the data file.

        Returns
        -------
        anndata.AnnData
            The loaded data.

        Raises
        ------
        ValueError
            If the file format is not supported.
        """
        raise NotImplementedError

    def load_csv(self, **kwargs) -> ad.AnnData:
        """Load data from CSV file.

        Parameters
        ----------
        **kwargs
            Additional arguments passed to pd.read_csv.

        Returns
        -------
        anndata.AnnData
            The loaded data.
        """
        raise NotImplementedError

    def load_h5ad(self) -> ad.AnnData:
        """Load data from h5ad file.

        Returns
        -------
        anndata.AnnData
            The loaded data.
        """
        raise NotImplementedError


def normalize_expression(
    adata: ad.AnnData,
    target_sum: float = 1e4,
    log_transform: bool = True,
) -> ad.AnnData:
    """Normalize gene expression data.

    Parameters
    ----------
    adata : anndata.AnnData
        The data to normalize.
    target_sum : float
        Target sum for normalization (default: 10,000 for TPM-like).
    log_transform : bool
        Whether to log-transform after normalization.

    Returns
    -------
    anndata.AnnData
        Normalized data.
    """
    raise NotImplementedError
