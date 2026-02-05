"""Data loading and preprocessing utilities.

This module handles loading scRNA-seq data from various formats
and preprocessing steps required for metabolic analysis.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import issparse

logger = logging.getLogger(__name__)


class DataLoader:
    """Load scRNA-seq data from various formats.

    Supports loading from CSV, TSV, h5ad (AnnData), MTX, and Seurat RDS formats.

    Parameters
    ----------
    filepath : str or Path
        Path to the data file.

    Attributes
    ----------
    adata : anndata.AnnData
        The loaded data as an AnnData object.

    Examples
    --------
    >>> from cellmetpro.core.preprocessing import DataLoader
    >>> loader = DataLoader("data/expression.h5ad")
    >>> adata = loader.load()
    >>> print(adata.shape)

    For Seurat objects (requires R and rpy2):
    >>> loader = DataLoader("data/seurat_object.rds")
    >>> adata = loader.load()
    """

    def __init__(self, filepath: str | Path) -> None:
        self.filepath = Path(filepath)
        self.adata: ad.AnnData | None = None

    def load(self) -> ad.AnnData:
        """Load the data file based on extension.

        Automatically detects file format from extension and loads
        appropriately.

        Returns
        -------
        anndata.AnnData
            The loaded data.

        Raises
        ------
        ValueError
            If the file format is not supported.
        FileNotFoundError
            If the file doesn't exist.
        """
        if not self.filepath.exists():
            raise FileNotFoundError(f"File not found: {self.filepath}")

        suffix = self.filepath.suffix.lower()
        logger.info(f"Loading data from {self.filepath} (format: {suffix})")

        if suffix == ".h5ad":
            self.adata = self.load_h5ad()
        elif suffix == ".csv":
            self.adata = self.load_csv()
        elif suffix in {".tsv", ".txt"}:
            self.adata = self.load_csv(sep="\t")
        elif suffix == ".mtx":
            self.adata = self.load_mtx()
        elif suffix == ".rds":
            self.adata = self.load_seurat()
        else:
            raise ValueError(
                f"Unsupported file format: {suffix}. "
                "Supported formats: .h5ad, .csv, .tsv, .txt, .mtx, .rds (Seurat)"
            )

        logger.info(
            f"Loaded data with {self.adata.n_obs} cells and {self.adata.n_vars} genes"
        )
        return self.adata

    def load_csv(self, sep: str = ",", **kwargs) -> ad.AnnData:
        """Load data from CSV/TSV file.

        Expects genes as rows and cells as columns.

        Parameters
        ----------
        sep : str, default=","
            Column separator.
        **kwargs
            Additional arguments passed to pd.read_csv.

        Returns
        -------
        anndata.AnnData
            The loaded data.
        """
        df = pd.read_csv(self.filepath, index_col=0, sep=sep, **kwargs)

        # Assume genes x cells format, transpose for AnnData
        # AnnData expects cells x genes
        adata = ad.AnnData(df.T)
        adata.var_names = df.index.astype(str)
        adata.obs_names = df.columns.astype(str)

        self.adata = adata
        return adata

    def load_h5ad(self) -> ad.AnnData:
        """Load data from h5ad file.

        Returns
        -------
        anndata.AnnData
            The loaded data.
        """
        self.adata = ad.read_h5ad(self.filepath)
        return self.adata

    def load_mtx(
        self,
        genes_file: str | Path | None = None,
        barcodes_file: str | Path | None = None,
    ) -> ad.AnnData:
        """Load data from Matrix Market (MTX) format.

        Looks for genes.tsv and barcodes.tsv in the same directory
        if not specified.

        Parameters
        ----------
        genes_file : str or Path, optional
            Path to genes file.
        barcodes_file : str or Path, optional
            Path to barcodes/cells file.

        Returns
        -------
        anndata.AnnData
            The loaded data.
        """
        from scipy.io import mmread

        # Load matrix
        matrix = mmread(self.filepath).T.tocsr()  # Transpose to cells x genes

        # Find gene and barcode files
        parent_dir = self.filepath.parent

        if genes_file is None:
            gene_file_names = [
                "genes.tsv",
                "features.tsv",
                "genes.tsv.gz",
                "features.tsv.gz",
            ]
            for name in gene_file_names:
                candidate = parent_dir / name
                if candidate.exists():
                    genes_file = candidate
                    break

        if barcodes_file is None:
            for name in ["barcodes.tsv", "cells.tsv", "barcodes.tsv.gz"]:
                candidate = parent_dir / name
                if candidate.exists():
                    barcodes_file = candidate
                    break

        # Load gene names
        if genes_file is not None:
            genes_df = pd.read_csv(genes_file, sep="\t", header=None)
            if genes_df.shape[1] == 1:
                gene_names = genes_df.iloc[:, 0].values
            else:
                gene_names = genes_df.iloc[:, 1].values
        else:
            gene_names = [f"gene_{i}" for i in range(matrix.shape[1])]

        # Load cell names
        if barcodes_file is not None:
            barcodes_df = pd.read_csv(barcodes_file, sep="\t", header=None)
            cell_names = barcodes_df.iloc[:, 0].values
        else:
            cell_names = [f"cell_{i}" for i in range(matrix.shape[0])]

        # Create AnnData
        self.adata = ad.AnnData(matrix)
        self.adata.var_names = pd.Index(gene_names).astype(str)
        self.adata.obs_names = pd.Index(cell_names).astype(str)

        return self.adata

    def load_seurat(
        self,
        assay: str | None = None,
        slot: str = "data",
    ) -> ad.AnnData:
        """Load data from Seurat RDS file.

        Requires R and rpy2 to be installed. The Seurat package must be
        installed in R.

        Parameters
        ----------
        assay : str, optional
            Seurat assay to extract (e.g., "RNA", "SCT"). If None, uses
            the default assay.
        slot : str, default="data"
            Data slot to extract: "counts", "data", or "scale.data".

        Returns
        -------
        anndata.AnnData
            The loaded data.

        Raises
        ------
        ImportError
            If rpy2 is not installed.
        RuntimeError
            If R or Seurat is not available.

        Notes
        -----
        For best compatibility, you can also export from R:

        .. code-block:: r

            # In R:
            library(Seurat)
            library(SeuratDisk)
            SaveH5Seurat(seurat_obj, "output.h5seurat")
            Convert("output.h5seurat", dest="h5ad")

        Examples
        --------
        >>> loader = DataLoader("seurat_object.rds")
        >>> adata = loader.load_seurat(assay="RNA", slot="data")
        """
        self.adata = load_seurat_rds(
            self.filepath,
            assay=assay,
            slot=slot,
        )
        return self.adata


def load_seurat_rds(
    filepath: str | Path,
    assay: str | None = None,
    slot: str = "data",
) -> ad.AnnData:
    """Load a Seurat object from RDS file.

    Converts a Seurat object to AnnData format using rpy2 to interface
    with R. Extracts the expression matrix and cell/gene metadata.

    Parameters
    ----------
    filepath : str or Path
        Path to the .rds file containing the Seurat object.
    assay : str, optional
        Seurat assay to extract (e.g., "RNA", "SCT"). If None, uses
        the default assay.
    slot : str, default="data"
        Data slot to extract:
        - "counts": Raw counts
        - "data": Normalized data (default)
        - "scale.data": Scaled data

    Returns
    -------
    anndata.AnnData
        The converted AnnData object containing:
        - X: Expression matrix from the specified slot
        - obs: Cell metadata from seurat_obj@meta.data
        - var: Gene metadata (names)
        - obsm: Embeddings (UMAP, PCA, etc.) if present

    Raises
    ------
    ImportError
        If rpy2 is not installed.
    RuntimeError
        If R or Seurat is not available, or the file is not a valid
        Seurat object.

    Examples
    --------
    >>> from cellmetpro.core.preprocessing import load_seurat_rds
    >>> adata = load_seurat_rds("pbmc3k.rds", assay="RNA", slot="data")
    >>> print(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    Notes
    -----
    Alternative: Export from R using SeuratDisk for better compatibility:

    .. code-block:: r

        # In R:
        library(Seurat)
        library(SeuratDisk)

        # Load your Seurat object
        seurat_obj <- readRDS("seurat_object.rds")

        # Convert to h5ad (AnnData format)
        SaveH5Seurat(seurat_obj, filename = "seurat.h5seurat")
        Convert("seurat.h5seurat", dest = "h5ad")

    Then load directly with AnnData:

    >>> import anndata as ad
    >>> adata = ad.read_h5ad("seurat.h5ad")
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.packages import importr
    except ImportError:
        raise ImportError(
            "rpy2 is required to load Seurat objects. Install with:\n"
            "  pip install rpy2\n\n"
            "You also need R and Seurat installed. Alternatively, export "
            "the Seurat object to h5ad format in R:\n"
            "  library(SeuratDisk)\n"
            "  SaveH5Seurat(obj, 'output.h5seurat')\n"
            "  Convert('output.h5seurat', dest='h5ad')"
        )

    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filepath}")

    logger.info(f"Loading Seurat object from {filepath}")

    # Activate automatic conversion
    numpy2ri.activate()
    pandas2ri.activate()

    try:
        # Import R packages
        base = importr("base")

        # Check if Seurat is installed
        try:
            importr("Seurat")
        except Exception:
            raise RuntimeError(
                "Seurat R package is not installed. Install in R with:\n"
                "  install.packages('Seurat')"
            )

        # Load the RDS file
        seurat_obj = base.readRDS(str(filepath))

        # Verify it's a Seurat object
        obj_class = ro.r("class")(seurat_obj)
        if "Seurat" not in list(obj_class):
            raise RuntimeError(
                f"File does not contain a Seurat object. Found: {list(obj_class)}"
            )

        # Get default assay if not specified
        if assay is None:
            assay = ro.r("DefaultAssay")(seurat_obj)[0]
            logger.info(f"Using default assay: {assay}")

        # Extract expression matrix
        # GetAssayData returns a sparse matrix
        expr_matrix = ro.r(
            f"as.matrix(GetAssayData({seurat_obj.r_repr()}, "
            f'assay="{assay}", slot="{slot}"))'
        )

        # Convert to numpy array
        expr_array = np.array(expr_matrix)

        # Get gene and cell names
        gene_names = list(ro.r("rownames")(expr_matrix))
        cell_names = list(ro.r("colnames")(expr_matrix))

        # Create AnnData (transpose: Seurat is genes x cells, AnnData is cells x genes)
        adata = ad.AnnData(expr_array.T)
        adata.var_names = pd.Index(gene_names).astype(str)
        adata.obs_names = pd.Index(cell_names).astype(str)

        # Extract cell metadata
        try:
            meta_data = ro.r("as.data.frame")(seurat_obj.rx2("meta.data"))
            meta_df = pd.DataFrame(meta_data)
            meta_df.index = adata.obs_names
            adata.obs = meta_df
            logger.info(f"Loaded {len(meta_df.columns)} metadata columns")
        except Exception as e:
            logger.warning(f"Could not extract metadata: {e}")

        # Extract embeddings (UMAP, PCA, etc.)
        try:
            reductions = ro.r("names")(seurat_obj.rx2("reductions"))
            if reductions is not ro.NULL:
                for red_name in list(reductions):
                    try:
                        red_data = ro.r(
                            f'Embeddings({seurat_obj.r_repr()}, reduction="{red_name}")'
                        )
                        red_array = np.array(red_data)
                        adata.obsm[f"X_{red_name}"] = red_array
                        logger.info(f"Loaded embedding: {red_name}")
                    except Exception:
                        pass
        except Exception as e:
            logger.warning(f"Could not extract embeddings: {e}")

        logger.info(
            f"Successfully loaded Seurat object: "
            f"{adata.n_obs} cells x {adata.n_vars} genes"
        )

        return adata

    finally:
        # Deactivate conversions
        numpy2ri.deactivate()
        pandas2ri.deactivate()


def normalize_expression(
    adata: ad.AnnData,
    target_sum: float = 1e4,
    log_transform: bool = True,
    copy: bool = True,
) -> ad.AnnData:
    """Normalize gene expression data.

    Performs library size normalization followed by optional log
    transformation.

    Parameters
    ----------
    adata : anndata.AnnData
        The data to normalize.
    target_sum : float, default=10000
        Target sum for normalization (default: 10,000 for CPM-like).
    log_transform : bool, default=True
        Whether to log-transform after normalization.
    copy : bool, default=True
        If True, return a copy. Otherwise, modify in place.

    Returns
    -------
    anndata.AnnData
        Normalized data.

    Examples
    --------
    >>> from cellmetpro.core.preprocessing import normalize_expression
    >>> adata_norm = normalize_expression(adata, target_sum=1e6)  # TPM
    """
    if copy:
        adata = adata.copy()

    # Get expression matrix
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X.copy()

    # Library size normalization
    library_sizes = X.sum(axis=1, keepdims=True)
    library_sizes[library_sizes == 0] = 1  # Avoid division by zero
    X = X / library_sizes * target_sum

    # Log transform
    if log_transform:
        X = np.log1p(X)

    adata.X = X

    # Store normalization info
    adata.uns["normalization"] = {
        "target_sum": target_sum,
        "log_transform": log_transform,
    }

    logger.info(f"Normalized expression (target_sum={target_sum}, log={log_transform})")

    return adata


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 200,
    max_genes: int | None = None,
    min_counts: int = 0,
    max_counts: int | None = None,
    copy: bool = True,
) -> ad.AnnData:
    """Filter cells based on quality metrics.

    Parameters
    ----------
    adata : anndata.AnnData
        The data to filter.
    min_genes : int, default=200
        Minimum genes expressed per cell.
    max_genes : int, optional
        Maximum genes expressed per cell.
    min_counts : int, default=0
        Minimum total counts per cell.
    max_counts : int, optional
        Maximum total counts per cell.
    copy : bool, default=True
        If True, return a copy.

    Returns
    -------
    anndata.AnnData
        Filtered data.
    """
    if copy:
        adata = adata.copy()

    # Calculate metrics
    if issparse(adata.X):
        n_genes = (adata.X > 0).sum(axis=1).A1
        total_counts = adata.X.sum(axis=1).A1
    else:
        n_genes = (adata.X > 0).sum(axis=1)
        total_counts = adata.X.sum(axis=1)

    # Apply filters
    keep = n_genes >= min_genes
    if max_genes is not None:
        keep &= n_genes <= max_genes
    if min_counts > 0:
        keep &= total_counts >= min_counts
    if max_counts is not None:
        keep &= total_counts <= max_counts

    n_removed = (~keep).sum()
    logger.info(f"Filtered {n_removed} cells ({keep.sum()} remaining)")

    return adata[keep].copy()


def filter_genes(
    adata: ad.AnnData,
    min_cells: int = 3,
    min_counts: int = 0,
    copy: bool = True,
) -> ad.AnnData:
    """Filter genes based on expression.

    Parameters
    ----------
    adata : anndata.AnnData
        The data to filter.
    min_cells : int, default=3
        Minimum cells expressing the gene.
    min_counts : int, default=0
        Minimum total counts for the gene.
    copy : bool, default=True
        If True, return a copy.

    Returns
    -------
    anndata.AnnData
        Filtered data.
    """
    if copy:
        adata = adata.copy()

    # Calculate metrics
    if issparse(adata.X):
        n_cells = (adata.X > 0).sum(axis=0).A1
        total_counts = adata.X.sum(axis=0).A1
    else:
        n_cells = (adata.X > 0).sum(axis=0)
        total_counts = adata.X.sum(axis=0)

    # Apply filters
    keep = n_cells >= min_cells
    if min_counts > 0:
        keep &= total_counts >= min_counts

    n_removed = (~keep).sum()
    logger.info(f"Filtered {n_removed} genes ({keep.sum()} remaining)")

    return adata[:, keep].copy()


def to_dataframe(
    adata: ad.AnnData,
    layer: str | None = None,
    genes_as_rows: bool = True,
) -> pd.DataFrame:
    """Convert AnnData to DataFrame.

    Parameters
    ----------
    adata : anndata.AnnData
        The data to convert.
    layer : str, optional
        Layer to use. If None, uses X.
    genes_as_rows : bool, default=True
        If True, returns genes x cells. Otherwise, cells x genes.

    Returns
    -------
    pd.DataFrame
        Expression data as DataFrame.
    """
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if issparse(X):
        X = X.toarray()

    if genes_as_rows:
        return pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)
    else:
        return pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
