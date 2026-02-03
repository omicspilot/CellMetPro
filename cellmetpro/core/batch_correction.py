"""Batch correction methods for multi-sample integration.

This module provides functions for correcting batch effects when
integrating metabolic profiles from multiple samples or experiments.

Methods available:
- Harmony: Fast integration using iterative clustering
- ComBat: Empirical Bayes batch correction
- Simple centering: Per-batch mean centering

References
----------
.. [1] Korsunsky et al. (2019). Fast, sensitive and accurate integration
       of single-cell data with Harmony. Nature Methods.
.. [2] Johnson et al. (2007). Adjusting batch effects in microarray
       expression data using empirical Bayes methods. Biostatistics.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    pass


def harmony_integrate(
    scores: pd.DataFrame,
    batch_labels: pd.Series | np.ndarray | list[str],
    n_components: int = 50,
    max_iter: int = 10,
    theta: float = 2.0,
    sigma: float = 0.1,
    random_state: int | None = None,
) -> pd.DataFrame:
    """Integrate batches using Harmony algorithm.

    Harmony iteratively clusters cells and adjusts the data to remove
    batch-specific variation while preserving biological variation.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    batch_labels : array-like
        Batch labels for each cell (length = n_cells).
    n_components : int, optional
        Number of PCA components to use. Default 50.
    max_iter : int, optional
        Maximum iterations for Harmony. Default 10.
    theta : float, optional
        Diversity clustering penalty. Higher values = more integration.
        Default 2.0.
    sigma : float, optional
        Width of soft kmeans clusters. Default 0.1.
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    pd.DataFrame
        Batch-corrected scores matrix.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> scores = pd.DataFrame(np.random.rand(100, 500))
    >>> batches = ['batch1'] * 250 + ['batch2'] * 250
    >>> corrected = harmony_integrate(scores, batches)

    Notes
    -----
    Requires the ``harmonypy`` package: ``pip install harmonypy``
    """
    try:
        import harmonypy as hm
    except ImportError as e:
        raise ImportError(
            "harmonypy is required for Harmony integration. "
            "Install it with: pip install harmonypy"
        ) from e

    # Ensure batch_labels is array-like
    if isinstance(batch_labels, pd.Series):
        batch_labels = batch_labels.values
    batch_labels = np.asarray(batch_labels)

    if len(batch_labels) != scores.shape[1]:
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must match "
            f"number of cells ({scores.shape[1]})"
        )

    # Validate we have multiple batches
    n_batches = len(np.unique(batch_labels))
    if n_batches < 2:
        import warnings

        warnings.warn(
            f"Only {n_batches} batch(es) found. Returning original data unchanged.",
            UserWarning,
            stacklevel=2,
        )
        return scores.copy()

    # Transpose: cells x reactions for PCA
    X = scores.T.values

    # Handle NaN values
    X = np.nan_to_num(X, nan=0.0)

    # PCA reduction
    from sklearn.decomposition import PCA

    n_pcs = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs, random_state=random_state)
    X_pca = pca.fit_transform(X)

    # Create metadata DataFrame for Harmony
    meta_data = pd.DataFrame({"batch": batch_labels})

    # Run Harmony
    ho = hm.run_harmony(
        X_pca,
        meta_data,
        "batch",
        max_iter_harmony=max_iter,
        theta=theta,
        sigma=sigma,
        random_state=random_state,
    )

    # Get corrected PCA embeddings
    X_corrected_pca = ho.Z_corr.T

    # Reconstruct full space using inverse PCA transform
    X_corrected = pca.inverse_transform(X_corrected_pca)

    # Return as DataFrame with original structure
    return pd.DataFrame(X_corrected.T, index=scores.index, columns=scores.columns)


def combat_correct(
    scores: pd.DataFrame,
    batch_labels: pd.Series | np.ndarray | list[str],
    covariates: pd.DataFrame | None = None,
    parametric: bool = True,
) -> pd.DataFrame:
    """Apply ComBat batch correction.

    ComBat uses an empirical Bayes framework to adjust for batch effects
    while preserving biological variation.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    batch_labels : array-like
        Batch labels for each cell.
    covariates : pd.DataFrame, optional
        Covariates to preserve (e.g., cell type). Rows are cells.
    parametric : bool, optional
        Use parametric (True) or non-parametric (False) estimation.
        Default True.

    Returns
    -------
    pd.DataFrame
        Batch-corrected scores matrix.

    Examples
    --------
    >>> corrected = combat_correct(scores, batch_labels)

    Notes
    -----
    This is a simplified implementation. For production use, consider
    using pycombat or inmoose packages.
    """
    # Ensure batch_labels is array
    if isinstance(batch_labels, pd.Series):
        batch_labels = batch_labels.values
    batch_labels = np.asarray(batch_labels)

    if len(batch_labels) != scores.shape[1]:
        raise ValueError(
            f"batch_labels length ({len(batch_labels)}) must match "
            f"number of cells ({scores.shape[1]})"
        )

    # Get unique batches
    batches = np.unique(batch_labels)
    n_batches = len(batches)

    if n_batches < 2:
        return scores.copy()

    # Work with transposed data: cells x reactions
    X = scores.T.values.copy()
    X = np.nan_to_num(X, nan=0.0)

    # Standardize
    grand_mean = X.mean(axis=0)
    grand_var = X.var(axis=0)
    grand_var[grand_var == 0] = 1  # Avoid division by zero

    # Compute batch-specific means and variances
    batch_means = np.zeros((n_batches, X.shape[1]))
    batch_vars = np.zeros((n_batches, X.shape[1]))
    batch_sizes = np.zeros(n_batches)

    for i, batch in enumerate(batches):
        mask = batch_labels == batch
        batch_sizes[i] = mask.sum()
        batch_means[i] = X[mask].mean(axis=0)
        batch_vars[i] = X[mask].var(axis=0)

    # Standardize each feature
    X_std = (X - grand_mean) / np.sqrt(grand_var)

    if parametric:
        # Parametric ComBat
        # Estimate gamma (location) and delta (scale) for each batch
        gamma_hat = batch_means - grand_mean
        delta_hat = batch_vars / grand_var
        delta_hat[delta_hat == 0] = 1

        # Empirical Bayes shrinkage (simplified)
        gamma_bar = gamma_hat.mean(axis=0)
        delta_bar = delta_hat.mean(axis=0)

        # Apply correction
        X_corrected = X_std.copy()
        for i, batch in enumerate(batches):
            mask = batch_labels == batch
            # Adjust for batch effect
            gamma_star = 0.5 * gamma_hat[i] + 0.5 * gamma_bar  # Shrinkage
            delta_star = 0.5 * delta_hat[i] + 0.5 * delta_bar

            X_corrected[mask] = (
                X_std[mask] - gamma_star / np.sqrt(grand_var)
            ) / np.sqrt(delta_star)

    else:
        # Non-parametric: simple mean/variance centering
        X_corrected = X_std.copy()
        for i, batch in enumerate(batches):
            mask = batch_labels == batch
            batch_data = X_std[mask]
            batch_mean = batch_data.mean(axis=0)
            batch_std = batch_data.std(axis=0)
            batch_std[batch_std == 0] = 1
            X_corrected[mask] = (batch_data - batch_mean) / batch_std

    # Rescale to original scale
    X_final = X_corrected * np.sqrt(grand_var) + grand_mean

    return pd.DataFrame(X_final.T, index=scores.index, columns=scores.columns)


def center_batches(
    scores: pd.DataFrame,
    batch_labels: pd.Series | np.ndarray | list[str],
) -> pd.DataFrame:
    """Simple batch centering by subtracting batch means.

    This is the simplest form of batch correction, which removes
    the mean of each batch to center all batches at zero.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    batch_labels : array-like
        Batch labels for each cell.

    Returns
    -------
    pd.DataFrame
        Centered scores matrix.

    Examples
    --------
    >>> centered = center_batches(scores, batch_labels)
    """
    if isinstance(batch_labels, pd.Series):
        batch_labels = batch_labels.values
    batch_labels = np.asarray(batch_labels)

    X = scores.values.copy()
    batches = np.unique(batch_labels)

    for batch in batches:
        mask = batch_labels == batch
        batch_mean = np.nanmean(X[:, mask], axis=1, keepdims=True)
        X[:, mask] = X[:, mask] - batch_mean

    # Add back grand mean
    grand_mean = np.nanmean(scores.values, axis=1, keepdims=True)
    X = X + grand_mean

    return pd.DataFrame(X, index=scores.index, columns=scores.columns)


def compute_integration_metrics(
    scores: pd.DataFrame,
    batch_labels: pd.Series | np.ndarray | list[str],
    cell_labels: pd.Series | np.ndarray | list[str] | None = None,
    n_neighbors: int = 15,
) -> dict[str, float]:
    """Compute integration quality metrics.

    Calculates metrics to assess how well batches are integrated:
    - Batch mixing score: How well batches mix in the embedding
    - Silhouette score: How well batches remain separated (lower = better mixing)
    - ARI with true labels: If cell labels provided, assess label preservation

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    batch_labels : array-like
        Batch labels for each cell.
    cell_labels : array-like, optional
        True cell type labels for assessing preservation.
    n_neighbors : int, optional
        Number of neighbors for local metrics. Default 15.

    Returns
    -------
    dict
        Dictionary containing:
        - 'batch_mixing': Batch mixing entropy (higher = better)
        - 'batch_silhouette': Batch silhouette score (lower = better mixing)
        - 'label_ari': ARI if cell_labels provided (higher = better preservation)

    Examples
    --------
    >>> metrics = compute_integration_metrics(scores, batch_labels)
    >>> print(f"Batch mixing: {metrics['batch_mixing']:.3f}")
    """
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import NearestNeighbors

    # Ensure arrays
    if isinstance(batch_labels, pd.Series):
        batch_labels = batch_labels.values
    batch_labels = np.asarray(batch_labels)

    # Transpose for cell-centric analysis
    X = scores.T.values
    X = np.nan_to_num(X, nan=0.0)

    # PCA reduction for efficiency
    from sklearn.decomposition import PCA

    n_pcs = min(50, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs)
    X_pca = pca.fit_transform(X)

    metrics: dict[str, float] = {}

    # Batch silhouette score (lower = better mixing)
    try:
        metrics["batch_silhouette"] = float(
            silhouette_score(
                X_pca, batch_labels, sample_size=min(5000, len(batch_labels))
            )
        )
    except ValueError:
        metrics["batch_silhouette"] = 0.0

    # Batch mixing entropy using k-NN
    try:
        nn = NearestNeighbors(n_neighbors=n_neighbors)
        nn.fit(X_pca)
        _, indices = nn.kneighbors(X_pca)

        batches = np.unique(batch_labels)
        n_batches = len(batches)

        # Compute local batch entropy for each cell
        entropies = []
        for i in range(len(X_pca)):
            neighbor_batches = batch_labels[indices[i]]
            _, counts = np.unique(neighbor_batches, return_counts=True)
            probs = counts / counts.sum()
            entropy = -np.sum(probs * np.log(probs + 1e-10))
            entropies.append(entropy)

        # Normalize by maximum entropy
        max_entropy = np.log(n_batches)
        metrics["batch_mixing"] = float(
            np.mean(entropies) / max_entropy if max_entropy > 0 else 0
        )
    except Exception:
        metrics["batch_mixing"] = 0.0

    # Label preservation (if cell labels provided)
    if cell_labels is not None:
        from sklearn.cluster import KMeans
        from sklearn.metrics import adjusted_rand_score

        if isinstance(cell_labels, pd.Series):
            cell_labels = cell_labels.values
        cell_labels = np.asarray(cell_labels)

        n_clusters = len(np.unique(cell_labels))
        kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
        predicted = kmeans.fit_predict(X_pca)

        metrics["label_ari"] = float(adjusted_rand_score(cell_labels, predicted))

    return metrics


def select_hvr(
    scores: pd.DataFrame,
    n_top: int = 2000,
    min_mean: float = 0.01,
    max_mean: float = 5.0,
    min_disp: float = 0.5,
) -> list[str]:
    """Select highly variable reactions for integration.

    Similar to highly variable gene selection in scRNA-seq,
    identifies reactions with high variance relative to their mean.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    n_top : int, optional
        Number of top HVRs to select. Default 2000.
    min_mean : float, optional
        Minimum mean expression. Default 0.01.
    max_mean : float, optional
        Maximum mean expression. Default 5.0.
    min_disp : float, optional
        Minimum normalized dispersion. Default 0.5.

    Returns
    -------
    list
        List of highly variable reaction names.

    Examples
    --------
    >>> hvr = select_hvr(scores, n_top=1000)
    >>> scores_filtered = scores.loc[hvr]
    """
    mean = scores.mean(axis=1)
    var = scores.var(axis=1)

    # Normalized dispersion - use absolute mean to handle negative values
    # This treats reactions with large magnitude means (+ or -) similarly
    mean_abs = np.abs(mean)
    mean_abs[mean_abs == 0] = 1e-10
    disp = var / mean_abs

    # Log transform for binning (use absolute mean for binning)
    mean_log = np.log1p(mean_abs)

    # Bin genes by mean expression
    n_bins = 20
    bins = pd.cut(mean_log, bins=n_bins, labels=False)

    # Normalize dispersion within each bin
    disp_norm = disp.copy()
    for b in range(n_bins):
        mask = bins == b
        if mask.sum() > 1:
            bin_mean = disp[mask].mean()
            bin_std = disp[mask].std()
            if bin_std > 0:
                disp_norm[mask] = (disp[mask] - bin_mean) / bin_std

    # Filter by mean and dispersion
    valid = (mean >= min_mean) & (mean <= max_mean) & (disp_norm >= min_disp)
    hvr_candidates = scores.index[valid]

    # Select top N by normalized dispersion
    hvr = list(disp_norm[hvr_candidates].nlargest(n_top).index)

    return hvr
