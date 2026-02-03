"""Metabolic trajectory analysis module.

This module provides functions for analyzing metabolic changes along
cellular trajectories, including:
- Pseudotime inference using diffusion maps
- Metabolic velocity estimation
- Branch point detection
- Trajectory-based differential analysis

References
----------
.. [1] Haghverdi et al. (2016). Diffusion pseudotime robustly
       reconstructs lineage branching. Nature Methods.
.. [2] La Manno et al. (2018). RNA velocity of single cells.
       Nature.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.neighbors import NearestNeighbors

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


def compute_pseudotime(
    scores: pd.DataFrame,
    root_cell: str | int | None = None,
    method: str = "dpt",
    n_neighbors: int = 15,
    n_components: int = 10,
    random_state: int | None = None,
) -> pd.Series:
    """Compute pseudotime ordering of cells based on metabolic profiles.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    root_cell : str or int, optional
        Root cell for pseudotime (cell name or index).
        If None, automatically selects based on extreme PC1 value.
    method : str, optional
        Pseudotime inference method:
        - "dpt": Diffusion pseudotime
        - "principal_curve": Principal curve projection
        - "correlation": Correlation-based ordering
    n_neighbors : int, optional
        Number of neighbors for graph construction. Default 15.
    n_components : int, optional
        Number of diffusion components. Default 10.
    random_state : int, optional
        Random seed for reproducibility.

    Returns
    -------
    pd.Series
        Pseudotime values for each cell, indexed by cell names.

    Examples
    --------
    >>> pseudotime = compute_pseudotime(scores, root_cell="cell_0")
    >>> # Order cells by pseudotime
    >>> ordered_cells = pseudotime.sort_values().index

    >>> # Auto-select root cell
    >>> pseudotime = compute_pseudotime(scores)
    """
    from sklearn.decomposition import PCA

    # Transpose: cells x reactions
    X = scores.T.values
    X = np.nan_to_num(X, nan=0.0)
    cell_names = scores.columns

    # PCA reduction
    n_pcs = min(n_components, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs, random_state=random_state)
    X_pca = pca.fit_transform(X)

    # Determine root cell
    if root_cell is None:
        # Use cell with most extreme PC1 value
        root_idx = int(np.argmin(X_pca[:, 0]))
    elif isinstance(root_cell, str):
        root_idx = int(np.where(cell_names == root_cell)[0][0])
    else:
        root_idx = int(root_cell)

    if method == "dpt":
        pseudotime_values = _compute_dpt(X_pca, root_idx, n_neighbors)
    elif method == "principal_curve":
        pseudotime_values = _compute_principal_curve_pseudotime(X_pca)
    elif method == "correlation":
        pseudotime_values = _compute_correlation_pseudotime(X_pca, root_idx)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Normalize to [0, 1]
    pt_min, pt_max = pseudotime_values.min(), pseudotime_values.max()
    if pt_max > pt_min:
        pseudotime_values = (pseudotime_values - pt_min) / (pt_max - pt_min)

    return pd.Series(pseudotime_values, index=cell_names, name="pseudotime")


def _compute_dpt(
    data: np.ndarray,
    root_idx: int,
    n_neighbors: int = 15,
    n_dcs: int = 10,
) -> np.ndarray:
    """Compute diffusion pseudotime using eigendecomposition.

    This implementation follows the approach from Haghverdi et al. (2016),
    computing diffusion pseudotime via diffusion components.

    Parameters
    ----------
    data : np.ndarray
        PCA-reduced data (cells x PCs).
    root_idx : int
        Index of root cell.
    n_neighbors : int
        Number of neighbors.
    n_dcs : int
        Number of diffusion components to use.

    Returns
    -------
    np.ndarray
        Pseudotime values.
    """
    n_cells = data.shape[0]

    # Build k-NN graph
    nn = NearestNeighbors(n_neighbors=min(n_neighbors, n_cells - 1))
    nn.fit(data)
    distances, indices = nn.kneighbors(data)

    # Compute adaptive bandwidth (local sigma)
    sigma = distances[:, -1]  # Distance to k-th neighbor
    sigma[sigma == 0] = 1e-10

    # Build symmetric affinity matrix using Gaussian kernel
    W = np.zeros((n_cells, n_cells))
    for i in range(n_cells):
        for j_idx, j in enumerate(indices[i]):
            d = distances[i, j_idx]
            # Use geometric mean of sigmas (adaptive bandwidth)
            s = np.sqrt(sigma[i] * sigma[j])
            W[i, j] = np.exp(-(d**2) / (2 * s**2))
            W[j, i] = W[i, j]

    # Normalize to get transition matrix (row-stochastic)
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    T = W / row_sums

    # Compute diffusion components via eigendecomposition
    # For symmetric normalization, use: T_sym = D^{-1/2} W D^{-1/2}
    d_sqrt_inv = 1.0 / np.sqrt(row_sums.flatten() + 1e-10)
    T_sym = d_sqrt_inv[:, None] * W * d_sqrt_inv[None, :]

    try:
        # Compute top eigenvectors
        from scipy.sparse.linalg import eigsh

        n_eigs = min(n_dcs + 1, n_cells - 2)
        eigenvalues, eigenvectors = eigsh(T_sym, k=n_eigs, which="LM")

        # Sort by eigenvalue (descending)
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]

        # Skip first eigenvector (trivial, all ones)
        # Diffusion components are scaled eigenvectors
        diffusion_components = eigenvectors[:, 1 : n_dcs + 1]

        # DPT = Euclidean distance in diffusion component space from root
        root_dc = diffusion_components[root_idx]
        pseudotime = np.sqrt(((diffusion_components - root_dc) ** 2).sum(axis=1))

    except Exception:
        # Fallback to power iteration method
        T_power = T.copy()
        for _ in range(5):  # More steps for better approximation
            T_power = T_power @ T
        pseudotime = 1 - T_power[root_idx, :]

    return np.asarray(pseudotime)


def _compute_principal_curve_pseudotime(data: np.ndarray) -> np.ndarray:
    """Compute pseudotime using principal curve projection.

    This uses an iterative projection approach that considers multiple
    PCs for better trajectory estimation.

    Parameters
    ----------
    data : np.ndarray
        PCA-reduced data.

    Returns
    -------
    np.ndarray
        Pseudotime values.
    """
    n_cells = data.shape[0]

    # Start with PC1 projection
    pc1 = data[:, 0]
    order = np.argsort(pc1)

    # Compute cumulative arc-length along the trajectory
    # Use multiple PCs for distance calculation
    n_pcs_use = min(3, data.shape[1])  # Use first 3 PCs
    data_subset = data[:, :n_pcs_use]

    # Order cells by PC1 and compute cumulative distance
    ordered_data = data_subset[order]
    distances = np.zeros(n_cells)
    for i in range(1, n_cells):
        distances[i] = distances[i - 1] + np.linalg.norm(
            ordered_data[i] - ordered_data[i - 1]
        )

    # Normalize to [0, 1]
    if distances[-1] > 0:
        distances = distances / distances[-1]

    # Map back to original cell order
    pseudotime = np.zeros(n_cells)
    pseudotime[order] = distances

    return pseudotime


def _compute_correlation_pseudotime(data: np.ndarray, root_idx: int) -> np.ndarray:
    """Compute pseudotime using correlation distances from root.

    Parameters
    ----------
    data : np.ndarray
        PCA-reduced data.
    root_idx : int
        Root cell index.

    Returns
    -------
    np.ndarray
        Pseudotime values.
    """
    root_profile = data[root_idx]

    # Compute correlation with root
    correlations = np.array(
        [np.corrcoef(root_profile, data[i])[0, 1] for i in range(len(data))]
    )

    # Convert to distance
    pseudotime = 1 - correlations

    return pseudotime


def compute_metabolic_velocity(
    scores: pd.DataFrame,
    pseudotime: pd.Series,
    window_size: int = 50,
    min_cells: int = 10,
) -> pd.DataFrame:
    """Compute metabolic velocity along pseudotime trajectory.

    Metabolic velocity represents the rate of change of reaction
    activity along the pseudotime axis.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    pseudotime : pd.Series
        Pseudotime values for cells.
    window_size : int, optional
        Number of cells in sliding window. Default 50.
    min_cells : int, optional
        Minimum cells required for velocity estimation. Default 10.

    Returns
    -------
    pd.DataFrame
        Velocity estimates for each reaction at each cell.
        Same shape as scores.

    Examples
    --------
    >>> pseudotime = compute_pseudotime(scores)
    >>> velocity = compute_metabolic_velocity(scores, pseudotime)
    >>> # High velocity = rapidly changing reaction
    >>> top_dynamic = velocity.abs().mean(axis=1).nlargest(20)
    """
    # Align cells
    common_cells = scores.columns.intersection(pseudotime.index)
    scores_aligned = scores[common_cells]
    pt_aligned = pseudotime[common_cells]

    # Sort by pseudotime
    order = pt_aligned.argsort()
    scores_sorted = scores_aligned.iloc[:, order]
    pt_sorted = pt_aligned.iloc[order]

    n_reactions, n_cells = scores_sorted.shape
    velocities = np.zeros((n_reactions, n_cells))

    half_window = window_size // 2

    for i in range(n_cells):
        # Define window
        start = max(0, i - half_window)
        end = min(n_cells, i + half_window + 1)

        if end - start < min_cells:
            velocities[:, i] = np.nan
            continue

        # Get window data
        pt_window = pt_sorted.iloc[start:end].values
        scores_window = scores_sorted.iloc[:, start:end].values

        # Compute velocity as slope of linear regression
        pt_centered = pt_window - pt_window.mean()
        pt_var = (pt_centered**2).sum()

        if pt_var > 0:
            for r in range(n_reactions):
                scores_centered = scores_window[r] - scores_window[r].mean()
                velocities[r, i] = (scores_centered * pt_centered).sum() / pt_var

    # Restore original order
    velocity_df = pd.DataFrame(
        velocities,
        index=scores_sorted.index,
        columns=scores_sorted.columns,
    )

    # Reorder to original cell order
    velocity_df = velocity_df[scores.columns.intersection(velocity_df.columns)]

    return velocity_df


def identify_branch_points(
    scores: pd.DataFrame,
    pseudotime: pd.Series,
    n_neighbors: int = 15,
    threshold: float = 0.3,
) -> pd.DataFrame:
    """Identify potential branch points in the trajectory.

    Branch points are locations where cells diverge into multiple
    lineages.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    pseudotime : pd.Series
        Pseudotime values for cells.
    n_neighbors : int, optional
        Number of neighbors for local analysis. Default 15.
    threshold : float, optional
        Threshold for branch point detection. Default 0.3.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - cell: Cell identifier
        - pseudotime: Pseudotime value
        - branch_score: Score indicating branch likelihood
        - is_branch: Boolean flag for branch points

    Examples
    --------
    >>> branches = identify_branch_points(scores, pseudotime)
    >>> branch_cells = branches[branches['is_branch']]['cell']
    """
    from sklearn.decomposition import PCA

    # Align and prepare data
    common_cells = scores.columns.intersection(pseudotime.index)
    X = scores[common_cells].T.values
    X = np.nan_to_num(X, nan=0.0)
    pt = pseudotime[common_cells].values
    cell_names = list(common_cells)

    # PCA reduction
    n_pcs = min(20, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_pcs)
    X_pca = pca.fit_transform(X)

    # Build k-NN graph
    nn = NearestNeighbors(n_neighbors=min(n_neighbors, len(X) - 1))
    nn.fit(X_pca)
    _, indices = nn.kneighbors(X_pca)

    # Compute branch scores
    branch_scores = np.zeros(len(X))

    for i in range(len(X)):
        neighbor_pts = pt[indices[i]]

        # Branch point: neighbors have divergent pseudotimes
        pt_std = np.std(neighbor_pts)

        # Also check for local density decrease
        distances = np.linalg.norm(X_pca[indices[i]] - X_pca[i], axis=1)
        avg_dist = np.mean(distances)

        # Combine metrics
        branch_scores[i] = pt_std * avg_dist

    # Normalize scores
    if branch_scores.max() > 0:
        branch_scores = branch_scores / branch_scores.max()

    # Identify branch points
    is_branch = branch_scores > threshold

    return pd.DataFrame(
        {
            "cell": cell_names,
            "pseudotime": pt,
            "branch_score": branch_scores,
            "is_branch": is_branch,
        }
    )


def trajectory_differential(
    scores: pd.DataFrame,
    pseudotime: pd.Series,
    n_bins: int = 10,
    method: str = "spearman",
) -> pd.DataFrame:
    """Identify reactions that change significantly along trajectory.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    pseudotime : pd.Series
        Pseudotime values for cells.
    n_bins : int, optional
        Number of bins for pseudotime (for trend analysis). Default 10.
    method : str, optional
        Correlation method: "spearman" or "pearson". Default "spearman".

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - reaction: Reaction ID
        - correlation: Correlation with pseudotime
        - pvalue: Statistical significance
        - trend: "increasing", "decreasing", or "non-monotonic"
        - early_mean: Mean in first pseudotime bin
        - late_mean: Mean in last pseudotime bin
        - log2fc: Log2 fold change (late/early)

    Examples
    --------
    >>> diff_results = trajectory_differential(scores, pseudotime)
    >>> increasing = diff_results[diff_results['trend'] == 'increasing']
    >>> decreasing = diff_results[diff_results['trend'] == 'decreasing']
    """
    # Align data
    common_cells = scores.columns.intersection(pseudotime.index)
    scores_aligned = scores[common_cells]
    pt_aligned = pseudotime[common_cells]

    results = []

    for rxn in scores_aligned.index:
        rxn_values = scores_aligned.loc[rxn].values
        pt_values = pt_aligned.values

        # Remove NaN values
        mask = ~np.isnan(rxn_values)
        if mask.sum() < 10:
            continue

        rxn_clean = rxn_values[mask]
        pt_clean = pt_values[mask]

        # Compute correlation
        if method == "spearman":
            corr, pval = stats.spearmanr(pt_clean, rxn_clean)
        else:
            corr, pval = stats.pearsonr(pt_clean, rxn_clean)

        # Bin cells by pseudotime
        bins = pd.cut(pt_clean, bins=n_bins, labels=False)

        # Early vs late comparison
        early_mask = bins == 0
        late_mask = bins == n_bins - 1

        early_mean = rxn_clean[early_mask].mean() if early_mask.sum() > 0 else np.nan
        late_mean = rxn_clean[late_mask].mean() if late_mask.sum() > 0 else np.nan

        # Log2 fold change
        if early_mean > 0 and late_mean > 0:
            log2fc = np.log2(late_mean / early_mean)
        else:
            log2fc = np.nan

        # Determine trend
        if pval < 0.05:
            if corr > 0.3:
                trend = "increasing"
            elif corr < -0.3:
                trend = "decreasing"
            else:
                trend = "non-monotonic"
        else:
            trend = "non-monotonic"

        results.append(
            {
                "reaction": rxn,
                "correlation": corr,
                "pvalue": pval,
                "trend": trend,
                "early_mean": early_mean,
                "late_mean": late_mean,
                "log2fc": log2fc,
            }
        )

    df = pd.DataFrame(results)

    # Adjust p-values (Benjamini-Hochberg)
    if len(df) > 0:
        from scipy.stats import rankdata

        n = len(df)
        ranked_pvals = rankdata(df["pvalue"])
        df["padj"] = df["pvalue"] * n / ranked_pvals
        df["padj"] = df["padj"].clip(upper=1.0)

    return df.sort_values("pvalue")


def fit_trajectory_genes(
    scores: pd.DataFrame,
    pseudotime: pd.Series,
    reactions: list[str] | None = None,
    n_knots: int = 5,
) -> dict[str, Any]:
    """Fit smooth trajectory curves for reactions.

    Uses spline fitting to model reaction activity along pseudotime.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).
    pseudotime : pd.Series
        Pseudotime values for cells.
    reactions : list[str], optional
        Reactions to fit. If None, fits top variable reactions.
    n_knots : int, optional
        Number of spline knots. Default 5.

    Returns
    -------
    dict
        Dictionary with:
        - 'fitted_values': DataFrame of fitted values
        - 'pseudotime_grid': Array of pseudotime evaluation points
        - 'r_squared': R-squared values for each reaction

    Examples
    --------
    >>> fit = fit_trajectory_genes(scores, pseudotime, n_knots=4)
    >>> fitted_values = fit['fitted_values']
    """
    from scipy.interpolate import UnivariateSpline

    # Align data
    common_cells = scores.columns.intersection(pseudotime.index)
    scores_aligned = scores[common_cells]
    pt_aligned = pseudotime[common_cells]

    # Select reactions
    if reactions is None:
        # Top variable reactions
        variance = scores_aligned.var(axis=1)
        reactions = variance.nlargest(50).index.tolist()

    # Sort by pseudotime
    order = pt_aligned.argsort()
    pt_sorted = pt_aligned.iloc[order].values
    scores_sorted = scores_aligned.iloc[:, order]

    # Create evaluation grid
    pt_grid = np.linspace(pt_sorted.min(), pt_sorted.max(), 100)

    fitted_values = {}
    r_squared = {}

    for rxn in reactions:
        y = scores_sorted.loc[rxn].values

        # Remove NaN
        mask = ~np.isnan(y)
        if mask.sum() < n_knots + 1:
            continue

        x_clean = pt_sorted[mask]
        y_clean = y[mask]

        try:
            # Fit smoothing spline
            spline = UnivariateSpline(x_clean, y_clean, k=3, s=len(x_clean))
            y_fitted = spline(pt_grid)
            fitted_values[rxn] = y_fitted

            # Compute R-squared
            y_pred = spline(x_clean)
            ss_res = np.sum((y_clean - y_pred) ** 2)
            ss_tot = np.sum((y_clean - y_clean.mean()) ** 2)
            r_squared[rxn] = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        except Exception as e:
            logger.warning(f"Failed to fit {rxn}: {e}")

    return {
        "fitted_values": pd.DataFrame(fitted_values, index=pt_grid),
        "pseudotime_grid": pt_grid,
        "r_squared": pd.Series(r_squared),
    }


def cluster_trajectory_patterns(
    trajectory_fit: dict[str, Any],
    n_clusters: int = 5,
) -> pd.Series:
    """Cluster reactions by their trajectory patterns.

    Parameters
    ----------
    trajectory_fit : dict
        Output from fit_trajectory_genes().
    n_clusters : int, optional
        Number of pattern clusters. Default 5.

    Returns
    -------
    pd.Series
        Cluster assignments for each reaction.

    Examples
    --------
    >>> fit = fit_trajectory_genes(scores, pseudotime)
    >>> clusters = cluster_trajectory_patterns(fit, n_clusters=4)
    >>> early_activated = clusters[clusters == 0].index
    """
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler

    fitted_values = trajectory_fit["fitted_values"]

    # Standardize patterns
    scaler = StandardScaler()
    X = scaler.fit_transform(fitted_values.T.values)

    # Cluster
    kmeans = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
    labels = kmeans.fit_predict(X)

    return pd.Series(labels, index=fitted_values.columns, name="pattern_cluster")
