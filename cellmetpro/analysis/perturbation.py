"""Perturbation prediction and simulation module.

This module provides functions for simulating and predicting the effects
of genetic and metabolic perturbations on cellular metabolism.

Features:
- Gene expression perturbation simulation
- Multi-gene knockout analysis
- Flux distribution comparison
- Essential reaction/gene identification
- Perturbation sensitivity analysis
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    import cobra

logger = logging.getLogger(__name__)


def simulate_expression_change(
    model: cobra.Model,
    gene_id: str,
    fold_change: float,
    method: str = "linear",
) -> pd.Series:
    """Simulate the effect of changing gene expression on fluxes.

    This function models the effect of gene expression changes by
    adjusting the flux bounds of reactions associated with the gene.

    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model.
    gene_id : str
        Gene ID to perturb.
    fold_change : float
        Expression fold change (e.g., 2.0 = 2x increase, 0.5 = 50% decrease).
    method : str, optional
        Method for mapping expression to flux:
        - "linear": Scale bounds linearly with expression
        - "log": Scale bounds with log(expression)
        - "threshold": Binary on/off based on fold_change

    Returns
    -------
    pd.Series
        Predicted flux distribution after perturbation.

    Examples
    --------
    >>> # Simulate 2x overexpression of a gene
    >>> fluxes = simulate_expression_change(model, "gene123", 2.0)

    >>> # Simulate 50% knockdown
    >>> fluxes = simulate_expression_change(model, "gene456", 0.5)

    Raises
    ------
    ValueError
        If gene_id is not found in the model.
    """

    model = model.copy()

    try:
        gene = model.genes.get_by_id(gene_id)
    except KeyError as e:
        raise ValueError(f"Gene '{gene_id}' not found in model") from e

    # Get reactions associated with this gene
    associated_reactions = list(gene.reactions)

    if not associated_reactions:
        logger.warning(f"Gene '{gene_id}' has no associated reactions")
        try:
            return model.optimize().fluxes
        except Exception:
            return pd.Series(0.0, index=[r.id for r in model.reactions])

    for rxn in associated_reactions:
        # Get original bounds
        lb, ub = rxn.lower_bound, rxn.upper_bound

        if method == "linear":
            # Scale bounds linearly
            scale = fold_change
        elif method == "log":
            # Logarithmic scaling (softer effect)
            scale = np.log2(fold_change + 1) / np.log2(2) if fold_change > 0 else 0
        elif method == "threshold":
            # Binary: if fold_change < 0.1, knock out; otherwise keep
            scale = 0.0 if fold_change < 0.1 else 1.0
        else:
            raise ValueError(f"Unknown method: {method}")

        # Apply scaling to bounds
        if scale == 0:
            rxn.knock_out()
        else:
            if lb < 0:  # Reversible or reverse-only
                rxn.lower_bound = lb * scale
            if ub > 0:  # Forward-capable
                rxn.upper_bound = ub * scale

    # Optimize and return fluxes
    try:
        solution = model.optimize()
        if solution.status == "optimal":
            return solution.fluxes
        else:
            logger.warning(f"Optimization status: {solution.status}")
            return pd.Series(0.0, index=[r.id for r in model.reactions])
    except Exception as e:
        logger.error(f"Optimization failed: {e}")
        return pd.Series(0.0, index=[r.id for r in model.reactions])


def multi_knockout(
    model: cobra.Model,
    genes: list[str],
    method: str = "simultaneous",
    return_intermediate: bool = False,
) -> pd.DataFrame | pd.Series:
    """Perform multi-gene knockout analysis.

    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model.
    genes : list[str]
        List of gene IDs to knock out.
    method : str, optional
        Knockout method:
        - "simultaneous": Knock out all genes at once
        - "sequential": Knock out genes one by one, keeping previous KOs
        - "individual": Knock out each gene independently (for comparison)
    return_intermediate : bool, optional
        If True and method="sequential", return intermediate states.

    Returns
    -------
    pd.DataFrame or pd.Series
        If method="individual" or return_intermediate=True: DataFrame with
        flux distributions for each knockout.
        Otherwise: Series with final flux distribution.

    Examples
    --------
    >>> # Simultaneous double knockout
    >>> fluxes = multi_knockout(model, ["gene1", "gene2"])

    >>> # Compare individual effects
    >>> all_fluxes = multi_knockout(model, ["gene1", "gene2"], method="individual")
    """

    if method == "simultaneous":
        model_copy = model.copy()
        for gene_id in genes:
            try:
                gene = model_copy.genes.get_by_id(gene_id)
                gene.knock_out()
            except KeyError:
                logger.warning(f"Gene '{gene_id}' not found, skipping")

        try:
            solution = model_copy.optimize()
            if solution.status == "optimal":
                return solution.fluxes
            else:
                return pd.Series(0.0, index=[r.id for r in model.reactions])
        except Exception as e:
            logger.error(f"Optimization failed: {e}")
            return pd.Series(0.0, index=[r.id for r in model.reactions])

    elif method == "sequential":
        results = {}
        model_copy = model.copy()

        for i, gene_id in enumerate(genes):
            try:
                gene = model_copy.genes.get_by_id(gene_id)
                gene.knock_out()
            except KeyError:
                logger.warning(f"Gene '{gene_id}' not found, skipping")
                continue

            try:
                solution = model_copy.optimize()
                if solution.status == "optimal":
                    results[f"KO_{i+1}_{gene_id}"] = solution.fluxes
                else:
                    results[f"KO_{i+1}_{gene_id}"] = pd.Series(
                        0.0, index=[r.id for r in model.reactions]
                    )
            except Exception:
                results[f"KO_{i+1}_{gene_id}"] = pd.Series(
                    0.0, index=[r.id for r in model.reactions]
                )

        if return_intermediate:
            return pd.DataFrame(results)
        else:
            # Return last state
            last_key = list(results.keys())[-1] if results else None
            if last_key:
                return results[last_key]
            return pd.Series(0.0, index=[r.id for r in model.reactions])

    elif method == "individual":
        results = {}

        # Wild-type reference
        try:
            wt_solution = model.optimize()
            if wt_solution.status == "optimal":
                results["wild_type"] = wt_solution.fluxes
        except Exception:
            pass

        for gene_id in genes:
            model_copy = model.copy()
            try:
                gene = model_copy.genes.get_by_id(gene_id)
                gene.knock_out()
                solution = model_copy.optimize()
                if solution.status == "optimal":
                    results[f"KO_{gene_id}"] = solution.fluxes
                else:
                    results[f"KO_{gene_id}"] = pd.Series(
                        0.0, index=[r.id for r in model.reactions]
                    )
            except KeyError:
                logger.warning(f"Gene '{gene_id}' not found")
            except Exception as e:
                logger.error(f"Failed for {gene_id}: {e}")

        return pd.DataFrame(results)

    else:
        raise ValueError(f"Unknown method: {method}")


def compare_flux_distributions(
    flux1: pd.Series,
    flux2: pd.Series,
    threshold: float = 0.01,
    method: str = "absolute",
) -> pd.DataFrame:
    """Compare two flux distributions.

    Parameters
    ----------
    flux1 : pd.Series
        First flux distribution (e.g., wild-type).
    flux2 : pd.Series
        Second flux distribution (e.g., knockout).
    threshold : float, optional
        Minimum flux magnitude to consider. Default 0.01.
    method : str, optional
        Comparison method:
        - "absolute": Absolute difference
        - "relative": Relative (fold) change
        - "both": Both metrics

    Returns
    -------
    pd.DataFrame
        DataFrame with columns for flux values and differences.

    Examples
    --------
    >>> wt_flux = fba.optimize()
    >>> ko_flux = fba.knockout("reaction1")
    >>> comparison = compare_flux_distributions(wt_flux, ko_flux)
    """
    # Align indices
    common_rxns = flux1.index.intersection(flux2.index)
    f1 = flux1[common_rxns]
    f2 = flux2[common_rxns]

    result = pd.DataFrame(
        {
            "flux_1": f1,
            "flux_2": f2,
        }
    )

    # Absolute difference
    result["abs_diff"] = f2 - f1

    if method in ["relative", "both"]:
        # Relative change (avoid division by zero)
        with np.errstate(divide="ignore", invalid="ignore"):
            fc = f2 / f1
            fc = np.where(np.isfinite(fc), fc, np.nan)
        result["fold_change"] = fc

        # Log2 fold change: handle flux direction changes properly
        # For positive fluxes: standard log2(f2/f1)
        # For sign changes: mark as inf or use signed magnitude
        log2_fc = np.full_like(fc, np.nan, dtype=float)

        # Both positive
        both_pos = (f1 > threshold) & (f2 > threshold)
        log2_fc = np.where(both_pos, np.log2(f2 / f1), log2_fc)

        # Both negative (use absolute values, negate result)
        both_neg = (f1 < -threshold) & (f2 < -threshold)
        log2_fc = np.where(both_neg, np.log2(np.abs(f2) / np.abs(f1)), log2_fc)

        # Sign change cases: use signed difference indicator
        sign_change = (f1 > threshold) & (f2 < -threshold) | (f1 < -threshold) & (
            f2 > threshold
        )
        log2_fc = np.where(sign_change, np.inf * np.sign(f2 - f1), log2_fc)

        result["log2_fc"] = log2_fc

    # Filter by threshold
    significant = (np.abs(f1) > threshold) | (np.abs(f2) > threshold)
    result["significant"] = significant

    # Sort by absolute difference
    result = result.sort_values("abs_diff", key=abs, ascending=False)

    return result


def identify_essential_genes(
    model: cobra.Model,
    genes: list[str] | None = None,
    threshold: float = 0.01,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Identify essential genes through systematic knockout.

    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model.
    genes : list[str], optional
        Genes to test. If None, tests all genes.
    threshold : float, optional
        Minimum fraction of wild-type growth for viability.
    n_jobs : int, optional
        Number of parallel jobs. Default 1 (serial).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - gene: Gene ID
        - essential: Boolean indicating essentiality
        - ko_growth: Growth rate after knockout
        - growth_ratio: Ratio of KO growth to WT growth

    Examples
    --------
    >>> essential_df = identify_essential_genes(model, threshold=0.05)
    >>> essential_genes = essential_df[essential_df['essential']]['gene'].tolist()
    """

    # Get wild-type growth
    wt_solution = model.optimize()
    if wt_solution.status != "optimal":
        raise ValueError("Wild-type optimization failed")
    wt_growth = wt_solution.objective_value

    if genes is None:
        genes = [g.id for g in model.genes]

    def test_gene(gene_id: str) -> dict[str, Any]:
        model_copy = model.copy()
        try:
            gene = model_copy.genes.get_by_id(gene_id)
            gene.knock_out()
            solution = model_copy.optimize()

            if solution.status == "optimal":
                ko_growth = solution.objective_value
                growth_ratio = ko_growth / wt_growth if wt_growth > 0 else 0
                is_essential = growth_ratio < threshold
            else:
                ko_growth = 0.0
                growth_ratio = 0.0
                is_essential = True

        except Exception:
            ko_growth = 0.0
            growth_ratio = 0.0
            is_essential = True

        return {
            "gene": gene_id,
            "essential": is_essential,
            "ko_growth": ko_growth,
            "growth_ratio": growth_ratio,
        }

    if n_jobs == 1:
        results = [test_gene(g) for g in genes]
    else:
        try:
            from joblib import Parallel, delayed

            results = Parallel(n_jobs=n_jobs)(delayed(test_gene)(g) for g in genes)
        except ImportError:
            logger.warning("joblib not available, running serially")
            results = [test_gene(g) for g in genes]

    return pd.DataFrame(results)


def sensitivity_analysis(
    model: cobra.Model,
    reactions: list[str] | None = None,
    perturbation_range: tuple[float, float] = (0.5, 2.0),
    n_steps: int = 10,
) -> pd.DataFrame:
    """Perform flux sensitivity analysis.

    Evaluates how objective value changes when reaction bounds
    are systematically varied.

    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model.
    reactions : list[str], optional
        Reactions to analyze. If None, analyzes top 20 by flux.
    perturbation_range : tuple, optional
        (min_scale, max_scale) for bound perturbation.
    n_steps : int, optional
        Number of perturbation steps.

    Returns
    -------
    pd.DataFrame
        Sensitivity coefficients and statistics for each reaction.

    Examples
    --------
    >>> sensitivity = sensitivity_analysis(model, n_steps=5)
    >>> most_sensitive = sensitivity.nlargest(10, 'sensitivity')
    """

    # Get baseline
    baseline = model.optimize()
    if baseline.status != "optimal":
        raise ValueError("Baseline optimization failed")
    baseline_obj = baseline.objective_value

    if reactions is None:
        # Select top reactions by flux
        reactions = baseline.fluxes.abs().nlargest(20).index.tolist()

    scales = np.linspace(perturbation_range[0], perturbation_range[1], n_steps)
    results = []

    for rxn_id in reactions:
        obj_values = []

        for scale in scales:
            model_copy = model.copy()
            rxn = model_copy.reactions.get_by_id(rxn_id)

            # Scale bounds
            orig_lb, orig_ub = rxn.lower_bound, rxn.upper_bound
            if orig_lb != 0:
                rxn.lower_bound = orig_lb * scale
            if orig_ub != 0:
                rxn.upper_bound = orig_ub * scale

            try:
                sol = model_copy.optimize()
                if sol.status == "optimal":
                    obj_values.append(sol.objective_value)
                else:
                    obj_values.append(np.nan)
            except Exception:
                obj_values.append(np.nan)

        obj_arr = np.array(obj_values)

        # Compute sensitivity (slope of objective vs scale)
        valid_mask = ~np.isnan(obj_arr)
        if valid_mask.sum() >= 2:
            slope, _ = np.polyfit(scales[valid_mask], obj_arr[valid_mask], 1)
            sensitivity = slope / baseline_obj if baseline_obj != 0 else 0
        else:
            sensitivity = np.nan

        results.append(
            {
                "reaction": rxn_id,
                "baseline_flux": baseline.fluxes[rxn_id],
                "sensitivity": sensitivity,
                "obj_min": np.nanmin(obj_arr),
                "obj_max": np.nanmax(obj_arr),
                "obj_range": np.nanmax(obj_arr) - np.nanmin(obj_arr),
            }
        )

    return pd.DataFrame(results).set_index("reaction")


def predict_synthetic_lethality(
    model: cobra.Model,
    genes: list[str],
    threshold: float = 0.01,
    max_combinations: int = 1000,
) -> pd.DataFrame:
    """Predict synthetic lethal gene pairs.

    Identifies gene pairs where individual knockouts are viable
    but the double knockout is lethal.

    Parameters
    ----------
    model : cobra.Model
        Genome-scale metabolic model.
    genes : list[str]
        Genes to test for synthetic lethality.
    threshold : float, optional
        Growth threshold for lethality.
    max_combinations : int, optional
        Maximum number of pairs to test.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - gene1, gene2: Gene pair
        - single_ko_1, single_ko_2: Single knockout growth ratios
        - double_ko: Double knockout growth ratio
        - synthetic_lethal: Boolean flag

    Examples
    --------
    >>> sl_pairs = predict_synthetic_lethality(model, gene_list)
    >>> lethal_pairs = sl_pairs[sl_pairs['synthetic_lethal']]
    """
    from itertools import combinations

    # Get wild-type growth
    wt_solution = model.optimize()
    if wt_solution.status != "optimal":
        raise ValueError("Wild-type optimization failed")
    wt_growth = wt_solution.objective_value

    # Test single knockouts first
    single_ko_growth = {}
    for gene_id in genes:
        model_copy = model.copy()
        try:
            gene = model_copy.genes.get_by_id(gene_id)
            gene.knock_out()
            sol = model_copy.optimize()
            if sol.status == "optimal":
                single_ko_growth[gene_id] = sol.objective_value / wt_growth
            else:
                single_ko_growth[gene_id] = 0.0
        except Exception:
            single_ko_growth[gene_id] = 0.0

    # Filter to non-essential genes only
    viable_genes = [g for g, ratio in single_ko_growth.items() if ratio >= threshold]

    # Test pairs
    pairs = list(combinations(viable_genes, 2))
    if len(pairs) > max_combinations:
        logger.warning(
            f"Limiting to {max_combinations} pairs (of {len(pairs)} possible)"
        )
        np.random.seed(42)
        indices = np.random.choice(len(pairs), max_combinations, replace=False)
        pairs = [pairs[i] for i in indices]

    results = []
    for gene1, gene2 in pairs:
        model_copy = model.copy()
        try:
            g1 = model_copy.genes.get_by_id(gene1)
            g2 = model_copy.genes.get_by_id(gene2)
            g1.knock_out()
            g2.knock_out()
            sol = model_copy.optimize()

            if sol.status == "optimal":
                double_ko_ratio = sol.objective_value / wt_growth
            else:
                double_ko_ratio = 0.0

        except Exception:
            double_ko_ratio = 0.0

        is_synthetic_lethal = (
            single_ko_growth[gene1] >= threshold
            and single_ko_growth[gene2] >= threshold
            and double_ko_ratio < threshold
        )

        results.append(
            {
                "gene1": gene1,
                "gene2": gene2,
                "single_ko_1": single_ko_growth[gene1],
                "single_ko_2": single_ko_growth[gene2],
                "double_ko": double_ko_ratio,
                "synthetic_lethal": is_synthetic_lethal,
            }
        )

    return pd.DataFrame(results)
