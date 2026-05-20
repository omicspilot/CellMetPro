"""Statistical comparison of metabolic states.

This module provides methods for identifying differentially active
reactions and pathways between cell populations or conditions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

if TYPE_CHECKING:
    pass


class DifferentialAnalysis:
    """Compare metabolic activity between groups.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    groups : pd.Series
        Group labels for each cell.

    Attributes
    ----------
    results : pd.DataFrame
        Differential analysis results with statistics.
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        groups: pd.Series,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.groups = groups
        self.results: pd.DataFrame | None = None

    def compare_groups(
        self,
        group1: str,
        group2: str,
        method: Literal["wilcoxon", "ttest", "mannwhitneyu"] = "wilcoxon",
    ) -> pd.DataFrame:
        """Compare two groups statistically.

        Parameters
        ----------
        group1 : str
            Name of first group.
        group2 : str
            Name of second group.
        method : str
            Statistical test to use.

        Returns
        -------
        pd.DataFrame
            Results with columns: reaction, group1_mean, group2_mean,
            log2fc, pvalue, padj_bh, padj_bonf.
        """

        valid_methods = {"wilcoxon", "ttest", "mannwhitneyu"}
        if method not in valid_methods:
            raise ValueError(
                f"Invalid method '{method}'. Must be one of: {valid_methods}"
            )

        # Ensure columns in reaction_scores match index in groups
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        # Filter cells by group
        cells1 = group_labels[group_labels == group1].index
        cells2 = group_labels[group_labels == group2].index

        results = []
        for reaction_id in scores.index:
            # reaction scores
            scores1 = scores.loc[reaction_id, cells1]
            scores2 = scores.loc[reaction_id, cells2]

            # means
            mean1 = np.mean(scores1)
            mean2 = np.mean(scores2)

            # log2 for change
            epsilon = 1e-9
            log2fc = np.log2((mean2 + epsilon) / (mean1 + epsilon))

            # run statistical test
            if method == "wilcoxon":
                stat, pval = stats.ranksums(scores1, scores2)
            elif method == "ttest":
                stat, pval = stats.ttest_ind(scores1, scores2)
            elif method == "mannwhitneyu":
                stat, pval = stats.mannwhitneyu(
                    scores1, scores2, alternative="two-sided"
                )

            results.append(
                {
                    "reaction": reaction_id,
                    "group1_mean": mean1,
                    "group2_mean": mean2,
                    "log2fc": log2fc,
                    "statistic": stat,
                    "pvalue": pval,
                }
            )

        # FDR correction
        results_df = pd.DataFrame(results)
        # Benjamini-Hochberg (less strict)
        _, padj_bh, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
        # Bonferroni (more strict)
        _, padj_bonf, _, _ = multipletests(results_df["pvalue"], method="bonferroni")

        results_df["padj_bh"] = padj_bh
        results_df["padj_bonf"] = padj_bonf

        # sort results based on p-val
        return results_df.sort_values("pvalue")

    def rank_reactions(
        self,
        group: str,
        n_top: int = 50,
    ) -> pd.DataFrame:
        """Rank reactions by activity in a group.

        Parameters
        ----------
        group : str
            Group to analyze.
        n_top : int
            Number of top reactions to return.

        Returns
        -------
        pd.DataFrame
            Top reactions ranked by mean activity (lowest penalty = highest activity).
        """
        # Get cells belonging to this group
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        group_cells = group_labels[group_labels == group].index

        # Compute statistics for each reaction
        results = []
        for reaction_id in scores.index:
            reaction_scores = scores.loc[reaction_id, group_cells]
            results.append(
                {
                    "reaction": reaction_id,
                    "mean_score": np.mean(reaction_scores),
                    "std_score": np.std(reaction_scores),
                    "min_score": np.min(reaction_scores),
                    "max_score": np.max(reaction_scores),
                    "n_cells": len(group_cells),
                }
            )

        # Sort by mean score (ascending = lowest penalty = highest activity)
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values("mean_score", ascending=True)

        # Add rank column
        results_df["rank"] = range(1, len(results_df) + 1)

        # Return top N
        return results_df.head(n_top).reset_index(drop=True)

    def compute_effect_size(
        self,
        group1: str,
        group2: str,
    ) -> pd.DataFrame:
        """Compute effect size (Cohen's d) between groups.

        Cohen's d measures the standardized difference between two means:
            d = (mean1 - mean2) / pooled_std

        Interpretation:
            |d| < 0.2  : negligible
            |d| 0.2-0.5: small
            |d| 0.5-0.8: medium
            |d| > 0.8  : large

        Parameters
        ----------
        group1 : str
            Name of first group.
        group2 : str
            Name of second group.

        Returns
        -------
        pd.DataFrame
            Effect sizes for each reaction with interpretation.
        """
        # Get common cells
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        # Filter cells by group
        cells1 = group_labels[group_labels == group1].index
        cells2 = group_labels[group_labels == group2].index

        n1 = len(cells1)
        n2 = len(cells2)

        results = []
        for reaction_id in scores.index:
            scores1 = scores.loc[reaction_id, cells1].values
            scores2 = scores.loc[reaction_id, cells2].values

            mean1 = np.mean(scores1)
            mean2 = np.mean(scores2)
            std1 = np.std(scores1, ddof=1)
            std2 = np.std(scores2, ddof=1)

            # Pooled standard deviation
            pooled_std = np.sqrt(
                ((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2)
            )

            # Cohen's d (avoid division by zero)
            if pooled_std > 0:
                cohens_d = (mean1 - mean2) / pooled_std
            else:
                cohens_d = 0.0

            # Interpretation
            abs_d = abs(cohens_d)
            if abs_d < 0.2:
                interpretation = "negligible"
            elif abs_d < 0.5:
                interpretation = "small"
            elif abs_d < 0.8:
                interpretation = "medium"
            else:
                interpretation = "large"

            results.append(
                {
                    "reaction": reaction_id,
                    "cohens_d": cohens_d,
                    "abs_cohens_d": abs_d,
                    "interpretation": interpretation,
                    "group1_mean": mean1,
                    "group2_mean": mean2,
                    "pooled_std": pooled_std,
                }
            )

        results_df = pd.DataFrame(results)
        return results_df.sort_values("abs_cohens_d", ascending=False)

    def compare_multiple_groups(
        self,
        groups: list[str] | None = None,
        method: Literal["anova", "kruskal"] = "kruskal",
    ) -> pd.DataFrame:
        """Compare metabolic activity across multiple groups.

        Performs ANOVA (parametric) or Kruskal-Wallis (non-parametric) tests
        for each reaction to identify those with significant differences
        across groups.

        Parameters
        ----------
        groups : list[str], optional
            Groups to compare. If None, uses all unique groups.
        method : str
            Statistical test: 'anova' for one-way ANOVA (parametric),
            'kruskal' for Kruskal-Wallis (non-parametric, recommended for
            non-normal distributions typical in single-cell data).

        Returns
        -------
        pd.DataFrame
            Results with columns: reaction, statistic, pvalue, padj_bh, padj_bonf,
            plus mean and std for each group.

        Notes
        -----
        The Kruskal-Wallis test is the non-parametric equivalent of one-way
        ANOVA and is recommended for single-cell data where normality
        assumptions may not hold.

        Example
        -------
        >>> da = DifferentialAnalysis(reaction_scores, cell_groups)
        >>> results = da.compare_multiple_groups(method="kruskal")
        >>> significant = results[results["padj_bh"] < 0.05]
        """
        valid_methods = {"anova", "kruskal"}
        if method not in valid_methods:
            raise ValueError(
                f"Invalid method '{method}'. Must be one of: {valid_methods}"
            )

        # Ensure columns in reaction_scores match index in groups
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]

        # Determine groups to compare
        if groups is None:
            groups = sorted(group_labels.unique())
        else:
            groups = [g for g in groups if g in group_labels.unique()]

        if len(groups) < 2:
            raise ValueError("Need at least 2 groups for comparison")

        # Get cells for each group
        group_cells = {g: group_labels[group_labels == g].index for g in groups}

        results = []
        for reaction_id in scores.index:
            # Get scores for each group
            group_scores = [
                scores.loc[reaction_id, cells].values for cells in group_cells.values()
            ]

            # Run statistical test
            if method == "anova":
                stat, pval = stats.f_oneway(*group_scores)
            else:  # kruskal
                stat, pval = stats.kruskal(*group_scores)

            # Compute group statistics
            row = {
                "reaction": reaction_id,
                "statistic": stat,
                "pvalue": pval,
                "n_groups": len(groups),
            }

            # Add mean and std for each group
            for g, cells in group_cells.items():
                g_scores = scores.loc[reaction_id, cells]
                row[f"{g}_mean"] = np.mean(g_scores)
                row[f"{g}_std"] = np.std(g_scores)
                row[f"{g}_n"] = len(cells)

            results.append(row)

        # FDR correction
        results_df = pd.DataFrame(results)
        _, padj_bh, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
        _, padj_bonf, _, _ = multipletests(results_df["pvalue"], method="bonferroni")

        results_df["padj_bh"] = padj_bh
        results_df["padj_bonf"] = padj_bonf

        return results_df.sort_values("pvalue")

    def posthoc_tests(
        self,
        reaction: str,
        groups: list[str] | None = None,
        method: Literal["tukey", "dunn", "conover"] = "dunn",
    ) -> pd.DataFrame:
        """Perform post-hoc pairwise comparisons after multi-group test.

        After finding significant differences with compare_multiple_groups(),
        use this method to identify which specific pairs of groups differ.

        Parameters
        ----------
        reaction : str
            Reaction ID to analyze.
        groups : list[str], optional
            Groups to compare. If None, uses all unique groups.
        method : str
            Post-hoc test method:
            - 'tukey': Tukey's HSD (for use after ANOVA, assumes normality)
            - 'dunn': Dunn's test with Bonferroni correction
              (for use after Kruskal-Wallis, non-parametric)
            - 'conover': Conover's test (for use after Kruskal-Wallis)

        Returns
        -------
        pd.DataFrame
            Pairwise comparisons with columns: group1, group2, statistic,
            pvalue, padj, significant.

        Example
        -------
        >>> da = DifferentialAnalysis(reaction_scores, cell_groups)
        >>> # First find reactions with significant multi-group differences
        >>> multi_results = da.compare_multiple_groups()
        >>> sig_reactions = multi_results[multi_results["padj_bh"] < 0.05]["reaction"]
        >>> # Then do post-hoc for specific reaction
        >>> posthoc = da.posthoc_tests("R1", method="dunn")
        """
        valid_methods = {"tukey", "dunn", "conover"}
        if method not in valid_methods:
            raise ValueError(
                f"Invalid method '{method}'. Must be one of: {valid_methods}"
            )

        if reaction not in self.reaction_scores.index:
            raise ValueError(f"Reaction {reaction} not found in data")

        # Ensure columns match
        common_cells = self.reaction_scores.columns.intersection(self.groups.index)
        scores = self.reaction_scores.loc[reaction, common_cells]
        group_labels = self.groups[common_cells]

        # Determine groups to compare
        if groups is None:
            groups = sorted(group_labels.unique())
        else:
            groups = [g for g in groups if g in group_labels.unique()]

        if len(groups) < 2:
            raise ValueError("Need at least 2 groups for comparison")

        # Get scores for each group
        group_data = {g: scores[group_labels == g].values for g in groups}

        if method == "tukey":
            return self._posthoc_tukey(group_data, groups)
        elif method == "dunn":
            return self._posthoc_dunn(group_data, groups)
        else:  # conover
            return self._posthoc_conover(group_data, groups)

    def _posthoc_tukey(
        self,
        group_data: dict[str, np.ndarray],
        groups: list[str],
    ) -> pd.DataFrame:
        """Tukey's HSD post-hoc test."""
        from statsmodels.stats.multicomp import pairwise_tukeyhsd

        # Prepare data for statsmodels
        all_values: list[float] = []
        all_groups: list[str] = []
        for g in groups:
            all_values.extend(group_data[g])
            all_groups.extend([g] * len(group_data[g]))

        # Run Tukey HSD
        tukey_result = pairwise_tukeyhsd(all_values, all_groups)

        # Extract results from summary table (more reliable than manual indexing)
        summary_data = tukey_result.summary().data[1:]  # Skip header row
        results = []
        for row in summary_data:
            group1, group2, meandiff, p_adj, lower, upper, reject = row
            results.append(
                {
                    "group1": group1,
                    "group2": group2,
                    "mean_diff": float(meandiff),
                    "pvalue": float(p_adj),
                    "significant": (
                        reject if isinstance(reject, bool) else reject == "True"
                    ),
                }
            )

        return pd.DataFrame(results)

    def _posthoc_dunn(
        self,
        group_data: dict[str, np.ndarray],
        groups: list[str],
    ) -> pd.DataFrame:
        """Dunn's test for pairwise comparisons after Kruskal-Wallis.

        Uses Bonferroni correction for multiple comparisons.
        """
        # Combine all data and compute ranks
        all_values = np.concatenate([group_data[g] for g in groups])
        ranks = stats.rankdata(all_values)

        # Map ranks back to groups
        group_ranks = {}
        idx = 0
        for g in groups:
            n = len(group_data[g])
            group_ranks[g] = ranks[idx : idx + n]
            idx += n

        N = len(all_values)  # Total sample size
        n_comparisons = len(groups) * (len(groups) - 1) // 2

        # Guard against edge case
        if N <= 1:
            raise ValueError("Need at least 2 total samples for Dunn's test")

        # Pre-compute tie correction term (shared across all pairs)
        _, counts = np.unique(all_values, return_counts=True)
        tie_sum = np.sum(counts**3 - counts)  # Sum of (t_i^3 - t_i) for tied groups

        results = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i + 1 :]:
                n1 = len(group_data[g1])
                n2 = len(group_data[g2])

                # Mean ranks
                mean_rank1 = np.mean(group_ranks[g1])
                mean_rank2 = np.mean(group_ranks[g2])

                # Standard error of difference with tie correction
                # Dunn (1964) formula: SE = sqrt(((N*(N+1)/12) - A) * (1/ni + 1/nj))
                # where A = tie_sum / (12*(N-1))
                variance_term = N * (N + 1) / 12
                if N > 1 and tie_sum > 0:
                    variance_term -= tie_sum / (12 * (N - 1))

                se = np.sqrt(max(0, variance_term) * (1 / n1 + 1 / n2))

                if se > 0:
                    z = (mean_rank1 - mean_rank2) / se
                    # Two-tailed p-value
                    pval = 2 * (1 - stats.norm.cdf(abs(z)))
                else:
                    z = 0.0
                    pval = 1.0

                # Bonferroni correction
                padj = min(pval * n_comparisons, 1.0)

                results.append(
                    {
                        "group1": g1,
                        "group2": g2,
                        "mean_rank1": mean_rank1,
                        "mean_rank2": mean_rank2,
                        "z_statistic": z,
                        "pvalue": pval,
                        "padj": padj,
                        "significant": padj < 0.05,
                    }
                )

        return pd.DataFrame(results).sort_values("pvalue")

    def _posthoc_conover(
        self,
        group_data: dict[str, np.ndarray],
        groups: list[str],
    ) -> pd.DataFrame:
        """Conover's test for pairwise comparisons after Kruskal-Wallis.

        Conover's test is more powerful than Dunn's test but assumes
        that the Kruskal-Wallis test was significant.
        """
        # First perform Kruskal-Wallis to get H statistic
        group_arrays = [group_data[g] for g in groups]
        H, _ = stats.kruskal(*group_arrays)

        # Combine all data and compute ranks
        all_values = np.concatenate(group_arrays)
        ranks = stats.rankdata(all_values)

        # Map ranks back to groups
        group_ranks = {}
        idx = 0
        for g in groups:
            n = len(group_data[g])
            group_ranks[g] = ranks[idx : idx + n]
            idx += n

        N = len(all_values)
        k = len(groups)
        n_comparisons = k * (k - 1) // 2

        # Guard against edge cases
        if N <= 1:
            raise ValueError("Need at least 2 total samples for Conover's test")
        if N <= k:
            raise ValueError(
                f"Need more samples than groups for Conover's test. "
                f"Got {N} samples and {k} groups."
            )

        # Calculate S^2 (variance of ranks)
        S2 = (1 / (N - 1)) * (np.sum(ranks**2) - N * ((N + 1) / 2) ** 2)

        # Calculate A term
        A = S2 * (N - 1 - H) / (N - k)

        results = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i + 1 :]:
                n1 = len(group_data[g1])
                n2 = len(group_data[g2])

                mean_rank1 = np.mean(group_ranks[g1])
                mean_rank2 = np.mean(group_ranks[g2])

                # Test statistic
                denominator = np.sqrt(A * (1 / n1 + 1 / n2))
                if denominator > 0:
                    t = (mean_rank1 - mean_rank2) / denominator
                    # Two-tailed p-value using t-distribution
                    df = N - k
                    pval = 2 * (1 - stats.t.cdf(abs(t), df))
                else:
                    t = 0.0
                    pval = 1.0

                # Bonferroni correction
                padj = min(pval * n_comparisons, 1.0)

                results.append(
                    {
                        "group1": g1,
                        "group2": g2,
                        "mean_rank1": mean_rank1,
                        "mean_rank2": mean_rank2,
                        "t_statistic": t,
                        "pvalue": pval,
                        "padj": padj,
                        "significant": padj < 0.05,
                    }
                )

        return pd.DataFrame(results).sort_values("pvalue")

    def all_pairwise_comparisons(
        self,
        groups: list[str] | None = None,
        method: Literal["wilcoxon", "ttest", "mannwhitneyu"] = "wilcoxon",
    ) -> pd.DataFrame:
        """Perform all pairwise comparisons between groups.

        This is useful when you want to see how each pair of groups
        differs for all reactions, rather than focusing on a single reaction.

        Parameters
        ----------
        groups : list[str], optional
            Groups to compare. If None, uses all unique groups.
        method : str
            Statistical test to use for pairwise comparisons.

        Returns
        -------
        pd.DataFrame
            Results for all pairwise comparisons with columns:
            comparison, reaction, group1_mean, group2_mean, log2fc, pvalue, padj.

        Example
        -------
        >>> da = DifferentialAnalysis(reaction_scores, cell_groups)
        >>> all_pairs = da.all_pairwise_comparisons()
        >>> # Filter for specific comparison
        >>> a_vs_b = all_pairs[all_pairs["comparison"] == "A_vs_B"]
        """
        # Determine groups to compare
        if groups is None:
            groups = sorted(self.groups.unique())
        else:
            groups = [g for g in groups if g in self.groups.unique()]

        if len(groups) < 2:
            raise ValueError("Need at least 2 groups for comparison")

        all_results = []
        for i, g1 in enumerate(groups):
            for g2 in groups[i + 1 :]:
                comparison_name = f"{g1}_vs_{g2}"
                result = self.compare_groups(g1, g2, method=method)
                result["comparison"] = comparison_name
                result["group1"] = g1
                result["group2"] = g2
                all_results.append(result)

        return pd.concat(all_results, ignore_index=True)


class PseudoBulkAnalysis:
    """Compare metabolic activity between groups using biological replicates.

    Aggregates reaction scores within each (sample, group) combination before
    running statistics, treating each biological replicate as an independent
    observation. This is the statistically correct approach when cells come
    from multiple donors, animals, or experimental batches — it avoids the
    pseudoreplication problem of treating cells as independent when they share
    a common source.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    groups : pd.Series
        Condition label for each cell (e.g. ``"control"``, ``"treatment"``).
    samples : pd.Series
        Biological replicate label for each cell (e.g. ``"donor1"``).
        Each sample must belong to exactly one group.

    Examples
    --------
    >>> pb = PseudoBulkAnalysis(reaction_scores, cell_groups, cell_samples)
    >>> result = pb.compare_groups("control", "treatment")
    >>> significant = result[result["padj_bh"] < 0.05]
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        groups: pd.Series,
        samples: pd.Series,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.groups = groups
        self.samples = samples
        self._pseudobulk: pd.DataFrame | None = None
        self._sample_groups: pd.Series | None = None

    def aggregate(
        self,
        method: Literal["mean", "sum"] = "mean",
    ) -> pd.DataFrame:
        """Aggregate cell-level scores to sample-level pseudo-bulk profiles.

        Parameters
        ----------
        method : str
            Aggregation function: ``"mean"`` (default) or ``"sum"``.

        Returns
        -------
        pd.DataFrame
            Pseudo-bulk scores (reactions x samples).

        Raises
        ------
        ValueError
            If any sample has cells assigned to more than one group.
        """
        common_cells = self.reaction_scores.columns.intersection(
            self.groups.index.intersection(self.samples.index)
        )
        scores = self.reaction_scores[common_cells]
        group_labels = self.groups[common_cells]
        sample_labels = self.samples[common_cells]

        # Verify each sample belongs to exactly one group
        sample_group_map: dict[str, str] = {}
        for cell in common_cells:
            s = sample_labels[cell]
            g = group_labels[cell]
            if s in sample_group_map and sample_group_map[s] != g:
                raise ValueError(
                    f"Sample '{s}' has cells in multiple groups "
                    f"({sample_group_map[s]!r} and {g!r}). "
                    "Each sample must belong to exactly one group."
                )
            sample_group_map[s] = g

        agg_func = np.mean if method == "mean" else np.sum
        pseudobulk = {
            sample: scores.loc[:, sample_labels[sample_labels == sample].index].apply(
                agg_func, axis=1
            )
            for sample in sample_labels.unique()
        }

        self._pseudobulk = pd.DataFrame(pseudobulk)
        self._sample_groups = pd.Series(sample_group_map)
        return self._pseudobulk

    def compare_groups(
        self,
        group1: str,
        group2: str,
        aggregate_method: Literal["mean", "sum"] = "mean",
        method: Literal["ttest", "wilcoxon", "mannwhitneyu"] = "ttest",
    ) -> pd.DataFrame:
        """Compare two groups using sample-level pseudo-bulk profiles.

        Parameters
        ----------
        group1 : str
            Name of the first group.
        group2 : str
            Name of the second group.
        aggregate_method : str
            How to aggregate cells within each sample.
        method : str
            Statistical test applied to the sample-level data.

        Returns
        -------
        pd.DataFrame
            Results with columns: reaction, group1_mean, group2_mean, log2fc,
            statistic, pvalue, padj_bh, padj_bonf, group1_n_samples,
            group2_n_samples.

        Raises
        ------
        ValueError
            If either group has fewer than 2 biological replicates.
        """
        valid_methods = {"ttest", "wilcoxon", "mannwhitneyu"}
        if method not in valid_methods:
            raise ValueError(
                f"Invalid method '{method}'. Must be one of: {valid_methods}"
            )

        pseudobulk = self.aggregate(method=aggregate_method)
        sample_groups = self._sample_groups  # set by aggregate()

        samples1 = sample_groups[sample_groups == group1].index
        samples2 = sample_groups[sample_groups == group2].index

        if len(samples1) < 2:
            raise ValueError(
                f"Group '{group1}' has {len(samples1)} sample(s). "
                "Need at least 2 biological replicates for pseudo-bulk analysis."
            )
        if len(samples2) < 2:
            raise ValueError(
                f"Group '{group2}' has {len(samples2)} sample(s). "
                "Need at least 2 biological replicates for pseudo-bulk analysis."
            )

        results = []
        for reaction_id in pseudobulk.index:
            vals1 = pseudobulk.loc[reaction_id, samples1].values
            vals2 = pseudobulk.loc[reaction_id, samples2].values

            mean1 = np.mean(vals1)
            mean2 = np.mean(vals2)
            epsilon = 1e-9
            log2fc = np.log2((mean2 + epsilon) / (mean1 + epsilon))

            if method == "ttest":
                stat, pval = stats.ttest_ind(vals1, vals2)
            elif method == "wilcoxon":
                stat, pval = stats.ranksums(vals1, vals2)
            else:  # mannwhitneyu
                stat, pval = stats.mannwhitneyu(vals1, vals2, alternative="two-sided")

            results.append(
                {
                    "reaction": reaction_id,
                    "group1_mean": mean1,
                    "group2_mean": mean2,
                    "log2fc": log2fc,
                    "statistic": stat,
                    "pvalue": pval,
                    "group1_n_samples": len(samples1),
                    "group2_n_samples": len(samples2),
                }
            )

        results_df = pd.DataFrame(results)
        _, padj_bh, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
        _, padj_bonf, _, _ = multipletests(results_df["pvalue"], method="bonferroni")
        results_df["padj_bh"] = padj_bh
        results_df["padj_bonf"] = padj_bonf

        return results_df.sort_values("pvalue")

    def compare_multiple_groups(
        self,
        groups: list[str] | None = None,
        aggregate_method: Literal["mean", "sum"] = "mean",
        method: Literal["anova", "kruskal"] = "anova",
    ) -> pd.DataFrame:
        """Compare metabolic activity across multiple groups.

        Parameters
        ----------
        groups : list[str], optional
            Groups to compare. If ``None``, uses all unique groups.
        aggregate_method : str
            How to aggregate cells within each sample.
        method : str
            Statistical test: ``"anova"`` (parametric) or ``"kruskal"``
            (non-parametric Kruskal-Wallis).

        Returns
        -------
        pd.DataFrame
            Results with per-group means and sample counts, statistic,
            pvalue, padj_bh, padj_bonf.

        Raises
        ------
        ValueError
            If any group has fewer than 2 biological replicates, or if
            fewer than 2 groups are available.
        """
        valid_methods = {"anova", "kruskal"}
        if method not in valid_methods:
            raise ValueError(
                f"Invalid method '{method}'. Must be one of: {valid_methods}"
            )

        pseudobulk = self.aggregate(method=aggregate_method)
        sample_groups = self._sample_groups

        if groups is None:
            groups = sorted(sample_groups.unique())
        else:
            groups = [g for g in groups if g in sample_groups.unique()]

        if len(groups) < 2:
            raise ValueError("Need at least 2 groups for comparison.")

        group_samples = {g: sample_groups[sample_groups == g].index for g in groups}

        for g, samps in group_samples.items():
            if len(samps) < 2:
                raise ValueError(
                    f"Group '{g}' has {len(samps)} sample(s). "
                    "Need at least 2 biological replicates."
                )

        results = []
        for reaction_id in pseudobulk.index:
            group_vals = [
                pseudobulk.loc[reaction_id, samps].values
                for samps in group_samples.values()
            ]

            if method == "anova":
                stat, pval = stats.f_oneway(*group_vals)
            else:
                stat, pval = stats.kruskal(*group_vals)

            row: dict = {
                "reaction": reaction_id,
                "statistic": stat,
                "pvalue": pval,
                "n_groups": len(groups),
            }
            for g, samps in group_samples.items():
                g_vals = pseudobulk.loc[reaction_id, samps]
                row[f"{g}_mean"] = np.mean(g_vals)
                row[f"{g}_std"] = np.std(g_vals)
                row[f"{g}_n_samples"] = len(samps)

            results.append(row)

        results_df = pd.DataFrame(results)
        _, padj_bh, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
        _, padj_bonf, _, _ = multipletests(results_df["pvalue"], method="bonferroni")
        results_df["padj_bh"] = padj_bh
        results_df["padj_bonf"] = padj_bonf

        return results_df.sort_values("pvalue")
