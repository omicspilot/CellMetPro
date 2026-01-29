"""Pathway-level aggregation and analysis.

This module aggregates reaction-level scores into pathway scores
and provides pathway enrichment analysis using Gene Ontology (GO) terms.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

if TYPE_CHECKING:
    import cobra


class PathwayAnalyzer:
    """Aggregate reactions into pathway-level scores.

    Parameters
    ----------
    reaction_scores : pd.DataFrame
        Reaction activity scores (reactions x cells).
    model : cobra.Model
        Metabolic model for pathway annotations.

    Attributes
    ----------
    pathway_scores : pd.DataFrame
        Pathway-level aggregated scores.
    pathway_mapping : dict
        Mapping of reactions to pathways.
    """

    def __init__(
        self,
        reaction_scores: pd.DataFrame,
        model: cobra.Model,
    ) -> None:
        self.reaction_scores = reaction_scores
        self.model = model
        self.pathway_scores: pd.DataFrame | None = None
        self.pathway_mapping: dict | None = None

    def get_pathway_mapping(self) -> dict[str, list[str]]:
        """Extract pathway-to-reaction mapping from model subsystems.

        Returns
        -------
        dict
            Mapping of pathway/subsystem names to reaction IDs.
        """
        pathway_mapping: dict[str, list[str]] = {}

        for reaction in self.model.reactions:
            subsystem = reaction.subsystem
            if subsystem:
                if subsystem not in pathway_mapping:
                    pathway_mapping[subsystem] = []
                pathway_mapping[subsystem].append(reaction.id)

        self.pathway_mapping = pathway_mapping
        return pathway_mapping

    def aggregate(
        self,
        method: Literal["mean", "median", "sum", "max"] = "mean",
    ) -> pd.DataFrame:
        """Aggregate reaction scores to pathway level.

        Parameters
        ----------
        method : str
            Aggregation method: 'mean', 'median', 'sum', or 'max'.

        Returns
        -------
        pd.DataFrame
            Pathway scores (pathways x cells).
        """
        if self.pathway_mapping is None:
            self.get_pathway_mapping()

        agg_funcs = {
            "mean": np.mean,
            "median": np.median,
            "sum": np.sum,
            "max": np.max,
        }
        agg_func = agg_funcs[method]

        results = {}
        for pathway, reactions in self.pathway_mapping.items():
            # Get reactions that exist in our scores
            available_reactions = [
                r for r in reactions if r in self.reaction_scores.index
            ]
            if available_reactions:
                pathway_data = self.reaction_scores.loc[available_reactions]
                results[pathway] = agg_func(pathway_data.values, axis=0)

        self.pathway_scores = pd.DataFrame(
            results,
            index=self.reaction_scores.columns
        ).T

        return self.pathway_scores

    def enrich(
        self,
        differential_results: pd.DataFrame,
        method: Literal["ora"] = "ora",
    ) -> pd.DataFrame:
        """Perform pathway enrichment analysis using subsystems.

        Parameters
        ----------
        differential_results : pd.DataFrame
            Results from differential analysis with 'reaction' column.
        method : str
            Enrichment method: 'ora' (Over-Representation Analysis).

        Returns
        -------
        pd.DataFrame
            Enrichment results with pathway names and statistics.
        """
        if self.pathway_mapping is None:
            self.get_pathway_mapping()

        # Get significant reactions (use padj_bh < 0.05 if available)
        if "padj_bh" in differential_results.columns:
            sig_reactions = set(
                differential_results[differential_results["padj_bh"] < 0.05]["reaction"]
            )
        else:
            sig_reactions = set(
                differential_results[differential_results["pvalue"] < 0.05]["reaction"]
            )

        # Background = all tested reactions
        background = set(differential_results["reaction"])

        return self._ora_enrichment(sig_reactions, background)

    def _ora_enrichment(
        self,
        sig_reactions: set,
        background: set,
    ) -> pd.DataFrame:
        """Over-Representation Analysis using Fisher's exact test."""
        results = []

        for pathway, reactions in self.pathway_mapping.items():
            pathway_reactions = set(reactions) & background

            if len(pathway_reactions) == 0:
                continue

            # Contingency table:
            #                    In Pathway  | Not in Pathway
            # Significant      |     a       |      b
            # Not Significant  |     c       |      d

            a = len(sig_reactions & pathway_reactions)
            b = len(sig_reactions - pathway_reactions)
            c = len(pathway_reactions - sig_reactions)
            d = len(background - sig_reactions - pathway_reactions)

            # Fisher's exact test
            _, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

            # Fold enrichment
            expected = len(sig_reactions) * len(pathway_reactions) / len(background)
            fold_enrichment = a / expected if expected > 0 else 0

            results.append({
                "pathway": pathway,
                "n_sig_in_pathway": a,
                "n_pathway_total": len(pathway_reactions),
                "n_sig_total": len(sig_reactions),
                "n_background": len(background),
                "fold_enrichment": fold_enrichment,
                "pvalue": pval,
            })

        results_df = pd.DataFrame(results)

        if len(results_df) > 0:
            _, padj, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
            results_df["padj"] = padj
            results_df = results_df.sort_values("pvalue")

        return results_df


class GOEnrichmentAnalyzer:
    """Perform Gene Ontology (GO) enrichment analysis.

    This class maps reactions to GO terms via their associated genes
    and performs over-representation analysis.

    Parameters
    ----------
    model : cobra.Model
        Metabolic model with gene-reaction associations.
    go_annotations : pd.DataFrame
        GO annotations with columns: 'gene_id', 'go_term', 'go_name', 'namespace'.

    Example
    -------
    >>> # Load GO annotations (from GAF file or custom mapping)
    >>> go_df = pd.DataFrame({
    ...     'gene_id': ['gene1', 'gene1', 'gene2'],
    ...     'go_term': ['GO:0006096', 'GO:0005737', 'GO:0006096'],
    ...     'go_name': ['glycolysis', 'cytoplasm', 'glycolysis'],
    ...     'namespace': ['BP', 'CC', 'BP']
    ... })
    >>> analyzer = GOEnrichmentAnalyzer(model, go_df)
    >>> enrichment = analyzer.enrich_reactions(significant_reactions)
    """

    def __init__(
        self,
        model: cobra.Model,
        go_annotations: pd.DataFrame,
    ) -> None:
        self.model = model
        self.go_annotations = go_annotations

        # Build gene-to-GO mapping
        self.gene_to_go = self._build_gene_to_go_mapping()

        # Build reaction-to-gene mapping from model
        self.reaction_to_genes = self._build_reaction_to_genes_mapping()

    def _build_gene_to_go_mapping(self) -> dict[str, set[str]]:
        """Build mapping from gene IDs to GO terms."""
        gene_to_go: dict[str, set[str]] = {}

        for _, row in self.go_annotations.iterrows():
            gene_id = row["gene_id"]
            go_term = row["go_term"]

            if gene_id not in gene_to_go:
                gene_to_go[gene_id] = set()
            gene_to_go[gene_id].add(go_term)

        return gene_to_go

    def _build_reaction_to_genes_mapping(self) -> dict[str, set[str]]:
        """Build mapping from reaction IDs to gene IDs."""
        reaction_to_genes: dict[str, set[str]] = {}

        for reaction in self.model.reactions:
            genes = {g.id for g in reaction.genes}
            if genes:
                reaction_to_genes[reaction.id] = genes

        return reaction_to_genes

    def get_reaction_go_terms(self, reaction_id: str) -> set[str]:
        """Get all GO terms associated with a reaction via its genes.

        Parameters
        ----------
        reaction_id : str
            Reaction ID.

        Returns
        -------
        set
            Set of GO term IDs.
        """
        go_terms: set[str] = set()

        genes = self.reaction_to_genes.get(reaction_id, set())
        for gene in genes:
            go_terms.update(self.gene_to_go.get(gene, set()))

        return go_terms

    def get_go_to_reactions_mapping(
        self,
        reactions: set[str] | None = None,
    ) -> dict[str, set[str]]:
        """Build mapping from GO terms to reactions.

        Parameters
        ----------
        reactions : set, optional
            Limit to these reactions. If None, use all model reactions.

        Returns
        -------
        dict
            Mapping of GO term IDs to reaction IDs.
        """
        if reactions is None:
            reactions = set(self.reaction_to_genes.keys())

        go_to_reactions: dict[str, set[str]] = {}

        for rxn_id in reactions:
            go_terms = self.get_reaction_go_terms(rxn_id)
            for go_term in go_terms:
                if go_term not in go_to_reactions:
                    go_to_reactions[go_term] = set()
                go_to_reactions[go_term].add(rxn_id)

        return go_to_reactions

    def enrich_reactions(
        self,
        significant_reactions: set[str],
        background_reactions: set[str] | None = None,
        min_genes: int = 3,
        namespace: str | None = None,
    ) -> pd.DataFrame:
        """Perform GO enrichment analysis on a set of reactions.

        Parameters
        ----------
        significant_reactions : set
            Set of significant reaction IDs.
        background_reactions : set, optional
            Background reaction set. If None, uses all model reactions
            with GO annotations.
        min_genes : int
            Minimum number of reactions in a GO term to consider.
        namespace : str, optional
            Filter by GO namespace: 'BP' (biological process),
            'MF' (molecular function), 'CC' (cellular component).
            If None, include all.

        Returns
        -------
        pd.DataFrame
            Enrichment results sorted by p-value.
        """
        # Set background
        if background_reactions is None:
            background_reactions = set(self.reaction_to_genes.keys())

        # Ensure significant reactions are subset of background
        significant_reactions = significant_reactions & background_reactions

        # Get GO-to-reactions mapping for background
        go_to_reactions = self.get_go_to_reactions_mapping(background_reactions)

        # Filter by namespace if specified
        if namespace:
            filtered = self.go_annotations[
                self.go_annotations["namespace"] == namespace
            ]
            namespace_terms = set(filtered["go_term"])
            go_to_reactions = {
                go: rxns for go, rxns in go_to_reactions.items()
                if go in namespace_terms
            }

        # Build GO term metadata
        go_metadata = self.go_annotations.drop_duplicates(
            subset=["go_term"]
        ).set_index("go_term")

        results = []
        for go_term, go_reactions in go_to_reactions.items():
            # Skip terms with too few reactions
            if len(go_reactions) < min_genes:
                continue

            # Fisher's exact test
            # Contingency table:
            #                    In GO term  | Not in GO term
            # Significant      |     a       |      b
            # Not Significant  |     c       |      d

            a = len(significant_reactions & go_reactions)
            b = len(significant_reactions - go_reactions)
            c = len(go_reactions - significant_reactions)
            d = len(background_reactions - significant_reactions - go_reactions)

            if a == 0:
                continue  # Skip if no significant reactions in this term

            _, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

            # Fold enrichment
            n_sig = len(significant_reactions)
            n_go = len(go_reactions)
            n_bg = len(background_reactions)
            expected = n_sig * n_go / n_bg
            fold_enrichment = a / expected if expected > 0 else 0

            # Get GO term name if available
            if go_term in go_metadata.index:
                go_name = go_metadata.loc[go_term, "go_name"]
                go_namespace = go_metadata.loc[go_term, "namespace"]
            else:
                go_name = ""
                go_namespace = ""

            results.append({
                "go_term": go_term,
                "go_name": go_name,
                "namespace": go_namespace,
                "n_sig_in_term": a,
                "n_term_total": len(go_reactions),
                "n_sig_total": len(significant_reactions),
                "n_background": len(background_reactions),
                "fold_enrichment": fold_enrichment,
                "pvalue": pval,
                "reactions": ",".join(sorted(significant_reactions & go_reactions)),
            })

        results_df = pd.DataFrame(results)

        if len(results_df) > 0:
            _, padj, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
            results_df["padj"] = padj
            results_df = results_df.sort_values("pvalue")

        return results_df.reset_index(drop=True)

    def enrich_from_differential(
        self,
        differential_results: pd.DataFrame,
        pvalue_threshold: float = 0.05,
        use_adjusted: bool = True,
        **kwargs,
    ) -> pd.DataFrame:
        """Perform GO enrichment using differential analysis results.

        Parameters
        ----------
        differential_results : pd.DataFrame
            Results from DifferentialAnalysis.compare_groups().
        pvalue_threshold : float
            P-value threshold for significance.
        use_adjusted : bool
            Use adjusted p-values (padj_bh) if True, else use raw pvalue.
        **kwargs
            Additional arguments passed to enrich_reactions().

        Returns
        -------
        pd.DataFrame
            GO enrichment results.
        """
        use_padj = use_adjusted and "padj_bh" in differential_results.columns
        pval_col = "padj_bh" if use_padj else "pvalue"

        sig_mask = differential_results[pval_col] < pvalue_threshold
        significant = set(differential_results[sig_mask]["reaction"])
        background = set(differential_results["reaction"])

        return self.enrich_reactions(significant, background, **kwargs)


def load_go_annotations_from_gaf(gaf_path: str) -> pd.DataFrame:
    """Load GO annotations from a GAF (Gene Association Format) file.

    Parameters
    ----------
    gaf_path : str
        Path to the GAF file.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: gene_id, go_term, go_name, namespace.

    Notes
    -----
    GAF format specification: http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
    """
    rows = []

    with open(gaf_path, "r") as f:
        for line in f:
            if line.startswith("!"):
                continue  # Skip header lines

            parts = line.strip().split("\t")
            if len(parts) < 15:
                continue

            gene_id = parts[2]  # DB Object Symbol
            go_term = parts[4]  # GO ID
            namespace_code = parts[8]  # Aspect (F, P, C)

            namespace_map = {"F": "MF", "P": "BP", "C": "CC"}
            namespace = namespace_map.get(namespace_code, namespace_code)

            rows.append({
                "gene_id": gene_id,
                "go_term": go_term,
                "go_name": "",  # GAF doesn't include names
                "namespace": namespace,
            })

    return pd.DataFrame(rows)


def create_go_annotations_from_dict(
    gene_to_go: dict[str, list[tuple[str, str, str]]],
) -> pd.DataFrame:
    """Create GO annotations DataFrame from a dictionary.

    Parameters
    ----------
    gene_to_go : dict
        Mapping of gene_id to list of (go_term, go_name, namespace) tuples.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: gene_id, go_term, go_name, namespace.

    Example
    -------
    >>> gene_to_go = {
    ...     'gene1': [
    ...         ('GO:0006096', 'glycolytic process', 'BP'),
    ...         ('GO:0005737', 'cytoplasm', 'CC'),
    ...     ],
    ...     'gene2': [
    ...         ('GO:0006096', 'glycolytic process', 'BP'),
    ...     ],
    ... }
    >>> go_df = create_go_annotations_from_dict(gene_to_go)
    """
    rows = []
    for gene_id, go_terms in gene_to_go.items():
        for go_term, go_name, namespace in go_terms:
            rows.append({
                "gene_id": gene_id,
                "go_term": go_term,
                "go_name": go_name,
                "namespace": namespace,
            })

    return pd.DataFrame(rows)


def subsystem_enrichment(
    significant_reactions: list[str] | set[str],
    background: list[str] | set[str],
    subsystem_mapping: dict[str, list[str]],
) -> pd.DataFrame:
    """Perform subsystem enrichment analysis using Fisher's exact test.

    This is a standalone function for enrichment analysis that doesn't
    require reaction scores or a full model object.

    Parameters
    ----------
    significant_reactions : list or set
        List of significant reaction IDs.
    background : list or set
        List of all reaction IDs (background set).
    subsystem_mapping : dict
        Mapping of subsystem names to lists of reaction IDs.

    Returns
    -------
    pd.DataFrame
        Enrichment results with columns: pathway, n_sig_in_pathway,
        n_pathway_total, n_sig_total, n_background, fold_enrichment,
        pvalue, padj.

    Example
    -------
    >>> sig_reactions = ["R1", "R2", "R3"]
    >>> background = ["R1", "R2", "R3", "R4", "R5", "R6"]
    >>> subsystem_map = {"Glycolysis": ["R1", "R2"], "TCA": ["R3", "R4"]}
    >>> results = subsystem_enrichment(sig_reactions, background, subsystem_map)
    """
    sig_reactions = set(significant_reactions)
    bg_reactions = set(background)

    results = []

    for pathway, reactions in subsystem_mapping.items():
        pathway_reactions = set(reactions) & bg_reactions

        if len(pathway_reactions) == 0:
            continue

        # Contingency table
        a = len(sig_reactions & pathway_reactions)
        b = len(sig_reactions - pathway_reactions)
        c = len(pathway_reactions - sig_reactions)
        d = len(bg_reactions - sig_reactions - pathway_reactions)

        # Fisher's exact test (one-sided, greater)
        _, pval = stats.fisher_exact([[a, b], [c, d]], alternative="greater")

        # Fold enrichment
        expected = len(sig_reactions) * len(pathway_reactions) / len(bg_reactions)
        fold_enrichment = a / expected if expected > 0 else 0

        results.append({
            "pathway": pathway,
            "n_sig_in_pathway": a,
            "n_pathway_total": len(pathway_reactions),
            "n_sig_total": len(sig_reactions),
            "n_background": len(bg_reactions),
            "fold_enrichment": fold_enrichment,
            "pvalue": pval,
        })

    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        _, padj, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
        results_df["padj"] = padj
        results_df = results_df.sort_values("pvalue")

    return results_df.reset_index(drop=True)
