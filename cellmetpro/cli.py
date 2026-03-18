"""Command-line interface for CellMetPro.

This module provides CLI commands for running metabolic analysis
pipelines from the command line.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from rich_argparse import RawDescriptionRichHelpFormatter, RichHelpFormatter

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("cellmetpro")


def _open_file(path: str) -> None:
    """Open a file with the OS default application."""
    import os
    import subprocess

    if sys.platform == "darwin":
        subprocess.run(["open", path], check=False)
    elif sys.platform.startswith("linux"):
        subprocess.run(["xdg-open", path], check=False)
    elif sys.platform == "win32":
        os.startfile(path)


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the CLI.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser.
    """
    # Configure rich-argparse styles
    RichHelpFormatter.styles["argparse.prog"] = "bold cyan"
    RichHelpFormatter.styles["argparse.args"] = "green"
    RichHelpFormatter.styles["argparse.groups"] = "bold magenta"
    RichHelpFormatter.styles["argparse.metavar"] = "yellow"

    parser = argparse.ArgumentParser(
        prog="cellmetpro",
        description="""
╔═══════════════════════════════════════════════════════════════════════╗
║                         CellMetPro:                                   ║
║         Cellular Metabolic Profiler for scRNA-seq data                ║
╚═══════════════════════════════════════════════════════════════════════╝

Analyze metabolic activity at single-cell resolution using the COMPASS
algorithm. Score reactions, identify metabolic heterogeneity, and
discover metabolic programs in your scRNA-seq data.
""",
        formatter_class=RawDescriptionRichHelpFormatter,
        epilog="""
Examples:
  cellmetpro run expression.h5ad -m human -o results/
  cellmetpro run data.csv -m /path/to/model.xml -o output/
  cellmetpro run seurat.rds -m human --assay RNA --slot data -o results/
  cellmetpro run large_data.h5ad --microcluster --cells-per-cluster 100
  cellmetpro differential scores.csv groups.csv --plot
  cellmetpro cluster scores.csv --method leiden --embedding umap
  cellmetpro report results/ -o my_report.html
  cellmetpro batch-correct scores.csv batches.csv --method harmony
  cellmetpro trajectory scores.csv --compute-velocity --differential
  cellmetpro dashboard results/
""",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Show version and exit",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output",
    )
    parser.add_argument(
        "-y",
        "--yes",
        action="store_true",
        help="Automatically answer yes to all prompts (e.g. model downloads)",
    )

    subparsers = parser.add_subparsers(
        dest="command", metavar="", help="Available commands"
    )

    # Run command
    run_parser = subparsers.add_parser(
        "run",
        help="Run COMPASS metabolic analysis pipeline",
        description=(
            "Run COMPASS algorithm to score metabolic reactions " "from scRNA-seq data"
        ),
        formatter_class=RichHelpFormatter,
    )
    run_parser.add_argument(
        "input",
        type=Path,
        help="Input file (h5ad, csv, tsv, mtx, or rds format)",
    )
    run_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results/)",
    )
    run_parser.add_argument(
        "-m",
        "--model",
        type=str,
        default="human",
        help="Metabolic model: 'human', 'mouse', 'recon2', or path to SBML/JSON file",
    )
    run_parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize expression data before analysis",
    )
    run_parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for normalization (default: 10000)",
    )

    # Seurat-specific options
    run_parser.add_argument(
        "--assay",
        type=str,
        default=None,
        help="Seurat assay to extract (default: DefaultAssay). Only for .rds files.",
    )
    run_parser.add_argument(
        "--slot",
        type=str,
        choices=["counts", "data", "scale.data"],
        default="data",
        help="Seurat slot to extract: counts, data, or scale.data (default: data)",
    )

    # COMPASS parameters
    run_parser.add_argument(
        "--beta",
        type=float,
        default=0.95,
        help="COMPASS beta parameter (default: 0.95)",
    )
    run_parser.add_argument(
        "--lambda",
        dest="lambda_penalty",
        type=float,
        default=0.0,
        help="Penalty smoothing lambda (0-1, default: 0)",
    )
    run_parser.add_argument(
        "--n-neighbors",
        type=int,
        default=30,
        help="Number of neighbors for penalty smoothing (default: 30)",
    )

    # Performance options
    run_parser.add_argument(
        "-j",
        "--n-processes",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1)",
    )
    run_parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="Number of cells per batch for memory efficiency (default: 100)",
    )
    run_parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable progress bars",
    )
    run_parser.add_argument(
        "--microcluster",
        action="store_true",
        help="Use microclustering for large datasets",
    )
    run_parser.add_argument(
        "--cells-per-cluster",
        type=int,
        default=100,
        help="Target cells per microcluster (default: 100)",
    )
    run_parser.add_argument(
        "--use-cache",
        action="store_true",
        help="Cache optimization results for faster reruns",
    )

    # Output options
    run_parser.add_argument(
        "--compute-exchange",
        action="store_true",
        help="Also compute uptake/secretion scores",
    )
    run_parser.add_argument(
        "--output-format",
        choices=["csv", "parquet", "h5ad"],
        default="csv",
        help="Output file format (default: csv)",
    )

    # Dashboard command
    dash_parser = subparsers.add_parser(
        "dashboard",
        help="Launch interactive Streamlit dashboard",
        formatter_class=RichHelpFormatter,
    )
    dash_parser.add_argument(
        "results",
        type=Path,
        help="Path to results directory from 'run' command",
    )
    dash_parser.add_argument(
        "-p",
        "--port",
        type=int,
        default=8501,
        help="Port for Streamlit server (default: 8501)",
    )

    # Info command
    info_parser = subparsers.add_parser(
        "info",
        help="Show information about a metabolic model",
        formatter_class=RichHelpFormatter,
    )
    info_parser.add_argument(
        "model",
        type=str,
        help="Model name or path to model file",
    )

    # Differential analysis command
    diff_parser = subparsers.add_parser(
        "differential",
        help="Run differential analysis between groups",
        description="Compare metabolic activity between cell groups",
        formatter_class=RichHelpFormatter,
    )
    diff_parser.add_argument(
        "scores",
        type=Path,
        help="Reaction scores file (CSV with reactions x cells)",
    )
    diff_parser.add_argument(
        "groups",
        type=Path,
        help="Group labels file (CSV with cell_id and group columns)",
    )
    diff_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/differential"),
        help="Output directory (default: results/differential/)",
    )
    diff_parser.add_argument(
        "--group1",
        type=str,
        help="First group name (for pairwise comparison)",
    )
    diff_parser.add_argument(
        "--group2",
        type=str,
        help="Second group name (for pairwise comparison)",
    )
    diff_parser.add_argument(
        "--method",
        choices=["wilcoxon", "ttest", "mannwhitneyu", "kruskal", "anova"],
        default="wilcoxon",
        help="Statistical test method (default: wilcoxon)",
    )
    diff_parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for significance (default: 0.05)",
    )
    diff_parser.add_argument(
        "--log2fc-threshold",
        type=float,
        default=0.5,
        help="Log2 fold change threshold (default: 0.5)",
    )
    diff_parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate volcano plot",
    )
    diff_parser.add_argument(
        "--interactive",
        action="store_true",
        help="Generate interactive HTML plots (requires --plot)",
    )

    # Clustering command
    cluster_parser = subparsers.add_parser(
        "cluster",
        help="Cluster cells based on metabolic profiles",
        description="Perform dimensionality reduction and clustering on metabolic data",
        formatter_class=RichHelpFormatter,
    )
    cluster_parser.add_argument(
        "scores",
        type=Path,
        help="Reaction scores file (CSV with reactions x cells)",
    )
    cluster_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/clustering"),
        help="Output directory (default: results/clustering/)",
    )
    cluster_parser.add_argument(
        "--n-clusters",
        type=int,
        help="Number of clusters (auto-detect if not specified)",
    )
    cluster_parser.add_argument(
        "--method",
        choices=["kmeans", "leiden", "louvain"],
        default="kmeans",
        help="Clustering method (default: kmeans)",
    )
    cluster_parser.add_argument(
        "--embedding",
        choices=["umap", "tsne", "pca"],
        default="umap",
        help="Embedding method for visualization (default: umap)",
    )
    cluster_parser.add_argument(
        "--n-pcs",
        type=int,
        default=50,
        help="Number of PCA components (default: 50)",
    )
    cluster_parser.add_argument(
        "--resolution",
        type=float,
        default=1.0,
        help="Resolution for Leiden/Louvain clustering (default: 1.0)",
    )
    cluster_parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate embedding plot",
    )
    cluster_parser.add_argument(
        "--interactive",
        action="store_true",
        help="Generate interactive HTML plots (requires --plot)",
    )

    # Pathway enrichment command
    pathway_parser = subparsers.add_parser(
        "pathway",
        help="Run pathway enrichment analysis",
        description=(
            "Perform GO term or subsystem enrichment analysis "
            "on significant reactions"
        ),
        formatter_class=RichHelpFormatter,
    )
    pathway_parser.add_argument(
        "reactions",
        type=Path,
        help="File with reaction list (one per line or CSV with 'reaction' column)",
    )
    pathway_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/pathway"),
        help="Output directory (default: results/pathway/)",
    )
    pathway_parser.add_argument(
        "--model",
        type=str,
        default="human",
        help="Metabolic model for annotation (default: human)",
    )
    pathway_parser.add_argument(
        "--background",
        type=Path,
        help="Background reaction list (default: all model reactions)",
    )
    pathway_parser.add_argument(
        "--method",
        choices=["go", "subsystem"],
        default="subsystem",
        help="Enrichment type: GO terms or subsystems (default: subsystem)",
    )
    pathway_parser.add_argument(
        "--namespace",
        choices=[
            "biological_process",
            "molecular_function",
            "cellular_component",
            "all",
        ],
        default="biological_process",
        help="GO namespace for GO enrichment (default: biological_process)",
    )
    pathway_parser.add_argument(
        "--go-annotations",
        type=Path,
        help="GO annotations file (GAF format) for GO enrichment",
    )
    pathway_parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold (default: 0.05)",
    )
    pathway_parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate enrichment dotplot",
    )
    pathway_parser.add_argument(
        "--interactive",
        action="store_true",
        help="Generate interactive HTML plots (requires --plot)",
    )

    # Report generation command
    report_parser = subparsers.add_parser(
        "report",
        help="Generate HTML report from analysis results",
        description="Generate a comprehensive HTML report summarizing analysis results",
        formatter_class=RichHelpFormatter,
    )
    report_parser.add_argument(
        "results",
        type=Path,
        help="Path to results directory from 'run' command",
    )
    report_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/report.html"),
        help="Output HTML file (default: results/report.html)",
    )
    report_parser.add_argument(
        "--n-top-diff",
        type=int,
        default=20,
        help="Number of top differential reactions to show (default: 20)",
    )
    report_parser.add_argument(
        "--no-figures",
        action="store_true",
        help="Exclude embedded figures from report",
    )

    # Batch correction command
    batch_parser = subparsers.add_parser(
        "batch-correct",
        help="Apply batch correction to reaction scores",
        description="Correct batch effects in multi-sample metabolic data",
        formatter_class=RichHelpFormatter,
    )
    batch_parser.add_argument(
        "scores",
        type=Path,
        help="Reaction scores file (CSV with reactions x cells)",
    )
    batch_parser.add_argument(
        "batches",
        type=Path,
        help="Batch labels file (CSV with cell_id and batch columns)",
    )
    batch_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/batch_corrected"),
        help="Output directory (default: results/batch_corrected/)",
    )
    batch_parser.add_argument(
        "--method",
        choices=["harmony", "combat", "center"],
        default="combat",
        help="Batch correction method (default: combat)",
    )
    batch_parser.add_argument(
        "--cell-labels",
        type=Path,
        help="Cell type labels file for integration quality metrics",
    )

    # Trajectory analysis command
    traj_parser = subparsers.add_parser(
        "trajectory",
        help="Perform trajectory/pseudotime analysis",
        description="Analyze metabolic changes along cellular trajectories",
        formatter_class=RichHelpFormatter,
    )
    traj_parser.add_argument(
        "scores",
        type=Path,
        help="Reaction scores file (CSV with reactions x cells)",
    )
    traj_parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("results/trajectory"),
        help="Output directory (default: results/trajectory/)",
    )
    traj_parser.add_argument(
        "--root-cell",
        type=str,
        help="Root cell ID for pseudotime (auto-detect if not specified)",
    )
    traj_parser.add_argument(
        "--method",
        choices=["dpt", "principal_curve", "correlation"],
        default="dpt",
        help="Pseudotime inference method (default: dpt)",
    )
    traj_parser.add_argument(
        "--compute-velocity",
        action="store_true",
        help="Also compute metabolic velocity",
    )
    traj_parser.add_argument(
        "--window-size",
        type=int,
        default=50,
        help="Window size for velocity computation (default: 50)",
    )
    traj_parser.add_argument(
        "--differential",
        action="store_true",
        help="Run trajectory differential analysis",
    )
    traj_parser.add_argument(
        "--plot",
        action="store_true",
        help="Generate trajectory plots",
    )

    return parser


def run_analysis(args: argparse.Namespace) -> int:
    """Run the COMPASS metabolic analysis pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """

    from cellmetpro.core import (
        CompassConfig,
        CompassScorer,
        DataLoader,
        microcluster,
        normalize_expression,
        to_dataframe,
        unpool_results,
    )
    from cellmetpro.models import load_gem

    logger.info("CellMetPro COMPASS Analysis")
    logger.info("=" * 50)
    logger.info(f"Input: {args.input}")
    logger.info(f"Model: {args.model}")
    logger.info(f"Output: {args.output}")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Load expression data
    logger.info("Loading expression data...")
    loader = DataLoader(args.input)

    # Handle Seurat files with special options
    if str(args.input).endswith(".rds"):
        logger.info("Detected Seurat RDS file")
        adata = loader.load_seurat(assay=args.assay, slot=args.slot)
    else:
        adata = loader.load()
    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    # Normalize if requested
    if args.normalize:
        logger.info("Normalizing expression data...")
        adata = normalize_expression(adata, target_sum=args.target_sum)

    # Convert to DataFrame (genes x cells)
    expression_df = to_dataframe(adata, genes_as_rows=True)

    # Load metabolic model
    logger.info(f"Loading metabolic model: {args.model}")
    model = load_gem(args.model, auto_confirm=getattr(args, "yes", False))
    logger.info(
        f"Model: {model.id} ({len(model.reactions)} reactions, "
        f"{len(model.metabolites)} metabolites, {len(model.genes)} genes)"
    )

    # Configure COMPASS
    config = CompassConfig(
        beta=args.beta,
        lambda_penalty=args.lambda_penalty,
        n_neighbors=args.n_neighbors,
        n_processes=args.n_processes,
        batch_size=args.batch_size,
        show_progress=not args.no_progress,
    )

    # Handle microclustering for large datasets
    microcluster_result = None
    if args.microcluster:
        from cellmetpro.core import MicroclusterConfig

        target = args.cells_per_cluster
        logger.info(f"Microclustering cells (target: {target} cells/cluster)...")
        mc_config = MicroclusterConfig(cells_per_cluster=args.cells_per_cluster)
        microcluster_result = microcluster(expression_df, mc_config)
        logger.info(f"Created {microcluster_result.n_clusters} microclusters")

        # Use pooled expression
        expression_df = microcluster_result.pooled_expression

    # Run COMPASS
    logger.info("Running COMPASS algorithm...")
    scorer = CompassScorer(model, expression_df, config)
    result = scorer.score(compute_exchange=args.compute_exchange)

    # Unpool results if microclustering was used
    if microcluster_result is not None:
        logger.info("Unpooling results to individual cells...")
        result.reaction_penalties = unpool_results(
            result.reaction_penalties, microcluster_result
        )
        result.reaction_scores = unpool_results(
            result.reaction_scores, microcluster_result
        )

    # Save results
    logger.info(f"Saving results to {args.output}/")

    if args.output_format == "csv":
        result.reaction_penalties.to_csv(args.output / "reaction_penalties.csv")
        result.reaction_scores.to_csv(args.output / "reaction_scores.csv")
        if result.uptake_scores is not None:
            result.uptake_scores.to_csv(args.output / "uptake_scores.csv")
        if result.secretion_scores is not None:
            result.secretion_scores.to_csv(args.output / "secretion_scores.csv")
    elif args.output_format == "parquet":
        result.reaction_penalties.to_parquet(args.output / "reaction_penalties.parquet")
        result.reaction_scores.to_parquet(args.output / "reaction_scores.parquet")
        if result.uptake_scores is not None:
            result.uptake_scores.to_parquet(args.output / "uptake_scores.parquet")
        if result.secretion_scores is not None:
            result.secretion_scores.to_parquet(args.output / "secretion_scores.parquet")
    elif args.output_format == "h5ad":
        import anndata as ad

        # Store results in AnnData format
        result_adata = ad.AnnData(result.reaction_scores.T)
        result_adata.layers["penalties"] = result.reaction_penalties.T.values
        result_adata.write(args.output / "compass_results.h5ad")

    # Save config
    import json

    config_dict = {
        "input": str(args.input),
        "model": args.model,
        "beta": config.beta,
        "lambda_penalty": config.lambda_penalty,
        "n_neighbors": config.n_neighbors,
        "n_processes": config.n_processes,
        "batch_size": config.batch_size,
        "microcluster": args.microcluster,
        "cells_per_cluster": args.cells_per_cluster if args.microcluster else None,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_reactions": len(result.reaction_scores),
    }
    with open(args.output / "config.json", "w") as f:
        json.dump(config_dict, f, indent=2)

    logger.info("Analysis complete!")
    logger.info(f"Results saved to: {args.output}/")
    logger.info(f"  - reaction_penalties.{args.output_format}")
    logger.info(f"  - reaction_scores.{args.output_format}")
    if args.compute_exchange:
        logger.info(f"  - uptake_scores.{args.output_format}")
        logger.info(f"  - secretion_scores.{args.output_format}")

    return 0


def run_dashboard(args: argparse.Namespace) -> int:
    """Launch the interactive dashboard.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import subprocess

    if not args.results.exists():
        logger.error(f"Results directory not found: {args.results}")
        return 1

    logger.info(f"Launching dashboard for results in {args.results}")
    logger.info(f"Running on port {args.port}")

    # Check if streamlit is available
    import importlib.util

    if importlib.util.find_spec("streamlit") is None:
        logger.error(
            "Streamlit not installed. Install with: pip install cellmetpro[dashboard]"
        )
        return 1

    # Get path to dashboard script
    from cellmetpro.visualization import dashboard

    dashboard_path = Path(dashboard.__file__)

    # Run streamlit
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(dashboard_path),
        "--server.port",
        str(args.port),
        "--",
        str(args.results),
    ]

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Dashboard failed: {e}")
        return 1
    except KeyboardInterrupt:
        logger.info("Dashboard stopped")

    return 0


def run_differential(args: argparse.Namespace) -> int:
    """Run differential analysis between groups.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import pandas as pd

    from cellmetpro.analysis.differential import DifferentialAnalysis

    logger.info("CellMetPro Differential Analysis")
    logger.info("=" * 50)

    # Load reaction scores
    logger.info(f"Loading reaction scores from {args.scores}")
    scores = pd.read_csv(args.scores, index_col=0)
    logger.info(f"Loaded {scores.shape[0]} reactions x {scores.shape[1]} cells")

    # Load group labels
    logger.info(f"Loading group labels from {args.groups}")
    sep = "\t" if str(args.groups).endswith((".tsv", ".txt")) else ","
    groups_df = pd.read_csv(args.groups, sep=sep)

    # Handle different group file formats
    group_col = "group" if "group" in groups_df.columns else groups_df.columns[1]
    if "cell_id" in groups_df.columns:
        groups = pd.Series(groups_df[group_col].values, index=groups_df["cell_id"])
    elif groups_df.shape[1] == 2:
        groups = pd.Series(groups_df.iloc[:, 1].values, index=groups_df.iloc[:, 0])
    else:
        # Assume first column is index
        groups = pd.Series(groups_df.iloc[:, 0].values, index=groups_df.index)

    # Normalize barcodes: strip 10x GEM well suffix (e.g. -1) if scores use bare barcodes
    if len(scores.columns.intersection(groups.index)) == 0:
        stripped = groups.index.str.replace(r"-\d+$", "", regex=True)
        if len(scores.columns.intersection(stripped)) > 0:
            groups.index = stripped

    # Check overlap
    common_cells = scores.columns.intersection(groups.index)
    logger.info(f"Found {len(common_cells)} cells in common")

    if len(common_cells) == 0:
        logger.error("No common cells between scores and groups")
        return 1

    unique_groups = groups[common_cells].unique()
    logger.info(f"Groups: {list(unique_groups)}")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Create differential analysis object
    da = DifferentialAnalysis(scores, groups)

    # Determine analysis type
    if args.group1 and args.group2:
        # Pairwise comparison
        logger.info(f"Running pairwise comparison: {args.group1} vs {args.group2}")
        result = da.compare_groups(args.group1, args.group2, method=args.method)
        output_file = args.output / f"differential_{args.group1}_vs_{args.group2}.csv"
    elif len(unique_groups) == 2:
        # Auto pairwise for 2 groups
        g1, g2 = sorted(unique_groups)
        logger.info(f"Running pairwise comparison: {g1} vs {g2}")
        result = da.compare_groups(g1, g2, method=args.method)
        output_file = args.output / f"differential_{g1}_vs_{g2}.csv"
    elif args.method in ["kruskal", "anova"]:
        # Multi-group comparison
        logger.info(f"Running multi-group comparison ({args.method})")
        result = da.compare_multiple_groups(method=args.method)
        output_file = args.output / "differential_multigroup.csv"
    else:
        # Default to all pairwise
        logger.info("Running all pairwise comparisons")
        result = da.all_pairwise_comparisons(method=args.method)
        output_file = args.output / "differential_all_pairwise.csv"

    # Filter significant results
    if "padj_bh" in result.columns:
        n_significant = (result["padj_bh"] < args.fdr_threshold).sum()
        fdr = args.fdr_threshold
        logger.info(f"Significant reactions (FDR < {fdr}): {n_significant}")

    # Save results
    result.to_csv(output_file, index=False)
    logger.info(f"Results saved to {output_file}")

    # Generate volcano plot if requested
    if args.plot and "log2fc" in result.columns and "padj_bh" in result.columns:
        logger.info("Generating volcano plot...")

        if args.interactive:
            from cellmetpro.visualization import plot_volcano_interactive

            plot_volcano_interactive(
                result,
                log2fc_threshold=args.log2fc_threshold,
                pvalue_threshold=args.fdr_threshold,
                save=str(args.output / "volcano_plot.html"),
            )
            logger.info(
                f"Interactive volcano plot saved to {args.output / 'volcano_plot.html'}"
            )
            _open_file(str(args.output / "volcano_plot.html"))
        else:
            import matplotlib

            matplotlib.use("Agg")
            from cellmetpro.visualization import plot_volcano

            plot_volcano(
                result,
                log2fc_threshold=args.log2fc_threshold,
                pvalue_threshold=args.fdr_threshold,
                save=str(args.output / "volcano_plot.png"),
            )
            logger.info(f"Volcano plot saved to {args.output / 'volcano_plot.png'}")
            _open_file(str(args.output / "volcano_plot.png"))

    logger.info("Differential analysis complete!")
    return 0


def run_cluster(args: argparse.Namespace) -> int:
    """Run clustering analysis on metabolic profiles.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import pandas as pd

    from cellmetpro.analysis.clustering import (
        MetabolicClustering,
        find_optimal_clusters,
    )

    logger.info("CellMetPro Clustering Analysis")
    logger.info("=" * 50)

    # Load reaction scores
    logger.info(f"Loading reaction scores from {args.scores}")
    scores = pd.read_csv(args.scores, index_col=0)
    logger.info(f"Loaded {scores.shape[0]} reactions x {scores.shape[1]} cells")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Create clustering object
    mc = MetabolicClustering(scores, n_clusters=args.n_clusters)

    # Compute PCA
    logger.info(f"Computing PCA with {args.n_pcs} components...")
    pca_result = mc.compute_pca(n_components=args.n_pcs)
    logger.info(f"PCA complete: {pca_result.shape}")

    # Auto-detect optimal clusters if not specified
    if args.n_clusters is None and args.method == "kmeans":
        logger.info("Finding optimal number of clusters...")
        n_clusters = find_optimal_clusters(
            pca_result, max_clusters=15, method="silhouette"
        )
        logger.info(f"Optimal clusters: {n_clusters}")
        mc.n_clusters = n_clusters

    # Compute embedding
    logger.info(f"Computing {args.embedding.upper()} embedding...")
    if args.embedding == "umap":
        try:
            embedding = mc.compute_umap()
        except ImportError:
            logger.warning("UMAP not available, falling back to t-SNE")
            embedding = mc.compute_tsne()
    elif args.embedding == "tsne":
        embedding = mc.compute_tsne()
    else:  # pca
        embedding = pca_result[:, :2]
        mc.embedding = embedding

    logger.info(f"Embedding complete: {embedding.shape}")

    # Run clustering
    logger.info(f"Running {args.method} clustering...")
    if args.method in ["leiden", "louvain"]:
        labels = mc.cluster(method=args.method, resolution=args.resolution)
    else:
        labels = mc.cluster(method=args.method)

    n_clusters = len(set(labels))
    logger.info(f"Found {n_clusters} clusters")

    # Export results
    result_df = mc.to_dataframe()
    result_df.to_csv(args.output / "clustering_results.csv", index=False)
    logger.info(f"Results saved to {args.output / 'clustering_results.csv'}")

    # Get cluster markers
    logger.info("Computing cluster markers...")
    markers = mc.get_cluster_markers(n_top=20)
    markers.to_csv(args.output / "cluster_markers.csv", index=False)
    logger.info(f"Markers saved to {args.output / 'cluster_markers.csv'}")

    # Generate plot if requested
    if args.plot:
        logger.info("Generating embedding plot...")

        if args.interactive:
            from cellmetpro.visualization import plot_embedding_interactive

            plot_embedding_interactive(
                embedding,
                color=labels.astype(str),
                labels=(
                    result_df["cell_id"].values
                    if "cell_id" in result_df.columns
                    else None
                ),
                title=f"Metabolic Clustering ({args.method})",
                xlabel=f"{args.embedding.upper()} 1",
                ylabel=f"{args.embedding.upper()} 2",
                save=str(args.output / "embedding_plot.html"),
            )
            logger.info(
                f"Interactive plot saved to {args.output / 'embedding_plot.html'}"
            )
            _open_file(str(args.output / "embedding_plot.html"))
        else:
            import matplotlib

            matplotlib.use("Agg")
            from cellmetpro.visualization import plot_embedding

            plot_embedding(
                embedding,
                color=labels.astype(str),
                title=f"Metabolic Clustering ({args.method})",
                xlabel=f"{args.embedding.upper()} 1",
                ylabel=f"{args.embedding.upper()} 2",
                save=str(args.output / "embedding_plot.png"),
            )
            logger.info(f"Plot saved to {args.output / 'embedding_plot.png'}")
            _open_file(str(args.output / "embedding_plot.png"))

    logger.info("Clustering analysis complete!")
    return 0


def run_pathway(args: argparse.Namespace) -> int:
    """Run pathway enrichment analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import pandas as pd

    from cellmetpro.analysis.pathway import (
        GOEnrichmentAnalyzer,
        load_go_annotations_from_gaf,
        subsystem_enrichment,
    )
    from cellmetpro.models import get_subsystem_reactions, load_gem

    logger.info("CellMetPro Pathway Enrichment Analysis")
    logger.info("=" * 50)

    # Load reaction list
    logger.info(f"Loading reactions from {args.reactions}")
    if args.reactions.suffix == ".csv":
        rxn_df = pd.read_csv(args.reactions)
        if "reaction" in rxn_df.columns:
            reactions = list(rxn_df["reaction"])
        else:
            reactions = list(rxn_df.iloc[:, 0])
    else:
        with open(args.reactions) as f:
            reactions = [line.strip() for line in f if line.strip()]

    logger.info(f"Loaded {len(reactions)} reactions")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Load metabolic model
    logger.info(f"Loading metabolic model: {args.model}")
    model = load_gem(args.model, auto_confirm=getattr(args, "yes", False))

    # Get background reactions
    if args.background:
        logger.info(f"Loading background from {args.background}")
        if args.background.suffix == ".csv":
            bg_df = pd.read_csv(args.background)
            if "reaction" in bg_df.columns:
                background = list(bg_df["reaction"])
            else:
                background = list(bg_df.iloc[:, 0])
        else:
            with open(args.background) as f:
                background = [line.strip() for line in f if line.strip()]
    else:
        background = [r.id for r in model.reactions]

    logger.info(f"Background: {len(background)} reactions")

    if args.method == "subsystem":
        # Subsystem enrichment
        logger.info("Running subsystem enrichment...")

        # Get subsystem mapping
        subsystem_mapping = get_subsystem_reactions(model)

        # Run enrichment using standalone function
        result = subsystem_enrichment(
            significant_reactions=reactions,
            background=background,
            subsystem_mapping=subsystem_mapping,
        )

        # Filter by FDR
        if len(result) > 0 and "padj" in result.columns:
            significant = result[result["padj"] < args.fdr_threshold]
            fdr = args.fdr_threshold
            logger.info(f"Significant subsystems (FDR < {fdr}): {len(significant)}")
        else:
            logger.info("No enriched subsystems found")

        # Save results
        result.to_csv(args.output / "subsystem_enrichment.csv", index=False)
        logger.info(f"Results saved to {args.output / 'subsystem_enrichment.csv'}")

    else:  # GO enrichment
        logger.info("Running GO term enrichment...")

        # Check if GO annotations file is provided
        if args.go_annotations is None:
            logger.warning(
                "GO enrichment requires GO annotations file. "
                "Please provide a GAF file using --go-annotations flag. "
                "Falling back to subsystem enrichment."
            )

            # Fallback to subsystem enrichment
            subsystem_mapping = get_subsystem_reactions(model)
            result = subsystem_enrichment(
                significant_reactions=reactions,
                background=background,
                subsystem_mapping=subsystem_mapping,
            )

            if len(result) > 0 and "padj" in result.columns:
                n_sig = (result["padj"] < args.fdr_threshold).sum()
                fdr = args.fdr_threshold
                logger.info(f"Significant subsystems (FDR < {fdr}): {n_sig}")

            result.to_csv(args.output / "subsystem_enrichment.csv", index=False)
            logger.info(f"Results saved to {args.output / 'subsystem_enrichment.csv'}")
        else:
            # Load GO annotations from GAF file
            logger.info(f"Loading GO annotations from {args.go_annotations}")
            go_annotations = load_go_annotations_from_gaf(str(args.go_annotations))
            logger.info(f"Loaded {len(go_annotations)} GO annotations")

            # Create GO enrichment analyzer
            analyzer = GOEnrichmentAnalyzer(model, go_annotations)

            # Map namespace argument to GO namespace codes
            namespace_map = {
                "biological_process": "BP",
                "molecular_function": "MF",
                "cellular_component": "CC",
                "all": None,
            }
            namespace = namespace_map.get(args.namespace)

            # Run GO enrichment
            result = analyzer.enrich_reactions(
                significant_reactions=set(reactions),
                background_reactions=set(background),
                namespace=namespace,
            )

            if len(result) > 0 and "padj" in result.columns:
                n_sig = (result["padj"] < args.fdr_threshold).sum()
                fdr = args.fdr_threshold
                logger.info(f"Significant GO terms (FDR < {fdr}): {n_sig}")
            else:
                logger.info("No enriched GO terms found")

            result.to_csv(args.output / "go_enrichment.csv", index=False)
            logger.info(f"Results saved to {args.output / 'go_enrichment.csv'}")

    # Generate plot if requested
    if args.plot and len(result) > 0:
        logger.info("Generating enrichment dotplot...")

        # Prepare data for plotting
        if args.method == "subsystem":
            plot_result = result.copy()
            has_go_name = "go_name" in plot_result.columns
            has_pathway = "pathway" in plot_result.columns
            if not has_go_name and has_pathway:
                plot_result["go_name"] = plot_result["pathway"]
                plot_result["go_term"] = plot_result["pathway"]
        else:
            plot_result = result

        if args.interactive:
            from cellmetpro.visualization import plot_enrichment_interactive

            plot_enrichment_interactive(
                plot_result,
                pvalue_threshold=args.fdr_threshold,
                title=f"{'Subsystem' if args.method == 'subsystem' else 'GO'} Enrichment",
                save=str(args.output / "enrichment_plot.html"),
            )
            logger.info(
                f"Interactive plot saved to {args.output / 'enrichment_plot.html'}"
            )
            _open_file(str(args.output / "enrichment_plot.html"))
        else:
            import matplotlib

            matplotlib.use("Agg")
            from cellmetpro.visualization import plot_enrichment_dotplot

            plot_enrichment_dotplot(
                plot_result,
                pvalue_threshold=args.fdr_threshold,
                title=f"{'Subsystem' if args.method == 'subsystem' else 'GO'} Enrichment",
                save=str(args.output / "enrichment_plot.png"),
            )
            logger.info(f"Plot saved to {args.output / 'enrichment_plot.png'}")
            _open_file(str(args.output / "enrichment_plot.png"))

    logger.info("Pathway enrichment analysis complete!")
    return 0


def run_report(args: argparse.Namespace) -> int:
    """Generate HTML report from analysis results.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    from cellmetpro.reporting import generate_html_report

    logger.info("CellMetPro Report Generation")
    logger.info("=" * 50)

    if not args.results.exists():
        logger.error(f"Results directory not found: {args.results}")
        return 1

    logger.info(f"Generating report from: {args.results}")

    output_path = args.output

    try:
        report_path = generate_html_report(
            args.results,
            output_path,
            n_top_diff=args.n_top_diff,
            include_figures=not args.no_figures,
        )
        logger.info(f"Report generated: {report_path}")
    except Exception as e:
        logger.error(f"Failed to generate report: {e}")
        return 1

    return 0


def run_batch_correct(args: argparse.Namespace) -> int:
    """Run batch correction on reaction scores.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import json

    import pandas as pd

    from cellmetpro.core.batch_correction import (
        center_batches,
        combat_correct,
        compute_integration_metrics,
        harmony_integrate,
    )

    logger.info("CellMetPro Batch Correction")
    logger.info("=" * 50)

    # Load reaction scores
    logger.info(f"Loading reaction scores from {args.scores}")
    scores = pd.read_csv(args.scores, index_col=0)
    logger.info(f"Loaded {scores.shape[0]} reactions x {scores.shape[1]} cells")

    # Load batch labels
    logger.info(f"Loading batch labels from {args.batches}")
    sep = "\t" if str(args.batches).endswith((".tsv", ".txt")) else ","
    batch_df = pd.read_csv(args.batches, sep=sep)

    if "cell_id" in batch_df.columns and "batch" in batch_df.columns:
        batch_labels = pd.Series(batch_df["batch"].values, index=batch_df["cell_id"])
    elif batch_df.shape[1] == 2:
        batch_labels = pd.Series(batch_df.iloc[:, 1].values, index=batch_df.iloc[:, 0])
    else:
        logger.error("Batch file must have cell_id and batch columns")
        return 1

    # Normalize barcodes: strip 10x GEM well suffix (e.g. -1) if scores use bare barcodes
    if len(scores.columns.intersection(batch_labels.index)) == 0:
        stripped = batch_labels.index.str.replace(r"-\d+$", "", regex=True)
        if len(scores.columns.intersection(stripped)) > 0:
            batch_labels.index = stripped

    # Align cells
    common_cells = scores.columns.intersection(batch_labels.index)
    logger.info(f"Found {len(common_cells)} cells in common")

    if len(common_cells) == 0:
        logger.error("No common cells between scores and batch labels")
        return 1

    scores_aligned = scores[common_cells]
    batch_aligned = batch_labels[common_cells]

    unique_batches = batch_aligned.unique()
    logger.info(f"Batches: {list(unique_batches)}")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Compute pre-correction metrics
    logger.info("Computing pre-correction integration metrics...")
    pre_metrics = compute_integration_metrics(scores_aligned, batch_aligned)
    logger.info(f"  Batch mixing: {pre_metrics['batch_mixing']:.3f}")
    logger.info(f"  Batch silhouette: {pre_metrics['batch_silhouette']:.3f}")

    # Run batch correction
    logger.info(f"Running {args.method} batch correction...")

    if args.method == "harmony":
        try:
            corrected = harmony_integrate(scores_aligned, batch_aligned)
        except ImportError:
            logger.error("harmonypy not installed. Install with: pip install harmonypy")
            return 1
    elif args.method == "combat":
        corrected = combat_correct(scores_aligned, batch_aligned)
    else:  # center
        corrected = center_batches(scores_aligned, batch_aligned)

    logger.info("Batch correction complete")

    # Compute post-correction metrics
    logger.info("Computing post-correction integration metrics...")
    post_metrics = compute_integration_metrics(corrected, batch_aligned)
    logger.info(f"  Batch mixing: {post_metrics['batch_mixing']:.3f}")
    logger.info(f"  Batch silhouette: {post_metrics['batch_silhouette']:.3f}")

    # Save results
    corrected.to_csv(args.output / "corrected_scores.csv")
    logger.info(f"Corrected scores saved to {args.output / 'corrected_scores.csv'}")

    # Save metrics
    metrics = {
        "method": args.method,
        "n_cells": len(common_cells),
        "n_batches": len(unique_batches),
        "batches": list(unique_batches),
        "pre_correction": pre_metrics,
        "post_correction": post_metrics,
    }
    with open(args.output / "batch_correction_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    logger.info(f"Metrics saved to {args.output / 'batch_correction_metrics.json'}")

    logger.info("Batch correction complete!")
    return 0


def run_trajectory(args: argparse.Namespace) -> int:
    """Run trajectory/pseudotime analysis.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    import pandas as pd

    from cellmetpro.analysis.trajectory import (
        compute_metabolic_velocity,
        compute_pseudotime,
        trajectory_differential,
    )

    logger.info("CellMetPro Trajectory Analysis")
    logger.info("=" * 50)

    # Load reaction scores
    logger.info(f"Loading reaction scores from {args.scores}")
    scores = pd.read_csv(args.scores, index_col=0)
    logger.info(f"Loaded {scores.shape[0]} reactions x {scores.shape[1]} cells")

    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)

    # Compute pseudotime
    logger.info(f"Computing pseudotime (method: {args.method})...")
    pseudotime = compute_pseudotime(
        scores,
        root_cell=args.root_cell,
        method=args.method,
    )
    logger.info(f"Pseudotime computed for {len(pseudotime)} cells")

    # Save pseudotime
    pt_df = pd.DataFrame({"cell_id": pseudotime.index, "pseudotime": pseudotime.values})
    pt_df.to_csv(args.output / "pseudotime.csv", index=False)
    logger.info(f"Pseudotime saved to {args.output / 'pseudotime.csv'}")

    # Compute velocity if requested
    if args.compute_velocity:
        logger.info(f"Computing metabolic velocity (window: {args.window_size})...")
        velocity = compute_metabolic_velocity(
            scores, pseudotime, window_size=args.window_size
        )
        velocity.to_csv(args.output / "metabolic_velocity.csv")
        logger.info(f"Velocity saved to {args.output / 'metabolic_velocity.csv'}")

        # Find top dynamic reactions
        dynamic_reactions = velocity.abs().mean(axis=1).nlargest(20)
        dynamic_df = pd.DataFrame(
            {
                "reaction": dynamic_reactions.index,
                "mean_abs_velocity": dynamic_reactions.values,
            }
        )
        dynamic_df.to_csv(args.output / "dynamic_reactions.csv", index=False)
        logger.info("Top 20 dynamic reactions saved")

    # Run trajectory differential if requested
    if args.differential:
        logger.info("Running trajectory differential analysis...")
        diff_results = trajectory_differential(scores, pseudotime)
        diff_results.to_csv(args.output / "trajectory_differential.csv", index=False)
        logger.info(
            f"Differential results saved to {args.output / 'trajectory_differential.csv'}"
        )

        # Report summary
        increasing = (diff_results["trend"] == "increasing").sum()
        decreasing = (diff_results["trend"] == "decreasing").sum()
        logger.info(f"  Increasing reactions: {increasing}")
        logger.info(f"  Decreasing reactions: {decreasing}")

    # Generate plots if requested
    if args.plot:
        logger.info("Generating trajectory plots...")

        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # Plot pseudotime distribution
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Histogram of pseudotime
        axes[0].hist(pseudotime.values, bins=50, edgecolor="black", alpha=0.7)
        axes[0].set_xlabel("Pseudotime")
        axes[0].set_ylabel("Number of cells")
        axes[0].set_title("Pseudotime Distribution")

        # Top reactions along pseudotime
        if args.differential:
            # Get top changing reactions
            top_increasing = diff_results[diff_results["trend"] == "increasing"].head(3)
            top_decreasing = diff_results[diff_results["trend"] == "decreasing"].head(3)

            order = pseudotime.argsort()
            pt_sorted = pseudotime.iloc[order].values

            for _, row in top_increasing.iterrows():
                rxn = row["reaction"]
                if rxn in scores.index:
                    y = scores.loc[rxn].iloc[order].values
                    axes[1].plot(pt_sorted, y, label=f"{rxn} (+)", alpha=0.7)

            for _, row in top_decreasing.iterrows():
                rxn = row["reaction"]
                if rxn in scores.index:
                    y = scores.loc[rxn].iloc[order].values
                    axes[1].plot(pt_sorted, y, label=f"{rxn} (-)", alpha=0.7, ls="--")

            axes[1].set_xlabel("Pseudotime")
            axes[1].set_ylabel("Reaction Score")
            axes[1].set_title("Top Trajectory-Associated Reactions")
            axes[1].legend(fontsize=8)
        else:
            # Just show variance along pseudotime
            axes[1].text(
                0.5,
                0.5,
                "Run with --differential\nto show trajectory genes",
                ha="center",
                va="center",
                transform=axes[1].transAxes,
            )

        plt.tight_layout()
        plt.savefig(args.output / "trajectory_plot.png", dpi=150)
        plt.close()
        logger.info(f"Plot saved to {args.output / 'trajectory_plot.png'}")
        _open_file(str(args.output / "trajectory_plot.png"))

    logger.info("Trajectory analysis complete!")
    return 0


def show_model_info(args: argparse.Namespace) -> int:
    """Show information about a metabolic model.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    from cellmetpro.models import get_subsystem_reactions, load_gem

    try:
        model = load_gem(args.model, auto_confirm=getattr(args, "yes", False))
    except (FileNotFoundError, ValueError) as e:
        logger.error(str(e))
        return 1

    print(f"\nModel Information: {model.id}")
    print("=" * 50)
    print(f"Name: {model.name or 'N/A'}")
    print(f"Reactions: {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes: {len(model.genes)}")
    print(f"Exchange reactions: {len(model.exchanges)}")

    # Count reactions with GPR rules
    n_with_gpr = sum(1 for r in model.reactions if r.gene_reaction_rule)
    print(f"Reactions with GPR rules: {n_with_gpr}")

    # Subsystems
    subsystems = get_subsystem_reactions(model)
    print(f"Subsystems: {len(subsystems)}")

    if subsystems:
        print("\nTop 10 subsystems by reaction count:")
        sorted_subsystems = sorted(
            subsystems.items(), key=lambda x: len(x[1]), reverse=True
        )[:10]
        for name, rxns in sorted_subsystems:
            print(f"  {name}: {len(rxns)} reactions")

    return 0


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the CLI.

    Parameters
    ----------
    argv : list[str], optional
        Command-line arguments. If None, uses sys.argv.

    Returns
    -------
    int
        Exit code.
    """
    import argcomplete

    parser = create_parser()
    argcomplete.autocomplete(parser)
    args = parser.parse_args(argv)

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    if args.version:
        from cellmetpro import __version__

        print(f"cellmetpro {__version__}")
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    try:
        if args.command == "run":
            return run_analysis(args)
        elif args.command == "dashboard":
            return run_dashboard(args)
        elif args.command == "info":
            return show_model_info(args)
        elif args.command == "differential":
            return run_differential(args)
        elif args.command == "cluster":
            return run_cluster(args)
        elif args.command == "pathway":
            return run_pathway(args)
        elif args.command == "report":
            return run_report(args)
        elif args.command == "batch-correct":
            return run_batch_correct(args)
        elif args.command == "trajectory":
            return run_trajectory(args)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        return 130
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
