"""Command-line interface for cellMetPro.

This module provides CLI commands for running metabolic analysis
pipelines from the command line.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def create_parser() -> argparse.ArgumentParser:
    """Create the argument parser for the CLI.

    Returns
    -------
    argparse.ArgumentParser
        Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        prog="cellmetpro",
        description="Cellular Metabolic Profiler - Analyze metabolic profiles from scRNA-seq data",
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Show version and exit",
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Run command
    run_parser = subparsers.add_parser(
        "run",
        help="Run metabolic analysis pipeline",
    )
    run_parser.add_argument(
        "input",
        type=Path,
        help="Input file (h5ad, csv, or tsv)",
    )
    run_parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("results"),
        help="Output directory (default: results)",
    )
    run_parser.add_argument(
        "-m", "--model",
        type=str,
        default="human",
        help="Metabolic model: 'human', 'mouse', or path to SBML file",
    )
    run_parser.add_argument(
        "--cluster-key",
        type=str,
        default="cluster",
        help="Column name for cluster labels in metadata",
    )

    # Dashboard command
    dash_parser = subparsers.add_parser(
        "dashboard",
        help="Launch interactive dashboard",
    )
    dash_parser.add_argument(
        "results",
        type=Path,
        help="Path to results directory from 'run' command",
    )
    dash_parser.add_argument(
        "-p", "--port",
        type=int,
        default=8501,
        help="Port for Streamlit server",
    )

    return parser


def run_analysis(args: argparse.Namespace) -> int:
    """Run the metabolic analysis pipeline.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    int
        Exit code (0 for success).
    """
    print(f"Running analysis on {args.input}")
    print(f"Using model: {args.model}")
    print(f"Output directory: {args.output}")
    # TODO: Implement analysis pipeline
    raise NotImplementedError("Analysis pipeline not yet implemented")


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
    print(f"Launching dashboard for results in {args.results}")
    print(f"Running on port {args.port}")
    # TODO: Implement dashboard launch
    raise NotImplementedError("Dashboard not yet implemented")


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
    parser = create_parser()
    args = parser.parse_args(argv)

    if args.version:
        from src import __version__
        print(f"cellmetpro {__version__}")
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "run":
        return run_analysis(args)
    elif args.command == "dashboard":
        return run_dashboard(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
