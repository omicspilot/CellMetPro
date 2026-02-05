"""Report generation for CellMetPro analysis results.

This module provides functions to generate HTML and PDF reports
summarizing metabolic analysis results.
"""

from __future__ import annotations

import base64
import json
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# HTML template for the report
HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CellMetPro Analysis Report</title>
    <style>
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --accent-color: #e74c3c;
            --background-color: #f8f9fa;
            --card-background: #ffffff;
            --text-color: #333333;
            --border-color: #dee2e6;
        }
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background-color: var(--background-color);
            color: var(--text-color);
            line-height: 1.6;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }
        header {
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 40px 20px;
            text-align: center;
            margin-bottom: 30px;
            border-radius: 8px;
        }
        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        header p {
            opacity: 0.9;
            font-size: 1.1em;
        }
        .card {
            background: var(--card-background);
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            overflow: hidden;
        }
        .card-header {
            background: var(--primary-color);
            color: white;
            padding: 15px 20px;
            font-size: 1.2em;
            font-weight: 600;
        }
        .card-body {
            padding: 20px;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }
        .stat-box {
            background: var(--background-color);
            border-radius: 8px;
            padding: 20px;
            text-align: center;
            border: 1px solid var(--border-color);
        }
        .stat-value {
            font-size: 2em;
            font-weight: bold;
            color: var(--secondary-color);
        }
        .stat-label {
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }
        th, td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }
        th {
            background: var(--background-color);
            font-weight: 600;
            color: var(--primary-color);
        }
        tr:hover {
            background: #f5f5f5;
        }
        .figure-container {
            margin: 20px 0;
            text-align: center;
        }
        .figure-container img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .figure-caption {
            margin-top: 10px;
            color: #666;
            font-style: italic;
        }
        .badge {
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: 500;
        }
        .badge-success { background: #d4edda; color: #155724; }
        .badge-warning { background: #fff3cd; color: #856404; }
        .badge-danger { background: #f8d7da; color: #721c24; }
        .badge-info { background: #d1ecf1; color: #0c5460; }
        footer {
            text-align: center;
            padding: 30px;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid var(--border-color);
            margin-top: 40px;
        }
        .toc {
            background: var(--background-color);
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 30px;
        }
        .toc h3 {
            margin-bottom: 15px;
            color: var(--primary-color);
        }
        .toc ul {
            list-style: none;
        }
        .toc li {
            padding: 5px 0;
        }
        .toc a {
            color: var(--secondary-color);
            text-decoration: none;
        }
        .toc a:hover {
            text-decoration: underline;
        }
        @media print {
            .card {
                break-inside: avoid;
            }
            header {
                background: var(--primary-color) !important;
                -webkit-print-color-adjust: exact;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>CellMetPro Analysis Report</h1>
            <p>Generated on {{ generation_date }}</p>
        </header>

        <div class="toc">
            <h3>Table of Contents</h3>
            <ul>
                <li><a href="#overview">1. Analysis Overview</a></li>
                <li><a href="#data-summary">2. Data Summary</a></li>
                {% if has_differential %}<li><a href="#differential">3. Differential Analysis</a></li>{% endif %}
                {% if has_enrichment %}<li><a href="#enrichment">4. Pathway Enrichment</a></li>{% endif %}
                {% if has_clustering %}<li><a href="#clustering">5. Clustering Results</a></li>{% endif %}
                {% if has_figures %}<li><a href="#figures">6. Visualizations</a></li>{% endif %}
            </ul>
        </div>

        <section id="overview" class="card">
            <div class="card-header">1. Analysis Overview</div>
            <div class="card-body">
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{{ n_cells }}</div>
                        <div class="stat-label">Cells Analyzed</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{{ n_reactions }}</div>
                        <div class="stat-label">Reactions Scored</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{{ n_genes }}</div>
                        <div class="stat-label">Genes Used</div>
                    </div>
                    <div class="stat-box">
                        <div class="stat-value">{{ model_name }}</div>
                        <div class="stat-label">Metabolic Model</div>
                    </div>
                </div>
                {% if config_table %}
                <h4>Configuration Parameters</h4>
                {{ config_table }}
                {% endif %}
            </div>
        </section>

        <section id="data-summary" class="card">
            <div class="card-header">2. Data Summary</div>
            <div class="card-body">
                <h4>Reaction Score Statistics</h4>
                {{ scores_summary_table }}

                {% if top_variable_reactions %}
                <h4 style="margin-top: 20px;">Top Variable Reactions</h4>
                {{ top_variable_reactions }}
                {% endif %}
            </div>
        </section>

        {% if has_differential %}
        <section id="differential" class="card">
            <div class="card-header">3. Differential Analysis</div>
            <div class="card-body">
                {% for comparison, diff_table in differential_tables.items() %}
                <h4>{{ comparison }}</h4>
                <p>Showing top {{ n_top_diff }} differentially active reactions (sorted by adjusted p-value)</p>
                {{ diff_table }}
                {% endfor %}
            </div>
        </section>
        {% endif %}

        {% if has_enrichment %}
        <section id="enrichment" class="card">
            <div class="card-header">4. Pathway Enrichment</div>
            <div class="card-body">
                {% for enrich_type, enrich_table in enrichment_tables.items() %}
                <h4>{{ enrich_type }} Enrichment</h4>
                {{ enrich_table }}
                {% endfor %}
            </div>
        </section>
        {% endif %}

        {% if has_clustering %}
        <section id="clustering" class="card">
            <div class="card-header">5. Clustering Results</div>
            <div class="card-body">
                <div class="stats-grid">
                    <div class="stat-box">
                        <div class="stat-value">{{ n_clusters }}</div>
                        <div class="stat-label">Clusters Identified</div>
                    </div>
                </div>
                <h4>Cluster Sizes</h4>
                {{ cluster_sizes_table }}
            </div>
        </section>
        {% endif %}

        {% if has_figures %}
        <section id="figures" class="card">
            <div class="card-header">6. Visualizations</div>
            <div class="card-body">
                {% for figure in figures %}
                <div class="figure-container">
                    <img src="data:image/png;base64,{{ figure.data }}" alt="{{ figure.caption }}">
                    <p class="figure-caption">{{ figure.caption }}</p>
                </div>
                {% endfor %}
            </div>
        </section>
        {% endif %}

        <footer>
            <p>Generated by CellMetPro v{{ version }} | {{ generation_date }}</p>
            <p>For more information, visit <a href="https://github.com/omicspilot/CellMetPro">CellMetPro Documentation</a></p>
        </footer>
    </div>
</body>
</html>
"""


def generate_summary_stats(scores: pd.DataFrame) -> dict[str, Any]:
    """Generate summary statistics for reaction scores.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores matrix (reactions x cells).

    Returns
    -------
    dict
        Dictionary containing summary statistics.
    """
    stats: dict[str, Any] = {
        "n_reactions": scores.shape[0],
        "n_cells": scores.shape[1],
        "mean_score": float(np.nanmean(scores.values)),
        "median_score": float(np.nanmedian(scores.values)),
        "std_score": float(np.nanstd(scores.values)),
        "min_score": float(np.nanmin(scores.values)),
        "max_score": float(np.nanmax(scores.values)),
        "pct_nonzero": float((scores.values != 0).sum() / scores.size * 100),
        "pct_nan": float(np.isnan(scores.values).sum() / scores.size * 100),
    }

    # Per-reaction stats
    reaction_means = scores.mean(axis=1)
    reaction_vars = scores.var(axis=1)

    stats["top_variable_reactions"] = reaction_vars.nlargest(10).index.tolist()
    stats["top_active_reactions"] = reaction_means.nlargest(10).index.tolist()

    return stats


def _df_to_html(df: pd.DataFrame, max_rows: int = 20) -> str:
    """Convert DataFrame to styled HTML table.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert.
    max_rows : int, optional
        Maximum number of rows to display.

    Returns
    -------
    str
        HTML table string.
    """
    # Make a copy to avoid modifying the original
    df = df.copy()

    if len(df) > max_rows:
        df = df.head(max_rows)

    # Format numeric columns
    for col in df.select_dtypes(include=[np.number]).columns:
        if df[col].dtype == float:
            col_max = df[col].abs().max()
            if pd.notna(col_max) and (col_max < 0.01 or col_max > 1000):
                df[col] = df[col].apply(lambda x: f"{x:.2e}" if pd.notna(x) else "")
            else:
                df[col] = df[col].apply(lambda x: f"{x:.4f}" if pd.notna(x) else "")

    return str(df.to_html(index=True, classes="data-table", border=0))


def _load_figure_as_base64(filepath: Path) -> str | None:
    """Load a figure file and convert to base64.

    Parameters
    ----------
    filepath : Path
        Path to the figure file.

    Returns
    -------
    str or None
        Base64-encoded image data, or None if loading fails.
    """
    try:
        with open(filepath, "rb") as f:
            return base64.b64encode(f.read()).decode("utf-8")
    except Exception:
        return None


class ReportGenerator:
    """Generate HTML reports for CellMetPro analysis results.

    Parameters
    ----------
    results_dir : str or Path
        Directory containing analysis results.

    Attributes
    ----------
    results_dir : Path
        Path to results directory.
    data : dict
        Loaded analysis data.

    Examples
    --------
    >>> generator = ReportGenerator("./results")
    >>> generator.generate("report.html")
    """

    def __init__(self, results_dir: str | Path):
        self.results_dir = Path(results_dir)
        self.data: dict[str, Any] = {}
        self._load_data()

    def _load_data(self) -> None:
        """Load all available data from results directory."""
        # Load scores
        scores_path = self.results_dir / "reaction_scores.csv"
        if scores_path.exists():
            self.data["scores"] = pd.read_csv(scores_path, index_col=0)

        # Load penalties
        penalties_path = self.results_dir / "reaction_penalties.csv"
        if penalties_path.exists():
            self.data["penalties"] = pd.read_csv(penalties_path, index_col=0)

        # Load config
        config_path = self.results_dir / "config.json"
        if config_path.exists():
            with open(config_path) as f:
                self.data["config"] = json.load(f)

        # Load differential results
        self.data["differential"] = {}
        for diff_file in self.results_dir.glob("differential_*.csv"):
            comparison = diff_file.stem.replace("differential_", "")
            self.data["differential"][comparison] = pd.read_csv(diff_file)

        # Load enrichment results
        self.data["enrichment"] = {}
        for enrich_file in self.results_dir.glob("*_enrichment.csv"):
            enrich_type = enrich_file.stem.replace("_enrichment", "")
            self.data["enrichment"][enrich_type] = pd.read_csv(enrich_file)

        # Load clustering results
        clustering_path = self.results_dir / "clustering_results.csv"
        if clustering_path.exists():
            self.data["clusters"] = pd.read_csv(clustering_path)

        # Find figures
        self.data["figures"] = []
        for ext in ["png", "jpg", "jpeg", "svg"]:
            self.data["figures"].extend(self.results_dir.glob(f"*.{ext}"))

    def generate(
        self,
        output_path: str | Path,
        title: str | None = None,
        n_top_diff: int = 20,
        include_figures: bool = True,
    ) -> Path:
        """Generate HTML report.

        Parameters
        ----------
        output_path : str or Path
            Path for output HTML file.
        title : str, optional
            Custom report title.
        n_top_diff : int, optional
            Number of top differential reactions to show.
        include_figures : bool, optional
            Whether to embed figures in the report.

        Returns
        -------
        Path
            Path to generated report.
        """
        output_path = Path(output_path)

        # Try to import jinja2, provide helpful error if not available
        try:
            from jinja2 import Template
        except ImportError:
            # Fall back to simple string formatting
            return self._generate_simple_report(output_path)

        # Prepare template variables
        template_vars: dict[str, Any] = {
            "generation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "version": "0.1.0",
            "n_cells": 0,
            "n_reactions": 0,
            "n_genes": 0,
            "model_name": "Unknown",
            "has_differential": bool(self.data.get("differential")),
            "has_enrichment": bool(self.data.get("enrichment")),
            "has_clustering": "clusters" in self.data,
            "has_figures": bool(self.data.get("figures")) and include_figures,
            "n_top_diff": n_top_diff,
        }

        # Add config info
        if "config" in self.data:
            config = self.data["config"]
            template_vars["n_cells"] = config.get("n_cells", 0)
            template_vars["n_genes"] = config.get("n_genes", 0)
            template_vars["n_reactions"] = config.get("n_reactions", 0)
            template_vars["model_name"] = config.get("model", "Unknown")

            # Config table
            config_df = pd.DataFrame(
                [(k, str(v)) for k, v in config.items()], columns=["Parameter", "Value"]
            )
            template_vars["config_table"] = _df_to_html(config_df)

        # Add scores summary
        if "scores" in self.data:
            scores = self.data["scores"]
            template_vars["n_reactions"] = scores.shape[0]
            template_vars["n_cells"] = scores.shape[1]

            stats = generate_summary_stats(scores)
            summary_df = pd.DataFrame(
                [
                    ("Mean", stats["mean_score"]),
                    ("Median", stats["median_score"]),
                    ("Std Dev", stats["std_score"]),
                    ("Min", stats["min_score"]),
                    ("Max", stats["max_score"]),
                    ("% Non-zero", stats["pct_nonzero"]),
                ],
                columns=["Statistic", "Value"],
            )
            template_vars["scores_summary_table"] = _df_to_html(summary_df)

            # Top variable reactions
            var_df = pd.DataFrame(
                {
                    "Reaction": stats["top_variable_reactions"],
                    "Variance": scores.var(axis=1)[stats["top_variable_reactions"]],
                    "Mean": scores.mean(axis=1)[stats["top_variable_reactions"]],
                }
            )
            template_vars["top_variable_reactions"] = _df_to_html(var_df)

        # Add differential tables
        if self.data.get("differential"):
            diff_tables = {}
            for comparison, df in self.data["differential"].items():
                # Sort by p-value and take top N
                if "padj" in df.columns:
                    df_sorted = df.nsmallest(n_top_diff, "padj")
                elif "padj_bh" in df.columns:
                    df_sorted = df.nsmallest(n_top_diff, "padj_bh")
                elif "pvalue" in df.columns:
                    df_sorted = df.nsmallest(n_top_diff, "pvalue")
                else:
                    df_sorted = df.head(n_top_diff)
                diff_tables[comparison] = _df_to_html(df_sorted)
            template_vars["differential_tables"] = diff_tables

        # Add enrichment tables
        if self.data.get("enrichment"):
            enrich_tables = {}
            for enrich_type, df in self.data["enrichment"].items():
                if "padj" in df.columns:
                    df_sorted = df.nsmallest(20, "padj")
                elif "pvalue" in df.columns:
                    df_sorted = df.nsmallest(20, "pvalue")
                else:
                    df_sorted = df.head(20)
                enrich_tables[enrich_type] = _df_to_html(df_sorted)
            template_vars["enrichment_tables"] = enrich_tables

        # Add clustering info
        if "clusters" in self.data:
            clusters = self.data["clusters"]
            if "cluster" in clusters.columns:
                cluster_counts = clusters["cluster"].value_counts().sort_index()
                template_vars["n_clusters"] = len(cluster_counts)
                sizes_df = pd.DataFrame(
                    {
                        "Cluster": cluster_counts.index,
                        "Size": cluster_counts.values,
                        "Percentage": (cluster_counts.values / len(clusters) * 100),
                    }
                )
                template_vars["cluster_sizes_table"] = _df_to_html(sizes_df)

        # Add figures
        if include_figures and self.data.get("figures"):
            figures = []
            for fig_path in self.data["figures"]:
                fig_data = _load_figure_as_base64(fig_path)
                if fig_data:
                    figures.append(
                        {
                            "data": fig_data,
                            "caption": fig_path.stem.replace("_", " ").title(),
                        }
                    )
            template_vars["figures"] = figures

        # Render template
        template = Template(HTML_TEMPLATE)
        html_content = template.render(**template_vars)

        # Write output
        with open(output_path, "w") as f:
            f.write(html_content)

        return output_path

    def _generate_simple_report(self, output_path: Path) -> Path:
        """Generate a simple HTML report without Jinja2.

        Parameters
        ----------
        output_path : Path
            Output file path.

        Returns
        -------
        Path
            Path to generated report.
        """
        date_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>CellMetPro Report</title>
    <style>
        body {{ font-family: sans-serif; margin: 40px; }}
        h1 {{ color: #2c3e50; }}
        table {{ border-collapse: collapse; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background: #f5f5f5; }}
    </style>
</head>
<body>
    <h1>CellMetPro Analysis Report</h1>
    <p>Generated: {date_str}</p>
"""

        if "scores" in self.data:
            scores = self.data["scores"]
            mean_score = np.nanmean(scores.values) if scores.size > 0 else 0.0
            html += f"""
    <h2>Data Summary</h2>
    <p>Reactions: {scores.shape[0]}</p>
    <p>Cells: {scores.shape[1]}</p>
    <p>Mean score: {mean_score:.4f}</p>
"""

        if "config" in self.data:
            html += "<h2>Configuration</h2><ul>"
            for k, v in self.data["config"].items():
                html += f"<li><strong>{k}:</strong> {v}</li>"
            html += "</ul>"

        html += """
    <p><em>Install jinja2 for enhanced report formatting: pip install jinja2</em></p>
</body>
</html>
"""

        with open(output_path, "w") as f:
            f.write(html)

        return output_path


def generate_html_report(
    results_dir: str | Path,
    output_path: str | Path | None = None,
    **kwargs: Any,
) -> Path:
    """Generate an HTML report from analysis results.

    This is a convenience function that creates a ReportGenerator
    and generates a report in one step.

    Parameters
    ----------
    results_dir : str or Path
        Directory containing analysis results.
    output_path : str or Path, optional
        Path for output HTML file. Defaults to results_dir/report.html.
    **kwargs
        Additional arguments passed to ReportGenerator.generate().

    Returns
    -------
    Path
        Path to generated report.

    Examples
    --------
    >>> generate_html_report("./results")
    PosixPath('./results/report.html')

    >>> generate_html_report("./results", "my_report.html", n_top_diff=50)
    PosixPath('my_report.html')
    """
    results_dir = Path(results_dir)

    if output_path is None:
        output_path = results_dir / "report.html"

    generator = ReportGenerator(results_dir)
    return generator.generate(output_path, **kwargs)
