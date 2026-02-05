"""Tests for the reporting module."""

import json

import numpy as np
import pandas as pd

from cellmetpro.reporting import (
    ReportGenerator,
    generate_html_report,
    generate_summary_stats,
)


class TestGenerateSummaryStats:
    """Tests for generate_summary_stats function."""

    def test_basic_stats(self):
        """Test basic summary statistics computation."""
        scores = pd.DataFrame(
            np.random.rand(20, 50),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(50)],
        )

        stats = generate_summary_stats(scores)

        assert "n_reactions" in stats
        assert "n_cells" in stats
        assert "mean_score" in stats
        assert "median_score" in stats
        assert stats["n_reactions"] == 20
        assert stats["n_cells"] == 50

    def test_with_nan_values(self):
        """Test stats with NaN values."""
        scores = pd.DataFrame(
            np.random.rand(10, 20),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(20)],
        )
        scores.iloc[0, 0] = np.nan
        scores.iloc[5, 10] = np.nan

        stats = generate_summary_stats(scores)

        assert stats["pct_nan"] > 0
        assert not np.isnan(stats["mean_score"])

    def test_top_variable_reactions(self):
        """Test identification of top variable reactions."""
        scores = pd.DataFrame(
            np.random.rand(20, 50),
            index=[f"RXN{i}" for i in range(20)],
            columns=[f"cell{i}" for i in range(50)],
        )
        # Make one reaction highly variable
        scores.loc["RXN0"] = np.linspace(0, 10, 50)

        stats = generate_summary_stats(scores)

        assert "top_variable_reactions" in stats
        assert "RXN0" in stats["top_variable_reactions"]


class TestReportGenerator:
    """Tests for ReportGenerator class."""

    def test_init_empty_directory(self, tmp_path):
        """Test initialization with empty directory."""
        generator = ReportGenerator(tmp_path)
        # Only differential, enrichment, and figures dicts/lists are initialized
        assert "scores" not in generator.data
        assert "config" not in generator.data
        assert "clusters" not in generator.data

    def test_load_scores(self, tmp_path):
        """Test loading reaction scores."""
        scores = pd.DataFrame(
            np.random.rand(10, 20),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(20)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        generator = ReportGenerator(tmp_path)

        assert "scores" in generator.data
        assert generator.data["scores"].shape == (10, 20)

    def test_load_config(self, tmp_path):
        """Test loading config file."""
        config = {"n_cells": 100, "model": "human"}
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        generator = ReportGenerator(tmp_path)

        assert "config" in generator.data
        assert generator.data["config"]["n_cells"] == 100

    def test_load_differential(self, tmp_path):
        """Test loading differential results."""
        diff = pd.DataFrame(
            {
                "reaction": ["RXN1", "RXN2"],
                "log2fc": [1.5, -0.8],
                "pvalue": [0.01, 0.05],
            }
        )
        diff.to_csv(tmp_path / "differential_A_vs_B.csv", index=False)

        generator = ReportGenerator(tmp_path)

        assert "differential" in generator.data
        assert "A_vs_B" in generator.data["differential"]

    def test_generate_report(self, tmp_path):
        """Test HTML report generation."""
        # Create test data
        scores = pd.DataFrame(
            np.random.rand(10, 20),
            index=[f"RXN{i}" for i in range(10)],
            columns=[f"cell{i}" for i in range(20)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        config = {"n_cells": 20, "n_genes": 100, "model": "human"}
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()
        with open(output_path) as f:
            html = f.read()
        assert "CellMetPro" in html
        assert "20" in html  # n_cells

    def test_generate_simple_report(self, tmp_path):
        """Test simple report generation (without jinja2)."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        generator = ReportGenerator(tmp_path)
        output_path = generator._generate_simple_report(tmp_path / "simple.html")

        assert output_path.exists()
        with open(output_path) as f:
            html = f.read()
        assert "CellMetPro" in html


class TestGenerateHtmlReport:
    """Tests for generate_html_report convenience function."""

    def test_convenience_function(self, tmp_path):
        """Test the convenience function."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        output_path = generate_html_report(tmp_path)

        assert output_path.exists()
        assert output_path.name == "report.html"

    def test_custom_output_path(self, tmp_path):
        """Test with custom output path."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        output_path = generate_html_report(tmp_path, tmp_path / "custom_report.html")

        assert output_path.name == "custom_report.html"
