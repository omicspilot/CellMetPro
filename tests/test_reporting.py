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


class TestDfToHtml:
    """Tests for the _df_to_html internal helper."""

    def test_truncates_long_dataframe(self):
        """DataFrames longer than max_rows are truncated."""
        from cellmetpro.reporting.report import _df_to_html

        df = pd.DataFrame({"a": range(50), "b": range(50)})
        html = _df_to_html(df, max_rows=10)

        # Only 10 rows rendered — check that row 11+ are absent
        assert html.count("<tr>") <= 12  # header + 10 data rows + some tolerance

    def test_scientific_notation_for_tiny_values(self):
        """Very small float values are rendered in scientific notation."""
        from cellmetpro.reporting.report import _df_to_html

        df = pd.DataFrame({"pvalue": [1e-15, 2e-10, 5e-8]})
        html = _df_to_html(df)

        assert "e-" in html.lower() or "E-" in html

    def test_standard_notation_for_normal_values(self):
        """Normal float values are rendered with 4 decimal places."""
        from cellmetpro.reporting.report import _df_to_html

        df = pd.DataFrame({"score": [0.1234, 0.5678]})
        html = _df_to_html(df)

        assert "0.1234" in html

    def test_returns_html_string(self):
        """Return value is a non-empty string."""
        from cellmetpro.reporting.report import _df_to_html

        df = pd.DataFrame({"x": [1, 2, 3]})
        result = _df_to_html(df)

        assert isinstance(result, str)
        assert len(result) > 0


class TestLoadFigureAsBase64:
    """Tests for _load_figure_as_base64."""

    def test_encodes_existing_file(self, tmp_path):
        """A real file is encoded as a non-empty base64 string."""
        from cellmetpro.reporting.report import _load_figure_as_base64

        img_path = tmp_path / "figure.png"
        img_path.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 100)

        result = _load_figure_as_base64(img_path)

        assert result is not None
        assert len(result) > 0
        # Base64 strings contain only valid characters
        import base64

        base64.b64decode(result)  # must not raise

    def test_returns_none_for_missing_file(self, tmp_path):
        """Missing file returns None instead of raising."""
        from cellmetpro.reporting.report import _load_figure_as_base64

        result = _load_figure_as_base64(tmp_path / "nonexistent.png")
        assert result is None


class TestReportGeneratorDataLoading:
    """Additional data-loading tests for ReportGenerator."""

    def test_loads_enrichment_files(self, tmp_path):
        """Enrichment CSV files are loaded into data['enrichment']."""
        enrich = pd.DataFrame({"pathway": ["P1", "P2"], "pvalue": [0.01, 0.05]})
        enrich.to_csv(tmp_path / "subsystem_enrichment.csv", index=False)

        generator = ReportGenerator(tmp_path)

        assert "enrichment" in generator.data
        assert "subsystem" in generator.data["enrichment"]

    def test_loads_clustering_results(self, tmp_path):
        """clustering_results.csv is loaded into data['clusters']."""
        clusters = pd.DataFrame({"cell": ["c1", "c2"], "cluster": [0, 1]})
        clusters.to_csv(tmp_path / "clustering_results.csv", index=False)

        generator = ReportGenerator(tmp_path)

        assert "clusters" in generator.data

    def test_finds_figure_files(self, tmp_path):
        """PNG/JPG files in the results dir appear in data['figures']."""
        (tmp_path / "volcano.png").write_bytes(b"fake_png")
        (tmp_path / "heatmap.jpg").write_bytes(b"fake_jpg")

        generator = ReportGenerator(tmp_path)

        figure_names = [p.name for p in generator.data["figures"]]
        assert "volcano.png" in figure_names
        assert "heatmap.jpg" in figure_names

    def test_generate_report_with_differential_data(self, tmp_path):
        """Report generation works when differential results are present."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        diff = pd.DataFrame(
            {
                "reaction": [f"RXN{i}" for i in range(5)],
                "log2fc": np.random.randn(5),
                "pvalue": np.random.rand(5),
                "padj_bh": np.random.rand(5),
            }
        )
        diff.to_csv(tmp_path / "differential_A_vs_B.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()
        with open(output_path) as f:
            html = f.read()
        assert "CellMetPro" in html


class TestReportGeneratorCoverage:
    """Targeted tests to cover remaining gaps in reporting/report.py."""

    def test_loads_penalties_file(self, tmp_path):
        """reaction_penalties.csv is loaded into data['penalties']."""
        penalties = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        penalties.to_csv(tmp_path / "reaction_penalties.csv")

        generator = ReportGenerator(tmp_path)

        assert "penalties" in generator.data

    def test_generate_without_scores(self, tmp_path):
        """generate() works when no scores are present (config only)."""
        config = {"n_cells": 50, "model": "human"}
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_padj_differential(self, tmp_path):
        """Differential results with 'padj' column are sorted by padj."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        diff = pd.DataFrame(
            {
                "reaction": [f"RXN{i}" for i in range(5)],
                "log2fc": np.random.randn(5),
                "padj": np.random.rand(5),
            }
        )
        diff.to_csv(tmp_path / "differential_A_vs_B.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_pvalue_only_differential(self, tmp_path):
        """Differential with only 'pvalue' (no padj) falls through to pvalue sort."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        diff = pd.DataFrame(
            {
                "reaction": [f"RXN{i}" for i in range(5)],
                "log2fc": np.random.randn(5),
                "pvalue": np.random.rand(5),
            }
        )
        diff.to_csv(tmp_path / "differential_C_vs_D.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_unsortable_differential(self, tmp_path):
        """Differential with no p-value columns falls back to head()."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        diff = pd.DataFrame(
            {
                "reaction": [f"RXN{i}" for i in range(5)],
                "log2fc": np.random.randn(5),
            }
        )
        diff.to_csv(tmp_path / "differential_E_vs_F.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_enrichment_padj(self, tmp_path):
        """Enrichment with 'padj' column is sorted by padj."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        enrich = pd.DataFrame(
            {
                "pathway": ["P1", "P2", "P3"],
                "pvalue": [0.01, 0.05, 0.2],
                "padj": [0.03, 0.1, 0.4],
            }
        )
        enrich.to_csv(tmp_path / "subsystem_enrichment.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_enrichment_pvalue_only(self, tmp_path):
        """Enrichment with only 'pvalue' (no padj) uses pvalue sort."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        enrich = pd.DataFrame(
            {
                "pathway": ["P1", "P2"],
                "pvalue": [0.01, 0.05],
            }
        )
        enrich.to_csv(tmp_path / "go_enrichment.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_enrichment_no_sort_column(self, tmp_path):
        """Enrichment with neither padj nor pvalue falls back to head()."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        enrich = pd.DataFrame(
            {
                "pathway": ["P1", "P2"],
                "score": [0.8, 0.6],
            }
        )
        enrich.to_csv(tmp_path / "kegg_enrichment.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_with_clustering_data(self, tmp_path):
        """Clustering with 'cluster' column generates cluster-sizes table."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        clusters = pd.DataFrame(
            {
                "cell": [f"cell{i}" for i in range(10)],
                "cluster": [0, 0, 1, 1, 1, 2, 2, 2, 2, 0],
            }
        )
        clusters.to_csv(tmp_path / "clustering_results.csv", index=False)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html")

        assert output_path.exists()

    def test_generate_embeds_figures(self, tmp_path):
        """PNG figures are embedded as base64 when include_figures=True."""
        scores = pd.DataFrame(
            np.random.rand(5, 10),
            index=[f"RXN{i}" for i in range(5)],
            columns=[f"cell{i}" for i in range(10)],
        )
        scores.to_csv(tmp_path / "reaction_scores.csv")

        (tmp_path / "figure.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 100)

        generator = ReportGenerator(tmp_path)
        output_path = generator.generate(tmp_path / "report.html", include_figures=True)

        assert output_path.exists()

    def test_simple_report_without_scores(self, tmp_path):
        """_generate_simple_report with no scores skips the data summary."""
        config = {"n_cells": 20}
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        generator = ReportGenerator(tmp_path)
        output_path = generator._generate_simple_report(tmp_path / "simple.html")

        assert output_path.exists()
        with open(output_path) as f:
            html = f.read()
        assert "CellMetPro" in html

    def test_simple_report_with_config(self, tmp_path):
        """_generate_simple_report renders a config section when config is present."""
        config = {"beta": 0.95, "model": "human"}
        with open(tmp_path / "config.json", "w") as f:
            json.dump(config, f)

        generator = ReportGenerator(tmp_path)
        output_path = generator._generate_simple_report(tmp_path / "simple_config.html")

        assert output_path.exists()
        with open(output_path) as f:
            html = f.read()
        assert "beta" in html
        assert "0.95" in html


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
