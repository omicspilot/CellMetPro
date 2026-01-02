"""Tests for CLI module."""

import json
import pytest
from pathlib import Path


class TestCLI:
    """Tests for command-line interface."""

    def test_parser_creation(self):
        """Test CLI parser creation."""
        from cellmetpro.cli import create_parser

        parser = create_parser()
        assert parser is not None
        assert parser.prog == "cellmetpro"

    def test_version_command(self, capsys):
        """Test --version command."""
        from cellmetpro.cli import main

        result = main(["--version"])

        assert result == 0
        captured = capsys.readouterr()
        assert "cellmetpro" in captured.out

    def test_help_command(self, capsys):
        """Test help output when no command given."""
        from cellmetpro.cli import main

        result = main([])

        assert result == 0

    def test_run_parser_arguments(self):
        """Test run command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        # Test default values
        args = parser.parse_args(["run", "test.h5ad"])
        assert args.input == Path("test.h5ad")
        assert args.output == Path("results")
        assert args.model == "human"
        assert args.beta == 0.95
        assert args.n_processes == 1

    def test_run_parser_custom_arguments(self):
        """Test run command with custom arguments."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args([
            "run", "data.csv",
            "-o", "output_dir",
            "-m", "mouse",
            "--beta", "0.9",
            "-j", "4",
            "--microcluster",
            "--cells-per-cluster", "50",
        ])

        assert args.input == Path("data.csv")
        assert args.output == Path("output_dir")
        assert args.model == "mouse"
        assert args.beta == 0.9
        assert args.n_processes == 4
        assert args.microcluster is True
        assert args.cells_per_cluster == 50

    def test_dashboard_parser_arguments(self):
        """Test dashboard command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(["dashboard", "results/", "-p", "8888"])

        assert args.results == Path("results/")
        assert args.port == 8888

    def test_info_parser_arguments(self):
        """Test info command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(["info", "human"])

        assert args.model == "human"

    def test_run_nonexistent_input(self, tmp_path, capsys):
        """Test run command with non-existent input file."""
        from cellmetpro.cli import main

        result = main(["run", str(tmp_path / "nonexistent.h5ad")])

        assert result == 1

    def test_run_analysis_creates_output(
        self, tmp_h5ad, tmp_model_json, tmp_path
    ):
        """Test that run command creates output files."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results"

        result = main([
            "run", str(tmp_h5ad),
            "-m", str(tmp_model_json),
            "-o", str(output_dir),
        ])

        assert result == 0
        assert output_dir.exists()
        assert (output_dir / "reaction_penalties.csv").exists()
        assert (output_dir / "reaction_scores.csv").exists()
        assert (output_dir / "config.json").exists()

        # Check config file content
        with open(output_dir / "config.json") as f:
            config = json.load(f)
        assert "beta" in config
        assert "n_cells" in config

    def test_run_analysis_with_normalization(
        self, tmp_h5ad, tmp_model_json, tmp_path
    ):
        """Test run command with normalization enabled."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_norm"

        result = main([
            "run", str(tmp_h5ad),
            "-m", str(tmp_model_json),
            "-o", str(output_dir),
            "--normalize",
            "--target-sum", "10000",
        ])

        assert result == 0
        assert output_dir.exists()

    def test_run_analysis_parquet_output(
        self, tmp_h5ad, tmp_model_json, tmp_path
    ):
        """Test run command with parquet output format."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_parquet"

        result = main([
            "run", str(tmp_h5ad),
            "-m", str(tmp_model_json),
            "-o", str(output_dir),
            "--output-format", "parquet",
        ])

        assert result == 0
        assert (output_dir / "reaction_penalties.parquet").exists()
        assert (output_dir / "reaction_scores.parquet").exists()

    def test_run_analysis_with_microclustering(
        self, tmp_h5ad, tmp_model_json, tmp_path
    ):
        """Test run command with microclustering."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_mc"

        result = main([
            "run", str(tmp_h5ad),
            "-m", str(tmp_model_json),
            "-o", str(output_dir),
            "--microcluster",
            "--cells-per-cluster", "2",
        ])

        assert result == 0
        assert output_dir.exists()

        # Check config records microclustering
        with open(output_dir / "config.json") as f:
            config = json.load(f)
        assert config["microcluster"] is True

    def test_info_with_model_file(self, tmp_model_json, capsys):
        """Test info command with model file."""
        from cellmetpro.cli import main

        result = main(["info", str(tmp_model_json)])

        assert result == 0
        captured = capsys.readouterr()
        assert "Reactions:" in captured.out
        assert "Metabolites:" in captured.out
        assert "Genes:" in captured.out

    def test_info_nonexistent_model(self, capsys):
        """Test info command with non-existent model."""
        from cellmetpro.cli import main

        result = main(["info", "nonexistent_model_xyz"])

        assert result == 1


class TestCLIVerbose:
    """Tests for CLI verbose mode."""

    def test_verbose_flag(self):
        """Test verbose flag is parsed."""
        from cellmetpro.cli import create_parser

        parser = create_parser()
        args = parser.parse_args(["-v", "--version"])

        assert args.verbose is True

    def test_verbose_with_error(self, tmp_path, capsys):
        """Test verbose mode shows traceback on error."""
        from cellmetpro.cli import main

        # This should fail and show traceback with -v
        result = main(["-v", "run", str(tmp_path / "nonexistent.h5ad")])

        assert result == 1
