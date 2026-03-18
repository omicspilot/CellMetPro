"""Tests for CLI module."""

import importlib.util
import json
from pathlib import Path

import pytest

_has_parquet = (
    importlib.util.find_spec("pyarrow") is not None
    or importlib.util.find_spec("fastparquet") is not None
)


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

        args = parser.parse_args(
            [
                "run",
                "data.csv",
                "-o",
                "output_dir",
                "-m",
                "mouse",
                "--beta",
                "0.9",
                "-j",
                "4",
                "--microcluster",
                "--cells-per-cluster",
                "50",
            ]
        )

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

    def test_run_analysis_creates_output(self, tmp_h5ad, tmp_model_json, tmp_path):
        """Test that run command creates output files."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results"

        result = main(
            [
                "run",
                str(tmp_h5ad),
                "-m",
                str(tmp_model_json),
                "-o",
                str(output_dir),
            ]
        )

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

    def test_run_analysis_with_normalization(self, tmp_h5ad, tmp_model_json, tmp_path):
        """Test run command with normalization enabled."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_norm"

        result = main(
            [
                "run",
                str(tmp_h5ad),
                "-m",
                str(tmp_model_json),
                "-o",
                str(output_dir),
                "--normalize",
                "--target-sum",
                "10000",
            ]
        )

        assert result == 0
        assert output_dir.exists()

    @pytest.mark.skipif(not _has_parquet, reason="pyarrow or fastparquet not installed")
    def test_run_analysis_parquet_output(self, tmp_h5ad, tmp_model_json, tmp_path):
        """Test run command with parquet output format."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_parquet"

        result = main(
            [
                "run",
                str(tmp_h5ad),
                "-m",
                str(tmp_model_json),
                "-o",
                str(output_dir),
                "--output-format",
                "parquet",
            ]
        )

        assert result == 0
        assert (output_dir / "reaction_penalties.parquet").exists()
        assert (output_dir / "reaction_scores.parquet").exists()

    def test_run_analysis_with_microclustering(
        self, tmp_h5ad, tmp_model_json, tmp_path
    ):
        """Test run command with microclustering."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "results_mc"

        result = main(
            [
                "run",
                str(tmp_h5ad),
                "-m",
                str(tmp_model_json),
                "-o",
                str(output_dir),
                "--microcluster",
                "--cells-per-cluster",
                "2",
            ]
        )

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
        args = parser.parse_args(["--verbose"])

        assert args.verbose is True

    def test_verbose_with_error(self, tmp_path, capsys):
        """Test verbose mode shows traceback on error."""
        from cellmetpro.cli import main

        # This should fail and show traceback with --verbose
        result = main(["--verbose", "run", str(tmp_path / "nonexistent.h5ad")])

        assert result == 1


class TestDifferentialCommand:
    """Tests for differential command."""

    def test_differential_parser_arguments(self):
        """Test differential command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "differential",
                "scores.csv",
                "groups.csv",
                "-o",
                "diff_results/",
                "--method",
                "wilcoxon",
                "--fdr-threshold",
                "0.01",
            ]
        )

        assert args.command == "differential"
        assert args.scores == Path("scores.csv")
        assert args.groups == Path("groups.csv")
        assert args.output == Path("diff_results/")
        assert args.method == "wilcoxon"
        assert args.fdr_threshold == 0.01

    def test_differential_pairwise_args(self):
        """Test differential command with pairwise groups."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "differential",
                "scores.csv",
                "groups.csv",
                "--group1",
                "control",
                "--group2",
                "treatment",
                "--plot",
            ]
        )

        assert args.group1 == "control"
        assert args.group2 == "treatment"
        assert args.plot is True

    def test_differential_multigroup_method(self):
        """Test differential command with multigroup methods."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        for method in ["kruskal", "anova"]:
            args = parser.parse_args(
                [
                    "differential",
                    "scores.csv",
                    "groups.csv",
                    "--method",
                    method,
                ]
            )
            assert args.method == method


class TestClusterCommand:
    """Tests for cluster command."""

    def test_cluster_parser_arguments(self):
        """Test cluster command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "cluster",
                "scores.csv",
                "-o",
                "cluster_results/",
                "--n-clusters",
                "5",
                "--method",
                "kmeans",
                "--embedding",
                "umap",
            ]
        )

        assert args.command == "cluster"
        assert args.scores == Path("scores.csv")
        assert args.n_clusters == 5
        assert args.method == "kmeans"
        assert args.embedding == "umap"

    def test_cluster_leiden_args(self):
        """Test cluster command with leiden method."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "cluster",
                "scores.csv",
                "--method",
                "leiden",
                "--resolution",
                "0.5",
                "--plot",
            ]
        )

        assert args.method == "leiden"
        assert args.resolution == 0.5
        assert args.plot is True

    def test_cluster_all_embeddings(self):
        """Test cluster command accepts all embedding types."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        for emb in ["umap", "tsne", "pca"]:
            args = parser.parse_args(
                [
                    "cluster",
                    "scores.csv",
                    "--embedding",
                    emb,
                ]
            )
            assert args.embedding == emb


class TestPathwayCommand:
    """Tests for pathway command."""

    def test_pathway_parser_arguments(self):
        """Test pathway command argument parsing."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "pathway",
                "reactions.txt",
                "-o",
                "pathway_results/",
                "--model",
                "human",
                "--method",
                "subsystem",
            ]
        )

        assert args.command == "pathway"
        assert args.reactions == Path("reactions.txt")
        assert args.model == "human"
        assert args.method == "subsystem"

    def test_pathway_go_args(self):
        """Test pathway command with GO enrichment."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "pathway",
                "reactions.csv",
                "--method",
                "go",
                "--namespace",
                "biological_process",
                "--fdr-threshold",
                "0.1",
                "--plot",
            ]
        )

        assert args.method == "go"
        assert args.namespace == "biological_process"
        assert args.fdr_threshold == 0.1
        assert args.plot is True

    def test_pathway_all_namespaces(self):
        """Test pathway command accepts all GO namespaces."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        namespaces = [
            "biological_process",
            "molecular_function",
            "cellular_component",
            "all",
        ]
        for ns in namespaces:
            args = parser.parse_args(
                [
                    "pathway",
                    "reactions.txt",
                    "--method",
                    "go",
                    "--namespace",
                    ns,
                ]
            )
            assert args.namespace == ns

    def test_pathway_go_annotations_arg(self):
        """Test pathway command with GO annotations file."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "pathway",
                "reactions.csv",
                "--method",
                "go",
                "--go-annotations",
                "annotations.gaf",
            ]
        )

        assert args.method == "go"
        assert args.go_annotations == Path("annotations.gaf")

    def test_pathway_go_annotations_default_none(self):
        """Test pathway command GO annotations default to None."""
        from cellmetpro.cli import create_parser

        parser = create_parser()

        args = parser.parse_args(
            [
                "pathway",
                "reactions.txt",
                "--method",
                "subsystem",
            ]
        )

        assert args.go_annotations is None


class TestDifferentialIntegration:
    """Integration tests for differential command."""

    @pytest.fixture
    def mock_scores_file(self, tmp_path):
        """Create mock reaction scores file."""
        import numpy as np
        import pandas as pd

        np.random.seed(42)
        reactions = [f"R{i}" for i in range(10)]
        cells = [f"cell_{i}" for i in range(20)]
        data = np.random.rand(10, 20)
        # Add signal
        data[0, :10] += 2.0  # R0 higher in first 10 cells
        df = pd.DataFrame(data, index=reactions, columns=cells)
        path = tmp_path / "scores.csv"
        df.to_csv(path)
        return path

    @pytest.fixture
    def mock_groups_file(self, tmp_path):
        """Create mock groups file."""
        import pandas as pd

        cells = [f"cell_{i}" for i in range(20)]
        groups = ["A"] * 10 + ["B"] * 10
        df = pd.DataFrame({"cell_id": cells, "group": groups})
        path = tmp_path / "groups.csv"
        df.to_csv(path, index=False)
        return path

    def test_differential_integration(
        self, mock_scores_file, mock_groups_file, tmp_path
    ):
        """Test differential command integration."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "diff_output"

        exit_code = main(
            [
                "differential",
                str(mock_scores_file),
                str(mock_groups_file),
                "-o",
                str(output_dir),
                "--method",
                "wilcoxon",
            ]
        )

        assert exit_code == 0
        assert output_dir.exists()
        # Check output file exists
        output_files = list(output_dir.glob("differential_*.csv"))
        assert len(output_files) > 0

    def test_differential_with_plot(self, mock_scores_file, mock_groups_file, tmp_path):
        """Test differential command with plot generation."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "diff_output"

        exit_code = main(
            [
                "differential",
                str(mock_scores_file),
                str(mock_groups_file),
                "-o",
                str(output_dir),
                "--plot",
            ]
        )

        assert exit_code == 0
        # Check plot file exists
        assert (output_dir / "volcano_plot.png").exists()

    def test_differential_missing_file(self, tmp_path):
        """Test differential command with missing file."""
        from cellmetpro.cli import main

        exit_code = main(
            [
                "differential",
                str(tmp_path / "nonexistent.csv"),
                str(tmp_path / "groups.csv"),
                "-o",
                str(tmp_path / "output"),
            ]
        )
        assert exit_code == 1


class TestClusterIntegration:
    """Integration tests for cluster command."""

    @pytest.fixture
    def mock_scores_file(self, tmp_path):
        """Create mock reaction scores file."""
        import numpy as np
        import pandas as pd

        np.random.seed(42)
        reactions = [f"R{i}" for i in range(10)]
        cells = [f"cell_{i}" for i in range(20)]
        data = np.random.rand(10, 20)
        df = pd.DataFrame(data, index=reactions, columns=cells)
        path = tmp_path / "scores.csv"
        df.to_csv(path)
        return path

    def test_cluster_integration(self, mock_scores_file, tmp_path):
        """Test cluster command integration."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "cluster_output"

        exit_code = main(
            [
                "cluster",
                str(mock_scores_file),
                "-o",
                str(output_dir),
                "--n-clusters",
                "2",
                "--method",
                "kmeans",
                "--embedding",
                "pca",  # Use PCA for speed
                "--n-pcs",
                "5",
            ]
        )

        assert exit_code == 0
        assert output_dir.exists()
        assert (output_dir / "clustering_results.csv").exists()
        assert (output_dir / "cluster_markers.csv").exists()

    def test_cluster_with_plot(self, mock_scores_file, tmp_path):
        """Test cluster command with plot generation."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "cluster_output"

        exit_code = main(
            [
                "cluster",
                str(mock_scores_file),
                "-o",
                str(output_dir),
                "--n-clusters",
                "2",
                "--embedding",
                "pca",
                "--n-pcs",
                "5",
                "--plot",
            ]
        )

        assert exit_code == 0
        assert (output_dir / "embedding_plot.png").exists()

    def test_cluster_missing_file(self, tmp_path):
        """Test cluster command with missing file."""
        from cellmetpro.cli import main

        exit_code = main(
            [
                "cluster",
                str(tmp_path / "nonexistent.csv"),
                "-o",
                str(tmp_path / "output"),
            ]
        )
        assert exit_code == 1


class TestPathwayIntegration:
    """Integration tests for pathway command."""

    @pytest.fixture
    def mock_reactions_file(self, tmp_path):
        """Create mock reactions file."""
        reactions = ["R1", "R2", "R3", "R4", "R5"]
        path = tmp_path / "reactions.txt"
        with open(path, "w") as f:
            for rxn in reactions:
                f.write(f"{rxn}\n")
        return path

    @pytest.fixture
    def mock_gaf_file(self, tmp_path):
        """Create mock GAF file with GO annotations.

        GAF 2.2 format columns (tab-separated):
        0: DB, 1: DB_Object_ID, 2: DB_Object_Symbol, 3: Qualifier,
        4: GO_ID, 5: DB:Reference, 6: Evidence_Code, 7: With/From,
        8: Aspect, 9: DB_Object_Name, 10: DB_Object_Synonym,
        11: DB_Object_Type, 12: Taxon, 13: Date, 14: Assigned_By
        """
        lines = [
            "!gaf-version: 2.2",
            "!generated-by: test",
            # Proper 15-column GAF format
            "\t".join(
                [
                    "UniProtKB",
                    "A0A001",
                    "gene1",
                    "",
                    "GO:0006096",
                    "PMID:123",
                    "IEA",
                    "",
                    "P",
                    "Gene 1",
                    "G1",
                    "protein",
                    "taxon:9606",
                    "20200101",
                    "Test",
                ]
            ),
            "\t".join(
                [
                    "UniProtKB",
                    "A0A001",
                    "gene1",
                    "",
                    "GO:0005737",
                    "PMID:123",
                    "IEA",
                    "",
                    "C",
                    "Gene 1",
                    "G1",
                    "protein",
                    "taxon:9606",
                    "20200101",
                    "Test",
                ]
            ),
            "\t".join(
                [
                    "UniProtKB",
                    "A0A002",
                    "gene2",
                    "",
                    "GO:0006096",
                    "PMID:456",
                    "IEA",
                    "",
                    "P",
                    "Gene 2",
                    "G2",
                    "protein",
                    "taxon:9606",
                    "20200101",
                    "Test",
                ]
            ),
            "\t".join(
                [
                    "UniProtKB",
                    "A0A003",
                    "gene3",
                    "",
                    "GO:0006099",
                    "PMID:789",
                    "IEA",
                    "",
                    "P",
                    "Gene 3",
                    "G3",
                    "protein",
                    "taxon:9606",
                    "20200101",
                    "Test",
                ]
            ),
        ]
        path = tmp_path / "annotations.gaf"
        with open(path, "w") as f:
            f.write("\n".join(lines))
        return path

    def test_pathway_subsystem_integration(
        self, mock_reactions_file, tmp_model_json, tmp_path
    ):
        """Test pathway command with subsystem enrichment."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "pathway_output"

        exit_code = main(
            [
                "pathway",
                str(mock_reactions_file),
                "-o",
                str(output_dir),
                "--model",
                str(tmp_model_json),
                "--method",
                "subsystem",
            ]
        )

        assert exit_code == 0
        assert output_dir.exists()
        assert (output_dir / "subsystem_enrichment.csv").exists()

    def test_pathway_go_without_annotations(
        self, mock_reactions_file, tmp_model_json, tmp_path
    ):
        """Test pathway GO command falls back without annotations file."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "pathway_output"

        exit_code = main(
            [
                "pathway",
                str(mock_reactions_file),
                "-o",
                str(output_dir),
                "--model",
                str(tmp_model_json),
                "--method",
                "go",
            ]
        )

        # Should succeed but fall back to subsystem
        assert exit_code == 0
        assert output_dir.exists()
        assert (output_dir / "subsystem_enrichment.csv").exists()

    def test_pathway_go_with_annotations(
        self, mock_reactions_file, mock_gaf_file, tmp_model_json, tmp_path
    ):
        """Test pathway GO command with annotations file."""
        from cellmetpro.cli import main

        output_dir = tmp_path / "pathway_output"

        exit_code = main(
            [
                "pathway",
                str(mock_reactions_file),
                "-o",
                str(output_dir),
                "--model",
                str(tmp_model_json),
                "--method",
                "go",
                "--go-annotations",
                str(mock_gaf_file),
            ]
        )

        assert exit_code == 0
        assert output_dir.exists()
        assert (output_dir / "go_enrichment.csv").exists()

    def test_pathway_missing_file(self, tmp_path):
        """Test pathway command with missing file."""
        from cellmetpro.cli import main

        exit_code = main(
            [
                "pathway",
                str(tmp_path / "nonexistent.txt"),
                "-o",
                str(tmp_path / "output"),
            ]
        )
        assert exit_code == 1
