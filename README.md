# CellMetPro

**Cellular Metabolic Profiler for scRNA-seq data**

Analyze metabolic activity at single-cell resolution using the COMPASS algorithm. Score reactions, identify metabolic heterogeneity, and discover metabolic programs in your scRNA-seq data.

> **Documentation**: For detailed tutorials and API reference, visit [omicspilot.com/projects/cellmetpro](https://omicspilot.com/projects/cellmetpro)

---

## Features

| Feature | Description |
|---------|-------------|
| **COMPASS Algorithm** | Score metabolic reactions from gene expression using genome-scale models |
| **Differential Analysis** | Compare metabolic activity between cell groups (Wilcoxon, t-test, ANOVA, Kruskal-Wallis) |
| **Pathway Enrichment** | GO term and subsystem enrichment analysis |
| **Metabolic Clustering** | PCA, UMAP, t-SNE embeddings with k-means, Leiden, Louvain clustering |
| **Visualization** | Volcano plots, heatmaps, dotplots, embedding plots |
| **CLI & Python API** | Full command-line interface and programmatic access |

---

## Installation

```bash
pip install cellmetpro
```

For development:
```bash
git clone https://github.com/omicspilot/CellMetPro.git
cd CellMetPro
pip install -e ".[dev]"
```

---

## Sample Data

CellMetPro includes sample datasets for testing:

```python
from cellmetpro.data import (
    load_sample_expression,
    load_sample_groups,
    load_sample_reaction_scores,
    create_sample_model,
)

# Load synthetic expression data (50 genes x 100 cells)
expression = load_sample_expression()
print(f"Expression: {expression.shape}")

# Load matching group annotations
groups = load_sample_groups()
print(f"Cell types: {groups['cell_type'].unique()}")

# Load pre-computed reaction scores for quick visualization
scores = load_sample_reaction_scores()

# Create a simple metabolic model for testing
model = create_sample_model()
print(f"Model reactions: {len(model.reactions)}")
```

The sample data includes:
- **Expression matrix**: 50 metabolic genes x 100 cells with 4 cell types (Proliferating, Quiescent, Hypoxic, Oxidative)
- **Group annotations**: Cell type and treatment labels
- **Reaction scores**: Pre-computed scores for differential analysis and visualization
- **Sample model**: Minimal glycolysis model with GPR rules

---

## Quick Start

### Command Line

```bash
# Run COMPASS analysis
cellmetpro run expression.h5ad -m human -o results/

# Differential analysis between groups
cellmetpro differential results/reaction_scores.csv groups.csv --plot

# Cluster cells by metabolic profile
cellmetpro cluster results/reaction_scores.csv --method leiden --embedding umap --plot

# Pathway enrichment
cellmetpro pathway significant_reactions.txt --method subsystem --plot
```

### Python API

```python
import cellmetpro as cmp

# Load data
loader = cmp.DataLoader("expression.h5ad")
adata = loader.load()

# Load metabolic model
model = cmp.load_gem("human")

# Run COMPASS
config = cmp.CompassConfig(beta=0.95, n_processes=4)
scorer = cmp.CompassScorer(model, adata, config)
result = scorer.score()

# Differential analysis
from cellmetpro.analysis import DifferentialAnalysis
da = DifferentialAnalysis(result.reaction_scores, cell_groups)
diff_results = da.compare_groups("control", "treatment")

# Visualize
from cellmetpro.visualization import plot_volcano
plot_volcano(diff_results, save="volcano.png")
```

---

## Supported Data Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| AnnData | `.h5ad` | Scanpy/AnnData objects |
| CSV | `.csv` | Comma-separated values |
| TSV | `.tsv` | Tab-separated values |

## Supported Models

| Model | Organism | Reactions | Genes |
|-------|----------|-----------|-------|
| `human` | Homo sapiens | ~13,000 | ~3,000 |
| `mouse` | Mus musculus | ~13,000 | ~3,000 |
| Custom | Any | User-defined | User-defined |

---

## CLI Commands

| Command | Description |
|---------|-------------|
| `cellmetpro run` | Run COMPASS metabolic analysis |
| `cellmetpro differential` | Compare groups statistically |
| `cellmetpro cluster` | Cluster cells by metabolic profile |
| `cellmetpro pathway` | Pathway enrichment analysis |
| `cellmetpro info` | Show model information |
| `cellmetpro dashboard` | Launch interactive dashboard |

Run `cellmetpro --help` or `cellmetpro <command> --help` for details.

---

## Analysis Modules

### Differential Analysis

```python
from cellmetpro.analysis import DifferentialAnalysis

da = DifferentialAnalysis(reaction_scores, groups)

# Pairwise comparison
results = da.compare_groups("A", "B", method="wilcoxon")

# Multi-group comparison
results = da.compare_multiple_groups(method="kruskal")

# Post-hoc tests
posthoc = da.posthoc_tests("reaction_id", method="dunn")

# Effect size
effect = da.compute_effect_size("A", "B")
```

### Clustering

```python
from cellmetpro.analysis import MetabolicClustering

mc = MetabolicClustering(reaction_scores, n_clusters=5)
mc.compute_pca(n_components=50)
mc.compute_umap()
labels = mc.cluster(method="leiden", resolution=1.0)
markers = mc.get_cluster_markers(n_top=20)
```

### Pathway Enrichment

```python
from cellmetpro.analysis import PathwayAnalyzer, GOEnrichmentAnalyzer

# Subsystem enrichment
pa = PathwayAnalyzer(subsystem_mapping)
results = pa.enrich(significant_reactions, background=all_reactions)

# GO enrichment
go = GOEnrichmentAnalyzer(model)
results = go.enrich_reactions(reactions, namespace="biological_process")
```

---

## Visualization

```python
from cellmetpro.visualization import (
    plot_volcano,
    plot_reaction_heatmap,
    plot_reaction_dotplot,
    plot_embedding,
    plot_enrichment_dotplot,
)

# Volcano plot
plot_volcano(diff_results, log2fc_threshold=0.5, pvalue_threshold=0.05)

# Heatmap with groups
plot_reaction_heatmap(scores, groups, reactions=top_reactions)

# Dotplot
plot_reaction_dotplot(scores, groups, reactions=markers)

# Embedding
plot_embedding(umap_coords, color=cluster_labels)

# Enrichment dotplot
plot_enrichment_dotplot(enrichment_results)
```

---

## Background

**COMPASS** (Characterizing Cell states through metabolic Profiling of the Transcriptome) integrates scRNA-seq data with Genome-Scale Metabolic Models (GEMs) to infer metabolic activity at single-cell resolution.

The algorithm:
1. Maps gene expression to reaction penalties
2. Optimizes flux through each reaction subject to stoichiometric constraints
3. Scores reactions based on consistency with expression data

---

## Citation

If you use CellMetPro in your research, please cite:

```
Wagner et al. (2021). Metabolic modeling of single Th17 cells reveals
regulators of autoimmunity. Cell, 184(16), 4168-4185.
```

---

## License

MIT License - see [LICENSE](LICENSE) for details.

---

## Links

- **Issues**: [github.com/omicspilot/CellMetPro/issues](https://github.com/omicspilot/CellMetPro/issues)
- [![CI](https://github.com/omicspilot/CellMetPro/actions/workflows/ci.yml/badge.svg)](https://github.com/omicspilot/CellMetPro/actions/workflows/ci.yml)
- [![codecov](https://codecov.io/gh/omicspilot/CellMetPro/branch/main/graph/badge.svg)](https://codecov.io/gh/omicspilot/CellMetPro)

- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
