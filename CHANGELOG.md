# Changelog

All notable changes to CellMetPro are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
This project adheres to [Semantic Versioning](https://semver.org/).

---

## [0.1.0] – 2026-03-18

First public release of CellMetPro — a Python package for analyzing and
visualizing cellular metabolic profiles from single-cell RNA-seq data using
flux balance analysis (FBA).

### Core engine

- **COMPASS algorithm** — cell-level reaction penalty and score computation
  from gene expression via GPR rules and FBA
- **FBA utilities** — cobra-based solver wrappers with multi-process support
- **Microclustering** — aggregate cells into metacells before scoring to reduce
  compute cost on large datasets
- **Batch correction** — mean-centering and ComBat integration across
  experimental batches; highly variable reaction (HVR) selection
- **Data loading** — CSV, TSV, h5ad (AnnData), MTX, and Seurat RDS formats

### Analysis modules

- **Differential analysis** — pairwise Wilcoxon / t-test / Mann-Whitney U
  with BH and Bonferroni correction; Cohen's d effect sizes; multi-group
  Kruskal-Wallis / ANOVA; Dunn, Tukey, and Conover post-hoc tests;
  all-pairwise comparisons
- **Clustering** — Leiden (via igraph/leidenalg), K-means, and hierarchical
  clustering on metabolic profiles; UMAP, t-SNE, and PCA embeddings; cluster
  marker identification
- **Pathway analysis** — reaction-to-pathway mapping; enrichment scoring;
  GO term enrichment with BH correction
- **Trajectory analysis** — pseudotime ordering; dynamic reaction detection;
  metabolic velocity computation
- **Perturbation analysis** — in-silico gene/reaction knock-out and
  over-expression simulations

### Visualization

- UMAP / t-SNE / PCA embeddings with continuous and categorical coloring;
  automatic outside-legend placement for >10 categories
- Reaction and pathway heatmaps with hierarchical clustering and z-score scaling
- Violin and box plots per reaction per group
- Dot plots for reaction activity and GO enrichment
- Volcano and MA plots for differential results
- PCA variance-explained bar chart
- Feature-on-embedding grid plots
- Interactive Plotly scatter plots and feature expression panels
- Streamlit dashboard (`cellmetpro dashboard`)

### CLI

- `cellmetpro run` — end-to-end scoring pipeline
- `cellmetpro differential` — group comparison with auto barcode normalization
  (handles 10x GEM well suffixes)
- `cellmetpro cluster` — Leiden / K-means clustering with UMAP embedding
- `cellmetpro pathway` — pathway enrichment analysis
- `cellmetpro trajectory` — pseudotime and dynamic reaction analysis
- `cellmetpro batch-correct` — batch effect correction
- `cellmetpro perturbation` — in-silico perturbation
- `cellmetpro report` — HTML summary report generation
- `cellmetpro dashboard` — interactive Streamlit dashboard
- `cellmetpro sample-data` — download built-in sample datasets
- Shell auto-completion via `argcomplete`
- `--plot` flag auto-opens generated figures with the OS default viewer

### Infrastructure

- CI pipeline (GitHub Actions) — test matrix: Ubuntu + macOS × Python
  3.10 / 3.11 / 3.12; black formatting check; ruff lint
- Automated PyPI publish on GitHub Release (OIDC trusted publishing)
- Codecov coverage reporting
- 474 unit and integration tests; 71 % line coverage

[0.1.0]: https://github.com/ndiayeoumar/CellMetPro/releases/tag/v0.1.0
