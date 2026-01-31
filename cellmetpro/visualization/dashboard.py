"""Streamlit dashboard for interactive exploration of CellMetPro results.

Launch with:
    cellmetpro dashboard results/

Or directly:
    streamlit run cellmetpro/visualization/dashboard.py -- results/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# Page config must be first Streamlit command
st.set_page_config(
    page_title="CellMetPro Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


def load_results(results_dir: Path) -> dict:
    """Load CellMetPro results from directory.

    Parameters
    ----------
    results_dir : Path
        Directory containing CellMetPro output files.

    Returns
    -------
    dict
        Dictionary with loaded data.
    """
    data: dict = {}

    # Load reaction scores
    scores_file = results_dir / "reaction_scores.csv"
    if scores_file.exists():
        data["scores"] = pd.read_csv(scores_file, index_col=0)

    # Load penalties
    penalties_file = results_dir / "reaction_penalties.csv"
    if penalties_file.exists():
        data["penalties"] = pd.read_csv(penalties_file, index_col=0)

    # Load config
    config_file = results_dir / "config.json"
    if config_file.exists():
        with open(config_file) as f:
            data["config"] = json.load(f)

    # Load differential results if present
    for diff_file in results_dir.glob("differential_*.csv"):
        if "differential" not in data:
            data["differential"] = {}
        name = diff_file.stem.replace("differential_", "")
        data["differential"][name] = pd.read_csv(diff_file)

    # Load clustering results if present
    cluster_file = results_dir / "clustering_results.csv"
    if cluster_file.exists():
        data["clusters"] = pd.read_csv(cluster_file)

    # Load enrichment results if present
    for enrich_file in results_dir.glob("*_enrichment.csv"):
        if "enrichment" not in data:
            data["enrichment"] = {}
        name = enrich_file.stem.replace("_enrichment", "")
        data["enrichment"][name] = pd.read_csv(enrich_file)

    return data


def compute_embedding(
    scores: pd.DataFrame, method: str = "umap", n_pcs: int = 50
) -> np.ndarray:
    """Compute dimensionality reduction embedding.

    Parameters
    ----------
    scores : pd.DataFrame
        Reaction scores (reactions x cells).
    method : str
        Embedding method: "umap", "tsne", or "pca".
    n_pcs : int
        Number of PCA components for preprocessing.

    Returns
    -------
    np.ndarray
        2D embedding coordinates.
    """
    # Transpose to cells x reactions
    X = scores.T.values

    # Handle missing values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # PCA preprocessing
    n_components = min(n_pcs, X.shape[0] - 1, X.shape[1])
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)

    if method == "pca":
        return np.asarray(X_pca[:, :2])
    elif method == "tsne":
        perplexity = min(30, X_pca.shape[0] - 1)
        tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
        return np.asarray(tsne.fit_transform(X_pca))
    else:  # umap
        try:
            import umap

            reducer = umap.UMAP(n_components=2, random_state=42)
            return np.asarray(reducer.fit_transform(X_pca))
        except ImportError:
            st.warning("UMAP not installed. Falling back to t-SNE.")
            perplexity = min(30, X_pca.shape[0] - 1)
            tsne = TSNE(n_components=2, perplexity=perplexity, random_state=42)
            return np.asarray(tsne.fit_transform(X_pca))


def page_overview(data: dict) -> None:
    """Render overview page."""
    st.header("Overview")

    if "config" in data:
        config = data["config"]
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Cells", config.get("n_cells", "N/A"))
        with col2:
            st.metric("Genes", config.get("n_genes", "N/A"))
        with col3:
            st.metric("Reactions", config.get("n_reactions", "N/A"))
        with col4:
            st.metric("Model", config.get("model", "N/A"))

        st.subheader("Analysis Parameters")
        st.json(config)

    if "scores" in data:
        scores = data["scores"]
        st.subheader("Reaction Scores Summary")

        col1, col2 = st.columns(2)

        with col1:
            st.write(
                f"**Dimensions:** {scores.shape[0]} reactions x {scores.shape[1]} cells"
            )

            # Score distribution
            flat_scores = scores.values.flatten()
            flat_scores = flat_scores[~np.isnan(flat_scores)]

            fig = px.histogram(
                x=flat_scores,
                nbins=50,
                title="Score Distribution",
                labels={"x": "Score", "y": "Count"},
            )
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            # Top variable reactions
            variance = scores.var(axis=1).sort_values(ascending=False)
            top_var = variance.head(20)

            fig = px.bar(
                x=top_var.values,
                y=top_var.index,
                orientation="h",
                title="Top 20 Variable Reactions",
                labels={"x": "Variance", "y": "Reaction"},
            )
            fig.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(fig, use_container_width=True)


def page_embedding(data: dict) -> None:
    """Render embedding visualization page."""
    st.header("Cell Embedding")

    if "scores" not in data:
        st.warning("No reaction scores loaded.")
        return

    scores = data["scores"]

    # Sidebar controls
    with st.sidebar:
        st.subheader("Embedding Settings")
        method = st.selectbox("Method", ["umap", "tsne", "pca"])
        n_pcs = st.slider("PCA components", 10, 100, 50)

        color_by = st.selectbox(
            "Color by",
            ["None", "Cluster", "Reaction Score"],
        )

        reaction = None
        if color_by == "Reaction Score":
            reaction = st.selectbox("Reaction", scores.index.tolist())

    # Compute embedding
    with st.spinner(f"Computing {method.upper()} embedding..."):
        embedding = compute_embedding(scores, method=method, n_pcs=n_pcs)

    # Prepare plot data
    plot_df = pd.DataFrame(
        {
            "x": embedding[:, 0],
            "y": embedding[:, 1],
            "cell": scores.columns,
        }
    )

    # Add color variable
    if color_by == "Cluster" and "clusters" in data:
        clusters = data["clusters"]
        if "cluster" in clusters.columns:
            plot_df["color"] = clusters["cluster"].astype(str).values
        else:
            plot_df["color"] = "0"
    elif color_by == "Reaction Score" and reaction is not None:
        plot_df["color"] = scores.loc[reaction].values
    else:
        plot_df["color"] = "All cells"

    # Create plot
    if color_by == "Reaction Score" and reaction is not None:
        fig = px.scatter(
            plot_df,
            x="x",
            y="y",
            color="color",
            hover_data=["cell"],
            title=f"{method.upper()} Embedding",
            labels={
                "x": f"{method.upper()} 1",
                "y": f"{method.upper()} 2",
                "color": reaction,
            },
            color_continuous_scale="Viridis",
        )
    else:
        fig = px.scatter(
            plot_df,
            x="x",
            y="y",
            color="color",
            hover_data=["cell"],
            title=f"{method.upper()} Embedding",
            labels={
                "x": f"{method.upper()} 1",
                "y": f"{method.upper()} 2",
                "color": "Group",
            },
        )

    fig.update_traces(marker={"size": 8, "opacity": 0.7})
    fig.update_layout(height=600)
    st.plotly_chart(fig, use_container_width=True)

    # Store embedding in session state for other pages
    st.session_state["embedding"] = embedding


def page_differential(data: dict) -> None:
    """Render differential analysis page."""
    st.header("Differential Analysis")

    if "scores" not in data:
        st.warning("No reaction scores loaded.")
        return

    # Check for pre-computed results
    if "differential" in data:
        st.subheader("Pre-computed Results")

        result_name = st.selectbox(
            "Select comparison", list(data["differential"].keys())
        )
        result = data["differential"][result_name]

        # Display volcano plot
        if "log2fc" in result.columns and "pvalue" in result.columns:
            col1, col2 = st.columns([3, 1])

            with col2:
                log2fc_thresh = st.slider("Log2FC threshold", 0.0, 2.0, 0.5)
                pval_thresh = st.slider("P-value threshold", 0.001, 0.1, 0.05)

            with col1:
                # Prepare volcano data
                volcano_df = result.copy()
                volcano_df["-log10(pvalue)"] = -np.log10(
                    volcano_df["pvalue"].clip(lower=1e-300)
                )

                # Classify points
                conditions = [
                    (volcano_df["log2fc"] >= log2fc_thresh)
                    & (volcano_df["pvalue"] <= pval_thresh),
                    (volcano_df["log2fc"] <= -log2fc_thresh)
                    & (volcano_df["pvalue"] <= pval_thresh),
                ]
                choices = ["Up", "Down"]
                volcano_df["Regulation"] = np.select(conditions, choices, default="NS")

                fig = px.scatter(
                    volcano_df,
                    x="log2fc",
                    y="-log10(pvalue)",
                    color="Regulation",
                    color_discrete_map={"Up": "red", "Down": "blue", "NS": "gray"},
                    hover_data=(
                        ["reaction"] if "reaction" in volcano_df.columns else None
                    ),
                    title="Volcano Plot",
                )
                fig.add_vline(x=log2fc_thresh, line_dash="dash", line_color="gray")
                fig.add_vline(x=-log2fc_thresh, line_dash="dash", line_color="gray")
                fig.add_hline(
                    y=-np.log10(pval_thresh), line_dash="dash", line_color="gray"
                )
                fig.update_layout(height=500)
                st.plotly_chart(fig, use_container_width=True)

        # Show table
        st.subheader("Results Table")
        st.dataframe(result, use_container_width=True)

        # Download button
        csv = result.to_csv(index=False)
        st.download_button(
            "Download CSV",
            csv,
            file_name=f"differential_{result_name}.csv",
            mime="text/csv",
        )

    else:
        st.info(
            "No pre-computed differential results found. "
            "Run `cellmetpro differential` to generate results."
        )

        # Simple on-the-fly comparison
        st.subheader("Quick Comparison (requires group annotations)")
        st.write(
            "To perform differential analysis, "
            "provide a groups file with cell annotations."
        )


def page_heatmap(data: dict) -> None:
    """Render heatmap visualization page."""
    st.header("Heatmap")

    if "scores" not in data:
        st.warning("No reaction scores loaded.")
        return

    scores = data["scores"]

    with st.sidebar:
        st.subheader("Heatmap Settings")

        # Reaction selection
        selection_method = st.radio("Reaction selection", ["Top variable", "Custom"])

        if selection_method == "Top variable":
            n_reactions = st.slider("Number of reactions", 10, 100, 30)
            variance = scores.var(axis=1).sort_values(ascending=False)
            selected_reactions = variance.head(n_reactions).index.tolist()
        else:
            selected_reactions = st.multiselect(
                "Select reactions",
                scores.index.tolist(),
                default=scores.index[:20].tolist(),
            )

        # Normalization
        normalize = st.checkbox("Z-score normalize", value=True)

        # Color scale
        color_scale = st.selectbox(
            "Color scale", ["RdBu_r", "Viridis", "Plasma", "Blues", "Reds"]
        )

    if not selected_reactions:
        st.warning("Please select at least one reaction.")
        return

    # Prepare data
    heatmap_data = scores.loc[selected_reactions]

    if normalize:
        # Z-score per reaction
        heatmap_data = (heatmap_data.T - heatmap_data.mean(axis=1)).T
        heatmap_data = (heatmap_data.T / heatmap_data.std(axis=1).replace(0, 1)).T

    # Create heatmap
    fig = px.imshow(
        heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        color_continuous_scale=color_scale,
        aspect="auto",
        title="Reaction Score Heatmap",
        labels={"color": "Z-score" if normalize else "Score"},
    )
    fig.update_layout(height=max(400, len(selected_reactions) * 15))
    st.plotly_chart(fig, use_container_width=True)


def page_enrichment(data: dict) -> None:
    """Render pathway enrichment page."""
    st.header("Pathway Enrichment")

    if "enrichment" not in data:
        st.info(
            "No enrichment results found. "
            "Run `cellmetpro pathway` to generate enrichment results."
        )
        return

    # Select enrichment result
    result_name = st.selectbox("Select enrichment", list(data["enrichment"].keys()))
    result = data["enrichment"][result_name]

    with st.sidebar:
        st.subheader("Filter Settings")
        if "padj" in result.columns:
            fdr_thresh = st.slider("FDR threshold", 0.001, 0.5, 0.05)
            filtered = result[result["padj"] <= fdr_thresh]
        else:
            filtered = result

        n_top = st.slider("Top pathways to show", 5, 50, 20)

    if len(filtered) == 0:
        st.warning("No significant pathways at this threshold.")
        return

    # Sort and limit
    if "padj" in filtered.columns:
        filtered = filtered.sort_values("padj").head(n_top)
    elif "pvalue" in filtered.columns:
        filtered = filtered.sort_values("pvalue").head(n_top)

    # Determine pathway name column
    name_col = None
    for col in ["pathway", "go_name", "term", "subsystem"]:
        if col in filtered.columns:
            name_col = col
            break

    if name_col is None:
        st.error("Could not find pathway name column.")
        return

    # Create dotplot
    if "padj" in filtered.columns and "odds_ratio" in filtered.columns:
        fig = px.scatter(
            filtered,
            x="odds_ratio",
            y=name_col,
            size="n_overlap" if "n_overlap" in filtered.columns else None,
            color=-np.log10(filtered["padj"]),
            color_continuous_scale="Reds",
            title="Enrichment Dotplot",
            labels={
                "odds_ratio": "Odds Ratio",
                name_col: "Pathway",
                "color": "-log10(FDR)",
            },
        )
        fig.update_layout(
            height=max(400, len(filtered) * 25),
            yaxis={"categoryorder": "total ascending"},
        )
        st.plotly_chart(fig, use_container_width=True)
    else:
        # Simple bar plot
        if "padj" in filtered.columns:
            filtered = filtered.copy()
            filtered["-log10(FDR)"] = -np.log10(filtered["padj"])
            fig = px.bar(
                filtered,
                x="-log10(FDR)",
                y=name_col,
                orientation="h",
                title="Pathway Enrichment",
            )
        else:
            fig = px.bar(
                filtered,
                x="pvalue" if "pvalue" in filtered.columns else filtered.columns[1],
                y=name_col,
                orientation="h",
                title="Pathway Enrichment",
            )
        fig.update_layout(
            height=max(400, len(filtered) * 25),
            yaxis={"categoryorder": "total ascending"},
        )
        st.plotly_chart(fig, use_container_width=True)

    # Show table
    st.subheader("Enrichment Results")
    st.dataframe(filtered, use_container_width=True)

    # Download
    csv = filtered.to_csv(index=False)
    st.download_button(
        "Download CSV",
        csv,
        file_name=f"enrichment_{result_name}_filtered.csv",
        mime="text/csv",
    )


def page_reaction_viewer(data: dict) -> None:
    """Render individual reaction viewer page."""
    st.header("Reaction Viewer")

    if "scores" not in data:
        st.warning("No reaction scores loaded.")
        return

    scores = data["scores"]

    # Reaction selection
    reaction = st.selectbox("Select reaction", scores.index.tolist())

    col1, col2 = st.columns(2)

    with col1:
        # Distribution plot
        values = scores.loc[reaction].values
        fig = px.histogram(
            x=values,
            nbins=30,
            title=f"Score Distribution: {reaction}",
            labels={"x": "Score", "y": "Count"},
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        # Stats
        st.subheader("Statistics")
        stats_df = pd.DataFrame(
            {
                "Statistic": ["Mean", "Median", "Std", "Min", "Max", "Non-zero %"],
                "Value": [
                    f"{np.nanmean(values):.4f}",
                    f"{np.nanmedian(values):.4f}",
                    f"{np.nanstd(values):.4f}",
                    f"{np.nanmin(values):.4f}",
                    f"{np.nanmax(values):.4f}",
                    f"{(values != 0).sum() / len(values) * 100:.1f}%",
                ],
            }
        )
        st.dataframe(stats_df, use_container_width=True, hide_index=True)

    # Show on embedding if available
    if "embedding" in st.session_state:
        st.subheader("Score on Embedding")
        embedding = st.session_state["embedding"]

        plot_df = pd.DataFrame(
            {
                "x": embedding[:, 0],
                "y": embedding[:, 1],
                "score": values,
                "cell": scores.columns,
            }
        )

        fig = px.scatter(
            plot_df,
            x="x",
            y="y",
            color="score",
            hover_data=["cell", "score"],
            title=f"{reaction} on Embedding",
            color_continuous_scale="Viridis",
        )
        fig.update_traces(marker={"size": 8, "opacity": 0.7})
        fig.update_layout(height=500)
        st.plotly_chart(fig, use_container_width=True)


def main() -> None:
    """Main dashboard entry point."""
    # Get results directory from command line
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[-1])
    else:
        results_dir = None

    # Sidebar
    st.sidebar.title("🧬 CellMetPro")
    st.sidebar.markdown("---")

    # Results directory input
    if results_dir is None or not results_dir.exists():
        st.sidebar.warning("No results directory specified.")
        results_path = st.sidebar.text_input("Results directory path", "results/")
        results_dir = Path(results_path)

    if not results_dir.exists():
        st.error(f"Results directory not found: {results_dir}")
        st.info(
            "Run `cellmetpro run` first to generate results, or specify a valid path."
        )
        return

    # Load data
    with st.spinner("Loading results..."):
        data = load_results(results_dir)

    if not data:
        st.error("No CellMetPro results found in the specified directory.")
        return

    st.sidebar.success(f"Loaded: {results_dir}")

    # Navigation
    st.sidebar.markdown("---")
    page = st.sidebar.radio(
        "Navigation",
        [
            "Overview",
            "Embedding",
            "Differential",
            "Heatmap",
            "Enrichment",
            "Reaction Viewer",
        ],
    )

    # Render selected page
    if page == "Overview":
        page_overview(data)
    elif page == "Embedding":
        page_embedding(data)
    elif page == "Differential":
        page_differential(data)
    elif page == "Heatmap":
        page_heatmap(data)
    elif page == "Enrichment":
        page_enrichment(data)
    elif page == "Reaction Viewer":
        page_reaction_viewer(data)

    # Footer
    st.sidebar.markdown("---")
    st.sidebar.markdown("""
        **CellMetPro Dashboard**

        [Documentation](https://omicspilot.com/projects/cellmetpro)
        """)


if __name__ == "__main__":
    main()
