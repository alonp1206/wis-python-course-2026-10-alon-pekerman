import contextlib
import io
import math
import os

import matplotlib

matplotlib.use("Agg")

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


def validate_count_matrix(counts_df):
    if counts_df.empty or counts_df.shape[1] == 0:
        raise ValueError("The count table must contain genes and samples.")
    if counts_df.columns.has_duplicates:
        raise ValueError("Sample identifiers in the count table must be unique.")
    if counts_df.index.isna().any() or (counts_df.index.astype(str).str.strip() == "").any():
        raise ValueError("Gene identifiers cannot be empty.")
    if counts_df.columns.isna().any() or (counts_df.columns.astype(str).str.strip() == "").any():
        raise ValueError("Sample identifiers cannot be empty.")

    try:
        numeric_counts = counts_df.apply(pd.to_numeric, errors="raise")
        values = numeric_counts.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Counts must be numeric raw read counts.") from exc
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Counts must be finite and non-negative.")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("DESeq2 requires non-negative integer raw counts.")
    numeric_counts = numeric_counts.astype(np.int64)
    if numeric_counts.index.has_duplicates:
        numeric_counts = numeric_counts.groupby(level=0, sort=False).sum()
    return numeric_counts


def validate_input_tables(counts_df, meta_df):
    counts_df = validate_count_matrix(counts_df)
    if meta_df.index.has_duplicates:
        raise ValueError("Sample identifiers in the metadata must be unique.")
    if meta_df.index.isna().any() or (meta_df.index.astype(str).str.strip() == "").any():
        raise ValueError("Metadata sample identifiers cannot be empty.")
    if not meta_df.columns.is_unique:
        raise ValueError("Metadata column names must be unique.")
    if meta_df.empty or meta_df.shape[1] == 0:
        raise ValueError("The metadata table must contain sample annotations.")
    return counts_df


def _log_cpm(counts_df):
    counts_df = validate_count_matrix(counts_df)
    library_sizes = counts_df.sum(axis=0)
    if (library_sizes <= 0).any():
        empty_samples = library_sizes.index[library_sizes <= 0].tolist()
        raise ValueError(f"Samples have zero library size: {', '.join(map(str, empty_samples))}")
    cpm = counts_df.div(library_sizes, axis=1) * 1_000_000
    return np.log2(cpm + 1)


def process_counts_and_filter(counts_df, min_count=10, min_samples=2):
    counts_df = validate_count_matrix(counts_df)
    if min_count < 0:
        raise ValueError("Minimum count must be zero or greater.")
    if min_samples < 1 or min_samples > counts_df.shape[1]:
        raise ValueError(f"Minimum samples must be between 1 and {counts_df.shape[1]}.")
    keep = counts_df.ge(min_count).sum(axis=1) >= min_samples
    return counts_df.loc[keep]


def filter_preview_curve(counts_df, min_samples=2, selected_min_count=10):
    return filter_preview_data(counts_df, min_samples, selected_min_count)["curve"]


def filter_preview_data(counts_df, min_samples=2, selected_min_count=10, bins=50):
    counts_df = validate_count_matrix(counts_df)
    if min_samples < 1 or min_samples > counts_df.shape[1]:
        raise ValueError(f"Minimum samples must be between 1 and {counts_df.shape[1]}.")
    if selected_min_count < 0:
        raise ValueError("Minimum count must be zero or greater.")
    if bins < 1:
        raise ValueError("Histogram bin count must be positive.")

    values = counts_df.to_numpy(dtype=np.int64)
    nth_largest = np.partition(values, values.shape[1] - min_samples, axis=1)[
        :, values.shape[1] - min_samples
    ]
    max_count = max(int(nth_largest.max(initial=0)), selected_min_count, 1)
    upper = max(100, max_count)
    bin_edges = np.expm1(np.linspace(0, np.log10(upper + 1), num=bins + 1))
    bin_edges[0] = 0
    bin_edges[-1] = np.nextafter(float(upper), np.inf)
    histogram_counts, _ = np.histogram(nth_largest, bins=bin_edges)
    total_genes = len(nth_largest)
    thresholds = np.unique(
        np.append(np.linspace(0, upper, num=101, dtype=int), selected_min_count)
    )
    retained = [int(np.count_nonzero(nth_largest >= threshold)) for threshold in thresholds]
    return {
        "histogram": [
            {
                "x": float(np.sqrt((left + 1) * (right + 1)) - 1),
                "y": float(count * 100 / total_genes),
                "bin_start": float(left),
                "bin_end": float(right),
            }
            for left, right, count in zip(
                bin_edges[:-1], bin_edges[1:], histogram_counts
            )
        ],
        "curve": [
            {
                "min_count": int(threshold),
                "genes_retained_percent": float(count * 100 / total_genes),
            }
            for threshold, count in zip(thresholds, retained)
        ],
        "max_count": upper,
        "selected_percent": float(
            np.count_nonzero(nth_largest >= selected_min_count) * 100 / total_genes
        ),
    }


def _top_variable_genes(expression_df, max_genes=1000):
    variances = expression_df.var(axis=1, ddof=0)
    variances = variances[variances > 0]
    if variances.empty:
        raise ValueError("No variable genes are available for this analysis.")
    return variances.nlargest(min(max_genes, len(variances))).index


def run_pca_analysis(
    counts_df,
    meta_df=None,
    color_by=None,
    output_dir="static",
):
    if counts_df.shape[1] < 2:
        raise ValueError("PCA requires at least two samples.")
    log_counts = _log_cpm(counts_df)
    if len(log_counts) < 2:
        raise ValueError("PCA requires at least two genes.")
    if (log_counts.var(axis=1, ddof=0) > 0).sum() < 2:
        raise ValueError("PCA requires at least two genes with variation across samples.")

    pca = PCA(n_components=2)
    components = pca.fit_transform(log_counts.T)
    fig, ax = plt.subplots(figsize=(10, 7))

    if color_by:
        if meta_df is None or color_by not in meta_df.columns:
            raise ValueError("Choose a valid metadata column for PCA coloring.")
        aligned_meta = meta_df.reindex(counts_df.columns)
        if aligned_meta[color_by].isna().any():
            raise ValueError(f"Metadata column '{color_by}' has missing sample values.")
        labels = aligned_meta[color_by].astype(str)
        categories = labels.drop_duplicates().tolist()
        palette = sns.color_palette("husl", n_colors=len(categories))
        for label, color in zip(categories, palette):
            mask = labels.to_numpy() == label
            ax.scatter(
                components[mask, 0],
                components[mask, 1],
                color=color,
                alpha=0.8,
                label=label,
            )
        ax.legend(title=color_by, bbox_to_anchor=(1.02, 1), loc="upper left")
    else:
        ax.scatter(components[:, 0], components[:, 1], color="steelblue", alpha=0.8)

    if len(counts_df.columns) <= 120:
        for i, sample in enumerate(counts_df.columns):
            ax.annotate(
                str(sample),
                (components[i, 0], components[i, 1]),
                fontsize=5,
                alpha=0.7,
            )

    ax.set_title(f"PCA of log2(CPM + 1), all {len(log_counts)} filtered genes")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "pca_plot.png"), dpi=160)
    plt.close(fig)
    return "pca_plot.png"


def run_deseq2(
    counts_df,
    meta_df,
    design_col,
    target_group,
    reference_group,
    lfc_threshold=1.0,
    padj_threshold=0.05,
    label_top_genes=False,
    output_dir="static",
    results_dir=None,
    comparison_mode="pairwise",
):
    counts_df = validate_count_matrix(counts_df)
    common_samples = [sample for sample in counts_df.columns if sample in meta_df.index]
    if len(common_samples) < 2:
        raise ValueError("DESeq2 requires at least two matched samples.")
    if design_col not in meta_df.columns:
        raise ValueError("Choose a valid metadata condition column.")
    if not math.isfinite(lfc_threshold) or lfc_threshold < 0:
        raise ValueError("The absolute log2 fold-change threshold must be zero or greater.")
    if not 0 < padj_threshold < 1:
        raise ValueError("The adjusted p-value threshold must be between 0 and 1.")

    counts_df = counts_df.loc[:, common_samples]
    metadata = meta_df.loc[common_samples, [design_col]].copy()
    if metadata[design_col].isna().any():
        raise ValueError(f"Metadata column '{design_col}' has missing values.")
    metadata[design_col] = metadata[design_col].astype(str)
    target_group = str(target_group)
    group_counts = metadata[design_col].value_counts()
    if target_group not in group_counts.index:
        raise ValueError(f"Group '{target_group}' is not present in the selected samples.")
    if comparison_mode not in {"pairwise", "one_vs_rest"}:
        raise ValueError("Choose pairwise or one-versus-rest comparison mode.")
    if group_counts[target_group] < 2:
        raise ValueError(
            f"Group '{target_group}' has fewer than two samples; DESeq2 requires replicates."
        )

    if comparison_mode == "one_vs_rest":
        other_groups = [
            str(group) for group in metadata[design_col].drop_duplicates()
            if str(group) != target_group
        ]
        if not other_groups:
            raise ValueError(
                "One-versus-rest comparison requires at least two condition levels."
            )
        pooled_reference_count = int(
            metadata[design_col].ne(target_group).sum()
        )
        if pooled_reference_count < 2:
            raise ValueError(
                "The pooled reference group has fewer than two samples; "
                "DESeq2 requires replicates."
            )
        reference_group = "all other groups"
        model_design_col = "__maven_contrast_group__"
        metadata = pd.DataFrame(
            {
                model_design_col: np.where(
                    metadata[design_col].eq(target_group),
                    "__maven_target__",
                    "__maven_rest__",
                )
            },
            index=metadata.index,
        )
        model_target_group = "__maven_target__"
        model_reference_group = "__maven_rest__"
        comparison_description = (
            f"{target_group} vs all other groups ({', '.join(other_groups)})"
        )
    else:
        reference_group = str(reference_group)
        if target_group == reference_group:
            raise ValueError("Target and reference groups must be different.")
        if reference_group not in group_counts.index:
            raise ValueError(
                f"Group '{reference_group}' is not present in the selected samples."
            )
        if group_counts[reference_group] < 2:
            raise ValueError(
                f"Group '{reference_group}' has fewer than two samples; "
                "DESeq2 requires replicates."
            )
        model_design_col = design_col
        model_target_group = target_group
        model_reference_group = reference_group
        comparison_description = f"{target_group} vs {reference_group}"

    os.makedirs(output_dir, exist_ok=True)
    results_dir = results_dir or output_dir
    os.makedirs(results_dir, exist_ok=True)
    n_cpus = min(4, os.cpu_count() or 1, max(1, counts_df.shape[0] // 5000))
    dds = DeseqDataSet(
        counts=counts_df.T,
        metadata=metadata,
        design_factors=model_design_col,
        n_cpus=n_cpus,
        # Cook's filtering remains enabled in DeseqStats; this avoids the
        # incompatible outlier-refit path in the pinned PyDESeq2/AnnData stack.
        refit_cooks=False,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=[model_design_col, model_target_group, model_reference_group],
        alpha=padj_threshold,
        n_cpus=n_cpus,
        quiet=True,
    )
    with contextlib.redirect_stdout(io.StringIO()):
        stats.summary()
    results = stats.results_df.copy()
    results["comparison"] = comparison_description
    results["comparison_mode"] = comparison_mode

    valid_padj = results["padj"].notna()
    significant = (
        valid_padj
        & (results["padj"] < padj_threshold)
        & (results["log2FoldChange"].abs() >= lfc_threshold)
    )
    results["is_significant"] = significant
    results["regulation"] = np.select(
        [
            significant & (results["log2FoldChange"] > 0),
            significant & (results["log2FoldChange"] < 0),
        ],
        ["Upregulated", "Downregulated"],
        default="Non-significant",
    )

    plot_results = results.loc[results["log2FoldChange"].notna()].copy()
    plot_results["plot_padj"] = (
        plot_results["padj"].fillna(1).clip(lower=np.finfo(float).tiny)
    )
    plot_results["neg_log10_padj"] = -np.log10(plot_results["plot_padj"])
    colors = {
        "Upregulated": "#d62728",
        "Downregulated": "#1f77b4",
        "Non-significant": "#9e9e9e",
    }

    fig, ax = plt.subplots(figsize=(9, 7))
    for category in ("Non-significant", "Downregulated", "Upregulated"):
        subset = plot_results[plot_results["regulation"] == category]
        ax.scatter(
            subset["log2FoldChange"],
            subset["neg_log10_padj"],
            s=12,
            alpha=0.65,
            color=colors[category],
            label=category,
            edgecolors="none",
        )
    ax.axhline(-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)
    ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=1)
    ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=1)
    ax.set_title(f"{comparison_description} ({design_col})")
    ax.set_xlabel("Log2 fold change")
    ax.set_ylabel("-Log10(adjusted p-value)")
    ax.legend()

    if label_top_genes:
        labeled = results.loc[significant].sort_values("padj").head(10)
        for gene, row in labeled.iterrows():
            ax.annotate(
                str(gene),
                (row["log2FoldChange"], -np.log10(max(row["padj"], np.finfo(float).tiny))),
                fontsize=7,
            )

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "volcano_plot.png"), dpi=160)
    plt.close(fig)
    results.to_csv(os.path.join(results_dir, "deseq_results.csv"), index_label="gene")
    return "volcano_plot.png", "deseq_results.csv", results.index[significant].tolist()


def run_kmeans_clustering(
    counts_df,
    n_clusters=3,
    significant_genes=None,
    output_dir="static",
    results_dir=None,
):
    os.makedirs(output_dir, exist_ok=True)
    results_dir = results_dir or output_dir
    os.makedirs(results_dir, exist_ok=True)
    log_counts = _log_cpm(counts_df)
    if significant_genes is None:
        selected_genes = _top_variable_genes(log_counts, max_genes=1000)
    else:
        available = log_counts.index.intersection(pd.Index(significant_genes))
        if len(available) == 0:
            raise ValueError("DESeq2 has no significant genes to cluster.")
        selected_genes = available
        if len(selected_genes) > 1000:
            selected_genes = _top_variable_genes(
                log_counts.loc[selected_genes], max_genes=1000
            )

    expression = log_counts.loc[selected_genes]
    standard_deviation = expression.std(axis=1, ddof=0)
    expression = expression.loc[standard_deviation > 0]
    if len(expression) < n_clusters:
        raise ValueError(
            f"K-Means needs at least {n_clusters} non-constant genes; "
            f"{len(expression)} are available."
        )
    scaled = expression.sub(expression.mean(axis=1), axis=0).div(
        expression.std(axis=1, ddof=0), axis=0
    )
    labels = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(
        scaled
    )
    assignments = pd.DataFrame(
        {"gene": scaled.index, "cluster": labels}
    ).sort_values(["cluster", "gene"], kind="stable")
    assignments.to_csv(
        os.path.join(results_dir, "cluster_assignments.csv"),
        index=False,
    )

    sorted_genes = assignments["gene"].tolist()
    cluster_by_gene = assignments.set_index("gene")["cluster"]
    heatmap = scaled.loc[sorted_genes]
    cluster_labels = cluster_by_gene.loc[sorted_genes].to_numpy()
    palette = sns.color_palette("husl", n_colors=n_clusters)
    cluster_colors = ListedColormap(palette)

    fig = plt.figure(figsize=(14, 8))
    grid = fig.add_gridspec(1, 2, width_ratios=[0.035, 1], wspace=0.04)
    annotation_ax = fig.add_subplot(grid[0, 0])
    heatmap_ax = fig.add_subplot(grid[0, 1])
    annotation_ax.imshow(
        cluster_labels[:, None],
        aspect="auto",
        interpolation="nearest",
        cmap=cluster_colors,
        vmin=-0.5,
        vmax=n_clusters - 0.5,
    )
    annotation_ax.set_axis_off()
    sns.heatmap(
        heatmap,
        ax=heatmap_ax,
        cmap="RdBu_r",
        center=0,
        yticklabels=False,
        xticklabels=True,
        cbar_kws={"label": "Row-wise z-score"},
    )
    heatmap_ax.set_title(f"K-Means on normalized expression (k={n_clusters})")
    heatmap_ax.set_xlabel("Samples")
    heatmap_ax.set_ylabel("Genes")
    legend = [
        Patch(facecolor=palette[i], label=f"Cluster {i}")
        for i in range(n_clusters)
    ]
    heatmap_ax.legend(handles=legend, title="Gene cluster", loc="upper left")
    fig.subplots_adjust(left=0.07, right=0.97, top=0.94, bottom=0.08, wspace=0.08)
    fig.savefig(os.path.join(output_dir, "kmeans_heatmap.png"), dpi=160)
    plt.close(fig)
    return "kmeans_heatmap.png", "cluster_assignments.csv"


def run_enrichment(
    results_df,
    gene_sets="KEGG_2021_Human",
    direction="combined",
    output_dir="static",
):
    if direction not in {"up", "down", "combined"}:
        raise ValueError("Enrichment direction must be up, down, or combined.")
    if "padj" not in results_df or "log2FoldChange" not in results_df:
        raise ValueError("DESeq2 results must include padj and log2FoldChange.")

    if "is_significant" not in results_df:
        raise ValueError("DESeq2 results must include the dual-threshold DEG classification.")
    significant = results_df.loc[
        results_df["padj"].notna() & results_df["is_significant"].astype(bool)
    ].copy()
    if direction == "up":
        significant = significant[significant["log2FoldChange"] > 0]
    elif direction == "down":
        significant = significant[significant["log2FoldChange"] < 0]
    significant = significant.sort_values("padj").head(500)
    gene_list = [
        str(gene)
        for gene in significant.index
        if str(gene).strip() and str(gene).lower() != "nan"
    ]
    if not gene_list:
        return None, []

    try:
        enrichment = gp.enrichr(
            gene_list=gene_list, gene_sets=gene_sets, outdir=None
        )
    except requests.exceptions.RequestException as exc:
        raise ValueError(f"Enrichr request failed: {exc}") from exc
    results = enrichment.results
    if results is None or results.empty:
        return None, []

    results = results.sort_values("Adjusted P-value").head(10).copy()
    results["-log10(P-value)"] = -np.log10(
        results["P-value"].clip(lower=np.finfo(float).tiny)
    )
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.barplot(data=results, x="-log10(P-value)", y="Term", color="steelblue", ax=ax)
    direction_title = {
        "up": "Upregulated",
        "down": "Downregulated",
        "combined": "Combined",
    }[direction]
    ax.set_title(f"{direction_title} enrichment ({gene_sets})")
    fig.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    fig.savefig(os.path.join(output_dir, "enrichment_plot.png"), dpi=160)
    plt.close(fig)
    table_data = results[
        ["Term", "Overlap", "P-value", "Adjusted P-value"]
    ].to_dict(orient="records")
    return "enrichment_plot.png", table_data
