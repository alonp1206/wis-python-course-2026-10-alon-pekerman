import csv
import os
import re
import uuid

import pandas as pd
from flask import (
    Flask,
    flash,
    jsonify,
    redirect,
    render_template,
    request,
    send_from_directory,
    session,
    url_for,
)

from pipeline import (
    filter_preview_data,
    process_counts_and_filter,
    run_deseq2,
    run_enrichment,
    run_kmeans_clustering,
    run_pca_analysis,
    validate_input_tables,
)


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
UPLOAD_ROOT = os.path.join(BASE_DIR, "uploads")
STATIC_ROOT = os.path.join(BASE_DIR, "static")
SECRET_KEY = os.environ.get("MAVEN_SECRET_KEY")
if not SECRET_KEY:
    raise RuntimeError(
        "Set the MAVEN_SECRET_KEY environment variable before starting the app."
    )

app = Flask(__name__, static_folder=STATIC_ROOT)
app.secret_key = SECRET_KEY
app.config["UPLOAD_FOLDER"] = UPLOAD_ROOT
app.config["MAX_CONTENT_LENGTH"] = 256 * 1024 * 1024
os.makedirs(UPLOAD_ROOT, exist_ok=True)
os.makedirs(STATIC_ROOT, exist_ok=True)

def _run_directory(run_id):
    if not run_id or not re.fullmatch(r"[a-f0-9]{32}", run_id):
        return None
    directory = os.path.join(UPLOAD_ROOT, run_id)
    return directory if os.path.isdir(directory) else None


def _read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as source:
        header = next(csv.reader(source), [])
    if len(header) < 2:
        raise ValueError(f"{os.path.basename(path)} must have an ID column and data columns.")
    if len(header) != len(set(header)):
        raise ValueError(f"{os.path.basename(path)} contains duplicate column names.")
    data = pd.read_csv(path, index_col=0, dtype={header[0]: str})
    data.index = pd.Index(
        [str(identifier) if pd.notna(identifier) else identifier for identifier in data.index]
    )
    data.columns = data.columns.astype(str)
    return data


def _load_tables(run_id):
    directory = _run_directory(run_id)
    if directory is None:
        raise ValueError("Analysis files are unavailable. Upload the input files again.")
    counts_df = _read_csv(os.path.join(directory, "counts.csv"))
    meta_df = _read_csv(os.path.join(directory, "metadata.csv"))
    original_gene_rows = len(counts_df)
    counts_df = validate_input_tables(counts_df, meta_df)

    count_samples = set(counts_df.columns)
    metadata_samples = set(meta_df.index)
    common_samples = [sample for sample in counts_df.columns if sample in metadata_samples]
    if not common_samples:
        raise ValueError("The count and metadata files have no matching sample IDs.")
    unmatched_counts = count_samples - metadata_samples
    unmatched_metadata = metadata_samples - count_samples
    warnings = []
    if unmatched_counts or unmatched_metadata:
        warnings.append(
            "Sample ID mismatch: "
            f"{len(unmatched_counts)} count-table samples and "
            f"{len(unmatched_metadata)} metadata samples have no match. "
            "Analyses will use matched samples only."
        )
    merged_gene_rows = original_gene_rows - len(counts_df)
    if merged_gene_rows:
        warnings.append(
            f"Merged {merged_gene_rows} duplicate gene row(s) by summing their raw "
            f"counts; {len(counts_df)} unique gene IDs remain."
        )
    return counts_df.loc[:, common_samples], meta_df.loc[common_samples], warnings


def _filtered_counts(run_id):
    directory = _run_directory(run_id)
    if directory is None:
        raise ValueError("Analysis files are unavailable. Upload the input files again.")
    filtered_path = os.path.join(directory, "filtered_counts.csv")
    if not os.path.isfile(filtered_path):
        raise ValueError("Filter the counts before running an analysis.")
    return _read_csv(filtered_path)


def _render_results(title, plot=None, message=None, download_artifact=None, table_data=None):
    run_id = session.get("run_id")
    plot_path = f"{run_id}/{plot}" if plot and run_id else None
    return render_template(
        "results.html",
        title=title,
        plot=plot_path,
        message=message,
        download_artifact=download_artifact,
        table_data=table_data or [],
    )


def _form_integer(name, default):
    value = request.form.get(name, str(default))
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name.replace('_', ' ').capitalize()} must be a whole number.") from exc


@app.route("/")
def index():
    return render_template("upload.html")


@app.route("/upload", methods=["POST"])
def upload_files():
    counts_file = request.files.get("counts_file")
    meta_file = request.files.get("meta_file")
    if not counts_file or not meta_file or not counts_file.filename or not meta_file.filename:
        flash("Choose both a count CSV and a metadata CSV.", "danger")
        return redirect(url_for("index"))
    if not counts_file.filename.lower().endswith(".csv") or not meta_file.filename.lower().endswith(".csv"):
        flash("Both uploaded files must use the .csv extension.", "danger")
        return redirect(url_for("index"))

    run_id = uuid.uuid4().hex
    run_dir = os.path.join(UPLOAD_ROOT, run_id)
    os.makedirs(run_dir)
    counts_path = os.path.join(run_dir, "counts.csv")
    metadata_path = os.path.join(run_dir, "metadata.csv")
    counts_file.save(counts_path)
    meta_file.save(metadata_path)

    try:
        counts_df, meta_df, warnings = _load_tables(run_id)
    except ValueError as exc:
        flash(str(exc), "danger")
        return redirect(url_for("index"))

    session.clear()
    session["run_id"] = run_id
    for warning in warnings:
        flash(warning, "warning")
    flash(
        f"Loaded {len(counts_df.columns)} matched samples and {len(counts_df)} genes.",
        "success",
    )
    return redirect(url_for("filter_page"))


@app.route("/filter", methods=["GET", "POST"])
def filter_page():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None:
        flash("Upload count and metadata files to start an analysis.", "warning")
        return redirect(url_for("index"))
    try:
        counts_df, _, warnings = _load_tables(run_id)
        if request.method == "POST":
            session.pop("filter_applied", None)
            session.pop("deseq_done", None)
            session.pop("cluster_file_ready", None)
            min_count = _form_integer("min_count", 10)
            min_samples = _form_integer("min_samples", 2)
            filtered_df = process_counts_and_filter(
                counts_df, min_count=min_count, min_samples=min_samples
            )
            filtered_df.to_csv(os.path.join(_run_directory(run_id), "filtered_counts.csv"))
            session["filter_applied"] = True
            session.pop("deseq_done", None)
            session.pop("cluster_file_ready", None)
            flash(
                f"Kept {len(filtered_df)} of {len(counts_df)} genes "
                f"(at least {min_count} reads in at least {min_samples} samples).",
                "success",
            )
            return redirect(url_for("dashboard_page"))

        default_samples = min(2, len(counts_df.columns))
        preview = filter_preview_data(counts_df, default_samples, 10)
        return render_template(
            "filter.html",
            sample_count=len(counts_df.columns),
            default_samples=default_samples,
            preview=preview,
            warnings=warnings,
        )
    except ValueError as exc:
        flash(str(exc), "danger")
        return redirect(url_for("filter_page") if request.method == "POST" else url_for("index"))


@app.route("/filter_preview")
def filter_preview():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None:
        return jsonify(error="Upload the input files again."), 400
    try:
        counts_df, _, _ = _load_tables(run_id)
        min_samples = int(request.args.get("min_samples", "2"))
        min_count = int(request.args.get("min_count", "10"))
        preview = filter_preview_data(counts_df, min_samples, min_count)
        return jsonify(preview)
    except ValueError as exc:
        return jsonify(error=str(exc)), 400


@app.route("/dashboard")
def dashboard_page():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None:
        flash("Upload count and metadata files to start an analysis.", "warning")
        return redirect(url_for("index"))
    if not session.get("filter_applied"):
        return redirect(url_for("filter_page"))
    try:
        _, meta_df, warnings = _load_tables(run_id)
    except ValueError as exc:
        flash(str(exc), "danger")
        return redirect(url_for("index"))
    columns = meta_df.columns.tolist()
    metadata_levels = {
        column: [str(value) for value in meta_df[column].dropna().drop_duplicates()]
        for column in columns
    }
    return render_template(
        "dashboard.html",
        columns=columns,
        metadata_levels=metadata_levels,
        default_color="Lineage" if "Lineage" in columns else "",
        warnings=warnings,
    )


@app.route("/run_pca")
def do_pca():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None or not session.get("filter_applied"):
        return redirect(url_for("filter_page"))
    try:
        counts_df = _filtered_counts(run_id)
        _, meta_df, _ = _load_tables(run_id)
        color_by = request.args.get("color_by") or None
        plot = run_pca_analysis(
            counts_df,
            meta_df=meta_df,
            color_by=color_by,
            output_dir=os.path.join(STATIC_ROOT, run_id),
        )
        return _render_results("PCA Results", plot=plot)
    except ValueError as exc:
        return _render_results("PCA Results", message=str(exc))


@app.route("/run_deseq", methods=["POST"])
def do_deseq():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None or not session.get("filter_applied"):
        return redirect(url_for("filter_page"))
    session.pop("deseq_done", None)
    session.pop("cluster_file_ready", None)
    try:
        counts_df = _filtered_counts(run_id)
        _, meta_df, _ = _load_tables(run_id)
        design_col = request.form.get("design_col", "")
        target_group = request.form.get("target_group", "")
        reference_group = request.form.get("reference_group", "")
        comparison_mode = request.form.get("comparison_mode", "pairwise")
        if comparison_mode not in {"pairwise", "one_vs_rest"}:
            raise ValueError("Choose pairwise or one-versus-rest comparison mode.")
        lfc_threshold = float(request.form.get("lfc_threshold", "1"))
        label_top_genes = request.form.get("label_top_genes") == "on"
        if comparison_mode == "one_vs_rest":
            reference_group = ""
        plot, _, significant_genes = run_deseq2(
            counts_df,
            meta_df,
            design_col,
            target_group,
            reference_group,
            comparison_mode=comparison_mode,
            lfc_threshold=lfc_threshold,
            label_top_genes=label_top_genes,
            output_dir=os.path.join(STATIC_ROOT, run_id),
            results_dir=_run_directory(run_id),
        )
        session["deseq_done"] = True
        session["cluster_file_ready"] = False
        if comparison_mode == "one_vs_rest":
            other_groups = [
                str(group)
                for group in meta_df[design_col].dropna().astype(str).drop_duplicates()
                if str(group) != target_group
            ]
            comparison_description = (
                f"{target_group} vs all other groups ({', '.join(other_groups)})"
            )
        else:
            comparison_description = f"{target_group} vs {reference_group}"
        return _render_results(
            "Differential Expression (DESeq2)",
            plot=plot,
            message=(
                f"Contrast: {comparison_description}. "
                f"{len(significant_genes)} genes meet padj < 0.05 and "
                f"|log2 fold change| >= {lfc_threshold:g}."
            ),
            download_artifact="deseq",
        )
    except (ValueError, TypeError) as exc:
        return _render_results("Differential Expression (DESeq2)", message=str(exc))


@app.route("/run_clustering", methods=["POST"])
def do_clustering():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None or not session.get("filter_applied"):
        return redirect(url_for("filter_page"))
    session["cluster_file_ready"] = False
    try:
        n_clusters = _form_integer("k_clusters", 3)
        if n_clusters < 2 or n_clusters > 10:
            raise ValueError("Choose between 2 and 10 clusters.")
        counts_df = _filtered_counts(run_id)
        significant_genes = None
        if session.get("deseq_done"):
            results = pd.read_csv(
                os.path.join(_run_directory(run_id), "deseq_results.csv"),
                index_col="gene",
            )
            significant_genes = results.index[
                results["is_significant"].astype(bool)
            ].tolist()
        plot, _ = run_kmeans_clustering(
            counts_df,
            n_clusters=n_clusters,
            significant_genes=significant_genes,
            output_dir=os.path.join(STATIC_ROOT, run_id),
            results_dir=_run_directory(run_id),
        )
        session["cluster_file_ready"] = True
        return _render_results(
            f"K-Means Clustering (k={n_clusters})",
            plot=plot,
            download_artifact="clusters",
        )
    except (ValueError, TypeError) as exc:
        return _render_results("K-Means Clustering", message=str(exc))


@app.route("/run_enrichment", methods=["POST"])
def do_enrichment():
    run_id = session.get("run_id")
    if _run_directory(run_id) is None or not session.get("filter_applied"):
        return redirect(url_for("filter_page"))
    if not session.get("deseq_done"):
        return _render_results(
            "Pathway Enrichment Analysis",
            message="Run DESeq2 first. Enrichment requires genes passing both DEG thresholds.",
        )
    try:
        results_df = pd.read_csv(
            os.path.join(_run_directory(run_id), "deseq_results.csv"),
            index_col="gene",
        )
        gene_sets = request.form.get("database", "KEGG_2021_Human")
        direction = request.form.get("direction", "combined")
        plot, table_data = run_enrichment(
            results_df,
            gene_sets=gene_sets,
            direction=direction,
            output_dir=os.path.join(STATIC_ROOT, run_id),
        )
        direction_label = {
            "up": "upregulated",
            "down": "downregulated",
            "combined": "combined",
        }.get(direction, direction)
        message = (
            f"Enrichment completed for {direction_label} genes."
            if plot
            else f"No qualifying genes or enriched pathways were found for {direction_label}."
        )
        return _render_results(
            "Pathway Enrichment Analysis",
            plot=plot,
            message=message,
            table_data=table_data,
        )
    except ValueError as exc:
        return _render_results("Pathway Enrichment Analysis", message=str(exc))


@app.route("/download/<artifact>")
def download_file(artifact):
    run_id = session.get("run_id")
    run_dir = _run_directory(run_id)
    allowed_files = {
        "deseq": ("deseq_results.csv", bool(session.get("deseq_done"))),
        "clusters": (
            "cluster_assignments.csv",
            bool(session.get("cluster_file_ready")),
        ),
    }
    if run_dir is None or artifact not in allowed_files:
        flash("That download is not available in this analysis.", "danger")
        return redirect(url_for("dashboard_page"))
    filename, ready = allowed_files[artifact]
    if not ready or not os.path.isfile(os.path.join(run_dir, filename)):
        flash("That download has not been generated yet.", "warning")
        return redirect(url_for("dashboard_page"))
    return send_from_directory(run_dir, filename, as_attachment=True)


if __name__ == "__main__":
    app.run(debug=False, port=5000)
