# MAVEN RNA-seq Analysis Web App — User Guide

MAVEN is a local Flask application for quality filtering, visualization,
differential-expression analysis, gene clustering, and pathway over-
representation analysis. It expects raw gene-level integer counts; CPM,
TPM, FPKM, and other normalized values are not suitable input for DESeq2.
The results are exploratory and should be interpreted with appropriate
biological and statistical review.

## Start the application

Use Python 3.11 in PowerShell from this project directory:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
$env:MAVEN_SECRET_KEY = (python -c "import secrets; print(secrets.token_hex(32))")
python app.py
```

Open <http://127.0.0.1:5000>. Keep the terminal open while using the app and
press **Ctrl+C** to stop Flask. Set `MAVEN_SECRET_KEY` to a stable, private
random value for each deployment; do not commit it. The app refuses to start
without this setting. The Flask development server is for local trusted use
only; do not expose it directly to the public internet.

## Input CSV files

Upload UTF-8 comma-separated files with one header row and no notes or blank
rows above it.

### Raw count table

- One or more rows per gene ID and one column per sample.
- The first column contains a gene identifier (for example, a symbol or
  Ensembl ID). Repeated identical IDs are combined into one gene by summing
  their raw counts for each sample; the first occurrence determines row order.
- Remaining columns contain finite, non-negative integer raw read counts.
- Sample column names must match metadata sample IDs exactly.
- Blank gene identifiers and duplicate or blank sample identifiers are
  rejected. Duplicate gene IDs are aggregated, not rejected.

```csv
gene_id,Control_1,Control_2,Treated_1,Treated_2
GeneA,120,98,240,221
GeneB,0,3,18,12
```

### Sample metadata

- One row per sample, with the first column containing unique sample IDs.
- Other columns contain sample annotations such as `Condition`, `Batch`,
  `CellType`, or `Lineage`.
- Values in the selected DESeq2 condition column must be present for all
  samples used in that contrast.

```csv
sample_id,Condition,Batch
Control_1,Control,Batch1
Control_2,Control,Batch2
Treated_1,Treated,Batch1
Treated_2,Treated,Batch2
```

If a sample appears in only one file, MAVEN warns and analyzes the matched
samples only. It rejects duplicate sample IDs, invalid count values, and files
with no matched samples. Identical gene IDs in multiple count-table rows are
merged by summing the raw counts across those rows before filtering and
analysis. Ensure repeated IDs genuinely refer to the same gene; rows for
different genes that happen to share a symbol will also be combined.

## Workflow

### 1. Upload and filter

Upload the count and metadata tables. The filter page displays an interactive
histogram and retained-gene curve. The histogram shows the distribution of
each gene's Nth-highest raw count across samples; the curve shows the number
of genes retained as the cutoff increases. Set:

- **Minimum raw reads per sample** (`min_count`; default `10`).
- **Minimum number of samples** (`N`; default `2`, capped at the available
  sample count).

A gene is retained if its raw count is at least `min_count` in at least `N`
matched samples. Drag the highlighted point along the retained-gene curve to
set the cutoff, or type the cutoff directly; changing either value updates
the plot and retained-gene count before applying the filter. Filtering uses
the matched sample columns. Revisit the filter step to change the thresholds;
doing so invalidates previous DESeq2 and clustering-download state.

### 2. PCA

Choose a metadata column to color samples (default: `Lineage`, when present)
or choose **No grouping**. PCA uses log2(CPM + 1), where each sample is
normalized by its total library size, and the top 1,000 variable genes (or
fewer if fewer are available). Points are labeled with sample IDs and a
legend identifies the selected groups. PCA is a descriptive visualization,
not a differential-expression test.

### 3. Differential expression with DESeq2

Choose one categorical metadata column, then select a **target condition
level**. Choose either a pairwise comparison with one reference level or
compare the target level against all other levels pooled together. For the
pooled comparison, all samples not in the target level are combined into one
reference group. The estimated log2 fold change is target relative to the
reference. Each resulting group must have at least two samples. This interface
fits a single-factor design; it does not adjust for additional covariates such
as batch.

Each selected group must have at least two matched samples. DEGs require both
`padj < 0.05` and `abs(log2FoldChange) >= threshold` (default threshold `1`,
editable on the dashboard). Missing adjusted p-values are treated as
non-significant, and zero adjusted p-values are safely bounded for plotting.
The volcano plot shows upregulated genes in red, downregulated genes in blue,
and other genes in gray, with effect-size and adjusted-p-value threshold
lines. Optionally label the ten most significant qualifying genes.
Cook's-distance filtering is retained; count refitting of outliers is disabled
for compatibility with the pinned PyDESeq2/AnnData versions.

Analysis forms display an indeterminate progress bar and a wait message while
uploads, filtering, PCA, DESeq2, clustering, or enrichment are running; actual
percent completion is not available.

Download the complete DESeq2 results CSV, which includes the selected
comparison, result statistics, and DEG classification columns.

### 4. K-means clustering

Before DESeq2 has been run, clustering uses up to 1,000 most-variable genes.
After a DESeq2 run, it uses the significant DEGs from that run. If that run
produced no DEGs, MAVEN reports that there are no genes to cluster instead of
silently switching to a different feature set.

Expression is transformed to log2(CPM + 1), then standardized row-wise across
samples to emphasize shared expression patterns. Constant genes are removed.
The heatmap uses a diverging scale centered at zero and a row annotation for
cluster membership. Choose `k` from 2 to 10; `k` cannot exceed the number of
usable genes. Download the gene-to-cluster CSV from the results page.

### 5. Directional pathway enrichment

Run DESeq2 first. MAVEN uses only DEGs passing both thresholds and sorts them
by adjusted p-value before selecting genes for Enrichr over-representation
analysis (up to 500 genes). Choose **Upregulated**, **Downregulated**, or
**Combined**. Each directional set is ranked and capped independently;
Combined uses the top 500 overall. Supported libraries are KEGG 2021 Human
and GO Biological Process 2021. This step requires internet access and
identifiers compatible with the selected library. It is over-representation
analysis, not a full ranked GSEA test.

## Sessions, files, and troubleshooting

Each upload receives a unique analysis directory, so concurrent analyses do
not overwrite one another. Download routes only serve the current session's
generated CSVs. Files remain on disk under `uploads/<analysis-id>/` and
`static/<analysis-id>/`; remove old analysis directories manually when they
are no longer needed. Do not delete files while an analysis is running.

- **Input errors:** Check the first-column IDs, unique sample identifiers, raw
  integer counts, CSV encoding, and matched sample names.
- **PCA errors:** At least two samples and two variable genes are needed.
  Remove zero-library samples and inspect the filtered data if CPM cannot be
  calculated.
- **DESeq2 errors:** Confirm both groups have at least two samples and that
  the selected condition values are non-empty. The single-factor model does
  not correct for batch or other covariates.
- **Clustering errors:** Ensure enough non-constant genes remain for the
  selected `k`. If DESeq2 has been run, its significant DEG set is used.
- **No enrichment result:** Confirm DESeq2 produced qualifying genes and
  that gene identifiers match the selected human gene-set library. Enrichr
  needs an internet connection.
- **Session unavailable:** Upload the files again. A new upload creates a
  new isolated analysis directory.
