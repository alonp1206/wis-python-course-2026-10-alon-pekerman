# MAVEN RNA-seq Analysis Web App — User Guide

This guide explains how to start and use the MAVEN web app, what the analysis
options do, and how to format the two CSV files it accepts. It describes the
current behavior of the app; it does not replace biological or statistical
review of the results.


## What the app does

MAVEN accepts a gene-by-sample read-count table and a sample metadata table.
It can filter low-count genes and then run:

- **PCA** to display samples in two dimensions after a `log2(count + 1)`
  transformation.
- **Differential expression (DESeq2)** using one metadata column as the
  condition/design factor. The app creates a volcano plot and a downloadable
  results CSV.
- **K-means clustering** of genes using their log-transformed counts, shown as
  a heatmap.
- **Pathway enrichment** using one of the Enrichr gene-set libraries offered
  in the page.

The web interface is implemented in `app.py`; filtering and analyses are
implemented in `pipeline.py`.

## Start the app

In PowerShell, change to the project folder and start Flask:

```powershell
Set-Location 'C:\Users\alon\Desktop\GIT\wis-python-course-2026-10-alon-pekerman\day04\MAVEN Project'
.\.venv\Scripts\python.exe app.py
```

Open <http://127.0.0.1:5000> in a browser. Keep the PowerShell window open while
using the app; press **Ctrl+C** there to stop the server. On a first launch,
scientific libraries may take time to initialize/build their font cache.

The `.venv` and dependency versions in `requirements.txt` are local to this
project. On a different computer, install a compatible Python (the prepared
environment uses Python 3.11), create a virtual environment, and install the
requirements before starting the app.

## Prepare the input CSV files

Both files must be ordinary comma-separated `.csv` files with a header row.
Save them as UTF-8 CSVs. Do not include an extra title line, notes, units, or
blank rows above the header.

### Read-count CSV

The count table must have:

- **One row per gene**.
- **One column per sample**.
- A first column containing unique gene identifiers (for example, gene
  symbols or Ensembl IDs). The app reads this first column as the row index;
  it should not be a numeric sample-count column.
- Numeric, non-negative count values in all remaining columns. Use raw integer
  read counts for DESeq2; do not use TPM, FPKM, CPM, normalized values, or
  percentages.
- Sample names in the header that match the sample identifiers in the
  metadata file exactly.

Example:

```csv
gene_id,Control_1,Control_2,Treated_1,Treated_2
GeneA,120,98,240,221
GeneB,0,3,18,12
GeneC,57,61,44,52
```

The first column is the gene ID. `Control_1`, `Control_2`, `Treated_1`, and
`Treated_2` are sample IDs; they must correspond to metadata rows below.
Avoid duplicate gene IDs or sample headers. Missing values and non-numeric
cells can cause analysis steps to fail.

### Metadata CSV

The metadata table must have:

- **One row per sample**.
- A first column containing the sample IDs. The app reads this column as the
  metadata index.
- One or more additional columns describing samples, such as `condition`,
  `batch`, or `tissue`.
- Sample IDs that match the read-count CSV column headers exactly (including
  spelling and capitalization).
- A condition/design column with a non-empty value for every sample when you
  run DESeq2.

Example:

```csv
sample_id,condition,batch
Control_1,Control,Batch1
Control_2,Control,Batch2
Treated_1,Treated,Batch1
Treated_2,Treated,Batch2
```

Here, `sample_id` matches the count-table headers, and `condition` and `batch`
are metadata columns. The app displays the columns after `sample_id` in the
DESeq2 selector. Select the column representing the comparison of interest.
For a basic two-group comparison, make sure the chosen column contains the
expected group labels for all samples and that each group has replicates.

**Important:** DESeq2 uses only sample IDs found in both files. The current app
does not warn if a sample ID is missing from either file; mismatches can
silently exclude samples. Check that the sets of count-table headers and
metadata sample IDs are identical before uploading.

## Use the web page

### 1. Upload

On the home page, choose the read-count CSV and metadata CSV, then select
**Upload & Continue**. Both files are required. The app saves the uploaded
files and moves to filtering.

### 2. Filter low-count genes

Enter **Minimum Total Count**. The app sums each gene's counts across all
samples and keeps genes whose sum is **greater than or equal to** the
threshold. The initial value is `10`. Use `0` to keep all rows under this
filter, or choose a threshold appropriate to your data and analysis.

Select **Filter Data & Run PCA** to apply the filter and continue to the
dashboard. Despite the button label, this step does not automatically run PCA;
it takes you to the dashboard, where you can choose an analysis.

### 3. Choose an analysis

On the dashboard, run one of the available analyses:

#### PCA & Transformation

Select **Run PCA**. The app applies `log2(count + 1)` and plots the samples on
the first two principal components. Each point is labeled with its sample ID.
This requires enough samples and genes with usable variation to calculate two
components.

#### Differential Expression (DESeq2)

Select a metadata column under **Condition Column**, then choose **Run DESeq2**.
The app aligns samples present in both tables, converts counts to integers,
runs PyDESeq2, and displays a volcano plot. It also writes a results CSV with
the statistics returned by PyDESeq2; use **Download CSV Data** to save it.
Genes with adjusted p-value (`padj`) below `0.05` are considered significant
for the app's later enrichment step.

Choose a meaningful categorical condition column and verify that it has
multiple groups and suitable replicates. The app reports an error on the
results page if the DESeq2 run fails.

#### K-Means Clustering

Enter the number of clusters (`k`, default `3`; the page allows `2` through
`10`) and choose **Run Clustering**. The app applies `log2(count + 1)`, clusters
genes by their sample profiles, and displays a heatmap. There must be at least
`k` genes remaining after filtering.

#### Pathway Enrichment

Choose **KEGG 2021 Human** or **GO Biological Process 2021**, then select **Run
Enrichment**. If DESeq2 has produced a non-empty significant-gene list in the
current browser session, that list is used. Otherwise, the app uses up to 500
genes with the highest variance across samples. Enrichr is queried through
GSEApy, so this step requires internet access. The page shows up to ten
returned terms and a plot when results are available.

### 4. Review results

Analysis results appear on a results page. Use **Back to Dashboard** to run
another analysis. The top navigation also provides **New Analysis (Upload)**.
Uploaded files and generated plots/CSVs are stored in the project's `uploads`
and `static` folders; output filenames are reused by later runs.

## Troubleshooting and limitations

- **Upload/filter/dashboard:** If a page reports missing session data, return
  to the home page and upload both files again.
- **DESeq2 error:** Check that counts are raw, numeric, non-negative integers;
  sample IDs match; the selected metadata column exists and has valid group
  labels; and enough samples/replicates are present.
- **PCA error:** PCA needs at least two samples and enough usable genes to
  calculate two components. Try checking the data and filtering threshold.
- **Clustering error:** Ensure at least `k` genes remain after filtering.
- **Enrichment error or empty result:** Check that gene identifiers are
  appropriate for the selected human gene-set library and that the app can
  access the internet.
- The app currently accepts uploads without validating their contents in
  advance. Use the CSV formats above to avoid downstream errors.
- Flask is configured with debug mode enabled and a development secret key.
  Run it only for local, trusted use; do not expose this development server to
  the public internet.