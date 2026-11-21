# MAVEN RNA-seq Analysis Web App

MAVEN is a local Flask web application for raw-count filtering, normalized
PCA, explicit-contrast DESeq2 analysis, expression-pattern clustering, and
directional Enrichr over-representation analysis.

## Quick start

From PowerShell in this directory:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
$env:MAVEN_SECRET_KEY = (python -c "import secrets; print(secrets.token_hex(32))")
python app.py
```

Open <http://127.0.0.1:5000>. The secret key is required and should remain
private. This app uses Flask's development server; keep it on a trusted local
machine.

Run the project workflow and input-validation tests from this directory:

```powershell
python -m unittest discover -s tests -v
```

Upload a gene-by-sample CSV of non-negative integer raw counts and a
sample-by-metadata CSV. Repeated identical gene IDs in the count file are
combined by summing their raw counts per sample. Filter genes by a minimum
read count in a minimum number of samples using the interactive histogram and
draggable retained-gene cutoff; the chart uses a logarithmic read-count axis
and percentages of genes. PCA includes all filtered genes. DESeq2 supports
pairwise and selected-level versus pooled-rest comparisons and uses up to four
CPU cores for larger fits while avoiding parallel startup overhead on small
inputs. For CSV formats, workflow details, thresholds, session storage, and
troubleshooting, see the [User Guide](./USER_GUIDE.md).

## Project files

- `app.py` — Flask routes, per-analysis session storage, and downloads.
- `pipeline.py` — filtering and analysis methods.
- `templates/` — upload, filter, dashboard, and result pages.
- `requirements.txt` — Python dependencies.

Uploaded data and generated files are stored in `uploads/` and `static/`,
which are intentionally excluded from version control.
