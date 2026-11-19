from flask import Flask, render_template, request, session, redirect, url_for
import os
from werkzeug.utils import secure_filename
import pandas as pd
from pipeline import process_counts_and_filter, run_pca_analysis, run_deseq2, run_kmeans_clustering, run_enrichment

app = Flask(__name__)
app.secret_key = 'super_secret_key_for_maven_project'
app.config['UPLOAD_FOLDER'] = 'uploads'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

@app.route('/')
def index():
    return render_template('upload.html')

@app.route('/upload', methods=['POST'])
def upload_files():
    counts_file = request.files.get('counts_file')
    meta_file = request.files.get('meta_file')
    if not counts_file or not meta_file: 
        return "Missing files", 400
    
    counts_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(counts_file.filename))
    meta_path = os.path.join(app.config['UPLOAD_FOLDER'], secure_filename(meta_file.filename))
    
    counts_file.save(counts_path)
    meta_file.save(meta_path)
    
    session['counts_path'] = counts_path
    session['meta_path'] = meta_path
    
    return redirect(url_for('filter_page'))

@app.route('/filter', methods=['GET', 'POST'])
def filter_page():
    if request.method == 'POST':
        min_count = int(request.form.get('min_count', 10))
        counts_path = session.get('counts_path')
        
        counts_df = pd.read_csv(counts_path, index_col=0)
        filtered_df = process_counts_and_filter(counts_df, min_count)
        
        filtered_path = os.path.join(app.config['UPLOAD_FOLDER'], 'filtered_counts.csv')
        filtered_df.to_csv(filtered_path)
        session['filtered_path'] = filtered_path
        
        return redirect(url_for('dashboard_page'))
        
    return render_template('filter.html')

@app.route('/dashboard')
def dashboard_page():
    meta_path = session.get('meta_path')
    if not meta_path: 
        return redirect(url_for('index'))
    
    meta_df = pd.read_csv(meta_path, index_col=0)
    columns = meta_df.columns.tolist()
    return render_template('dashboard.html', columns=columns)

@app.route('/run_pca')
def do_pca():
    filtered_path = session.get('filtered_path')
    df = pd.read_csv(filtered_path, index_col=0)
    plot = run_pca_analysis(df)
    return render_template('results.html', title="PCA Results", plot=plot)

@app.route('/run_deseq', methods=['POST'])
def do_deseq():
    design_col = request.form.get('design_col')
    counts_path = session.get('filtered_path')
    meta_path = session.get('meta_path')
    
    counts_df = pd.read_csv(counts_path, index_col=0)
    meta_df = pd.read_csv(meta_path, index_col=0)
    
    try:
        plot, csv_file, sig_genes = run_deseq2(counts_df, meta_df, design_col)
        session['sig_genes'] = sig_genes
        msg = f"DESeq2 completed successfully. {len(sig_genes)} significant genes found."
    except Exception as e:
        plot = None
        csv_file = None
        msg = f"Error running DESeq2: {str(e)}"
        
    return render_template('results.html', title="Differential Expression (DESeq2)", plot=plot, message=msg, download_file=csv_file)

@app.route('/run_clustering', methods=['POST'])
def do_clustering():
    k = int(request.form.get('k_clusters', 3))
    counts_path = session.get('filtered_path')
    df = pd.read_csv(counts_path, index_col=0)
    
    plot = run_kmeans_clustering(df, n_clusters=k)
    return render_template('results.html', title=f"K-Means Clustering (k={k})", plot=plot)

@app.route('/run_enrichment', methods=['POST'])
def do_enrichment():
    db = request.form.get('database')
    sig_genes = session.get('sig_genes', [])
    
    if not sig_genes:
        counts_path = session.get('filtered_path')
        df = pd.read_csv(counts_path, index_col=0)
        # Fallback: top 500 variable genes if DESeq2 wasn't run
        sig_genes = df.var(axis=1).nlargest(500).index.tolist()
    
    try:
        plot, table_data = run_enrichment(sig_genes, gene_sets=db)
        if plot:
            msg = f"Enrichment analysis completed using {db}."
        else:
            msg = "No significant enrichment pathways found for these genes."
    except Exception as e:
        plot = None
        table_data = []
        msg = f"Error running Enrichment: {str(e)}"
        
    return render_template('results.html', title="Pathway Enrichment Analysis", plot=plot, message=msg, table_data=table_data)

if __name__ == '__main__':
    app.run(debug=True, port=5000)