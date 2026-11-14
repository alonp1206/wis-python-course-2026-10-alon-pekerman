import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp

def process_counts_and_filter(counts_df, min_count=10):
    return counts_df[counts_df.sum(axis=1) >= min_count]

def run_pca_analysis(counts_df, output_dir='static'):
    os.makedirs(output_dir, exist_ok=True)
    log_counts = np.log2(counts_df + 1)
    pca = PCA(n_components=2)
    components = pca.fit_transform(log_counts.T)
    
    plt.figure(figsize=(8, 6))
    plt.scatter(components[:, 0], components[:, 1], c='blue', alpha=0.7)
    plt.title('PCA of Gene Expression')
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    
    for i, sample in enumerate(counts_df.columns):
        plt.annotate(sample, (components[i, 0], components[i, 1]))
        
    plot_path = os.path.join(output_dir, 'pca_plot.png')
    plt.savefig(plot_path)
    plt.close()
    return 'pca_plot.png'

def run_deseq2(counts_df, meta_df, design_col, output_dir='static'):
    os.makedirs(output_dir, exist_ok=True)
    
    # Clean up types and indices
    counts_df = counts_df.astype(int)
    common_samples = counts_df.columns.intersection(meta_df.index)
    counts_df = counts_df[common_samples].T
    meta_df = meta_df.loc[common_samples]
    
    # Run PyDESeq2
    dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design_factors=design_col)
    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()
    res_df = stat_res.results_df
    
    # Generate Volcano Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=res_df['log2FoldChange'], y=-np.log10(res_df['padj']), alpha=0.6)
    plt.axhline(-np.log10(0.05), color='red', linestyle='--')
    plt.title(f'Volcano Plot ({design_col})')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(padj)')
    
    volcano_path = os.path.join(output_dir, 'volcano_plot.png')
    plt.savefig(volcano_path)
    plt.close()
    
    # Save CSV
    csv_name = 'deseq_results.csv'
    csv_path = os.path.join(output_dir, csv_name)
    res_df.to_csv(csv_path)
    
    # Extract significant genes (padj < 0.05)
    sig_genes = res_df[res_df['padj'] < 0.05].index.tolist()
    return 'volcano_plot.png', csv_name, sig_genes

def run_kmeans_clustering(counts_df, n_clusters=3, output_dir='static'):
    os.makedirs(output_dir, exist_ok=True)
    log_counts = np.log2(counts_df + 1)
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(log_counts)
    log_counts['Cluster'] = clusters
    log_counts = log_counts.sort_values('Cluster')
    log_counts.drop('Cluster', axis=1, inplace=True)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(log_counts, cmap="viridis", yticklabels=False, xticklabels=True)
    plt.title(f'K-Means Clustering (k={n_clusters})')
    
    heatmap_path = os.path.join(output_dir, 'kmeans_heatmap.png')
    plt.tight_layout()
    plt.savefig(heatmap_path)
    plt.close()
    
    return 'kmeans_heatmap.png'

def run_enrichment(gene_list, gene_sets='KEGG_2021_Human', output_dir='static'):
    os.makedirs(output_dir, exist_ok=True)
    gene_list = [str(g) for g in gene_list if str(g) != 'nan']
    
    # Run GSEAPy Enrichr
    enr = gp.enrichr(gene_list=gene_list[:500], gene_sets=gene_sets, outdir=None)
    res_df = enr.results
    
    if res_df.empty:
        return None, []
        
    top_res = res_df.head(10).copy()
    top_res['-log10(P-value)'] = -np.log10(top_res['P-value'])
    
    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_res, x='-log10(P-value)', y='Term', color='steelblue')
    plt.title(f'Enrichment Analysis ({gene_sets})')
    plt.tight_layout()
    
    enr_plot_path = os.path.join(output_dir, 'enrichment_plot.png')
    plt.savefig(enr_plot_path)
    plt.close()
    
    table_data = top_res[['Term', 'Overlap', 'P-value', 'Adjusted P-value']].to_dict(orient='records')
    return 'enrichment_plot.png', table_data
