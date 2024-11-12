import os
import numpy as np
import scanpy as sc
import requests
import yaml
from bs4 import BeautifulSoup
from pdb import set_trace
import multiprocessing
import anndata
import matplotlib.pyplot as plt
# from mygene import MyGeneInfo
import seaborn as sns
import pandas as pd

sc.set_figure_params(figsize=(6, 6), frameon=False)
sc.settings.n_jobs=2

# mg = MyGeneInfo()

# Number of cores
num_cores = multiprocessing.cpu_count()
print("Number of cores: ", num_cores)

root_dir = '/grand/GeomicVar/rbr7/Geneformer/'
dirs = ["HNSCC_data"]


for dir in dirs:
    print("-----------------------Entering directory ---------------------: ", dir)
    adata = sc.read_h5ad(os.path.join(root_dir, dir, "hnscc_mod_B_filtered.h5ad"))

    # keep only the filtered rows
    # adata = adata[adata.obs['filter_pass'] == 1]
    print(adata)
    print(" +++++++++ Cell-level data (the indices are the bar-codes for each cell) +++++++++++++++")
    print(adata.obs.head())
    print(" +++++++++ Gene-level data (the indices represent all genes in the dataset) ++++++++++++ ")
    print(adata.var.head())
    total_read_count = np.sum(adata.X.A[:, :], axis=1)
    # total_read_count = np.sum(adata.X.A[:, :], axis=1)
    adata.obs['total_read_count'] = total_read_count
    adata.var['mitochondrial_count'] = adata.var['Gene name'].str.startswith('MT-')
    # pct_mito = mito_counts / total_counts


    sc.pp.calculate_qc_metrics(adata, qc_vars=['mitochondrial_count'], percent_top=None, log1p=False, inplace=True)
    # sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True)  
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_read_count', 'pct_counts_mitochondrial_count'], multi_panel=True)  
    plt.savefig(str(dir) + '-violin_plot_filtered.png', dpi=300, bbox_inches='tight')

    # sns.jointplot(data=adata.obs, x="log1p_total_counts", y="log1p_n_genes_by_counts", kind="hex",)
    # plt.savefig(str(dir) + '-jointplot_filtered.png', dpi=300, bbox_inches='tight')

    # sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
    # plt.savefig(str(dir) + '-total_counts_vs_n_genes_by_counts_filtered.png', dpi=300, bbox_inches='tight')

    # ax = sns.histplot(adata.obs["pct_counts_mt"])
    # plt.savefig(str(dir) + '-mito_hist_filtered.png', dpi=300, bbox_inches='tight')

    
    # PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # Scatter plot with PC 
    # sc.pl.pca(adata, color='CST3')
    # plt.savefig(str(dir) + '-PC_filtered.png', dpi=300, bbox_inches='tight')


    # calculate a nearest neighbor graph between cells
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
    # perform UMAP on the data
    sc.tl.louvain(adata, resolution = 0.25)

    # Visualizes the UMAP embedding produced by sc.tl.umap.
    sc.pl.umap(adata,color=['cell_type'])
    plt.savefig(str(dir) + '-umap_celltype_filtered.png', dpi=300, bbox_inches='tight')

    sc.pl.umap(adata, color=['disease'])
    plt.savefig(str(dir) + '-umap_disease_filtered.png', dpi=300, bbox_inches='tight')


#set_trace()

