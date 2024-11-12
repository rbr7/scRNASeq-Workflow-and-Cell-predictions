import os
import numpy as np
import scanpy as sc
import requests
import yaml
from bs4 import BeautifulSoup
from pdb import set_trace
from joblib import Parallel, delayed
import multiprocessing
import matplotlib.pyplot as plt

# Number of cores 
num_cores = multiprocessing.cpu_count()
print("NUmber of cores: ", num_cores)

root_dir = '/gpfs/alpine/med112/proj-shared/tnandi/Geneformer/cancer_cellxgene_datasets/'
dirs = os.listdir(root_dir)
# dirs = ['glioblastoma', 'follicular_lymphoma', 'hepatoblastoma']
dirs = ['hepatoblastoma']
# dirs = ['B-cell_non-Hodgkin_lymphoma', 'pilocytic_astrocytoma']

with open('./ensembl_genes/geneid_proteincoding_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_to_search = f.read().splitlines()

# loop through all datasets
for dir in dirs:
    if os.path.isdir(dir): # exclude all files listed by os.listdir
        print("Entering directory ", dir)
        adata = sc.read_h5ad(os.path.join(root_dir, dir, "local.h5ad"))
        total_read_count = adata.obs['nCount_RNA']
        # Calculate total mitochondrial read count for each cell
        pct_mt_read_count = adata.obs['percent.mt']

        gene_counts = (adata.X > 0).sum(axis=1)
        plt.bar(gene_counts)
        plt.show()
        plt.savefig('hepatoblastoma.gene_count.png')

        set_trace()
        mito_genes = adata.var['feature_name'][adata.var['feature_name'].str.startswith('MI')]
        mito_ids = np.array(mito_genes.index)
        mito_counts = np.sum(adata[:, mito_ids].X, axis=1)

        # Get the mean and std of read count across the samples
        mean_read_count         = np.mean(total_read_count) 
        std_read_count          = np.std(total_read_count) 
        mean_mt_read_count      = np.mean(pct_mt_read_count)
        std_mt_read_count       = np.std(pct_mt_read_count)

    
        # Create a new column "filter_pass" and set 1 for passed samples 
        adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) <= 3 * std_read_count) & (np.abs(pct_mt_read_count - mean_mt_read_count) <= 3 * std_mt_read_count) & (total_read_count > 7), 1, 0)

        # Loop over all genes and create a new column called "gene_type" in adata.var
        gene_types = []
        genes =  adata.var.index
        print("looping over ", len(genes), " genes")
        # for gene in genes:
        #       print("gene: ", gene)  
        #       # If gene is in ./ensembl_genes/geneid_proteincoding_mirna.csv, set gene_types=1, otherwise set it to 0
        #       if gene in gene_ids_to_search:
        #           gene_types.append(1)
        #       else:
        #           gene_types.append(0)


        # function to check each gene
        def check_gene(gene):
            print(gene)
            if gene in gene_ids_to_search:
                return 1
            else:
                return 0

        # parallel loop
        gene_types = Parallel(n_jobs=num_cores)(delayed(check_gene)(gene) for gene in genes)

        adata.var['proteincoding_or_mirna_gene'] = gene_types


        # Save the modified AnnData object back to a new h5ad file
        adata.write(os.path.join(root_dir, dir, "local.with_filter.h5ad"))

        # Write as loom file
        adata.write_loom(os.path.join(root_dir, dir, "local.with_filter.loom"))
    

        # set_trace()


# Load the h5ad file
# adata = sc.read_h5ad("pilocytic_astrocytoma/local.h5ad")


