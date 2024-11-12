import os
import numpy as np
import scanpy as sc
import requests
import yaml
from bs4 import BeautifulSoup
from pdb import set_trace
from joblib import Parallel, delayed
import multiprocessing

# Number of cores 
num_cores = multiprocessing.cpu_count()
print("NUmber of cores: ", num_cores)

root_dir = '/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/'
dirs = os.listdir(root_dir)
# dirs = ['HLCA'] # has issues with memory allotment


with open('./ensembl_genes/geneid_proteincoding_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_coding_mirna = f.read().splitlines()

with open('./ensembl_genes/geneid_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_mirna = f.read().splitlines()


# loop through all datasets
for dir in dirs:
    if os.path.isdir(dir): # exclude all files listed by os.listdir
        print("-----------------------Entering directory ---------------------: ", dir)
        adata = sc.read_h5ad(os.path.join(root_dir, dir, "local.h5ad"))
        # Rename TCR/BCR for follicular_lymphoma to avoid this error: KeyError: 'Attribute name cannot contain slash (/)' 
        if dir == 'follicular_lymphoma':
            adata.obs.rename(columns={'TCR/BCR':'TCRBCR'}, inplace=True)
        # total_read_count = adata.obs['nCount_RNA']
        total_read_count = np.sum(adata.raw.X.A[:, :], axis=1)   # alternative way to access the raw counts
        mito_genes = adata.var['feature_name'][adata.var['feature_name'].str.startswith('MIR')]
        print("Mitochondrial genes in this dataset: ", list(mito_genes), len(mito_genes))
        mito_ids = np.array(mito_genes.index)
        mito_indices = [adata.raw.var_names.get_loc(id) for id in mito_ids]
        mito_counts = np.sum(adata.raw.X[:, mito_indices], axis=1).flatten().tolist()[0] # converts to a list
        min_read_count         = np.min(total_read_count) 
        mean_read_count         = np.mean(total_read_count) 
        std_read_count          = np.std(total_read_count) 
        mean_mt_read_count      = np.mean(mito_counts)
        std_mt_read_count       = np.std(mito_counts)

        # Create a gene attribute "ensembl_id" [Ensembl ID for each gene]
        adata.var["ensembl_id"] = np.array(adata.var.index)
        # Create a cell attribute "n_counts" [total read counts in that cell]
        adata.obs["n_counts"] = total_read_count
        # Create a cell attribute "filter_pass" and set 1 for passed samples 
        adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) >= 3 * std_read_count) | (np.abs(mito_counts - mean_mt_read_count) >= 3 * std_mt_read_count) | (total_read_count < 7), 0, 1)

        gene_types = []
        genes =  adata.var.index
        print("Number of cells: ", adata.raw.X.shape[0])
        print("Number of cells after filtering: ", adata.obs['filter_pass'].value_counts()[1])
        print("Number of genes: ", len(genes))
        print("Minimum number of read in any cell: ", min_read_count)
        print("Mean number of mitochondrial genes per cell: ", mean_mt_read_count)
        #print("Max number of mitochondrial genes in any cell: ", np.max(mito_counts))
        # set_trace()


        # function to check each gene [sets 1 to coding or miRNA genes, and 0 for others]
        def check_gene(gene):
            print(gene)
            if gene in gene_ids_coding_mirna:
                return 1
            else:
                return 0

        # parallel loop
        # gene_types = Parallel(n_jobs=num_cores)(delayed(check_gene)(gene) for gene in genes)

        # adata.var['proteincoding_or_mirna_gene'] = gene_types # accordingly change obtain_nonzero_median_digests.py [coding_miRNA_loc = np.where((data.ra.proteincoding_or_mirna_gene == 1))[0] instead of np.where((data.ra.gene_type == "protein_coding") | (data.ra.gene_type == "miRNA"))[0]]

        # Save the modified AnnData object back to a new h5ad file
        # adata.write(os.path.join(root_dir, dir, "local.with_filter.h5ad"))
        adata.write_loom(os.path.join(root_dir, dir, f"local.with_filter.{dir}.loom"))





