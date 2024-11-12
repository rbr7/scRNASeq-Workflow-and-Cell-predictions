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

# root_dir = '/gpfs/alpine/med112/proj-shared/tnandi/Geneformer/cancer_cellxgene_datasets/'
root_dir = '/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/'
dirs = os.listdir(root_dir)
# dirs = ['glioblastoma', 'follicular_lymphoma', 'hepatoblastoma']
# dirs = ['hepatoblastoma']
# dirs = ['human_breast_all_cells']
# dirs = ['pilocytic_astrocytoma']
# dirs = ['B-cell_non-Hodgkin_lymphoma', 'pilocytic_astrocytoma', 'follicular_lymphoma', 'gastric']
# dirs = ['follicular_lymphoma'] # follicular_lymphoma has issues (fixed) with KeyError: 'Attribute name cannot contain slash (/)'
# dirs = ['gastric'] # has issues (fixed by moving to Python 3.7 from 3.8 on Polaris): UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 2: ordinal not in range(128)
# dirs = ['lung_cancer']
# dirs = ['ovarian_cancer'] # has this error (fixed by moving to Python 3.7 on Polaris) ValueError: Argument must be a list, tuple, numpy matrix, numpy ndarray or sparse matrix.
# dirs = ['glioblastoma'] # has issues (fixed by moving to Python 3.7 from 3.8 on Polaris): UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 2: ordinal not in range(128)
dirs = ['HLCA'] # has issues with memory allotment

with open('./ensembl_genes/geneid_proteincoding_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_coding_mirna = f.read().splitlines()

with open('./ensembl_genes/geneid_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_mirna = f.read().splitlines()


# loop through all datasets
for dir in dirs:
    if os.path.isdir(dir): # exclude all files listed by os.listdir
        print("Entering directory ", dir)
        adata = sc.read_h5ad(os.path.join(root_dir, dir, "local.h5ad"))
        # Rename TCR/BCR for follicular_lymphoma to avoid this error: KeyError: 'Attribute name cannot contain slash (/)' 
        # adata.obs.rename(columns={'TCR/BCR':'TCRBCR'}, inplace=True)
        # total_read_count = adata.obs['nCount_RNA']
        total_read_count = np.sum(adata.raw.X.A[:, :], axis=1)   # alternative way to access the raw counts
        # set_trace()
        # Calculate total mitochondrial read count for each cell
        # pct_mt_read_count = adata.obs['percent.mt']
        mito_genes = adata.var['feature_name'][adata.var['feature_name'].str.startswith('MIR')]
        mito_ids = np.array(mito_genes.index)
        mito_indices = [adata.raw.var_names.get_loc(id) for id in mito_ids]
        mito_counts = np.sum(adata.raw.X[:, mito_indices], axis=1).flatten().tolist()[0] # converts to a list
        pct_mito = (mito_counts / total_read_count) * 100

        # for df in [adata.obs, adata.var]:
        #    for col in df.columns:
        #        if df[col].dtype == object:
        #            print(col, df[col].str.encode(encoding='ascii', errors='ignore').str.contains(r'[^\x00-\x7F]'))
        # set_trace()
        # mito_counts = np.sum(adata[:, mito_ids].X, axis=1)
        # mito_counts = np.sum(adata.raw.X[:, mito_ids], axis=1) # have to use the raw counts
        # set_trace()
        # cell_0_expr = adata.X.A[0, :] # Gets the dense representation of the expression level for all genes for cell 0
        # Get the mean and std of read count across the samples
        mean_read_count         = np.mean(total_read_count) 
        std_read_count          = np.std(total_read_count) 
        mean_mt_read_count      = np.mean(mito_counts)
        std_mt_read_count       = np.std(mito_counts)
        # set_trace()

        # Create a new column "filter_pass" and set 1 for passed samples 
        adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) <= 3 * std_read_count) & (np.abs(mito_counts - mean_mt_read_count) <= 3 * std_mt_read_count) & (total_read_count > 7), 1, 0)

        # Loop over all genes and create a new column called "gene_type" in adata.var
        gene_types = []
        genes =  adata.var.index
        print("looping over ", len(genes), " genes")


        # function to check each gene [sets 1 to coding or miRNA genes, and 0 for others]
        def check_gene(gene):
            print(gene)
            if gene in gene_ids_coding_mirna:
                return 1
            else:
                return 0

        # parallel loop
        gene_types = Parallel(n_jobs=num_cores)(delayed(check_gene)(gene) for gene in genes)

        adata.var['proteincoding_or_mirna_gene'] = gene_types # accordingly change obtain_nonzero_median_digests.py [coding_miRNA_loc = np.where((data.ra.proteincoding_or_mirna_gene == 1))[0] instead of np.where((data.ra.gene_type == "protein_coding") | (data.ra.gene_type == "miRNA"))[0]]

        # Save the modified AnnData object back to a new h5ad file
        adata.write(os.path.join(root_dir, dir, "local.with_filter.h5ad"))

        # Write as loom file

        # adata.write_loom(os.path.join(root_dir, dir, "local.with_filter.loom"))
        # to avoid UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 2: ordinal not in range(128)
        #loom_file_path = os.path.join(root_dir, dir, "local.with_filter.loom")
        # Specify the encoding when writing the Loom file
        #with open(loom_file_path, "w", encoding="utf-8") as loom_file:
        #    adata.write_loom(loom_file)

        # set_trace()


# Load the h5ad file
# adata = sc.read_h5ad("pilocytic_astrocytoma/local.h5ad")


