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
root_dir = '/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/'
dirs = os.listdir(root_dir)
dirs = ['pilocytic_astrocytoma']
dirs = ['B-cell_non-Hodgkin_lymphoma'] #, 'pilocytic_astrocytoma', 'follicular_lymphoma', 'gastric']
dirs = ['gastric'] # has issues (fixed by moving to Python 3.7 from 3.8 on Polaris): UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 2: ordinal not in range(128)
dirs = ['follicular_lymphoma'] # follicular_lymphoma has issues (fixed) with KeyError: 'Attribute name cannot contain slash (/)'
# dirs = ['human_breast_all_cells']
# dirs = ['ovarian_cancer'] # has this error (fixed by moving to Python 3.7 on Polaris) ValueError: Argument must be a list, tuple, numpy matrix, numpy ndarray or sparse matrix.
dirs = ['glioblastoma'] # has issues (fixed by moving to Python 3.7 from 3.8 on Polaris): UnicodeDecodeError: 'ascii' codec can't decode byte 0xc3 in position 2: ordinal not in range(128)
# dirs = ['lung_cancer']
dirs = ['HLCA_5'] # HLCA has issues with memory allocation. So it has been split into 5 h5ad files

# dirs = ['pilocytic_astrocytoma', 'hepatoblastoma', 'B-cell_non-Hodgkin_lymphoma', 'follicular_lymphoma', 'gastric', 'human_breast_all_cells', 'ovarian_cancer', 'glioblastoma', 'lung_cancer', 'HLCA'] # HLCA has issues with memory allotment (can split it into 2)

# dirs = ['ovarian_cancer'] #, 'glioblastoma', 'lung_cancer', 'HLCA']
dirs = ['breast_cancer_1', 'breast_cancer_2', 'breast_cancer_3'] # breast cancer has issues with memory allocation. So it has been split into 5 h5ad files

dirs = ["B-cell_non-Hodgkin_lymphoma", "follicular_lymphoma", "hepatoblastoma", "lung_adenocarcinoma", "pilocytic_astrocytoma", "breast_cancer_1", "breast_cancer_2", "breast_cancer_3", "glioblastoma_1", "glioblastoma_2", "glioblastoma_3", "glioblastoma_4", "HLCA_1", "HLCA_2", "HLCA_3", "HLCA_4", "HLCA_5", "HLCA_6", "HLCA_7", "HLCA_8"]

dirs = ["lung_cancer_1", "lung_cancer_2", "lung_cancer_3", "lung_cancer_4", "lung_cancer_5", "ovarian_cancer_1", "ovarian_cancer_2", "ovarian_cancer_3", "ovarian_cancer_4"]
dirs = ['hepatoblastoma']


# loop through all datasets
for dir in dirs:
    if os.path.isdir(dir): # exclude all files listed by os.listdir
        print("-----------------------Entering directory ---------------------: ", dir)
        adata = sc.read_h5ad(os.path.join(root_dir, dir, "local.h5ad"))
        print("cell and gene attributes: ", adata)
        # print("gene attributes: ", adata.var)
        # Rename TCR/BCR for follicular_lymphoma to avoid this error: KeyError: 'Attribute name cannot contain slash (/)' 
        if dir == 'follicular_lymphoma':
            adata.obs.rename(columns={'TCR/BCR':'TCRBCR'}, inplace=True)
        # total_read_count = adata.obs['nCount_RNA']
        total_read_count = np.sum(adata.raw.X.A[:, :], axis=1)   # alternative way to access the raw counts
        
        mito_genes = adata.var['feature_name'][adata.var['feature_name'].str.startswith('MT-')]
        mito_mask = adata.raw.var_names.isin(mito_genes.index)
        print("Mitochondrial genes in this dataset: ", list(mito_genes), len(mito_genes))
        mito_counts = np.sum(adata.raw.X[:, mito_mask], axis=1).A1 # flatten to an 1D array


        # Get the mean and std of read count across the samples
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
        if mito_counts.sum() == 0: # lung cancer dataset have mitochondrial genes already removed
            adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) >= 3 * std_read_count) | (total_read_count < 7), 0, 1)
        else:
            adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) >= 3 * std_read_count) | (np.abs(mito_counts - mean_mt_read_count) >= 3 * std_mt_read_count) | (total_read_count < 7), 0, 1)
        set_trace()

        gene_types = []
        genes =  adata.var.index
        # set_trace()
        print("Number of cells: ", adata.raw.X.shape[0])
        print("Number of cells after filtering: ", adata.obs['filter_pass'].value_counts()[1])
        print("Number of genes: ", len(genes))
        print("Minimum number of reads in any cell: ", min_read_count)
        print("Mean number of mitochondrial genes per cell: ", mean_mt_read_count)
        #print("Max number of mitochondrial genes in any cell: ", np.max(mito_counts))
        # set_trace()



        # Save the modified AnnData object back to a new h5ad file
        adata.write(os.path.join(root_dir, dir,  f"local.with_filter.{dir}.h5ad"))
        # adata.write_loom(os.path.join(root_dir, dir, f"local.with_filter.{dir}.loom"))





