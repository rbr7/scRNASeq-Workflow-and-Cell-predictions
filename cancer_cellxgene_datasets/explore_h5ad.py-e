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
print("Number of cores: ", num_cores)

# root_dir = '/gpfs/alpine/med112/proj-shared/tnandi/Geneformer/cancer_cellxgene_datasets/'
root_dir = '/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/'
dirs = os.listdir(root_dir)
dirs = ['pilocytic_astrocytoma', 'hepatoblastoma', 'B-cell_non-Hodgkin_lymphoma', 'follicular_lymphoma', 'gastric', 'human_breast_all_cells', 'ovarian_cancer', 'glioblastoma', 'lung_cancer', 'HLCA'] # HLCA has issues with memory allotment

with open('./ensembl_genes/geneid_proteincoding_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_coding_mirna = f.read().splitlines()

with open('./ensembl_genes/geneid_mirna.csv') as f: # The genes listed in this file are either protein coding or miRNA
    gene_ids_mirna = f.read().splitlines()


# loop through all datasets
for dir in dirs:
    if os.path.isdir(dir): # exclude all files listed by os.listdir
        print("-----------------------------------------------------------------------")
        print("Entering directory ", dir)
        adata = sc.read_h5ad(os.path.join(root_dir, dir, "local.h5ad"))
        # Rename TCR/BCR for follicular_lymphoma to avoid this error: KeyError: 'Attribute name cannot contain slash (/)' 
        # adata.obs.rename(columns={'TCR/BCR':'TCRBCR'}, inplace=True)
        # total_read_count = adata.obs['nCount_RNA']
        print(adata)
        set_trace()





