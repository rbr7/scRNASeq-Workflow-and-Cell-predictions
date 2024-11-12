import anndata 
import numpy as np
import os
from pdb import set_trace


root_dir = '/grand/GeomicVar/tarak/Geneformer/HNSCC_data'

adata = anndata.read_loom('hnscc_mod_B.loom')
total_read_count = adata.obs['n_counts'].values
mito_genes = adata.var['Gene name'][adata.var['Gene name'].str.startswith('MT-')]
mito_mask = (adata.var['Gene name'].index).isin(mito_genes.index)
print("Mitochondrial genes in this dataset: ", list(mito_genes), len(mito_genes))
mito_counts = np.sum(adata.X[:, mito_mask], axis=1).A1 # flatten to an 1D array


# Get the mean and std of read count across the samples
min_read_count         = np.min(total_read_count)
mean_read_count         = np.mean(total_read_count)
std_read_count          = np.std(total_read_count)
mean_mt_read_count      = np.mean(mito_counts)
std_mt_read_count       = np.std(mito_counts)
# set_trace()
# Create a gene attribute "ensembl_id" [Ensembl ID for each gene]
# adata.var["ensembl_id"] = np.array(adata.var.index)
# Create a cell attribute "n_counts" [total read counts in that cell]
# adata.obs["n_counts"] = total_read_count
# Create a cell attribute "filter_pass" and set 1 for passed samples
if mito_counts.sum() == 0: # lung cancer dataset have mitochondrial genes already removed
    adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) >= 3 * std_read_count) | (total_read_count < 7), 0, 1)
else:
    adata.obs["filter_pass"] = np.where((np.abs(total_read_count - mean_read_count) >= 3 * std_read_count) | (np.abs(mito_counts - mean_mt_read_count) >= 3 * std_mt_read_count) | (total_read_count < 7), 0, 1)
# Rename the existing 'cell_type' column to 'cell_type_ID'
adata.obs.rename(columns={'cell_type': 'cell_type_ID'}, inplace=True)

# Create a new 'cell_type' column based on the values in 'cell_type_ID'
cell_type_mapping = {
    3: 'CAFs',
    8: 'B cells',
    1: 'Epithelial',
    10: 'TAMs',
    5: 'Endothelial',
    7: 'T cells',
    4: 'Myofib',
    2: 'Salivary',
    6: 'NK',
    9: 'plasma cells',
    11: 'pDC',
    12: 'mature DC',
    13: 'Mast cells',
}

adata.obs['cell_type'] = adata.obs['cell_type_ID'].map(cell_type_mapping)


genes =  adata.var['Gene name']
# set_trace()
print("Number of cells: ", adata.X.shape[0])
print("Number of cells after filtering: ", adata.obs['filter_pass'].value_counts()[1])
print("Number of genes: ", len(genes))
print("Minimum number of reads in any cell: ", min_read_count)
print("Mean number of mitochondrial genes per cell: ", mean_mt_read_count)

# set_trace()
# Save the modified AnnData object back to a new loom file
adata.write_loom(os.path.join(root_dir, f"hnscc_mod_B_filtered.loom"))
# Save the modified AnnData object back to a new h5ad file
adata.write(os.path.join(root_dir, f"hnscc_mod_B_filtered.h5ad"))



set_trace()
