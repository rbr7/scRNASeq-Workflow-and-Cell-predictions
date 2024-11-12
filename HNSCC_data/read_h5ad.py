import scanpy as sc
from pdb import set_trace

# Load the h5ad file
# adata = sc.read_h5ad("gastric/local.with_filter.h5ad")
# adata = sc.read_h5ad("lung_adenocarcinoma/local.h5ad")
adata = sc.read_h5ad("hnscc_mod.h5ad")
set_trace()

print("basic information about the dataset")
print(adata)
print("---------------------------------------------")
print("Access the raw counts matrix")
print(adata.raw)
print("---------------------------------------------")

print("Print the variable (gene) names")
print(adata.var_names)
print("---------------------------------------------")

print("print the observation (cell) names")
print(adata.obs_names)
print("---------------------------------------------")

print("Access the cell annotations")
print(adata.obs)
print("---------------------------------------------")

print("Access the gene annotations")
print(adata.var)
print("---------------------------------------------")

set_trace()

# Perform clustering
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# Visualize the clustering results
sc.pl.umap(adata, color="leiden")

