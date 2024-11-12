import anndata
from pdb import set_trace
import pandas as df
import biomart
import scanpy as sc
import h5py
import numpy as np
import loompy as lp
import pybiomart
import pandas as pd

# adata = anndata.read_loom("./glioblastoma_2/local.with_filter.glioblastoma_2.loom")
# adata = anndata.read_loom("/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/glioblastoma/local.with_filter.glioblastoma.loom")
adata = anndata.read_loom("/lus/grand/projects/GeomicVar/rbr7/Geneformer/HNSCC_data/test.loom")

# Map gene names to Ensembl IDs
symbols = adata.var_names

dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

# Fetch all Ensembl gene IDs and gene names
attributes = ['ensembl_gene_id', 'external_gene_name']

print("querying database to map gene symbols and ensembl ids")
ensembl_df = dataset.query(attributes=attributes)

print(ensembl_df.columns)

# Filter the DataFrame to include only rows where the external_gene_name is in the 'symbols' list
print("filtering df to keep mapping of only the genes symbols present in this database")
filtered_df = ensembl_df[ensembl_df['Gene name'].isin(symbols)]

# set the 'Gene name' the column as index
filtered_genename_index_df = filtered_df.set_index('Gene name')
# a single gene symbol may map to multiple ensembl IDs
# for now, keep the first entry for such instances and remove others
# unique_df = filtered_df.drop_duplicates(subset='Gene name')



# Extract the list of 'Gene name' values from unique_df
# gene_names_to_keep = unique_df['Gene name'].tolist()

# Create a dictionary mapping gene names to column indices
# gene_to_index = {gene_name: index for index, gene_name in enumerate(adata.var_names)}

# Create a list of integer indices for the selected columns
# selected_indices = [gene_to_index[gene_name] for gene_name in gene_names_to_keep]



# create a df that'll be converted to adata.var attribute (that is a pandas df, and in this case will only consist of the gene symbol as the index and the ensembl ID as the column)
common_elements = set(filtered_df['Gene name']) & set(adata.var_names)
common_elements_list = list(common_elements)
gene_names_to_keep = common_elements_list
gene_to_index = {gene_name: index for index, gene_name in enumerate(adata.var_names)}
selected_indices = [gene_to_index[gene_name] for gene_name in gene_names_to_keep]

# Create a subset of the raw matrix using the selected_indices
bdata = adata.layers['counts'][:, selected_indices]
dense_bdata = bdata.toarray()

# create an anndata object from dense_bdata
cdata = anndata.AnnData(X=dense_bdata)

# create the obs attributes
cdata.obs['n_counts'] = dense_bdata.sum(axis=1)
cdata.obs['P_mid'] = adata.obs['P_Mid'].values
cdata.obs['cell_type'] = adata.obs['cell_type'].values
cdata.obs['origin'] = adata.obs['origin'].values
cdata.obs['patientID'] = adata.obs['patientID'].values
cdata.obs['sampleID'] = adata.obs['sampleID'].values
cdata.obs['percent.mt'] = adata.obs['percent.mt'].values

# Create a DataFrame for the var attribute with gene symbols as index and gene ensembl ID as the columns

filtered_df_sorted = filtered_df.sort_values(by='Gene name', ascending=True)
filtered_df_sorted_2 = filtered_df_sorted.sort_values(by='Gene stable ID', ascending=True)
filtered_df_no_duplicates = filtered_df_sorted_2.drop_duplicates(subset='Gene name', keep='first')
filtered_df_arranged = filtered_df_no_duplicates.set_index('Gene name').loc[gene_names_to_keep] # re-order to maintain consistency with the adata matrix subsetting
# symbols = adata.var_names[selected_indices]
# gene_ensembl_ids = filtered_df_no_duplicates[filtered_df_no_duplicates['Gene name'].isin(gene_names_to_keep)]['Gene stable ID']



# var_df = pd.DataFrame({'ensembl_id': filtered_df_arranged['Gene stable ID']}, index=filtered_df_arranged.index)
var_df = pd.DataFrame({'ensembl_id': filtered_df_arranged['Gene stable ID']})

# Add the var attribute to the AnnData object
cdata.var = var_df


cdata.write('hnscc_mod_B.h5ad')
cdata.write_loom('hnscc_mod_B.loom')


# gene_df = pd.DataFrame(index=adata.var_names, columns=filtered_genename_index_df.loc[list(common_elements)])
# gene_df['ensembl_id'] = adata.var_names
# adata.var = gene_df


set_trace()


