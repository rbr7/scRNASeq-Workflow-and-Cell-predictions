import biomart 
import scanpy as sc
import h5py
import numpy as np
import loompy as lp
from pdb import set_trace
import pybiomart


data = lp.connect("/lus/grand/projects/GeomicVar/tarak/Geneformer/HNSCC_data/test.loom")

print("Data shape:", data.shape)
print("Available layers:", data.layers.keys())
print("Row attributes:", data.ra.keys())
print("Column attributes:", data.ca.keys())

print(data.attrs.keys())

# Map gene names to Ensembl IDs
symbols = data.ra['Gene']

dataset = pybiomart.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

# Fetch all Ensembl gene IDs and gene names
attributes = ['ensembl_gene_id', 'external_gene_name']
ensembl_df = dataset.query(attributes=attributes)

print(ensembl_df.columns)

# Filter the DataFrame to include only rows where the external_gene_name is in the 'symbols' list
filtered_df = ensembl_df[ensembl_df['Gene name'].isin(symbols)]

# a single gene symbol may map to multiple ensembl IDs
# for now, keep the first entry for such instances and remove others
unique_df = filtered_df.drop_duplicates(subset='Gene name')


# Extract the list of 'Gene name' values from unique_df
gene_names_to_keep = unique_df['Gene name'].tolist()

# Create a dictionary mapping gene names to column indices
gene_to_index = {gene_name: index for index, gene_name in enumerate(data.ra['Gene'])}

# Create a list of integer indices for the selected columns
selected_indices = [gene_to_index[gene_name] for gene_name in gene_names_to_keep]

set_trace()

# Sort the selected_indices in ascending order
selected_indices.sort()

# Use the sorted selected_indices to filter the Loom dataset
print("Creating filtered dataset")
filtered_data = data[:, selected_indices]

new_loom_file_path = "/lus/grand/projects/GeomicVar/tarak/Geneformer/HNSCC_data/test_mod.loom"

# Create a dictionary to map Ensembl IDs to gene symbols from ensembl_df
# ensembl_to_gene_symbol = dict(zip(ensembl_df['ensembl_gene_id'], ensembl_df['external_gene_name']))
# ensembl_to_gene_symbol = dict(zip(filtered_df['Gene stable ID'], filtered_df['Gene name']))
gene_symbol_to_ensembl = dict(zip(filtered_df['Gene name'], filtered_df['Gene stable ID']))


# Rename the column attribute 'nCount_RNA' to 'n_counts'
# data.ca['n_counts'] = data.ca['nCount_RNA']
# del data.ca['nCount_RNA']

print("Creating the updated loom file")
# Create a new Loom file with the filtered_data, specifying layers, row_attrs, and col_attrs
with lp.create(new_loom_file_path, layers={'': filtered_data, 'counts': filtered_data}, row_attrs={'ensembl_id': [gene_symbol_to_ensembl.get(gene_symbol, '') for gene_symbol in data.ra['Gene']]} if "Gene" in data.ra else {}, col_attrs={attr: data.ca[attr] for attr in data.ca}) as new_data:

    # Set the row attribute 'ensembl_id' based on the mapping
    new_data.ra['ensembl_id'] = [gene_symbol_to_ensembl.get(gene_symbol, '') for gene_symbol in data.ra['Gene']]
    
    # Copy the remaining row attributes from the original data to the new data
    for attr_name in data.ra.keys():
        if attr_name != 'Gene':
            new_data.ra[attr_name] = data.ra[attr_name]

    # Copy the modified column attributes from the original data to the new data
    for attr_name, attr_values in data.ca.items():
        new_data.ca[attr_name] = attr_values

# Confirm that the new Loom file has been created with the desired attributes and data
print("New Loom file created:", new_loom_file_path)










set_trace()





