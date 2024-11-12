import anndata
from pdb import set_trace
# code to split the HLCA h5ad file as filtering the single file leads to memory issues

# Load the original h5ad file
original_data = anndata.read_h5ad("HLCA/local.h5ad")


n_cells_subset1 = 1000000

# Divide the data (e.g., split cells into two subsets)
subset1 = original_data[:n_cells_subset1]
subset2 = original_data[n_cells_subset1:]

set_trace()

# Save the subsets as separate h5ad files
subset1.write_h5ad("HLCA/local.subset1.h5ad")
subset2.write_h5ad("HLCA/local.subset2.h5ad")

