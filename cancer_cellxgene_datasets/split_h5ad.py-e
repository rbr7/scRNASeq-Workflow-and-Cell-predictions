import anndata
import os
from pdb import set_trace

# Load the original h5ad file
cancer_type = "glioblastoma"
num_splits = 4

original_data = anndata.read_h5ad(f"{cancer_type}/local.h5ad")

total_cells = len(original_data)
cells_per_subset = total_cells // num_splits

# set_trace()

for i in range(num_splits):
    print(f"Creating split {i+1} of {num_splits} for {cancer_type}")
    # Calculate start and end indices for the current subset
    start_idx = i * cells_per_subset
    end_idx = (i + 1) * cells_per_subset if i < (num_splits - 1) else total_cells

    # Create a subset of the data
    subset = original_data[start_idx:end_idx]

    # Create a directory for the subset if it doesn't exist
    subset_folder = f"{cancer_type}_{i + 1}"
    os.makedirs(subset_folder, exist_ok=True)

    # Save the subset as 'local.h5ad' inside the subset folder
    subset.write_h5ad(os.path.join(subset_folder, "local.h5ad"))
