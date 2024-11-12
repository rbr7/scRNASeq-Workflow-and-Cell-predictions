import os
import shutil

# Base directory containing the subdirectories with loom files
base_directory = '/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets'
subdirs = ['B-cell_non-Hodgkin_lymphoma', 'follicular_lymphoma' , 'hepatoblastoma', 'lung_adenocarcinoma', 'pilocytic_astrocytoma', 'ovarian_cancer_1', 'ovarian_cancer_2', 'ovarian_cancer_3', 'ovarian_cancer_4', 'lung_cancer_1', 'lung_cancer_2', 'lung_cancer_3', 'lung_cancer_4', 'lung_cancer_5', 'glioblastoma_1', 'glioblastoma_2', 'glioblastoma_3', 'glioblastoma_4', 'HLCA_1', 'HLCA_2', 'HLCA_3', 'HLCA_4', 'HLCA_5', 'HLCA_6', 'HLCA_7', 'HLCA_8', 'breast_cancer_1', 'breast_cancer_2', 'breast_cancer_3', 'gastric']

print("NUmber of subdirs: ", len(subdirs))

# Create the target directory if it doesn't exist
target_dir = 'loom_files'
if not os.path.exists(target_dir):
    os.mkdir(target_dir)

# Function to copy loom files from a directory to the target directory
def copy_loom_files(source_dir, target_dir):
    for root, _, files in os.walk(source_dir):
        for file in files:
            if file.endswith('.loom'):
                source_path = os.path.join(root, file)
                target_path = os.path.join(target_dir, file)
                shutil.copy(source_path, target_path)
                print(f"Copied {source_path} to {target_path}")

# Iterate through the list of subdirectories and copy loom files
for subdir in subdirs:
    full_subdir_path = os.path.join(base_directory, subdir)
    if os.path.exists(full_subdir_path):
        copy_loom_files(full_subdir_path, target_dir)
    else:
        print(f"Subdirectory '{subdir}' not found.")

print("Loom files copying complete.")

