# This script randomly shuffles the rows of the combined dataset to reduce training bias

from datasets import load_dataset, load_from_disk, concatenate_datasets
import os
import random
from pdb import set_trace
import multiprocessing
from functools import partial
import numpy as np

# load data 
# dataset = "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset"

# loop over all .dataset dirs 

directory_path = '/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/'
items = os.listdir(directory_path)
dataset_directories = [item for item in items if os.path.isdir(os.path.join(directory_path, item)) and item.endswith('.dataset')]

# dataset = "/lus/grand/projects/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/B-cell_non-Hodgkin_lymphoma.dataset"

for dataset in sorted(dataset_directories):
    print("Loading ", dataset)
    dataset = directory_path + str(dataset)
    dataset = load_from_disk(dataset)
    num_cells = len(dataset)
    print("Number of cells: ", num_cells)
    df = dataset.to_pandas()
    tissues = df['tissue']
    diseases = df['disease']
    print("tissues: ", sorted(tissues.unique()))
    print("diseases: ", sorted(diseases.unique()))


set_trace()
