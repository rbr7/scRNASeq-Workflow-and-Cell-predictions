# This script randomly shuffles the rows of the combined dataset to reduce training bias

from datasets import load_dataset, load_from_disk, concatenate_datasets
import os
import random
from pdb import set_trace
import multiprocessing
from functools import partial
import numpy as np

# load data 
# dataset = "/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset"
dataset = "/lus/grand/projects/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/B-cell_non-Hodgkin_lymphoma.dataset"

print("Loading ", dataset)

dataset = load_from_disk(dataset)

num_cells = len(dataset)

df = dataset.to_pandas()

set_trace()
