# This script randomly shuffles the rows of the combined dataset to reduce training bias

from datasets import load_dataset, load_from_disk, concatenate_datasets
import os
import random
from pdb import set_trace

dataset = "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset"

print("Loading ", dataset)
dataset = load_from_disk(dataset)
# set_trace()

num_rows = len(dataset)

dataset = dataset.shuffle(seed=42)
dataset = dataset.flatten_indices()

#set_trace()

dataset.save_to_disk("/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined_6M_shuffled.dataset")
