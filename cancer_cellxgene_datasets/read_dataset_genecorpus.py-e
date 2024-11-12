# This script randomly shuffles the rows of the combined dataset to reduce training bias

from datasets import load_dataset, load_from_disk, concatenate_datasets
import os
import random
from pdb import set_trace
import multiprocessing
from functools import partial
import numpy as np

# load data 
dataset = "/lus/grand/projects/GeomicVar/tarak/Geneformer/Genecorpus-30M/genecorpus_30M_2048.dataset"

print("Loading ", dataset)

dataset = load_from_disk(dataset)

num_cells = len(dataset)

print("Converting to pandas df")
df = dataset.to_pandas()

set_trace()
