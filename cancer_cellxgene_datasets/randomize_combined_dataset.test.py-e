# This script randomly shuffles the rows of the combined dataset to reduce training bias

from datasets import load_dataset, load_from_disk, concatenate_datasets
import os
import random
from pdb import set_trace
import multiprocessing
from functools import partial
import numpy as np

# load data and add index
dataset = "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset"
print("Loading ", dataset)
dataset = load_from_disk(dataset)
num_cells = len(dataset)
cell_index = np.arange(num_cells)
dataset = dataset.add_column('index', cell_index)
random_id = [random.random() for i in range(num_cells)]
dataset = dataset.add_column('random_id', random_id)
dataset = dataset.sort('random_id')
dataset = dataset.remove_columns('random_id')
#dataset = dataset.sort('index')
#set_trace()

#num_rows = len(dataset)
#randomized_indices = random.sample(range(num_rows), num_rows)
#randomized_dataset = dataset.select(randomized_indices)
#randomized_dataset.save_to_disk("/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined_randomized_cancer_3M.dataset")

def parallel_sort(partition):
  # Sort each partition by index  
  return partition.sort('index')

if __name__ == '__main__':
  # Split into partitions
  num_partitions = 10
  partitions = np.array_split(dataset.to_numpy(), num_partitions) 

  pool = multiprocessing.Pool(num_partitions)
  # Parallel sort
  print("Parallel sort begins ")
  sorted_partitions = pool.map(parallel_sort, partitions)  

  # Concatenate partitions
  print("Concatenate partitions ")
  dataset = ConcatenateDataset(sorted_partitions)

  pool.close()
  pool.join()
