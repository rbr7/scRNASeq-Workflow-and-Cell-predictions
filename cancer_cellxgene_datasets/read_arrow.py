from datasets import Dataset
import pandas as pd
from pdb import set_trace

# dataset = Dataset.from_file("./human_dcm_hcm_nf.dataset/dataset.arrow")
dataset = Dataset.from_file("/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset/data-00000-of-00037.arrow")
# dataset = Dataset.from_file("/grand/RNAtoImage/rbr7_files/Geneformer/Genecorpus-30M/genecorpus_30M_2048.dataset/dataset.arrow")
set_trace()

# Access dataset metadata
print(dataset.info)

# Convert to pandas
df = dataset.to_pandas()

set_trace()

