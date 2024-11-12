import scanpy as sc
import h5py
import numpy as np
import loompy as lp
from pdb import set_trace

# data = lp.connect("HLCA_1/local.with_filter.HLCA_1.loom")
# data = lp.connect("glioblastoma_3/local.with_filter.glioblastoma_3.loom")
# data = lp.connect("/lus/grand/projects/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/glioblastoma/local.with_filter.glioblastoma.loom")
data = lp.connect("/lus/grand/projects/GeomicVar/tarak/Geneformer/HNSCC_data/hnscc/hnscc_mod_B_filtered.loom")
# data = lp.connect("hnscc_mod_B.loom")

print("Data shape:", data.shape)
print("Available layers:", data.layers.keys())
print("Row attributes:", data.ra.keys())
print("Column attributes:", data.ca.keys())

print(data.attrs.keys())

set_trace()

