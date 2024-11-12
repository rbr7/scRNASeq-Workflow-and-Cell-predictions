import loompy
import scanpy as sc
import pandas
import numpy
import scipy
from pdb import set_trace

adata = sc.read_loom('test.loom')

set_trace()

# Move embeddings info to the right place and right format
x = pandas.Series.to_numpy(adata.obs['_X'])
y = pandas.Series.to_numpy(adata.obs['_Y'])
xy = numpy.stack((x,y)).transpose().reshape(-1,2)
adata.obsm['X_test'] = xy

# Only include necessary metadata:
adata.obs['Clusters'] = pandas.Categorical(adata.obs['Clusters'])
adata.obs = adata.obs[{'Clusters','Age','Sex'}]

# Change the matrix format
adata.X = scipy.sparse.csc_matrix(adata.X)

# Make variable and observation names unique
adata.var_names_make_unique()
adata.obs_names_make_unique()

# Write h5ad file
data.write('filename.h5ad')
