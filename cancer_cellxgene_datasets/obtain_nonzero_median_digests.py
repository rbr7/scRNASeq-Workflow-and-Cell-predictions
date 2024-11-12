import os
import numpy as np
import loompy as lp
import pandas as pd
import crick
import pickle
import math
from tqdm.notebook import tqdm
from pdb import set_trace

# input_file = "hepatoblastoma.loom"
input_file = "local.with_filter.loom"
current_database = "hepatoblastoma"

rootdir = f"./{current_database}/data/"
output_file = input_file.replace(".loom", ".gene_median_digest_dict.pickle")
outdir = rootdir.replace("/data/", "/tdigest/")

# Create output directory if it doesn't exist 
os.makedirs(outdir, exist_ok=True)

with lp.connect(f"{rootdir}{input_file}") as data:
    # define coordinates of protein-coding or miRNA genes
    # set_trace()
    # coding_miRNA_loc = np.where((data.ra.gene_type == "protein_coding") | (data.ra.gene_type == "miRNA"))[0]
    coding_miRNA_loc = np.where((data.ra.proteincoding_or_mirna_gene == 1))[0]
    coding_miRNA_genes = data.ra["feature_id"][coding_miRNA_loc]

    # initiate tdigests
    median_digests = [crick.tdigest.TDigest() for _ in range(len(coding_miRNA_loc))]

    # initiate progress meters
    progress = tqdm(total=len(coding_miRNA_loc))
    last_view_row = 0
    progress.update(0)

    for (ix, selection, view) in data.scan(items=coding_miRNA_loc, axis=0):
        print("ix: ", ix)
        # define coordinates of cells passing filter
        filter_passed_loc = np.where(view.ca.filter_pass == 1)[0]
        subview = view.view[:, filter_passed_loc]
        # normalize by total counts per cell and multiply by 10,000 to allocate bits to precision
        subview_norm_array = subview[:,:]/subview.ca.nCount_RNA*10_000
         # if integer, convert to float to prevent error with filling with nan
        if np.issubdtype(subview_norm_array.dtype, np.integer):
            subview_norm_array = subview_norm_array.astype(np.float32)
        # mask zeroes from distribution tdigest by filling with nan
        nonzero_data = np.ma.masked_equal(subview_norm_array, 0.0).filled(np.nan)
        # update tdigests
        [median_digests[i+last_view_row].update(nonzero_data[i,:]) for i in range(nonzero_data.shape[0])]
        # update progress meters
        progress.update(view.shape[0])
        last_view_row = last_view_row + view.shape[0]

median_digest_dict = dict(zip(coding_miRNA_genes, median_digests))
with open(f"{outdir}{output_file}", "wb") as fp:
    pickle.dump(median_digest_dict, fp)
