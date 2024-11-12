from geneformer import TranscriptomeTokenizer



## tk = TranscriptomeTokenizer(nproc=64)

# tokenize the 6M cell cancer corpus
# tk = TranscriptomeTokenizer(custom_attr_name_dict={"cell_type": "cell_type", "disease": "disease", "tissue":"tissue"}, nproc=64)
# tk.tokenize_data(f"/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/loom_files", "/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus", "cancer")

# tokenize the HNSCC dataset
tk = TranscriptomeTokenizer(custom_attr_name_dict={"cell_type": "cell_type", "origin": "origin", "P_mid":"P_mid"}, nproc=64)
tk.tokenize_data(f"/grand/GeomicVar/rbr7/Geneformer/HNSCC_data", "/grand/GeomicVar/rbr7/Geneformer/HNSCC_data/HNSCC", "hsncc")
