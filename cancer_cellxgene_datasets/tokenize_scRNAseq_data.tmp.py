from geneformer import TranscriptomeTokenizer

cancer_types = ["B-cell_non-Hodgkin_lymphoma", "follicular_lymphoma", "hepatoblastoma", "", "", "", "", "", "", ""]


for cancer_type in cancer_types:
    print("Tokenizing ", cancer_type)
    tk = TranscriptomeTokenizer(custom_attr_name_dict={"cell_type": "cell_type", "disease": "disease", "tissue":"tissue"}, nproc=64)
    # tk = TranscriptomeTokenizer(nproc=64)
    tk.tokenize_data(f"/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/{cancer_type}/", "/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus", cancer_type)
