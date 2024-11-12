from geneformer import TranscriptomeTokenizer

# cancer_types = ['B-cell_non-Hodgkin_lymphoma', 'follicular_lymphoma' , 'hepatoblastoma', 'lung_adenocarcinoma', 'pilocytic_astrocytoma', 'ovarian_cancer_1', 'ovarian_cancer_2', 'ovarian_cancer_3', 'ovarian_cancer_4', 'lung_cancer_1']

cancer_types = ['hnscc'] 


# compulsory entries in the loom file are:
# 1. row (gene) attribute: "ensembl_id"; Ensembl ID for each gene
# 2. col (cell) attribute: "n_counts"; total read counts in that cell
for cancer_type in cancer_types:
    # print("Tokenizing ", cancer_type)
    tk = TranscriptomeTokenizer(custom_attr_name_dict={"P_mid":"P_mid", "cell_type": "cell_type", "origin":"origin", "patientID":"patientID", "sampleID":"sampleID"}, nproc=64)
    # tk = TranscriptomeTokenizer(nproc=64)
    # tk.tokenize_data("loom_data_directory", "output_directory", "output_prefix")
    tk.tokenize_data(f"/grand/GeomicVar/rbr7/Geneformer/HNSCC_data/{cancer_type}/", "/grand/GeomicVar/rbr7/Geneformer/HNSCC_data/hnscc_corpus", cancer_type)
