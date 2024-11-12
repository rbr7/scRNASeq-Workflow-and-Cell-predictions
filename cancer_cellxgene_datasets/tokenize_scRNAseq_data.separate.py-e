from geneformer import TranscriptomeTokenizer

# cancer_types = ['B-cell_non-Hodgkin_lymphoma', 'follicular_lymphoma' , 'hepatoblastoma', 'lung_adenocarcinoma', 'pilocytic_astrocytoma', 'ovarian_cancer_1', 'ovarian_cancer_2', 'ovarian_cancer_3', 'ovarian_cancer_4', 'lung_cancer_1']

# cancer_types = ['lung_cancer_2', 'lung_cancer_3', 'lung_cancer_4', 'lung_cancer_5', 'glioblastoma_1', 'glioblastoma_2'] #

# cancer_types = ['glioblastoma_3', 'glioblastoma_4'] # has issues at 
cancer_types = ['glioblastoma_3'] 


for cancer_type in cancer_types:
    # print("Tokenizing ", cancer_type)
    tk = TranscriptomeTokenizer(custom_attr_name_dict={"cell_type": "cell_type", "disease": "disease", "tissue":"tissue"}, nproc=64)
    # tk = TranscriptomeTokenizer(nproc=64)
    # tk.tokenize_data("loom_data_directory", "output_directory", "output_prefix")
    tk.tokenize_data(f"/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/{cancer_type}/", "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus", cancer_type)
