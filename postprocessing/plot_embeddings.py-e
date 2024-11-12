from geneformer import EmbExtractor

# initiate EmbExtractor
embex = EmbExtractor(model_type="Pretrained",
                     num_classes=0, # 0 for pretrained model, and non-zero for fine tuned gene or cell classifier model
                     # filter_data={"cell_type":["fibroblast","T cell"]}, # embeddings will be extracted only for these cells [e.g., for certain cell types]
                     # filter_data={"disease":["normal"], "tissue":["lung", "breast", "brain", "bone marrow"]}, # embeddings will be extracted only for these cells [e.g., for certain cell types]
                     # filter_data={"tissue":["liver"], "disease":["normal", "blastoma"]}, # embeddings will be extracted only for these cells [e.g., for certain cell types]
                     # filter_data={"disease":["normal", "glioblastoma", "lung adenocarcinoma", "malignant ovarian serous tumor", "pulmonary fibrosis", "squamous cell lung carcinoma", "breast cancer", "follicular lymphoma"]}, 
                     # filter_data={"tissue":["lung"], "disease":["normal", "lung adenocarcinoma"]},
                     # filter_data={"tissue": ["brain"], "disease":["normal", "glioblastoma"]},
                     filter_data={"tissue": ["breast"]},
                     # filter_data={"tissue": ["liver", "brain"], "disease":["normal", "blastoma", "glioblastoma"]},
                     max_ncells=20000, # Maximum number of (randomly chosen) cells to extract embeddings from
                     emb_mode="cell",
                     emb_layer=-1,  # [0, -1] for last layer or 2nd to last layer
                     # emb_label=["disease","cell_type"], # associated with the labels for each cell in the csv file
                     emb_label=["cell_type"], # associated with the labels for each cell in the csv file
                     labels_to_plot=["cell_type"],  # used for plotting and needs to be a subset of 'emb_label'
                     forward_batch_size=50,
                     nproc=16)

# extracts embedding from input data
# example dataset: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/example_input_files/cell_classification/disease_classification/human_dcm_hcm_nf.dataset

'''
embs = embex.extract_embs("../fine_tuned_models/geneformer-6L-30M_CellClassifier_cardiomyopathies_220224",
                          "path/to/input_data/",
                          "path/to/output_directory/",
                          "output_prefix")
'''
# embs = embex.extract_embs("/lus/grand/projects/GeomicVar/tarak/Geneformer/cancer_train/models/230910_102418_geneformer_30M_L6_emb256_SL2048_E10_B12_LR0.001_LSlinear_WU10000_Oadamw_DS8/checkpoint-439712/", "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined_6M_shuffled.dataset/", "./output_trained","breast_cell_type")
embs = embex.extract_embs("/lus/grand/projects/GeomicVar/tarak/Geneformer/cancer_train/models/230904_161950_geneformer_30M_L6_emb256_SL2048_E3_B12_LR0.001_LSlinear_WU10000_Oadamw_DS4/checkpoint-10", "/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined_6M_shuffled.dataset/", "./output_untrained","breast_cell_type")


# plot UMAP of cell embeddings
# note: scanpy umap necessarily saves figs to figures directory
'''
embex.plot_embs(embs=embs,
                plot_style="umap",
                output_directory="path/to/output_directory/",
                output_prefix="emb_plot")
'''
embex.plot_embs(embs=embs,
                plot_style="umap",
                # output_directory="./output",
                output_directory="./output_untrained",
                output_prefix="untrained",
                kwargs_dict={"palette": "Set1", "size": 200})
