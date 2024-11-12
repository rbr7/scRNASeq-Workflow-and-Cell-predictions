from datasets import load_dataset, load_from_disk, concatenate_datasets
import os

dataset_dirs = os.listdir("/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/")
dataset_paths = [os.path.join("/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/", d) for d in dataset_dirs if d.endswith(".dataset")]

datasets = []
for dataset_path in dataset_paths:
    print("Loading ", dataset_path)
    dataset = load_from_disk(dataset_path)
    datasets.append(dataset)

combined = concatenate_datasets(datasets) 

combined.save_to_disk("/grand/GeomicVar/rbr7/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset")
