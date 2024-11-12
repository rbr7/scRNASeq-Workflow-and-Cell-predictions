from datasets import load_from_disk
from pdb import set_trace
import pickle

# Write out the token length (number of genes used) for each cell
print("Loading from combined.dataset")
# cancer_corpus = load_from_disk("/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined.dataset")
cancer_corpus = load_from_disk("/grand/GeomicVar/tarak/Geneformer/cancer_cellxgene_datasets/cancer_corpus/combined_6M_shuffled.dataset")
print("Getting token length")
cancer_corpus_lengths = cancer_corpus["length"]

print("Total number of cells = ", len(cancer_corpus_lengths))

with open('cancercorpus_6M_2048_lengths.pkl', 'wb') as file:
  pickle.dump(cancer_corpus_lengths, file)

set_trace()
