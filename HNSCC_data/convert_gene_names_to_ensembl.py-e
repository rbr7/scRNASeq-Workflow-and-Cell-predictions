from Bio import SeqIO
from biomart import BiomartServer

server = BiomartServer("http://www.ensembl.org/biomart")

dataset = server.datasets["hsapiens_gene_ensembl"]

gene_names = ["SAMD11", "NOC2L", "KLHL17", "PLEKHN1", "PERM1"]  # Replace with your gene names

query = dataset.query(attributes=["ensembl_gene_id"], filters={"external_gene_name": gene_names})

results = list(query.results)

gene_name_to_id = {}
for result in results:
    gene_name = result["external_gene_name"]
    ensembl_id = result["ensembl_gene_id"]
    gene_name_to_id[gene_name] = ensembl_id


gene_name = "SAMD11"  # Replace with the gene name you want to look up
ensembl_id = gene_name_to_id.get(gene_name, "Gene not found in Ensembl")
print(f"{gene_name}: {ensembl_id}")



