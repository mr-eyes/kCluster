import sys
from Bio import SeqIO

if len(sys.argv) < 3:
    exit("Please pass [1]: file_path [2] RNA_type\npython filter_multifasta.py <fasta_file> <rna_type(s)>")

rna_types = []
if len(sys.argv) == 3:
    rna_types = sys.argv[2]

elif len(sys.argv) > 3:
    rna_types = sys.argv[2:]
    
if "," in rna_types:
    rna_types = rna_types.split(",")
else:
    rna_types = [rna_types]

input_file = sys.argv[1]

output_file = "-".join(rna_types) + "_" + input_file

fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

with open(output_file,"w") as handle:
    for fasta in fasta_sequences:
        description = fasta.description        
        for _type in rna_types:
            if _type in str(description):
                SeqIO.write(fasta, handle, "fasta")
