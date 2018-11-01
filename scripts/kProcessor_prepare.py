from Bio import SeqIO
import sys
import os
input_file = sys.argv[1]
output_file = "min" + "_" + os.path.basename(input_file)
names_file = "min" + "_" + os.path.basename(input_file) + ".names"
fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
names = open(names_file,"w")
with open(output_file,"w") as handle:
    for fasta in fasta_sequences:
        description = fasta.description.split("|")[0]
        names.write(description + "\t" + description + "\n")
        fasta.description = fasta.id = fasta.name = description
        SeqIO.write(fasta, handle, "fasta")
names.close()
