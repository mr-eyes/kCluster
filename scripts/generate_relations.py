"""
generate_relations.py 

input: 
    - MAP index file
    - Names MAP file

output:
    - TSV files holds the information of shared kmers count between each two transcripts.
"""

from __future__ import division
from collections import Counter
import itertools
import re
import gc
import tqdm
import json
import sys
import os
from Bio import SeqIO
    
names_map_file = ""
map_index_file = ""
fasta_file = ""
output_file = ""
kmer_size = 25

if len(sys.argv) < 5:
    exit("Please pass <map_index_file> <names_map_file> <fasta_file> <kmer_size>")

else:
    map_index_file = sys.argv[1]
    names_map_file = sys.argv[2]
    fasta_file = sys.argv[3]
    kmer_size = int(sys.argv[4])

if len(sys.argv) == 6:
        output_file = sys.argv[5]

else:
    output_file = os.path.basename(map_index_file).split(".")[0]


print ("Reading names...")
id_to_name = {}
name_to_id = {}
with open(names_map_file) as namesMap:
    for name in namesMap:
        _id = int(re.findall(r'\t(\d+)', name)[0])
        _name = re.findall(r'(.*)\t', name)[0]
        id_to_name[_id] = _name
        name_to_id[_name] = _id

print ("Calculating Number of kmers...")
seq_to_kmersNo = {}
for seq_record in SeqIO.parse(fasta_file, "fasta"):
    seq_to_kmersNo[name_to_id[seq_record.id]] = len(seq_record) - kmer_size + 1


print ("Reading Colors & Groups...")
groups = {}
colors = []
with open(map_index_file) as MAP:
    for line in MAP:
        if ":" in line:
            colors.append(int(re.findall(r':(\d+)',line)[0]))
        elif "-" in line:
            groups[int(re.findall(r'(\d+)-', line)[0])] = [int(i) for i in re.findall(r"(\d+),", line)]

colors = Counter(colors)
print ("Done Counting Colors...")

edges = {}
nodes = set()

print ("Processing groups...")


for color, tr_ids in tqdm.tqdm(groups.items()):
    color_count = colors[color]

    if len(tr_ids) == 1:
        nodes.add(tr_ids[0])
        continue
    
    for tr in tr_ids:
        nodes.add(tr)

    for combination in itertools.combinations(tr_ids,2):
        _seq1 = combination[0]
        _seq2 = combination[1]

        if _seq1 in edges:
            if _seq2 in edges[_seq1]:
                edges[_seq1][_seq2] += color_count
            else:
                edges[_seq1][_seq2] = color_count

        else:
            edges[_seq1] = {_seq2: color_count}
        
        # Write Reverse Edge
        # if _seq2 in edges:
        #     if _seq1 in edges[_seq2]:
        #         edges[_seq2][_seq1] += color_count
        #     else:
        #         edges[_seq2][_seq1] = color_count

        # else:
        #     edges[_seq2] = {_seq1: color_count}



print ("Writing TSV file ...")

tsv = open(output_file + ".tsv", "w")
tsv.write("seq_1\tseq_2\tshared_kmers\tnorm%\n")

for _1st, info in tqdm.tqdm(edges.items()):
    for _2nd, _no_shared_kmers in info.items():
        _smallest_kmers_no = min(seq_to_kmersNo[_1st], seq_to_kmersNo[_2nd])
        _similarity = _no_shared_kmers / _smallest_kmers_no  # Normalized Weight
        _similarity *= 100

        l = "%d\t%d\t%d\t%.2f\n" % (_1st, _2nd, _no_shared_kmers, _similarity)
        tsv.write(l)

tsv.close()
