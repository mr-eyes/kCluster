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
import gzip
import sys
import os
from Bio import SeqIO
<<<<<<< HEAD

=======
    
>>>>>>> f7c5e66c613b4e3f531813996de5a8248429b7fc
names_map_file = ""
map_index_file = ""
#fasta_file = ""
output_file = ""
#kmer_size = 25

if len(sys.argv) < 3:
    exit("Please pass <map_index_file> <names_map_file> <new:namesFile>")

else:
    map_index_file = sys.argv[1]
    names_map_file = sys.argv[2]
    # fasta_file = sys.argv[3]
    # kmer_size = int(sys.argv[4]) No need after the last edit

if len(sys.argv) == 4:
        output_file = sys.argv[3]

else:
    output_file = os.path.basename(map_index_file).split(".")[0]

"""
The following edit has been added for the m:1 relation ship of the namesMap file. {multiple transcripts : single gene}
"""

print ("Parsing namesMap...")

# This is gene to id, and id to gene

id_to_name = {}
name_to_id = {}

with open(names_map_file) as namesMap:
    for name in namesMap:
        _id = int(re.findall(r'\t(\d+)', name)[0])
        _name = re.findall(r'(.*)\t', name)[0]
        id_to_name[_id] = _name
        name_to_id[_name] = _id


print ("Reading Colors & Groups...")
groups = {}
colors = []
<<<<<<< HEAD

if ".gz" in map_index_file:
    with gzip.open(map_index_file, 'rt') as MAP:
        for line in MAP:
            if ":" in line:
                colors.append(int(line.strip().split(":")[-1]))
            elif "-" in line:
                line = line[0:-2].strip().split("-")
                group = int(line[0])
                values = list(map(int, line[1].split(",")))
                groups[group] = values


else:
    with open(map_index_file, 'r') as MAP:
        for line in MAP:
            if ":" in line:
                colors.append(int(line.strip().split(":")[-1]))
            elif "-" in line:
                line = line[0:-2].strip().split("-")
                group = int(line[0])
                values = list(map(int, line[1].split(",")))
                groups[group] = values
=======
with open(map_index_file) as MAP:
    for line in MAP:
        if ":" in line:
            colors.append(int(line.strip().split(":")[-1]))
        elif "-" in line:
            line = line[0:-2].strip().split("-")
            group = int(line[0])
            values = map(int, line[1].split(","))
            groups[group] = values

>>>>>>> f7c5e66c613b4e3f531813996de5a8248429b7fc

print ("Counting Colors...")
colors = Counter(colors)

print ("Calculating Number of kmers...")
readID_to_kmersNo = {}
for k, v in groups.items():
<<<<<<< HEAD
    colors_count = colors[k]
    for read_id in v:
        if read_id not in readID_to_kmersNo:
            readID_to_kmersNo[read_id] = colors_count
        else:
            readID_to_kmersNo[read_id] += colors_count




#print ("Collecting Garbage")
gc.collect()
=======
    _count = colors[k]
    for read_id in v:
        if read_id not in readID_to_kmersNo:
            readID_to_kmersNo[read_id] = _count
        else:
            readID_to_kmersNo[read_id] += _count
>>>>>>> f7c5e66c613b4e3f531813996de5a8248429b7fc

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

    for combination in itertools.combinations(tr_ids, 2):
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
        #
        # else:
        #     edges[_seq2] = {_seq1: color_count}


<<<<<<< HEAD
# del colors
# del groups
# gc.collect()

=======
>>>>>>> f7c5e66c613b4e3f531813996de5a8248429b7fc

print ("Writing TSV file ...")

tsv = open(output_file + ".tsv", "w")
#tsv.write("seq_1\tseq_2\tdist%\n")

for _1st, info in tqdm.tqdm(edges.items()):
    for _2nd, _no_shared_kmers in info.items():
        _smallest_kmers_no = min(readID_to_kmersNo[_1st], readID_to_kmersNo[_2nd])
        _similarity = _no_shared_kmers / _smallest_kmers_no  # Normalized Weight
        _similarity *= 100
        #_dissimilarity = 100 - _similarity

        #l = "%d\t%d\t%d\t%.2f\n" % (_1st, _2nd, _no_shared_kmers, _dissimilarity)
        l = "%d,%d,%.2f\n" % (_1st, _2nd, _similarity)
        tsv.write(l)

tsv.close()
