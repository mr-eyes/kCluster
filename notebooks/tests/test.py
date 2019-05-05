"""
Simulate the virtualQs manually
"""

from Bio import SeqIO
import json
from itertools import combinations

# CONSTANTS
fasta_file_path = "seq.fa"
KMER_SIZE = 25
MIN_Q = 1
MAX_Q = 25
STEP_Q = 1


def read_seqs(fasta_file):
    seqs = {}
    records = SeqIO.parse(fasta_file, "fasta")
    for seq in records:
        seqs[seq.id] = str(seq.seq)

    return seqs


def seq_to_kmers(seq, Q):
    kmers = []

    for i in range(len(seq) - Q + 1):
        #print(seq[i:i+Q])
        kmers.append(seq[i:i+Q])

    return kmers


def common_kmers(kmers1, kmers2):
    assert(len(kmers1) == len(kmers2))
    count = 0
    for i in range(len(kmers1)):
        if kmers1[i] == kmers2[i]:
            count += 1

    return(count)


seqs = read_seqs(fasta_file_path)

virtualQs_kmers = {}

for seq_id, seq in seqs.items():
    virtualQs_kmers[seq_id] = dict()
    for Q in range(MIN_Q, MAX_Q+1, STEP_Q):
        virtualQs_kmers[seq_id][Q] = seq_to_kmers(seq, Q)


# Empty Dictionary to hold the number of shared kmers
virtualQs_common = {}

all_seqs_ids = combinations(seqs.keys(), 2)

for seq_comb in all_seqs_ids:
    key = ":".join(list(seq_comb))

    seq_id_1 = seq_comb[0]
    seq_id_2 = seq_comb[1]
    virtualQs_common[key] = dict()

    for Q in range(MIN_Q, MAX_Q+1, STEP_Q):
        kmers1 = virtualQs_kmers[seq_id_1][Q]
        kmers2 = virtualQs_kmers[seq_id_2][Q]
        virtualQs_common[key][f"Q{Q:02d}"] = common_kmers(kmers1, kmers2)


print(json.dumps(virtualQs_common, sort_keys=True, indent=4, separators=(',', ': ')))