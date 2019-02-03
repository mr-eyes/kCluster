from __future__ import division
import sys
sys.path.insert(0, '../kProcessor')
import kProcessor as KP
from collections import Counter, defaultdict
from Bio import SeqIO
import re


# ----------------------------------------------
#                  Methods                     #
#-----------------------------------------------

def find(x):
    l = leaders[x]
    if l is not None:
        l = find(l)
        leaders[x] = l
        return l
    return x


def union(x, y):
    lx, ly = find(x), find(y)
    if lx != ly:
        leaders[lx] = ly

# ----------------------------------------------
#                  Constants                   #
#-----------------------------------------------

help_message = "<query_fasta> <ref_fasta> <kprocessor_index> <kprocessor_index_namesmap> <output_prefix> <op:kClusters_file> <op:threshold %>"

query_read_fa = ""
index_fasta = ""
index_file = ""
index_namesmap = ""
output_prefix = ""
clstr_file_path = ""
query_threshold = 0.0

if len(sys.argv) < 6:
    exit("run: python assign_query_read.py " + help_message)

else:
    query_read_fa = sys.argv[1]
    index_fasta = sys.argv[2]
    index_file = sys.argv[3].replace(".map", "")
    index_namesmap = sys.argv[4]
    output_prefix = sys.argv[5]

    if len(sys.argv) == 8:
        clstr_file_path = sys.argv[6]
        query_threshold = float(sys.argv[7])

# ----------------------------------------------
#                  Reading KCLusters           #
#-----------------------------------------------


transcript_to_cluster = {}

with open(clstr_file_path, "r") as clusters:
    next(clusters)
    for line in clusters:
        line = line.replace("\n", "")
        fields = line.split("\t")
        cluster_id = int(fields[0])
        _transcripts = fields[1].split(",")
        for tr in _transcripts:    
            transcript_to_cluster[tr] = cluster_id


# ----------------------------------------------
#                  Read NamesMap file          #
#-----------------------------------------------

id_to_name = {}
name_to_id = {}
with open(index_namesmap) as namesMap:
    for name in namesMap:
        _id = int(re.findall(r'\t(\d+)', name)[0])
        _name = re.findall(r'(.*)\t', name)[0]
        id_to_name[_id] = _name
        name_to_id[_name] = _id

# ----------------------------------------------
#                  Load Index                   #
#-----------------------------------------------

kFrameMQF = KP.kDataFrame.load(index_file, "MAP")
kSize = int(kFrameMQF.getkSize())


# ----------------------------------------------
#            Number of kmer per transcript     #
#-----------------------------------------------

seq_to_kmersNo = {}
for seq_record in SeqIO.parse(index_fasta, "fasta"):
    seq_to_kmersNo[name_to_id[seq_record.id]] = len(seq_record) - kSize + 1

# ----------------------------------------------
#                  Query all kmers             #
#-----------------------------------------------

readID_to_colors = {}
query_name_to_id = {}
query_id_to_name = {}
query_kmers_count = {}

ids = 1
for record in SeqIO.parse(query_read_fa, "fasta"):
    read_header = record.id
    read_seq = record.seq
    query_name_to_id[read_header] = ids
    query_id_to_name[ids] = read_header
    ids += 1
    all_colors = []
    query_kmers_count[read_header] = len(read_seq) - kSize + 1
    for i in range(len(read_seq) - kSize + 1):
        kmer = str(read_seq[i:i+kSize])
        colors = list(kFrameMQF.getColors(kmer))
        for col in colors:
            tr_header = id_to_name[col]
            all_colors.append(tr_header)
    
    readID_to_colors[read_header] = dict(Counter(all_colors))

# ----------------------------------------------
#                  Pairwise matrix             #
#-----------------------------------------------
pairwise_result = {}

for query_name, matches in readID_to_colors.items():
    pairwise_result[query_name] = {}
    for match, shared_count in matches.items():
        match_kmers = seq_to_kmersNo[name_to_id[match]]
        min_kmers_count = min(match_kmers, query_kmers_count)
        normalized_value = shared_count / min_kmers_count
        normalized_value *= 100
        pairwise_result[query_name][match] = [shared_count, normalized_value]


# ----------------------------------------------
#                 Writing Pairwise             #
#-----------------------------------------------

pairwise_tsv = open(output_prefix + "_pairwise_matrix.tsv" , "w")
pairwise_tsv.write("query\tref_seq\tkmers\tnorm\n")
for q_record, ref_records in pairwise_result.items():
    for ref_record, all_count in ref_records.items():
        seq1 = q_record
        seq2 = ref_record
        kCount = all_count[0]  # type: int
        norm = all_count[1]
        l = "%s\t%s\t%d\t%.2f\n" % (seq1, seq2, kCount, norm)
        pairwise_tsv.write(l)

pairwise_tsv.close()


# ----------------------------------------------
#            Decide to continue                #
#-----------------------------------------------

if len(sys.argv) != 8:
    exit()

# ----------------------------------------------
#            Clustering pairwise_result        #
#-----------------------------------------------

source = []
target = []
for q_record, ref_records in pairwise_result.items():
    for ref_record, all_count in ref_records.items():
        seq1 = q_record
        seq2 = ref_record
        norm = all_count[1]
        if norm < query_threshold:
                continue

        source.append(seq1)
        target.append(seq2)
    
leaders = defaultdict(lambda: None)
for i in range(len(source)):
    union(source[i], target[i])

components = defaultdict(set)

for x in leaders:
    components[find(x)].add(x)

clusters = {}
cluster_id = 0

for comp in components.values():
    clusters[cluster_id] = comp
    cluster_id += 1


# ----------------------------------------------
#           clusters reads to cluster ID       #
#-----------------------------------------------


new_read_to_cluster = {}
for cluster_id, cluster_reads in clusters.items():
    for read_id in cluster_reads:
        new_read_to_cluster[read_id] = cluster_id


# ----------------------------------------------
#                Recommend Best Cluster        #
#-----------------------------------------------

best_cluster_id = 0
best_cluster_norm = 0

for cluster_id, cluster_reads in clusters.items():
    query_seqs = cluster_reads.intersection(set(source))
    ref_seqs = cluster_reads - query_seqs
    total_norms = 0
    print "D: Cluster ID: ", cluster_id
    print "Cluster Reads:", str(cluster_reads)
    print "query seqs: ", str(query_seqs)
    print "ref_seqs", str(ref_seqs)
    print "-----------------------------"

    for q_seq in query_seqs:
        for rf_seq in ref_seqs:
            print "total_norms, q_seq:",q_seq," | rf_seq:",rf_seq
            total_norms += pairwise_result[q_seq][rf_seq][1]
            
    
    if total_norms > best_cluster_norm:
        best_cluster_norm = total_norms
        best_cluster_id = cluster_id

best_cluster_seqs = clusters[best_cluster_id].intersection(target)
opposite_cluster = []

for best_seq in best_cluster_seqs:
    opposite_cluster.append(transcript_to_cluster[best_seq])

best_opposite_cluster = max(opposite_cluster, key=opposite_cluster.count)

    
# ----------------------------------------------
#                  Writing Clusters            #
#-----------------------------------------------

print ("Writing results...")
res = open(output_prefix + "_clusters.tsv", "w")
res.write("cluster_id\ttranscripts_ids\n")

for cluster_id, cluster_reads in clusters.items():
    res.write(str(cluster_id) + '\t' + ','.join(cluster_reads) + "\n")

res.close()


# ----------------------------------------------
#              Printing Summary Report         #
#-----------------------------------------------

"""
Best query Cluster
Best reference Cluster
"""

report = ("At cut-off threshold %.2f Best query Cluster is cluster:%d & Best matched reference cluster:%d") % (query_threshold, best_cluster_id, best_opposite_cluster)
report += ("with best cluster")
print (report)
