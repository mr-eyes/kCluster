"""
A script to append Cluster ID to fasta headers
"""


import re
import sys

fasta_file = ""
clstr_file = ""
threshold = ""

if len(sys.argv) < 4:
    exit("run: python <fasta_file> <clstr_file> <threshold%>")

else:
    fasta_file = sys.argv[1]
    clstr_file = sys.argv[2]
    threshold = sys.argv[3]


clstr_file = open(clstr_file, "r")
clstr_data = clstr_file.read()
clstr_file.close()

rep = {"\t": ",", "at +/": "", "at -/": "","...": ",", "nt": "", "%": "", " ": ""}
rep = dict((re.escape(k), v) for k, v in rep.iteritems())
pattern = re.compile("|".join(rep.keys()))
clstr_data = pattern.sub(lambda m: rep[re.escape(m.group(0))], clstr_data)
all_clusters = clstr_data.split(">Cluster")

clusters_transcripts_ids = {}
transcript_to_cluster = {}

for i in range(1, len(all_clusters), 1):
    cluster = all_clusters[i]
    cluster = cluster.split("\n")
    cluster_id = int(cluster[0])

    for item in cluster[1:-1]:
        item = item.replace(">", "").split(",")
        transcript_id = item[2].split("|")[0]
        gene_id = item[2].split("|")[1]
        
        transcript_to_cluster[transcript_id] = cluster_id
        
    
        if cluster_id in clusters_transcripts_ids:
            clusters_transcripts_ids[cluster_id].append(transcript_id)

        else:
            clusters_transcripts_ids[cluster_id] = [transcript_id]


new_fasta = open("clstr" + threshold + "_" + fasta_file, 'w')

with open (fasta_file, 'r') as fasta:
    for line in fasta:
        if line[0] == ">":
            line = line[1:-1]
            transcript_id = line.split("|")[0]
            line = ">" + line + "CLSTR_" + str(transcript_to_cluster[transcript_id]) + "|\n"
            new_fasta.write(line)
        else:
            new_fasta.write(line)

new_fasta.close()
