"""
Output
+--------+----------+--------+----------+
| kCl ID | kCl Type | CDH ID | CDH Type |
+--------+----------+--------+----------+
|   1    |    CC    |   1    |    CM    |
+--------+----------+--------+----------+
|   1    |    CC    |   2    |    CM    |
+--------+----------+--------+----------+
|   2    |    IC    |   3    |    CM    |
+--------+----------+--------+----------+
"""


import re
from tqdm import tqdm
import sys

if len(sys.argv) < 7:
    exit("run: python generate_main_tsv.py <cdhit*clstr> <cdhit_assessement*tsv> <kCluster_clusters*tsv> <kCluster_assessement*tsv> <fasta_file> <output_file>")
else:
    cdhit_clstrs_file = sys.argv[1]
    cdhit_assessement_file = sys.argv[2]
    kCluster_clusters_file = sys.argv[3]
    kCluster_assessement_file = sys.argv[4]
    fasta_file = sys.argv[5]
    output_file = sys.argv[6]

# Reading CDHIT Assessment
cdhit_types = {}
with open(cdhit_assessement_file, 'r') as cdhit:
    next(cdhit) #skip header
    for cline in cdhit:
        cline = cline.split()
        cluster_id = cline[0]
        cluster_type = cline[1][0] + cline[2][0]
        cdhit_types[cluster_id] = cluster_type


# Reading CDHIT transcripts to clusterID
clstr_file = open(cdhit_clstrs_file, "r")
clstr_data = clstr_file.read()
clstr_file.close()

rep = {"\t": ",", "at +/": "", "at -/": "","...": ",", "nt": "", "%": "", " ": ""}
rep = dict((re.escape(k), v) for k, v in rep.iteritems())
pattern = re.compile("|".join(rep.keys()))
clstr_data = pattern.sub(lambda m: rep[re.escape(m.group(0))], clstr_data)
all_clusters = clstr_data.split(">Cluster")

cdhit_tr_cluster = {}

for i in tqdm(range(1, len(all_clusters), 1)):
    cluster = all_clusters[i]
    cluster = cluster.split("\n")
    cluster_id = cluster[0]

    for item in cluster[1:-1]:
        item = item.replace(">", "").split(",")
        transcript_id = item[2].split("|")[0]
        gene_id = item[2].split("|")[1]
        cdhit_tr_cluster[transcript_id] = cluster_id


# Reading KCluster Assessment
kCluster_types = {}
with open(kCluster_assessement_file, 'r') as kCluster:
    next(kCluster)  # skip header
    for cline in kCluster:
        cline = cline.split()
        cluster_id = cline[0]
        cluster_type = cline[1][0] + cline[2][0]
        kCluster_types[cluster_id] = cluster_type


# Reading kCluster transcripts to clusterID
kCluster_tr_cluster = {}
with open(kCluster_clusters_file, "r") as kCluster:
    next(kCluster)
    for line in kCluster:
        cline = line.split()
        cluster_id = cline[0]
        for tr in cline[1].split(","):
            kCluster_tr_cluster[tr] = cluster_id

unique_result = set()

with open(fasta_file, 'r') as fasta, open("total_" + output_file , 'w') as total_bio_assess:
    total_bio_assess.write("kCl_id\tkCl_type\tCDH_id\tCDH_type\n")
    for line in fasta:
        if line[0] == ">":
            tr_id = line.split("|")[0][1:]
            CDH_ID = cdhit_tr_cluster[tr_id]
            CDH_type = cdhit_types[CDH_ID]
            KCL_ID = kCluster_tr_cluster[tr_id]
            KCL_type = kCluster_types[KCL_ID]
            _record = "%s\t%s\t%s\t%s\n" % (KCL_ID, KCL_type, CDH_ID, CDH_type)
            unique_result.add(_record)
            total_bio_assess.write(_record)

with open("uniq_" + output_file, 'w') as bio_assess:
    bio_assess.write("kCl_id\tkCl_type\tCDH_id\tCDH_type\n")
    for res in sorted(unique_result):
         bio_assess.write(res)
