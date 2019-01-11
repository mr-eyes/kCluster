from prettytable import PrettyTable as PT
import sys
import re


"""
Parsing CDHIT CLSTR FILE
"""
CDH_seqs = {}
# Reading CDHIT transcripts to clusterID
clstr_file = open("cdhit/cdhit_pc2_transcripts_80.fa.clstr", "r")
clstr_data = clstr_file.read()
clstr_file.close()

rep = {"\t": ",", "at +/": "", "at -/": "",
       "...": ",", "nt": "", "%": "", " ": ""}
rep = dict((re.escape(k), v) for k, v in rep.iteritems())
pattern = re.compile("|".join(rep.keys()))
clstr_data = pattern.sub(lambda m: rep[re.escape(m.group(0))], clstr_data)
all_clusters = clstr_data.split(">Cluster")

cdhit_tr_cluster = {}
for i in range(1, len(all_clusters), 1):
    cluster = all_clusters[i]
    cluster = cluster.split("\n")
    cluster_id = cluster[0]
    for item in cluster[1:-1]:
        if cluster_id in CDH_seqs:
            CDH_seqs[cluster_id] += 1
        else:
            CDH_seqs[cluster_id] = 1


"""
Parsing kCluster Clusters file
"""
kCl_seqs = {}
with open("clusters_c66.0_66%_clusters.tsv", "r") as kCL:
    next(kCL)
    for line in kCL:
        cline = line.split()
        cluster_id = cline[0]
        no_of_seqs = len(cline[1].split(","))
        kCl_seqs[cluster_id] = no_of_seqs


"""
Constructing kCL to CDH
"""
kCl_to_CDH = {"CC": dict(), "IC": dict(), "CM": dict(), "IM": dict()}
with open("uniq_bio_assess.tsv", 'r') as tsv:
    next(tsv)
    for line in tsv:
        line = line.split()
        kCl_ID = line[0]
        kCl_type = line[1]
        CDH_ID = line[2]
        CDH_type = line[3]

        if kCl_ID in kCl_to_CDH[kCl_type]:
            kCl_to_CDH[kCl_type][kCl_ID].append(CDH_ID)
        else:
            kCl_to_CDH[kCl_type][kCl_ID] = [CDH_ID]


"""
CDHIT To kCluster
"""
CDH_to_kCl = {"CC": dict(), "IC": dict(), "CM": dict(), "IM": dict()}
with open("uniq_bio_assess.tsv", 'r') as tsv:
    next(tsv)
    for line in tsv:
        line = line.split()
        kCl_ID = line[0]
        kCl_type = line[1]
        CDH_ID = line[2]
        CDH_type = line[3]

        if CDH_ID in CDH_to_kCl[kCl_type]:
            CDH_to_kCl[CDH_type][CDH_ID].append(kCl_ID)
        else:
            CDH_to_kCl[CDH_type][CDH_ID] = [kCl_ID]


"""
CDHIT -> kCl
"""
summary = PT()
summary.field_names = ["kCl Type", "kCluster",
                       "CDHIT", "kCluster Seqs", "CDHIT Seqs"]
summary_total = ["Total", 0, 0, 0, 0]

for TYPE in ["CC", "IC", "IM", "CM"]:
    values_len = 0
    all_vals = set()
    keys_len = 0
    kCL_total_seqs = 0
    CDH_total_seqs = 0
    for key, val in kCl_to_CDH[TYPE].iteritems():
        keys_len += 1
        kCL_total_seqs += kCl_seqs[key]
        for _CDH_cluster in val:
            all_vals.add(_CDH_cluster)

    for _CDH_cluster in all_vals:
        values_len += 1
        CDH_total_seqs += CDH_seqs[_CDH_cluster]

    summary.add_row([TYPE, keys_len, values_len,
                     kCL_total_seqs, CDH_total_seqs])
    summary_total[1] += keys_len
    summary_total[2] += values_len
    summary_total[3] += kCL_total_seqs
    summary_total[4] += CDH_total_seqs

summary.add_row(["---", "---", "---", "---", "---"])
summary.add_row(summary_total)


"""
kCl -> CDHIT
"""
summary2 = PT()
summary2.field_names = ["CDHIT Type", "CDHIT",
                        "kCluster", "CDHIT Seqs", "kCluster Seqs"]

summary2_total = ["Total", 0, 0, 0, 0]
for TYPE in ["CC", "IC", "IM", "CM"]:
    values_len = 0
    all_vals = set()
    keys_len = 0
    all_keys = set()
    kCL_total_seqs = 0
    CDH_total_seqs = 0

    for key, val in CDH_to_kCl[TYPE].iteritems():
        all_keys.add(key)

        for _kCl_cluster in val:
            all_vals.add(_kCl_cluster)

    for _kCl_cluster in all_vals:
        values_len += 1
        kCL_total_seqs += kCl_seqs[_kCl_cluster]

    for key in all_keys:
        keys_len += 1
        CDH_total_seqs += CDH_seqs[key]

    summary2.add_row([TYPE, keys_len, values_len,
                      CDH_total_seqs, kCL_total_seqs])
    summary2_total[1] += keys_len
    summary2_total[2] += values_len
    summary2_total[3] += CDH_total_seqs
    summary2_total[4] += kCL_total_seqs

summary2.add_row(["---", "---", "---", "---", "---"])
summary2.add_row(summary2_total)


"""
Printing Results
"""
print "CDHIT -> kCl"
print summary
print "\n", ("~" * 65), "\n"
print "kCl -> CDHIT"
print summary2
