from prettytable import PrettyTable as PT
import re
from pprint import pprint, pformat
import sys

if len(sys.argv) < 5:
    exit("run: python generate_main_tsv.py <cdhit*clstr> <kCluster_clusters*tsv> <unique_bio_assess*tsv> <total_bio_assess*tsv>")
else:
    cdhit_clstrs_file = sys.argv[1]
    kCluster_clusters_file = sys.argv[2]
    unique_bio_assess = sys.argv[3]
    inunique_bio_assess = sys.argv[4]



CDH_seqs = {}
# Reading CDHIT transcripts to clusterID
clstr_file = open(cdhit_clstrs_file, "r")
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
with open(kCluster_clusters_file, "r") as kCL:
    next(kCL)
    for line in kCL:
        cline = line.split()
        cluster_id = cline[0]
        no_of_seqs = len(cline[1].split(","))
        kCl_seqs[cluster_id] = no_of_seqs



"""
---------------------------------------------------------------------
"""


empty_d = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}
kCl_to_CDH = {"CC": dict(empty_d), "IC": dict(
    empty_d), "CM": dict(empty_d), "IM": dict(empty_d)}
with open(inunique_bio_assess, 'r') as tsv:
    next(tsv)
    for line in tsv:
        line = line.split()
        kCl_ID = line[0]
        kCl_type = line[1]
        CDH_ID = line[2]
        CDH_type = line[3]
        kCl_to_CDH[kCl_type][CDH_type] += 1

pr = PT()
pr.field_names = ["kCl/CDH", "CC", "IC", "CM", "IM", "Total"]
total_clusters = 0
TYPES = ["CC", "IC", "CM", "IM"]

cols = [0, 0, 0, 0]
total_sum = 0
for TYPE in TYPES:
    row = []
    for T in TYPES:
        total_clusters += kCl_to_CDH[TYPE][T]
        row.append(kCl_to_CDH[TYPE][T])
    pr.add_row([TYPE]+row+[str(sum(row))])
    total_sum += sum(row)
    for i in range(len(row)):
        cols[i] += row[i]
s_cols = map(str, cols)
pr.add_row(["---", "---", "---", "---", "---", "---"])
pr.add_row(["Total"] + s_cols + [str(total_sum)])

print pr



"""
---------------------------------------------------------------------
"""



empty_d = {"CC":0, "IC":0, "CM":0, "IM":0}
pairwise_count = {"CC":dict(empty_d), "IC":dict(empty_d), "CM":dict(empty_d), "IM":dict(empty_d)}
with open(unique_bio_assess, 'r') as tsv:
    next(tsv)
    for line in tsv:
        line = line.split()
        kCl_ID = line[0]
        kCl_type = line[1]
        CDH_ID = line[2]
        CDH_type = line[3]
        pairwise_count[kCl_type][CDH_type] += 1


"""
---------------------------------------------------------------------
"""





"""

"""
kCl_to_CDH = {"CC": dict(), "IC": dict(), "CM": dict(), "IM": dict()}
CDH_to_kCl = {"CC": dict(), "IC": dict(), "CM": dict(), "IM": dict()}
CDH_to_type = {}
kCl_to_type = {}

with open(unique_bio_assess, 'r') as tsv:
    next(tsv)
    for line in tsv:
        line = line.split()
        kCl_ID = line[0]
        kCl_type = line[1]
        CDH_ID = line[2]
        CDH_type = line[3]
        CDH_to_type[CDH_ID] = CDH_type
        kCl_to_type[kCl_ID] = kCl_type

        # kCl To CDHIT
        if kCl_ID in kCl_to_CDH[kCl_type]:
            kCl_to_CDH[kCl_type][kCl_ID].append(CDH_ID)
        else:
            kCl_to_CDH[kCl_type][kCl_ID] = [CDH_ID]

        # CDHIT To kCl
        if CDH_ID in CDH_to_kCl[CDH_type]:
            CDH_to_kCl[CDH_type][CDH_ID].append(kCl_ID)
        else:
            CDH_to_kCl[CDH_type][CDH_ID] = [kCl_ID]


"""
---------------------------------------------------------------------
"""

empty_dic1 = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}
kCl_to_CDH_clstrs = {"CC": dict(empty_dic1), "IC": dict(empty_dic1), "CM": dict(empty_dic1), "IM": dict(empty_dic1)}
KCL = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}
uq_CDH = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}

for TYPE in ["CC","IC","IM","CM"]:
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
        _cdh_type = CDH_to_type[_CDH_cluster]
        kCl_to_CDH_clstrs[TYPE][_cdh_type] += 1
        values_len += 1
    
        uq_CDH[TYPE] += 1
    
    KCL[TYPE] = keys_len


CDH_to_kCl_clstrs = {"CC": dict(empty_dic1), "IC": dict(empty_dic1), "CM": dict(empty_dic1), "IM": dict(empty_dic1)}
CDH = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}
uq_KCL = {"CC": 0, "IC": 0, "CM": 0, "IM": 0}

for TYPE in ["CC","IC","IM","CM"]:
    values_len = 0
    all_vals = set()
    keys_len = 0
    all_keys = set()
    kCL_total_seqs = 0
    CDH_total_seqs = 0
    
    for key, val in CDH_to_kCl[TYPE].iteritems():
        all_keys.add(key)
        keys_len += 1
        
        for _kCl_cluster in val:
            all_vals.add(_kCl_cluster)
    
    for _kCl_cluster in all_vals:
        values_len += 1
        _kCl_type = kCl_to_type[_kCl_cluster]
        CDH_to_kCl_clstrs[TYPE][_kCl_type] += 1
        kCL_total_seqs += kCl_seqs[_kCl_cluster]
        uq_KCL[TYPE] += 1
    
    
    CDH[TYPE] = keys_len
    

result = PT()
head = []

for TYPE in ["CC", "IC", "IM", "CM"]:
    head.append(TYPE+ " " + str(CDH[TYPE]))

result.field_names = ["kCl/CDH"] + head + ["Total","uq_Total"]

col_total = {"CC": 0, "IC": 0, "CM": 0, "IM": 0, "Total": 0}


for k, v in pairwise_count.iteritems():
    total = sum(v.values())
    str_v = map(str, v.values())
    col_total["Total"] += total
    result.add_row([str(KCL[k]) + " " + k] + str_v + [str(total)] + [str(uq_CDH[k])])


#Total Row
for i in pairwise_count.values():
    for k, v in i.iteritems():
        col_total[k] += v

str_col_total = []
for TYPE in ["CC", "IC", "IM", "CM", "Total"]:
    str_col_total.append(str(col_total[TYPE]))

str_col_total.append(str(sum(uq_CDH.values()))) # Add cdhit uq total


result.add_row(["----", "----", "----", "----", "----", "----", "----"])
result.add_row(["Total"] + str_col_total)

# kCl uq_total
uq_total_str = []
for TYPE in ["CC", "IC", "IM", "CM"]:
    uq_total_str.append(str(uq_KCL[TYPE]))

result.add_row(["uq_Total"] + uq_total_str + ["----", "----"])

print result
