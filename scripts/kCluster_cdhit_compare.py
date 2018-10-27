from __future__ import division
import tqdm
import re
import sys
import json
import math


def transcripts_to_clusters(transcipts_ids):
    clusters = {}  # {cluster:counted}
    for transcript in transcipts_ids:
        cluster = transcript_cluster[transcript]

        if cluster in clusters:
            clusters[cluster] += 1
        else:
            clusters[cluster] = 1

    return clusters


def transcripts_to_genes(transcipts_ids):
    genes = {}  # {cluster:counted}
    for transcript in transcipts_ids:
        gene = transcript_gene[transcript]

        if gene in genes:
            genes[gene] += 1
        else:
            genes[gene] = 1

    return genes


def how_many_clusters(cluster):
    return len(transcripts_to_clusters(cluster))


def how_many_genes(cluster):
    return len(transcripts_to_genes(cluster))


def how_many_complete_clusters(cluster):
    clusters = transcripts_to_clusters(cluster)
    complete = 0
    for key, value in clusters.iteritems():
        if len(cluster_transcripts[key]) == value:
            complete += 1

    return complete


def how_many_complete_genes(cluster):
    genes = transcripts_to_genes(cluster)
    complete = 0
    for key, value in genes.iteritems():
        if len(gene_transcripts[key]) == value:
            complete += 1

    return complete


def Q1(cluster):
    cluster_genes = transcripts_to_genes(cluster)
    for key, value in cluster_genes.iteritems():
        if len(gene_transcripts[key]) != value:
            return False

    return True


def Q2(cluster):
    if len(transcripts_to_clusters(cluster)) > 1:
        return True
    else:
        return False


def build_stats(cluster_type, no_clusters, no_genes, no_complete_genes, no_complete_clusters):
    stats["clusters"][cluster_type].append(no_clusters)
    stats["genes"][cluster_type].append(no_genes)
    stats["complete-genes"][cluster_type].append(no_complete_genes)
    stats["complete-clusters"][cluster_type].append(no_complete_clusters)


def _mean(lst):
    return round(sum(lst) / len(lst), 2)


def _std(lst):
    mean = sum(lst) / len(lst)   # mean
    var = sum(pow(x-mean, 2) for x in lst) / len(lst)  # variance
    std = math.sqrt(var)  # standard deviation
    return round(std, 2)


def _median(lst):
    n = len(lst)
    if n < 1:
            return None
    if n % 2 == 1:
            return sorted(lst)[n//2]
    else:
            return sum(sorted(lst)[n//2-1:n//2+1])/2.0

fasta_file_path = ""
clstr_file_path = ""
output_file = ""

if len(sys.argv) < 3:
    sys.exit(
        "Kindly pass positional arguments, ex: python clusters_assessment.py [fasta_file] [clstr_file]")

else:
    fasta_file_path = sys.argv[1]
    clstr_file_path = sys.argv[2]
    output_file = sys.argv[3]

_complete_mixed = 0
_complete_clean = 0
_incomplete_mixed = 0
_incomplete_clean = 0

stats = {"clusters": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "genes": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "complete-genes": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []},
         "complete-clusters": {"_complete_mixed": [], "_complete_clean": [], "_incomplete_mixed": [], "_incomplete_clean": []}
         }


cluster_transcripts = {}  # cluster_0 : gene1,gene2,...
transcript_cluster = {}  # gene1:cluster1, gene2:cluster1, gene3:cluster3

gene_transcripts = {}  # gene_id: transcript1,transcript2,.....
transcript_gene = {}   # transcript1: gene9, transcript2: gene1, ....

print "Reading Fasta..."
with open(fasta_file_path) as fa:
    for line in fa:
        if line[0] != ">":
            continue
        fields = line.split("|")
        transcript_id = fields[0][1:]
        gene = fields[-2]
        cluster = fields[-2]
        transcript_cluster[transcript_id] = cluster
        transcript_gene[transcript_id] = gene

        if gene in gene_transcripts:
            gene_transcripts[gene].append(transcript_id)
        else:
            gene_transcripts[gene] = [transcript_id]

        if cluster in cluster_transcripts:
            cluster_transcripts[cluster].append(transcript_id)
        else:
            cluster_transcripts[cluster] = [transcript_id]


clusters_transcripts_ids = {}

with open(clstr_file_path, "r") as clusters:
    next(clusters)
    for line in clusters:
        fields = line.split("\t")
        cluster_id = int(fields[0])
        _transcripts = fields[1].split(",")
        transcripts = []
        for tr in _transcripts:
            transcripts.append(tr.replace("\n", ""))

        clusters_transcripts_ids[cluster_id] = transcripts

res = open(output_file, "w")
res.write("cluster_id\tQ1\tQ2\tclusters\tcomplete_clusters\tkCluster\tcdhit\n")


for cluster_id, transcripts_ids in sorted(clusters_transcripts_ids.iteritems()):
    q1 = Q1(transcripts_ids)
    q2 = Q2(transcripts_ids)
    no_clusters = how_many_clusters(transcripts_ids)
    no_genes = how_many_genes(transcripts_ids)
    no_complete_clusters = how_many_complete_clusters(transcripts_ids)
    no_complete_genes = how_many_complete_genes(transcripts_ids)

    ans1 = ""
    ans2 = ""

    kCluster_clusters = ""
    cdhit_clusters = ""

    if q1 == True and q2 == True:
        ans1, ans2 = "Complete", "Mixed"
        _complete_mixed += 1
        cdhit_clusters = ",".join(transcripts_to_clusters(transcripts_ids))
        build_stats("_complete_mixed", no_clusters, no_genes,
                    no_complete_genes, no_complete_clusters)
    if q1 == True and q2 == False:
        ans1, ans2 = "Complete", "Clean"
        cdhit_clusters = ",".join(transcripts_to_clusters(transcripts_ids))
        _complete_clean += 1
        build_stats("_complete_clean", no_clusters, no_genes,
                    no_complete_genes, no_complete_clusters)
    if q1 == False and q2 == True:
        ans1, ans2 = "InComplete", "Mixed"
        cdhit_clusters = ",".join(transcripts_to_clusters(transcripts_ids))
        _incomplete_mixed += 1
        build_stats("_incomplete_mixed", no_clusters, no_genes,
                    no_complete_genes, no_complete_clusters)
    if q1 == False and q2 == False:
        ans1, ans2 = "InComplete", "Clean"
        cdhit_clusters = ",".join(transcripts_to_clusters(transcripts_ids))
        _incomplete_clean += 1
        build_stats("_incomplete_clean", no_clusters, no_genes,
                    no_complete_genes, no_complete_clusters)

    line = str(cluster_id) + "\t" + ans1 + "\t" + ans2 + "\t" + str(no_clusters) + "\t" + str(no_complete_clusters) + "\t" + str(cluster_id) + "\t" + cdhit_clusters + "\t" + "\n"
    res.write(line)

res.close()

# Writing summary file of counts ________________________________

summary = open(output_file.split(".")[0] + "_summary.txt", "w")
summary.write(
    ("%d Complete Mixed Components | [_complete_mixed]\n") % (_complete_mixed))
summary.write(
    ("%d Complete Clean Components | [_complete_clean]\n") % (_complete_clean))
summary.write(("%d Incomplete Mixed  Components | [_incomplete_mixed]\n") % (
    _incomplete_mixed))
summary.write(("%d Incomplete Clean Components | [_incomplete_clean]\n") % (
    _incomplete_clean))
summary.close()

# Writing statistics json file ________________________________

# EXIT
exit()

json_output = {}
for cluster_type in ["_incomplete_mixed", "_incomplete_clean", "_complete_mixed", "_complete_clean"]:
    result = {
        'mean': {
            'no_genes': _mean(stats["genes"][cluster_type]),
            'complete_genes':  _mean(stats["complete-genes"][cluster_type]),
            'no_clusters':  _mean(stats["clusters"][cluster_type]),
            'complete_clusters':  _mean(stats["complete-clusters"][cluster_type])
        },
        'std': {
            'no_genes': _std(stats["genes"][cluster_type]),
            'complete_genes':  _std(stats["complete-genes"][cluster_type]),
            'no_clusters':  _std(stats["clusters"][cluster_type]),
            'complete_clusters':  _std(stats["complete-clusters"][cluster_type])
        },
        'min': {"no_genes": min(stats["genes"][cluster_type]),
                "complete_genes": min(stats["complete-genes"][cluster_type]),
                "no_clusters": min(stats["clusters"][cluster_type]),
                "complete_clusters": min(stats["complete-clusters"][cluster_type])},
        'max': {"no_genes": max(stats["genes"][cluster_type]),
                "complete_genes": max(stats["complete-genes"][cluster_type]),
                "no_clusters": max(stats["clusters"][cluster_type]),
                "complete_clusters": max(stats["complete-clusters"][cluster_type])},
        'median': {"no_genes": _median(stats["genes"][cluster_type]),
                "complete_genes": _median(stats["complete-genes"][cluster_type]),
                "no_clusters": _median(stats["clusters"][cluster_type]),
                "complete_clusters": _median(stats["complete-clusters"][cluster_type])}
    }
    json_output[cluster_type] = result


json_file = open(output_file.split(".")[0] + "_stats.json", "w")
json_file.write(json.dumps(json_output, sort_keys=True,
                           indent=4, separators=(',', ': ')))
json_file.close()
