from collections import Counter
import itertools
import re
import gc
import tqdm
import json
import sys
import os

def ids_to_names(cluster):
    new_cluster = []

    for id in cluster:
        new_cluster.append(names_map[id])
    
    return new_cluster
    

def get_all_connected_groups(graph):
    already_seen = set()
    result = []
    for node in graph:
        #print (node)
        if node not in already_seen:
            connected_group, already_seen = get_connected_group(node, already_seen)
            result.append(connected_group)
    return result


def get_connected_group(node, already_seen):
        result = []
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            already_seen.add(node)
            nodes = nodes or graph[node] - already_seen
            result.append(node)
        return result, already_seen


cut_off_threshold = 0.0
relations_file = ""
names_map_file = ""
output_file = ""

if len(sys.argv) < 3:
    exit("Please pass <relations_tsv_file> <names_map_file> <op: cuttof_threshold %> <op: output_file_name> <op: output_file_name>")

else:
    relations_file = sys.argv[1]
    names_map_file = sys.argv[2]
    
if len(sys.argv) >= 4:
    cut_off_threshold = float(sys.argv[3])
    
else:
    cut_off_threshold = 0.00000000001

if len(sys.argv) == 5:
        output_file = sys.argv[4]
        
else:
    output_file = os.path.basename(names_map_file).split(".")[0]




print ("Reading names...")
names_map = {}
with open(names_map_file) as namesMap:
    for name in namesMap:
        names_map[int(re.findall(r'\t(\d+)', name)[0])] = re.findall(r'(.*)\t', name)[0]


graph = {}
print ("Building the graph...")


with open(relations_file) as REL:
    next(REL)
    for line in REL:
        info = re.findall(r'(\d+(?:\.\d+)?)' , line)
        _seq1 = int(info[0])
        _seq2 = int(info[1])
        _norm = float(info[3])

        if _norm < cut_off_threshold:
            continue
        
        if _seq1 in graph:
            graph[_seq1].add(_seq2)
        else:
            graph[_seq1] = {_seq2}
        
        if _seq2 in graph:
            graph[_seq2].add(_seq1)
        else:
            graph[_seq2] = {_seq1}

            
for i in range(1, len(names_map), 1):
    if i not in graph:
        graph[i] = set()

print ("Clustering...")
components = get_all_connected_groups(graph)

print ("Writing results...")
res = open("clusters_c" + str(cut_off_threshold) + "_" + output_file + ".tsv", "w")
res.write("cluster_id\ttranscripts_ids\n")
cluster_id = 0
for i in range(len(components)):
    component = components[i]
    res.write(str(cluster_id) + '\t' + ','.join(ids_to_names(component)) + "\n")
    cluster_id += 1

res.close()

print "Number of clusters: %d" % (cluster_id)
