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


cut_off_threshold = 1
names_map_file = ""
map_index_file = ""
output_file = ""

if len(sys.argv) < 3:
    exit("Please pass <map_index_file> <names_map_file> <cuttof_threshold>")

else:
    map_index_file = sys.argv[1]
    names_map_file = sys.argv[2]
    cut_off_threshold = int(sys.argv[3])

if len(sys.argv) == 5:
        output_file = sys.argv[4]

else:
    output_file = os.path.basename(map_index_file).split(".")[0]


print ("Reading names...")
names_map = {}
with open(names_map_file) as namesMap:
    for name in namesMap:
        names_map[int(re.findall(r'\t(\d+)', name)[0])] = re.findall(r'(.*)\t', name)[0]

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
#print ("Collecting Garbage")
gc.collect()

graph = {}

print ("Processing groups...")
for color, tr_ids in tqdm.tqdm(groups.items()):
    color_count = colors[color]
    #print color_count

    if len(tr_ids) == 1:
        continue

    if color_count < cut_off_threshold:
        continue

    for combination in itertools.combinations(tr_ids,2):
        if combination[0] in graph:
            graph[combination[0]].add(combination[1])
        else:
            graph[combination[0]] = {combination[1]}
        
        if combination[1] in graph:
            graph[combination[1]].add(combination[0])
        else:
            graph[combination[1]] = {combination[0]}
            
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
    #if len(component) > 0:
    res.write(str(cluster_id) + '\t' + ','.join(ids_to_names(component)) + "\n")
    cluster_id += 1

res.close()

print "Number of clusters: %d" % (cluster_id)
