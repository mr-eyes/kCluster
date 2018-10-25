from collections import Counter
import itertools
import re
import gc
import tqdm
import json
import sys

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


names_map_file = ""
map_index_file = ""
output_file = ""

if len(sys.argv) < 3:
    exit("Please pass <map_index_file> <names_map_file>")

else:
    map_index_file = sys.argv[1]
    names_map_file = sys.argv[2]

if len(sys.argv) == 4:
        output_file = sys.argv[3]

else:
    output_file = map_index_file.split(".")[0]


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
print ("Collecting Garbage")
gc.collect()

graph = {}

print ("Processing groups...")
for color, tr_ids in tqdm.tqdm(groups.items()):
    color_count = colors[color]
    if len(tr_ids) == 1:
        #graph[tr_ids[0]] = tr_ids[0]
        continue

    for combination in itertools.combinations(tr_ids,2):
        if combination[0] in graph:
            graph[combination[0]].add(combination[1])        
        else:
            graph[combination[0]] = {combination[1]}
            
for i in range(len(names_map)):
    if i not in graph:
        graph[i] = set()

print ("Clustering...")
components = get_all_connected_groups(graph)

print ("Writing results...")
res = open("clusters_" + output_file + ".tsv", "w")
res.write("cluster_id\ttranscripts_ids\n")
cluster_id = 0
for i in tqdm.tqdm(range(len(components))):
    component = components[i]
    if len(component) > 1:
        res.write(str(cluster_id) + '\t' + ','.join(ids_to_names(component)) + "\n")
        cluster_id += 1

res.close()
