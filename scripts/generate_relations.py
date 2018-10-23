"""
generate_relations.py 

input: 
    - MAP index file
    - Names MAP file

output:
    - JSON & TSV files holds the information of shared kmers count between each two transcripts.
"""


from collections import Counter
import itertools
import re
import gc
import tqdm
import json
import sys

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


names_map = []
with open(names_map_file) as namesMap:
    for name in namesMap:
        names_map.append(int(re.findall(r'\t(\d+)', name)[0]))

print ("Done reading names...")

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

edges = {}

print ("Processing groups...")
for color, tr_ids in tqdm.tqdm(groups.items()):
    color_count = colors[color]
    if len(tr_ids) == 1:
        continue

    for combination in itertools.combinations(tr_ids,2):
        if combination[0] in edges:
            if combination[1] in edges[combination[0]]:
                edges[combination[0]][combination[1]] += color_count
            else:
                edges[combination[0]][combination[1]] = color_count
        
        else:
            edges[combination[0]] = {combination[1]: color_count}

del colors
del groups
gc.collect()

print ("Writing JSON file ...")
json_file = open(output_file + ".json", "w")
json_file.write(json.dumps(edges, sort_keys=True,
                           indent=4, separators=(',', ': ')))
json_file.close()


print ("Writing TSV file ...")

tsv = open(output_file + ".tsv", "w")

for _1st, info in tqdm.tqdm(edges.items()):
    for _2nd, _weight in info.items():
        l = "%d\t%d\t%d\n" % (_1st, _2nd, _weight)
        tsv.write(l)

tsv.close()
