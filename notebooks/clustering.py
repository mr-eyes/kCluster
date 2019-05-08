from collections import Counter, defaultdict
import itertools
import re
import gc
import tqdm
import json
import sys
import os


class kClusters:

    names_map = {}
    source = []
    target = []
    components = defaultdict(set)

    def __init__(self, pairwise_file, names_map_file, cut_off_threshold = 0.0):
        
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file = pairwise_file

        with open(names_map_file) as namesMap:
            next(namesMap)
            for name in namesMap:
                name = name.strip().split()
                self.names_map[int(name[0])] = name[1]


    def ids_to_names(self, cluster):
        new_cluster = []
        for id in cluster:
            new_cluster.append(self.names_map[id])

        return new_cluster


    def build_graph(self):
        print ("Building the graph...")

        with open(self.pairwise_file) as PRW:
            next(PRW)
            for line in PRW:
                info = re.findall(r'(\d+(?:\.\d+)?)', line)
                _seq1 = int(info[0])
                _seq2 = int(info[1])
                _norm = float(info[-1])

                if _norm < self.cut_off_threshold:
                    continue

                self.source.append(_seq1)
                self.target.append(_seq2)

        for i in range(1, len(self.names_map) + 1, 1):
                self.source.append(i)
                self.target.append(i)

    def clustering(self):
        registers = defaultdict(lambda: None)

        def find(x):
            l = registers[x]
            if l is not None:
                l = find(l)
                registers[x] = l
                return l
            return x

        def union(x, y):
            lx, ly = find(x), find(y)
            if lx != ly:
                registers[lx] = ly

        
        for i in range(len(self.source)):
            union(self.source.pop(), self.target.pop())


        for x in registers:
            self.components[find(x)].add(x)


    def export_kCluster(self):
        kCluster_file_name = f"kCluster_{self.cut_off_threshold:.2f}%_"
        kCluster_file_name += os.path.basename(self.pairwise_file).split(".")[0]
        kCluster_file_name += ".tsv"
        
        with open(kCluster_file_name, 'w') as kClusters:
            kClusters.write("kClust_id\tseqs_ids\n")
            for cluster_id, (k, v) in enumerate(self.components.items()):
                #kClusters.write(f"{cluster_id + 1}\t{','.join(self.ids_to_names(v))}\n")
                kClusters.write(f"{cluster_id}\t{','.join(list(map(str,v)))}\n")  

        print(f"Total Number Of Clusters: {cluster_id}")




cut_off_threshold = 0.0
relations_file = ""
names_map_file = ""
output_file = ""

if len(sys.argv) < 3:
    exit("Please pass <relations_tsv_file> <names_map_file> <op: cuttof_threshold %> <op: output_file_name>")

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

kCl = kClusters(relations_file, names_map_file, cut_off_threshold)
kCl.build_graph()
kCl.clustering()
kCl.export_kCluster()