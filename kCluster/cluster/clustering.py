from __future__ import division
from collections import defaultdict
import itertools
import sys
import os
import sqlite3
import click

class kClusters:

    source = []
    target = []
    names_map = dict()
    components = defaultdict(set)

    def __init__(self, sqlite_file, userQs, cut_off_threshold):

        try:
            self.conn = sqlite3.connect(sqlite_file)
        
        except sqlite3.Error as err:
            print(err)
            print(f"couldn't connect to {sqlite_file}", file = sys.stderr)
            sys.exit(1)

        # Check user selected Qs
        sqliteQs = [int(Q.split("_")[1]) for Q in self.sqlite_getQs()]
        
        if userQs == [-1]:
            self.userQs = sqliteQs
        else:
            commonQs = set(userQs).intersection(set(sqliteQs))
            ln_commonQs = len(commonQs)
            ln_userQs = len(userQs)
            
            if not ln_commonQs:
                print("invalid virtualQs range, none of them found in the sqlite database", file = sys.stderr)
                sys.exit(1)
            
            elif ln_commonQs == ln_userQs:
                self.userQs = userQs
            
            else:
                unfound = set(userQs) - set(sqliteQs)
                print(f"The following Qs {unfound} couldn't be found in the sqlite database, processing the others..", file = sys.stderr)
                self.userQs = list(commonQs)
            
        
        self.cut_off_threshold = cut_off_threshold
        self.sqlite_file = sqlite_file
        self.kSize = self.conn.execute("SELECT value from meta_info WHERE key='kSize'").fetchone()[0]
        self.sqlite_get_namesmap()

    def ids_to_names(self, cluster):
        new_cluster = []
        for id in cluster:
            new_cluster.append(self.names_map[id])

        return new_cluster

    def sqlite_getQs(self):
        """get all the Qs values in the database
        Returns:
            set of table Qs values.
        """

        gold_names = {'ID', 'seq1', 'seq2', 'min_kmers'}
        cursor = self.conn.execute('select * from virtualQs')
        cols_names = set(map(lambda x: x[0], cursor.description))
        if len(gold_names.intersection(cols_names)) != 4:
            return False
        else:
            return cols_names - gold_names

    def sqlite_get_namesmap(self):
        cursor = self.conn.execute("SELECT * FROM namesmap")
        cursor = cursor.fetchall()
        for row in cursor:
            self.names_map[row[1]] = row[2]

    def build_graph(self):
        user_Qs = self.userQs
        user_Qs.sort(reverse = True)
        indeces = {i:idx+3 for idx,i in enumerate(user_Qs)}
        query_Qs = ", ".join([f"Q_{Q}" for Q in user_Qs])

        curs = self.conn.execute(f'select seq1, seq2, min_kmers, {query_Qs} from virtualQs')

        for row in curs:
            seq1 = row[0]
            seq2 = row[1]
            min_kmers = row[2]
            similarity = 0.0
            Q_prev = 0
            
            for Q_val, idx in indeces.items():
                Q_curr = row[idx]
                sim = ( (Q_curr - Q_prev) * (Q_val / self.kSize) ) / min_kmers
                similarity += sim
                Q_prev = Q_curr

            if similarity < self.cut_off_threshold:
                continue
            
            self.source.append(seq1)
            self.target.append(seq2)


        for i in range(1, len(self.names_map) + 1, 1):
                if i not in self.source and i not in self.target:
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
        kCluster_file_name = f"clusters_{self.cut_off_threshold:.2f}%_"
        kCluster_file_name += os.path.basename(self.sqlite_file).split(".")[0]
        kCluster_file_name += ".tsv"
        
        with open(kCluster_file_name, 'w') as kClusters:
            kClusters.write("kClust_id\tseqs_ids\n")
            for cluster_id, (k, v) in enumerate(self.components.items(), 1):
                kClusters.write(f"{cluster_id}\t{','.join(self.ids_to_names(v))}\n")

        print(f"Total Number Of Clusters: {cluster_id}", file = sys.stderr)


@click.command()
@click.option('-m','--min-q', required=False, type=int, default = None, help="minimum virtualQ")
@click.option('-M','--max-q', required=False, type=int, default = None, help="maximum virtualQ")
@click.option('-s','--step-q', required=False, type=int, default = None, help="virtualQs range step")
@click.option('-c','--cutoff', required=False, type=click.FloatRange(0, 1, clamp=False) , default = 0.0, show_default=True, help="cluster sequences with (similarity > cutoff)")
@click.option('-d', '--db', required=True, type=click.Path(exists=True), help="sqlite database file")
def main(min_q, max_q, step_q, db, cutoff):
    """Sequences clustering regarding user-selected virtualQs."""

    scanQs = map(bool, [min_q, max_q, step_q])

    if True in scanQs and False in scanQs:
        print("Please complete the virtualQs range", file = sys.stderr)
        exit(1)
    
    elif not(min_q and max_q and step_q):
        print("processing all virtualQs in the sqlite database", file = sys.stderr)
        userQs = [-1]
    
    else:
        userQs = [Q for Q in range(min_q, max_q + 1, step_q)]
    
    kCl = kClusters(db, userQs, cutoff)
    kCl.build_graph()
    kCl.clustering()
    kCl.export_kCluster()


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter