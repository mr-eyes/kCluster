from __future__ import division
from collections import defaultdict
import itertools
import sys
import os
import sqlite3
import click
from src.click_context import cli
import glob

class kClusters:

    source = []
    target = []
    names_map = dict()
    components = defaultdict(set)

    def sqlite_constructor(self, sqlite_file, userQs):
        try:
            self.conn = sqlite3.connect(sqlite_file)

        except sqlite3.Error as err:
            self.Logger.ERROR(f"couldn't connect to {sqlite_file}\n{err}")

        # Check user selected Qs
        sqliteQs = [int(Q.split("_")[1]) for Q in self.sqlite_getQs()]

        if userQs == [-1]:
            self.userQs = sqliteQs
        else:
            commonQs = set(userQs).intersection(set(sqliteQs))
            ln_commonQs = len(commonQs)
            ln_userQs = len(userQs)

            if not ln_commonQs:
                self.Logger.ERROR("invalid virtualQs range, none of them found in the sqlite database")

            elif ln_commonQs == ln_userQs:
                self.userQs = userQs

            else:
                unfound = set(userQs) - set(sqliteQs)
                self.Logger.WARNING(
                    f"The following Qs '{unfound}' couldn't be found in the sqlite database, processing the others..")
                self.userQs = list(commonQs)


        self.sqlite_file = sqlite_file
        self.kSize = self.conn.execute("SELECT value from meta_info WHERE key='kSize'").fetchone()[0]
        self.sqlite_get_namesmap()

    def tsv_constructor(self, tsv_file, userQs):
        tsvQs = set()
        userQs = set(userQs)
        with open(tsv_file, 'r') as pairwise_tsv:
            header_line = next(pairwise_tsv).strip().split()[3:]
            tsvQs = {int(Q.split("_")[1]) for Q in header_line}
            commonQs = set(userQs).intersection(set(tsvQs))
            ln_commonQs = len(commonQs)
            ln_userQs = len(userQs)

            if not ln_commonQs:
                self.Logger.ERROR("invalid virtualQs range, none of them found in the pairwise TSV file")

            elif ln_commonQs == ln_userQs:
                self.userQs = userQs

            else:
                unfound = userQs - tsvQs
                self.Logger.WARNING(f"The following Qs '{unfound}' couldn't be found in the pairwise TSV, processing the others..")
                self.userQs = list(commonQs)

            self.kSize = max(tsvQs)
            self.tsv_get_namesmap()




    def __init__(self,logger_obj, pairwise_file_type, pairwise_file, userQs, cut_off_threshold):
        self.Logger = logger_obj
        self.cut_off_threshold = cut_off_threshold
        self.pairwise_file_type = pairwise_file_type
        self.file_name = pairwise_file

        if pairwise_file_type == "sqlite":
            self.sqlite_constructor(pairwise_file, userQs)
        elif pairwise_file_type == "tsv":
            self.tsv_constructor(pairwise_file, userQs)

    def ids_to_names(self, cluster):
        new_cluster = []
        for id in cluster:
            new_cluster.append(self.names_map[int(id)])

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

    def tsv_get_namesmap(self):
        if len(os.path.dirname(self.file_name)) > 1:
            namesMap_file = glob.glob(os.path.dirname(self.file_name) + "/*map")[0].split(".")[0] + ".namesMap"
        else:
            namesMap_file = glob.glob("*map")[0].split(".")[0] + ".namesMap"

        with open(namesMap_file, 'r') as namesMap:
            next(namesMap) #skip the header
            for row in namesMap:
                row = row.strip().split()
                self.names_map[int(row[0])] = row[1]


    def build_graph(self):
        if self.pairwise_file_type == "tsv":
            self.tsv_build_graph()
        elif self.pairwise_file_type == "sqlite":
            self.sqlite_build_graph()

    def tsv_build_graph(self):
        user_Qs = list(self.userQs)
        user_Qs.sort(reverse=False)
        indeces = {i: idx + 3 for idx, i in enumerate(user_Qs)}
        _temp = indeces.copy()
        for q, idx in _temp.items():
            if int(q) < 5:
                del indeces[q]

        inv_indeces = {v: k for k, v in indeces.items()}
        indeces_values = sorted(inv_indeces.keys(), reverse= True)

        with open(self.file_name , 'r') as pairwise_tsv:
            next(pairwise_tsv) #skip header
            for row in pairwise_tsv:
                row = row.strip().split()
                seq1 = int(row[0])
                seq2 = int(row[1])
                min_kmers = int(row[2])
                similarity = 0.0
                Q_prev = 0

                for idx in indeces_values:
                    Q_val = inv_indeces[idx]
                    Q_curr = int(row[idx])
                    sim = ((Q_curr - Q_prev) * (Q_val / self.kSize)) / min_kmers
                    sim = abs(sim)
                    #print(f"TSV: sim({seq1},{seq2}) = (({Q_curr} - {Q_prev}) * ({Q_val} / {self.kSize})) / {min_kmers} = {sim}")
                    similarity += abs(sim)
                    Q_prev = Q_curr

                if similarity < self.cut_off_threshold:
                    continue

                self.source.append(seq1)
                self.target.append(seq2)

            for i in range(1, len(self.names_map) + 1, 1):
                self.source.append(i)
                self.target.append(i)

    def sqlite_build_graph(self):
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
                sim = ((Q_curr - Q_prev) * (Q_val / self.kSize)) / min_kmers
                sim = abs(sim)
                #print(f"SQLITE: sim({seq1},{seq2}) = (({Q_curr} - {Q_prev}) * ({Q_val} / {self.kSize})) / {min_kmers} = {sim}")
                similarity += sim
                Q_prev = Q_curr

            if similarity < self.cut_off_threshold:
                continue

            self.source.append(seq1)
            self.target.append(seq2)


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
        kCluster_file_name = f"clusters_{self.cut_off_threshold:.2f}%_"
        kCluster_file_name += os.path.basename(self.file_name).split(".")[0]
        kCluster_file_name += ".tsv"

        with open(kCluster_file_name, 'w') as kClusters:
            kClusters.write("kClust_id\tseqs_ids\n")
            for cluster_id, (k, v) in enumerate(self.components.items(), 1):
                kClusters.write(f"{cluster_id}\t{','.join(self.ids_to_names(v))}\n")

        print(f"Total Number Of Clusters: {cluster_id}", file = sys.stderr)


@cli.command(name = "cluster", help_priority=3)
@click.option('-m','--min-q', required=False, type=int, default = None, help="minimum virtualQ")
@click.option('-M','--max-q', required=False, type=int, default = None, help="maximum virtualQ")
@click.option('-s','--step-q', required=False, type=int, default = None, help="virtualQs range step")
@click.option('--qs', required=False, type=click.STRING, default = None, help="comma separated virtualQs, ex: '25,30,31'")
@click.option('-c','--cutoff', required=False, type=click.FloatRange(0, 1, clamp=False), default = 0.0, show_default=True, help="cluster sequences with (similarity > cutoff)")
@click.option('-d', '--db', required=False, type=click.Path(exists=True), default=None,  help="sqlite database file")
@click.option('-t', '--tsv', required=False, type=click.Path(exists=True), default=None, help="sqlite database file")
@click.pass_context
def main(ctx, min_q, max_q, step_q, qs, db, tsv, cutoff):
    """Sequences clustering regarding user-selected virtualQs."""
    scanQs = map(bool, [min_q, max_q, step_q])

    if True in scanQs and False in scanQs:
        ctx.obj.ERROR("Please complete the virtualQs range")

    elif len(qs):
        userQs = list(map(int, qs.split(",")))

    elif not(min_q and max_q and step_q):
        ctx.obj.WARNING("processing all virtualQs in the sqlite database")
        userQs = [-1]

    else:
        userQs = [Q for Q in range(min_q, max_q + 1, step_q)]

    pairwise_file = ""

    if db == None and tsv:
        file_type = "tsv"
        pairwise_file = tsv
    elif tsv == None and db:
        file_type = "sqlite"
        pairwise_file = db
    elif tsv and db:
        ctx.obj.ERROR("Can't select both sqlite and tsv.")
    elif tsv == None and db == None:
        ctx.obj.ERROR("You must select pairwise matrix for clustering.")


    # kCl = kClusters(logger_obj= ctx.obj, sqlite_file = db, userQs = userQs, cut_off_threshold = cutoff)
    kCl = kClusters(logger_obj=ctx.obj, pairwise_file_type=file_type, pairwise_file = pairwise_file, userQs = userQs, cut_off_threshold = cutoff)
    kCl.build_graph()
    kCl.clustering()
    kCl.export_kCluster()