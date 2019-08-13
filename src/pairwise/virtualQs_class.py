#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import kProcessor as kp
from itertools import combinations
import json
import sys
import os
import pickle
import sqlite3


class virtualQs:
    """Holds the superColors and superColorsCount tables."""

    __kSize = None
    overwrite = False

    def __init__(self, logger_obj, index_prefix: str):
        """VirtualQs class constructor.

        Args:
            index_file_name (string): coloredKDataFrame index file.
            output_prefix (string): prefix for the output files(s)

        """
        self.Logger = logger_obj
        self.kf = kp.kDataFrame.load(index_prefix)
        self.__kSize = self.kf.ksize()

        self.index_prefix = index_prefix

        if self.__kSize is None:
            print("error loading the index", file=sys.stderr)
            sys.exit(1)

        # Convert colors to IDs
        self.color_to_ids = {}
        with open(index_prefix + "colors.intvectors", 'r') as colors:
            next(colors)  # skip the first line (Number of colors)
            for line in colors:
                values = list(map(int, line.strip().split()))
                self.color_to_ids[values[0]] = values[2:]

    def __mask(self, Q):
        """create a bit mask given kmer size and Q value."""
        return (~(-1 << Q * 2)) << (self.__kSize * 2 - Q * 2)

    def set_params(self, minQ: int, maxQ: int, stepQ: int):
        """virtualQs parameters setting.

        Args:
            minQ (int): minimum virtual Q (>= 1).
            maxQ (int): minimum virtual Q (<= kmer size).
            stepQ (int): virtual Q step (< maxQ)

        """

        if maxQ > self.__kSize or maxQ == -1:
            print(
                "[INFO] auto reinitializing Q with kSize %d" % (self.__kSize), file=sys.stderr)
            self.__maxQ = self.__kSize

        elif maxQ is 0:
            self.__maxQ = self.__kSize

        else:
            self.__maxQ = maxQ

        if (minQ < 5):
            print("[WARNING] minQ shouldn't be less than 5, auto reinitializing minQ to 5", file=sys.stderr)
            self.__minQ = 5
        elif minQ > self.__maxQ:
            print("[WARNING] minQ shouldn't exceed the maxQ, auto reinitializing minQ to maxQ", file=sys.stderr)
            self.__minQ = self.__maxQ
        else:
            self.__minQ = minQ

        if (stepQ < 1):
            print("auto resetting Q step to 2", file=sys.stderr)
            self.__minQ = 1
        else:
            self.__stepQ = stepQ

    def prepare(self):
        """
        prepare masks and initialize data structure for the generating the virtualQs
        """

        self.superColors = {}
        self.temp_superColors = {}
        self.superColorsCount = {}
        self.edges = {}
        self.seq_to_kmers_no = {}
        self.__masks = {}

        # Determine minQ and maxQ and get list of masks & superColorsDIct initialization
        # for Q in range(maxQ, minQ-1, -stepQ):
        for Q in self.mainQs:
            self.__masks[Q] = self.__mask(Q)
            self.superColors[Q] = {}
            self.superColorsCount[Q] = {}
            self.temp_superColors[Q] = []

        """
        Must consider Q=kSize for calculation of actual kmer size.
        Adding it manually and creating a mask for it. 
        """

        if self.kSize not in self.masks:
            self.__masks[self.kSize] = self.__mask(self.kSize)
            self.superColors[self.kSize] = {}
            self.superColorsCount[self.kSize] = {}
            self.temp_superColors[self.kSize] = []

    def set_mainQs(self, Qs):
        """
        set the main Qs values that will be processed.
        """
        self.mainQs = [int(Q.split("_")[1]) for Q in Qs]

    def calculate_kmers_number(self):
        """
        Calculating number of kmers of each sequences
        """

        for super_color, colors in self.superColors[self.kSize].items():
            tr_ids = list({i for c in colors for i in self.color_to_ids[c]})
            color_count = self.superColorsCount[self.kSize][super_color]

            # For loop to calculate number of kmers per seq_id
            for tr_id in tr_ids:
                if tr_id not in self.seq_to_kmers_no:
                    self.seq_to_kmers_no[tr_id] = color_count
                else:
                    self.seq_to_kmers_no[tr_id] += color_count

    def pairwise(self):
        """
        pairwise similarity matrix construction for all main Qs
        """

        ## TODO This function can be parallelized

        # Calculating number of kmers per each sequence
        self.calculate_kmers_number()

        # Calculating pairwise distance
        for Q in self.mainQs:
            for color, colors in self.superColors[Q].items():
                tr_ids = list({i for c in colors for i in self.color_to_ids[c]})
                color_count = self.superColorsCount[Q][color]

                for combination in combinations(tr_ids, 2):
                    _seq1 = combination[0]
                    _seq2 = combination[1]

                    seq_pair = tuple(sorted((_seq1, _seq2)))
                    if seq_pair in self.edges:
                        self.edges[seq_pair][Q] += color_count

                    else:
                        self.edges[seq_pair] = dict()

                        for _Q in self.mainQs:
                            self.edges[seq_pair][_Q] = 0

                        self.edges[seq_pair][Q] = color_count

    def sqlite_table_exists(self, table_name):
        """Check if the sqlite table exists
        Args:
            table_name (str): table name to be checked.
        Returns:
            boolean value, True for existance.
        """

        res = self.conn.execute(f"SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{table_name}'")
        return bool(res.fetchone()[0] == 1)

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

    def sqlite_getOldQs(self, duplicateQs):
        """
        fetch the Qs names for previously cached Qs and duplicate Qs that will be replaced.

        Args:
            duplicateQs (list): list of duplicate Qs between user predefined Qs and Qs stored in the database.

        Returns:
            List of new Qs values.

        """
        sqliteQs = self.sqlite_getQs()
        cachedQs = dict()
        result = dict()

        cached_before = [Q.replace("cached_", "") for Q in sqliteQs if "cached" in Q]

        for Q in cached_before:
            _Q = Q.split("_")
            key = f"{_Q[0]}_{_Q[1]}"

            if key not in duplicateQs:
                continue

            count = int(_Q[-1])

            if key not in cachedQs:
                cachedQs[key] = count
            else:
                if count > cachedQs[key]:
                    cachedQs[key] = count

        for Q in duplicateQs:
            if Q in cachedQs:
                result[Q] = cachedQs[Q] + 1
            else:
                result[Q] = 0

        return result

    def _sqlite_insert(self):
        for seq_pair, Qs in self.edges.items():
            seq1 = seq_pair[0]
            seq2 = seq_pair[1]
            min_kmers = min(self.seq_to_kmers_no[seq1], self.seq_to_kmers_no[seq2])
            Qs_cols = ", ".join([f"Q_{Q}" for Q in Qs])
            Qs_vals = ", ".join([f"{Q}" for Q in Qs.values()])
            self.conn.execute(
                f"INSERT INTO virtualQs (seq1, seq2, min_kmers, {Qs_cols}) VALUES ({seq1}, {seq2}, {min_kmers}, {Qs_vals})")

    def _sqlite_update(self):
        for seq_pair, Qs in self.edges.items():
            seq1 = seq_pair[0]
            seq2 = seq_pair[1]
            new_values = [f"Q_{Q}={val}" for Q, val in Qs.items()]
            new_values = ", ".join(new_values)
            self.conn.execute(f"UPDATE virtualQs SET {new_values} WHERE seq1={seq1} AND seq2={seq2}")

    def sqlite_dropQs(self, duplicateQs):
        # Until now, sqlite does not support dropping columns.
        # So, I cache the Qs by changing their names.

        cachedQs = self.sqlite_getOldQs(duplicateQs)
        for Q in duplicateQs:
            new_value = f"cached_{Q}_{cachedQs[Q]}"
            self.conn.execute(f"ALTER TABLE virtualQs RENAME COLUMN {Q} TO {new_value};")

        self.conn.commit()

    def sqlite_createQs(self, Qs):
        # Qs is a list of {Qn,Qn+1, ...}

        _Qs = [int(Q.split("_")[1]) for Q in Qs]
        _Qs = sorted(_Qs)

        for Q in _Qs:
            self.conn.execute(f"ALTER TABLE virtualQs ADD COLUMN Q_{Q} INT DEFAULT 0 NOT NULL;")

        self.conn.commit()

    def sqlite_create_tables(self):
        """Create virtualQs and meta_info tables"""

        # Create virtualQs table
        self.conn.execute('''CREATE TABLE virtualQs
        (ID INTEGER PRIMARY KEY AUTOINCREMENT,
         seq1            INT     NOT NULL,
         seq2            INT     NOT NULL,
         min_kmers       INT     NOT NULL,
         UNIQUE (seq1, seq2));''')

        # Create meta information table
        self.conn.execute('''CREATE TABLE meta_info
                            (ID INTEGER PRIMARY KEY AUTOINCREMENT,
                            key              TEXT     NOT NULL,
                            value            INT     NOT NULL);''')

        # Just insert the kSize
        self.conn.execute(f"INSERT OR IGNORE INTO meta_info (key, value) VALUES ('kSize', {self.kSize})")

    def sqlite_import_namesmap(self):
        """Import kProcessor namesmap file in the sqlite database"""

        self.conn.execute('''CREATE TABLE namesmap
                            (ID INTEGER PRIMARY KEY AUTOINCREMENT,
                             seq_id     INT     NOT NULL,
                             seq_name   TEXT    NOT NULL);''')

        with open(self.index_prefix + ".namesMap", 'r') as namesmap:
            next(namesmap)
            for record in namesmap:
                record = record.strip().split(" ")
                seq_id = record[0]
                seq_name = record[1]

                self.conn.execute(f"INSERT INTO namesmap (seq_id, seq_name) VALUES ({seq_id}, '{seq_name}')")

    def sqlite_initiate(self, prefix, force_write=False, backup=False):
        """Initialize the Sqlite database for the processing.
        Args:
            prefix (str): database file name prefix.
            force_write (bool): boolean flag to allow or disallow force write.
            backup (bool): boolean flag to allow ir disallow backing up overwritten Qs.
        """
        self.force_write = force_write
        self.backup = backup

        database_path = prefix + "_kCluster.sqlite"
        db_exist = os.path.isfile(database_path)
        self.conn = sqlite3.connect(database_path)

        # Get user Qs
        self.cursor = self.conn.cursor()
        userQs = {f"Q_{i}" for i in range(self._virtualQs__minQ, self._virtualQs__maxQ + 1, self._virtualQs__stepQ)}

        # Database file exists
        if db_exist:

            # If the table does not exist
            if not self.sqlite_table_exists("virtualQs"):

                # Setting the main insertion function to insert
                self.sqlite_insert = self._sqlite_insert

                # IF table does not exist, create it.
                self.sqlite_create_tables()

                # Import namesMap file
                self.sqlite_import_namesmap()


            # Table exists, validate the columns
            else:
                # Scan the table
                sqliteQs = self.sqlite_getQs()

                # There are no Qs.
                if not sqliteQs:
                    print("sqlite db file is corrupted.", file=sys.stderr)
                    sys.exit(1)

                else:
                    # Qs found

                    # Setting the main insertion function to update
                    self.sqlite_insert = self._sqlite_update

                    # User has enabled force write to over-write the Qs.
                    if self.force_write:
                        missingQs = userQs - sqliteQs
                        duplicated = userQs.intersection(sqliteQs)
                        self.sqlite_dropQs(duplicated)
                        self.sqlite_createQs(userQs)
                        self.set_mainQs(userQs)
                        if len(duplicated):
                            self.overwrite = True


                    else:
                        missingQs = userQs - sqliteQs
                        self.sqlite_createQs(missingQs)
                        self.set_mainQs(missingQs)


        # Database does not exist at all
        else:
            # print("database does not exist, creating..")

            # setting the main insertion function to insert
            self.sqlite_insert = self._sqlite_insert
            # IF table does not exist, create it.

            self.sqlite_create_tables()

            # Import namesMap file
            self.sqlite_import_namesmap()

            self.sqlite_createQs(userQs)
            self.set_mainQs(userQs)

        self.conn.commit()
        self.prepare()

    def export_pairwise(self):
        """pairwise similarity matrix exporting as TSV file

        Args:
            prefix (str): exported file name prefix.
            Q (int): Q value for the pairwise matrix construction.

        """

        self.sqlite_insert()

        if not self.backup and self.overwrite:
            Qs = [Q for Q in self.sqlite_getQs() if Q[0] == "Q"]
            gold_names = ['ID', 'seq1', 'seq2', 'min_kmers']
            preserved_Qs = gold_names + Qs
            Qs_columns = ", ".join(preserved_Qs)

            self.conn.execute(f"CREATE TABLE virtualQs_backup AS SELECT {Qs_columns} FROM virtualQs;")
            self.conn.execute("DROP TABLE virtualQs;")
            self.conn.execute("ALTER TABLE virtualQs_backup RENAME TO virtualQs;")

        self.conn.commit()

    def export_superColors(self, prefix, Q, method="json"):
        """superColors table exporting

        Retrieves Q value needs to be exported and output file format

        Args:
            prefix: exported file name prefix.
            Q: Q value to be extracted from the superColors tables.
            method: specify the output file format pickle or json.
        """

        if Q not in self.superColors and Q not in self.superColorsCount:
            print("virtualQ: {} does not exist".format(Q), file=sys.stderr)
            sys.exit(1)

        if method == "pickle":
            suffix = ".pickle"
        elif method == "json":
            suffix = ".json"
        else:
            print("export only in [pickle,json]", file=sys.stderr)
            sys.exit(1)

        virtualQs_file_name = prefix + "_K" + str(self.kSize) + "_Q" + str(Q) + suffix
        virtualQs_count_file_name = prefix + "_K" + str(self.kSize) + "_Q" + str(Q) + "_counts" + suffix

        if method == "pickle":
            with open(virtualQs_file_name, "wb") as f:
                pickle.dump(self.superColors[Q], f, pickle.HIGHEST_PROTOCOL)

            with open(virtualQs_count_file_name, "wb") as f:
                pickle.dump(
                    self.superColorsCount[Q], f, pickle.HIGHEST_PROTOCOL)

        elif method == "json":
            with open(virtualQs_file_name, "w") as f:
                f.write(json.dumps(
                    self.superColors[Q], sort_keys=True, indent=4, separators=(',', ': ')))

            with open(virtualQs_count_file_name, "w") as f:
                f.write(json.dumps(
                    self.superColorsCount[Q], sort_keys=True, indent=4, separators=(',', ': ')))

    # Take list of colors, sort it, and create a return a hash value
    @staticmethod
    def create_super_color(colors):
        return hash(tuple(sorted(set(colors))))

    @property
    def get_params(self):
        return {"minQ": self.__minQ, "maxQ": self.__maxQ, "stepQ": self.__stepQ}

    @property
    def kSize(self):
        return (self.__kSize)

    @property
    def masks(self):
        return self.__masks

    @staticmethod
    def int_to_str(kmer, kSize):
        _map = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
        kmer_str = ""
        for i in range(kSize, 0, -1):
            base = (kmer >> (i * 2 - 2)) & 3
            ch = _map[base]
            kmer_str += ch

        return kmer_str
