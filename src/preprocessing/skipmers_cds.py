#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import click
from src.click_context import cli


class Preprocess:
    namesfile = dict()
    part_to_len = dict()
    temp_headers = dict()
    diff_threshold = 10
    cds_file = ""

    def __init__(self, fasta_file, threshold):
        self.cds_file = fasta_file
        self.diff_threshold = threshold

    def analyse_line(self, line):
        splitted2 = line.split()
        _name, _type, _len, _range = splitted2[0], splitted2[1], int(splitted2[2].split(":")[-1]), splitted2[3]
        _part = splitted2[0].split(".")[-1].replace("p", "")
        _group = "".join(splitted2[0].split(".")[:-1])
        return {"name": _name, "part": _part, "len": _len, "group": _group}

    def select_parts(self):

        first_part = next(iter(self.part_to_len.keys()))
        parts = [first_part]
        max_len = self.part_to_len.pop(first_part)

        for _part, _len in self.part_to_len.items():
            if max_len - _len < self.diff_threshold:
                parts.append(_part)
                max_len = _len
            else:
                break

        names = [self.temp_headers[part] for part in parts]
        return names

    def parse(self):
        with open(self.cds_file, 'r') as cds:
            line = next(cds).strip()
            splitted = self.analyse_line(line)
            prev_group = splitted["group"]
            self.part_to_len[splitted["part"]] = splitted["len"]
            self.temp_headers[splitted["part"]] = line

            for line in cds:
                line = line.strip()
                if line[0] != ">":
                    continue

                splitted = self.analyse_line(line)
                if splitted["group"] == prev_group:
                    self.part_to_len[splitted["part"]] = splitted["len"]
                    self.temp_headers[splitted["part"]] = line

                else:
                    parts = self.select_parts()
                    self.namesfile[prev_group] = parts
                    prev_group = splitted["group"]
                    self.part_to_len.clear()
                    self.part_to_len[splitted["part"]] = splitted["len"]
                    self.temp_headers[splitted["part"]] = line

            else:
                parts = self.select_parts()
                self.namesfile[prev_group] = parts
                self.part_to_len.clear()
                self.part_to_len[splitted["part"]] = splitted["len"]
                self.temp_headers.clear()


@cli.command(name="preprocess_cds", help_priority=9)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-t', '--diff-threshold', "diff_threshold", required=False, default=10, type=click.INT,
              help="minimum length difference")
@click.option('-l', '--cds-len', "cds_length", required=False, default=50, type=click.INT, help="Minimum CDS length")
@click.pass_context
def preprocess_cds(ctx, fasta_file, diff_threshold, cds_length):
    '''Preprocess protein coding transcript to extract CDS'''

    output_names_file = fasta_file + ".names"

    CDS = Preprocess(fasta_file, diff_threshold)
    CDS.parse()

    with open(output_names_file, 'w') as output:
        for groupName, headers in CDS.namesfile.items():
            for header in headers:
                output.write(f"{header}\t{groupName}\n")