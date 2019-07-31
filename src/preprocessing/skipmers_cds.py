#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import click
import subprocess
import errno
import glob
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
            if max_len - _len < max_len * (self.diff_threshold / 100):
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

    @staticmethod
    def tool_exist(name):
        try:
            devnull = open(os.devnull)
            subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
        except OSError as e:
            if e.errno == errno.ENOENT:
                return False
        return True

    @staticmethod
    def run_command(command):
        try:
            devnull = open(os.devnull)
            subprocess.Popen(command, stdout=devnull, stderr=devnull).communicate()
        except OSError as e:
            if e.errno == errno.ENOENT:
                return False
        return True


@cli.command(name="preprocess_cds", help_priority=9)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-t', '--diff-threshold', "diff_threshold", required=False, default=10, type=click.INT,
              help="minimum length difference %")
@click.option('-l', '--cds-len', "cds_length", required=False, default=50, type=click.INT, help="Minimum CDS length")
@click.option('-s', '--strand', 'strand', is_flag=True, required=False, help="Examine only the top strand")
@click.pass_context
def preprocess_cds(ctx, fasta_file, diff_threshold, cds_length, strand):
    '''Preprocess protein coding transcript to extract CDS'''

    if not Preprocess.tool_exist("TransDecoder.LongOrfs"):
        ctx.obj.ERROR("TransDecoder is not installed")
    else:
        ctx.obj.INFO("Processing..")

    """Construct Transdecoder Commands"""

    commands = list()
    cmd = f"TransDecoder.LongOrfs -m {cds_length}"

    if strand:
        cmd += " -S"

    cmd += f" -t {fasta_file}"
    cmd = cmd.split()
    Preprocess.run_command(cmd)

    cmd = "rm -rf".split() + glob.glob("*__checkpoints*")
    Preprocess.run_command(cmd)

    cmd = "rm -rf".split() + glob.glob("*cmds")
    Preprocess.run_command(cmd)


    cmd = "mkdir -p transdecoder_cds".split()
    Preprocess.run_command(cmd)


    cmd = f"mv {os.path.basename(fasta_file)}.transdecoder_dir/longest_orfs.cds ./transdecoder_cds/".split()
    Preprocess.run_command(cmd)


    new_fasta_file = f"./transdecoder_cds/cds_{os.path.basename(fasta_file)}"
    cmd = f"mv ./transdecoder_cds/longest_orfs.cds {new_fasta_file}".split()
    Preprocess.run_command(cmd)


    cmd = f"rm -rf".split() + glob.glob("*transdecoder_dir*")
    Preprocess.run_command(cmd)



    """END Transdecoder Commands Construction"""

    CDS = Preprocess(new_fasta_file, diff_threshold)
    CDS.parse()
    ctx.obj.INFO("Writing the names file...")
    output_names_file = "transdecoder_cds/cds_" + fasta_file + ".names"
    with open(output_names_file, 'w') as output:
        for groupName, headers in CDS.namesfile.items():
            for header in headers:
                output.write(f"{header[1:]}\t{groupName[1:]}\n")

    ctx.obj.SUCCESS("Completed..")
