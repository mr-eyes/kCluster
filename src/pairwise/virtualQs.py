#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
import os
import click
import subprocess
import errno
from src.click_context import cli
from src.pairwise.virtualQs_class import  virtualQs

def prepare_params(ksize: int, minQ: int, maxQ: int, stepQ: int):
    """Verify virtualQs

    Args:
        ksize (int) : index kmersize
        minQ (int): minimum virtual Q (>= 1).
        maxQ (int): minimum virtual Q (<= kmer size).
        stepQ (int): virtual Q step (< maxQ)

    """

    __maxQ, __minQ, __stepQ = 0, 0, 0

    if maxQ > ksize or maxQ == -1:
        print("[INFO] auto reinitializing Q with kSize %d" % (ksize), file=sys.stderr)
        __maxQ = ksize

    elif maxQ is 0:
        __maxQ = ksize

    else:
        __maxQ = maxQ

    if (minQ < 5):
        print("[WARNING] minQ shouldn't be less than 5, auto reinitializing minQ to 5", file=sys.stderr)
        __minQ = 5
    elif minQ > __maxQ:
        print("[WARNING] minQ shouldn't exceed the maxQ, auto reinitializing minQ to maxQ", file=sys.stderr)
        __minQ = __maxQ
    else:
        __minQ = minQ

    if (stepQ < 1):
        print("auto resetting Q step to 2", file=sys.stderr)
        __minQ = 1
    else:
        __stepQ = stepQ

    all_Qs = []
    for i in range(__minQ, __maxQ + 1, __stepQ):
        all_Qs.append(i)

    return ",".join(map(str, all_Qs))

def run_command(command):
    try:
        devnull = open(os.devnull)
        subprocess.Popen(command, stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

@cli.command(name = "pairwise", help_priority=2)
@click.option('-m','--min-q', 'min_q', required=False, type=click.INT, default = 5, show_default=True, help="minimum virtualQ")
@click.option('-M','--max-q', 'max_q', required=False, type=click.INT, default = -1, help="maximum virtualQ")
@click.option('-s','--step-q', 'step_q', required=False, type=click.INT, default = 2, show_default=True,  help="virtualQs range step")
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="kProcessor index file prefix")
@click.option('-o', '--output-prefix', required=False, type=click.STRING, default=None, help="virtualQs output file(s) prefix")
@click.option('--force','force_write', is_flag=True, help="Overwrite the already proessed virtualQs")
@click.option('--backup', is_flag=True, help="Back up old virtualQs")
@click.option('--export-colors', required=False, type=click.Choice(['json', 'pickle']), default=None, help="export supercolors data [debugging purposes]")
@click.pass_context
def main(ctx, min_q, max_q, step_q, index_prefix, output_prefix, force_write, backup, export_colors):
    """
    Generating pairwise  matrices for single/multiple virtualQs.
    """
    if not os.path.isfile(index_prefix + ".map") and not os.path.isfile(index_prefix + ".mqf"):
        print(f"Index prefix {index_prefix} Does not exist!", file = sys.stderr)
        sys.exit(1)

    if not output_prefix:
        output_prefix = os.path.basename(index_prefix)

    ksize = int()

    with open(index_prefix + ".extra") as extra:
        ksize = int(next(extra).strip())

    Qs_list_str = prepare_params(ksize, min_q, max_q, step_q)

    

    ctx.obj.SUCCESS("All completed...")