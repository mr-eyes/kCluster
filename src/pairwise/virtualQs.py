#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
import os
import click
from src.click_context import cli
from src.pairwise.virtualQs_class import  virtualQs

# TODO use logging instead of normal prints.
# TODO check best practices for sqlite for enhancements
# TODO refactor the whole code to multiple classes
# TODO add validation callbacks for click


@cli.command(name = "pairwise", help_priority=2)
@click.option('-m','--min-q', 'min_q', required=False, type=click.INT, default = 5, show_default=True, help="minimum virtualQ")
@click.option('-M','--max-q', 'max_q', required=False, type=click.INT, default = -1, help="maximum virtualQ")
@click.option('-s','--step-q', 'step_q', required=False, type=click.INT, default = 2, show_default=True,  help="virtualQs range step")
@click.option('-i', '--index-prefix', required=True, type=click.STRING, help="kProcessor index file prefix")
@click.option('-o', '--output-prefix', required=False, type=click.STRING, default=None, help="virtualQs output file(s) prefix")
@click.option('--force','force_write', is_flag=True, help="Overwrite the already proessed virtualQs")
@click.option('--backup', is_flag=True, help="Back up old virtualQs")
# @click.option('--simple', 'simple_output', is_flag=True, required=False, help="export in a tsv output [seq1,seq2,shared] no virtualQs")
@click.option('--export-colors', required=False, type=click.Choice(['json', 'pickle']), default=None, help="export supercolors data [debugging purposes]")
@click.pass_context
def main(ctx, min_q, max_q, step_q, index_prefix, output_prefix, force_write, backup, export_colors):
    """
    Generating pairwise  matrices for single/multiple virtualQs.
    """

    for suffix in [".map", ".mqf", ".phmap"]:
        if os.path.isfile(index_prefix + suffix):
            break
    else:
        print(f"Index prefix {index_prefix} Does not exist!", file=sys.stderr)
        sys.exit(1)

    if not output_prefix:
        output_prefix = os.path.basename(index_prefix)

    ctx.obj.INFO("Loading the index...")
    VQ = virtualQs(logger_obj = ctx.obj, index_prefix=index_prefix)
    VQ.set_params(minQ=min_q, maxQ=max_q, stepQ=step_q)
    VQ.sqlite_initiate(output_prefix, force_write, backup)

    it = VQ.kf.begin()
    prev_kmer = it.getHashedKmer()
    prev_kmer_color = it.getCount()

    ctx.obj.INFO("Processing...")
    # Iterate over all kmers.
    while it != VQ.kf.end():
        it.next()
        curr_kmer = it.getHashedKmer()
        curr_kmer_color = it.getCount()

        # Apply XOR to kmer1 and kmer2 (single time per iteration)
        xor = prev_kmer ^ curr_kmer

        # Apply all masks with all Qs
        for Q, MASK in VQ.masks.items():

            # True if there's match, False if not
            matched = not bool(xor & MASK)

            if matched:
                VQ.temp_superColors[Q] += [prev_kmer_color, curr_kmer_color]

            else:
                VQ.temp_superColors[Q].append(prev_kmer_color)
                super_color_id = VQ.create_super_color(VQ.temp_superColors[Q])

                # print("Matching Q%d %s & %s | FALSE | [prevC:%d, currC=%d]" % (Q, VQ.int_to_str(
                #     prev_kmer, VQ.kSize), VQ.int_to_str(curr_kmer, VQ.kSize),  prev_kmer_color, curr_kmer_color))


                # Check if the superColor already exist
                # If yes: increment the count to one
                # If No:  Insert the new superColor and set the count to 1
                if super_color_id not in VQ.superColors[Q]:
                    VQ.superColors[Q][super_color_id] = list(set(VQ.temp_superColors[Q]))
                    VQ.superColorsCount[Q][super_color_id] = 1

                else:
                    # IF the supercolor already exist, just increment it
                    VQ.superColorsCount[Q][super_color_id] += 1


                VQ.temp_superColors[Q] = [curr_kmer_color]

        prev_kmer = curr_kmer
        prev_kmer_color = curr_kmer_color

    # If the last iteration got a match, push it to the superColors

    for Q, colors in VQ.temp_superColors.items():
        colors.remove(curr_kmer_color)
        if not len(colors):
            continue
        super_color_id = VQ.create_super_color(colors)

        if super_color_id not in VQ.superColors[Q]:

            VQ.superColors[Q][super_color_id] = list(set(colors))
            VQ.superColorsCount[Q][super_color_id] = 1

        else:
            # IF the supercolor already exist, just increment it
            VQ.superColorsCount[Q][super_color_id] += 1


    # Save all Qs to files.
    if export_colors:
        _dir_name = os.path.dirname(output_prefix)
        if not os.path.isdir(_dir_name):
            os.mkdir(_dir_name)

        for Q in VQ.mainQs:
            VQ.export_superColors(output_prefix, Q, export_colors)

    ctx.obj.INFO("Constructing the pairwise matrix...")
    # Construct pairwise matrices    
    VQ.pairwise()

    ctx.obj.INFO("Saving results...")
    # export pairwise matrices
    VQ.export_pairwise()
    ctx.obj.SUCCESS("All completed...")