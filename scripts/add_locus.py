#!/usr/bin/python
# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse
import os


def Main():

    intervals = []
    genes = []
    chromosomes = []

    with open(gtf_file_path, 'r') as gtf:
        for line in gtf:
            if '\tgene\t' not in line:
                continue

            fields = line.split()

            interval = (int(fields[3]), int(fields[4]))
            gene_id = (fields[9])[1:-2]
            chromosome = fields[0]

            intervals.append(interval)
            genes.append(gene_id)
            chromosomes.append(chromosome)

    loci = {}
    temp_locus = set()
    locus_index = 0
    largest_end = intervals[0][1]

    for v in range(1, len(intervals), 1):
        _prev_interval = intervals[v - 1]
        _curr_interval = intervals[v]
        _prev_chr = chromosomes[v - 1]
        _curr_chr = chromosomes[v]
        _prev_gene_id = genes[v - 1]
        _curr_gene_id = genes[v]

        if _prev_interval[1] > largest_end:
            largest_end = _prev_interval[1]

        if _curr_chr != _prev_chr:
            temp_locus.add(_prev_gene_id)

            for gene in temp_locus:
                loci[gene] = locus_index

            locus_index += 1
            temp_locus = set()
            largest_end = _curr_interval[1]

        elif _curr_interval[0] < largest_end:
            temp_locus.add(_prev_gene_id)

        else:

            temp_locus.add(_prev_gene_id)
            for gene in temp_locus:
                loci[gene] = locus_index

            locus_index += 1
            temp_locus = set()
            largest_end = _curr_interval[1]
            temp_locus.add(_curr_gene_id)

        if v == len(intervals) - 1:
            temp_locus.add(_curr_gene_id)
            for gene in temp_locus:
                loci[gene] = locus_index
            locus_index += 1

    with open(output_fasta, 'w') as output_handle:
        sequences = SeqIO.parse(open(original_fasta_file), 'fasta')
        for seq in sequences:
            gene_id = seq.id.split('|')[1]
            locus = 'Locus_' + str(loci[gene_id])
            seq.id = seq.description = seq.name + locus + '|'
            SeqIO.write(seq, output_handle, 'fasta')


if __name__ == '__main__':
    original_fasta_file = ''
    gtf_file_path = ''
    output_fasta = ''

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-i', action='store', dest='original_fasta_file', help='The Original fasta file')

    parser.add_argument('-g', action='store', dest='gtf_file',
                        help='Set GTF File path')

    parser.add_argument('-o', action='store', dest='output_fasta_file',
                        help='New fasta file name')

    args = parser.parse_args()

    if args.gtf_file:
        if os.path.isfile(args.gtf_file):
            gtf_file_path = args.gtf_file
        else:
            exit('File: ' + args.gtf_file + 'Does not exist!')

        if args.original_fasta_file:
            if os.path.isfile(args.original_fasta_file):
                original_fasta_file = args.original_fasta_file
            else:
                exit('File: ' + args.original_fasta_file
                     + 'Does not exist!')

        if args.output_fasta_file:
            output_fasta = args.output_fasta_file
        else:
            output_fasta = 'new_' + original_fasta_file
    else:

        parser.print_help()
        exit(0)

    Main()
