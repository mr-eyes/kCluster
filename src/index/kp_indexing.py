import sys
import os
import click
from src.click_context import cli
#from src.lib.custom_logger import Logger

try:
    import kProcessor as kp
except ImportError:
    click.secho("kProcessor package could not found.", fg="red", bold=True , file = sys.stderr) 

class Index:

    def __init__(self,logger_obj, fasta_file, names_file):
        self.Logger = logger_obj
        self.fasta_file = fasta_file
        self.names_file = names_file

    def validate_names(self):
        '''validate names file for indexing'''

        self.Logger.INFO("validating names file..")

        with open(self.names_file) as names:
            for i, line in enumerate(names, 1):
                if len(line.strip().split("\t")) != 2:
                    self.Logger.ERROR(f"invalid names line detected at L{i}: '{line.strip()}'")

    def index(self, mode, params, noncanonical = False):
        """
        peform indexing with given kSize
        """

        self.Logger.INFO(f"Indexing by {mode} with params: {params}")

        try:
            # self.idx = kp.kDataFrameMQF(params["k_size"], 28, 3)
            self.idx = kp.kDataFrameMAP(params["k_size"])

            if noncanonical:
                self.Logger.INFO("nonCanonical mode..")
                kp.kmerDecoder_setHashing(self.idx, 2, False)

            self.ckf = kp.index(self.idx, {"mode": 1}, self.fasta_file, 100, self.names_file)
            self.Logger.SUCCESS("Indexing Completed")
        except Exception as e:
            print(e)
            self.Logger.ERROR("Indexing failed")

    def write_to_disk(self, output_prefix):
        """save index file to disk"""

        try:
            self.ckf.save(output_prefix)
        except:
            self.Logger.ERROR("saving index to disk failed")


def validate_kSize(ctx, param, value):
    if not value % 2:
        raise click.BadParameter(f"kmer size: {value} is even, please enter an odd value.")
    return value


def validate_ORF(ctx, param, value):
    if value not in [0, 1, 2, 3]:
        raise click.BadParameter(f"Please select ORF 1, 2 or 3")
    return value


"""
KMERS
"""


@cli.command(name="index_kmers", help_priority=1)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "kSize", callback=validate_kSize, required=True,
              type=click.IntRange(7, 31, clamp=False), help="kmer size")
@click.option('--noncanonical', is_flag=True)
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def kmers(ctx, fasta_file, names_file, kSize, output_prefix, noncanonical):
    '''FASTA file indexing by Kmers'''

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    mode = "kmers"
    params = {"k_size": kSize}

    idx = Index(logger_obj=ctx.obj, fasta_file=fasta_file, names_file=names_file)
    idx.validate_names()
    idx.index(mode, params, noncanonical)
    idx.write_to_disk(output_prefix)


"""
SKIPMERS
"""


@cli.command(name="index_skipmers", help_priority=9)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "skipmers_kSize", required=True, type=click.INT, help="kmer size")
@click.option('-m', '--cycle-bases', "skipmers_m", required=True, type=click.INT, help="used bases per cycle")
@click.option('-n', '--cycle-length', "skipmers_n", required=True, type=click.INT, help="kmer size(cycle length")
@click.option('--orf', "orf", required=False, type=click.INT, callback=validate_ORF, default=0,
              help="select ORF <1,2,3>")
@click.option('-o', '--output', "output_prefix", required=False, default=None, help="index output file prefix")
@click.pass_context
def skipmers(ctx, fasta_file, names_file, skipmers_kSize, skipmers_m, skipmers_n, orf, output_prefix):
    '''FASTA file indexing by Skipmers'''

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    mode = "skipmers"
    params = {"k_size": skipmers_kSize, "m": skipmers_m, "n": skipmers_n, "orf": orf}

    idx = Index(logger_obj=ctx.obj, fasta_file=fasta_file, names_file=names_file)
    idx.validate_names()
    idx.index(mode, params)
    idx.write_to_disk(output_prefix)


""" Return back later"""
# @cli.command(name = "index", help_priority=1)
# @click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
# @click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
# @click.option('-k', '--kmer-size', "kSize", callback=validate_kSize, required=True, type=click.IntRange(15, 31, clamp=False), help = "kmer size" )
# @click.option('-q', '--mqf-q', "mqf_q", required=True, type=click.INT, default=27 , help = "MQF Q Value" )
# @click.option('-o', '--output', "output_prefix", required=False, default=None, help = "index output file prefix")
# @click.pass_context
# def main(ctx, fasta_file, names_file, kSize, mqf_q, output_prefix):
#     '''FASTA file indexing'''
#
#     if not output_prefix:
#         output_prefix = os.path.basename(fasta_file)
#         output_prefix = os.path.splitext(output_prefix)[0]
#         output_prefix = "idx" + "_" + output_prefix
#
#     idx = Index(logger_obj = ctx.obj,fasta_file =  fasta_file,names_file = names_file)
#     idx.validate_names()
#     idx.index(kSize, mqf_q)
#     idx.write_to_disk(output_prefix)