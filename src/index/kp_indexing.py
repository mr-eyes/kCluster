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

    def index(self, kSize, mqf_q):
        """
        peform indexing with given kSize
        """

        self.Logger.INFO(f"kSize:{kSize}, Q:{mqf_q}")
        self.Logger.INFO("Indexing..")

        try:
            KD = kp.initialize_kmerDecoder(self.fasta_file, 1000, "kmers" , {"k_size" : kSize})
            self.idx = kp.kDataFrameMAP(kSize)
            kp.index(KD, self.names_file, self.idx)
            self.Logger.SUCCESS("Indexing Completed")
        except Exception as e:
            print(e)
            self.Logger.ERROR("Indexing failed")

    def write_to_disk(self, output_prefix):
        '''save index file to disk'''

        try:
            self.idx.save(output_prefix)
        except:
            self.Logger.ERROR("saving index to disk failed")


def validate_kSize(ctx, param, value):
    if not value % 2:
        raise click.BadParameter(f"kmer size: {value} is even, please enter an odd value.")
    return value

@cli.command(name = "index", help_priority=1)
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "kSize", callback=validate_kSize, required=True, type=click.IntRange(15, 31, clamp=False), help = "kmer size" )
@click.option('-q', '--mqf-q', "mqf_q", required=True, type=click.INT, default=27 , help = "MQF Q Value" )
@click.option('-o', '--output', "output_prefix", required=False, default=None, help = "index output file prefix")
@click.pass_context
def main(ctx, fasta_file, names_file, kSize, mqf_q, output_prefix):
    '''FASTA file indexing'''  

    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    idx = Index(logger_obj = ctx.obj,fasta_file =  fasta_file,names_file = names_file)
    idx.validate_names()
    idx.index(kSize, mqf_q)
    idx.write_to_disk(output_prefix)