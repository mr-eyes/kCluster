import click
import sys
import os

try:
    from kProcessor import index as kp_index
except ImportError:
    click.secho("kProcessor package could not found.", fg="red", bold=True , file = sys.stderr) 

class Index:

    def __init__(self, fasta_file, names_file):
        self.fasta_file = fasta_file
        self.names_file = names_file

    def validate_names(self):
        '''validate names file for indexing'''

        click.secho("[INFO] validating names file..", fg="green", bold=True, file = sys.stderr)
        with open(self.names_file) as names:
            for i, line in enumerate(names, 1):
                if len(line.strip().split("\t")) != 2:
                    click.secho(f"[ERROR] invalid names line detected at L{i}: '{line.strip()}'", fg="red", bold=True, file = sys.stderr)
                    sys.exit(1)

    def index(self, kSize):
        """
        peform indexing with given kSize
        """

        click.secho("[INFO] Indexing..", fg="green", bold=True, file = sys.stderr)

        try:
            self.idx = kp_index(self.fasta_file, self.names_file, kSize)
        except:
            click.secho("[ERROR] Indexing failed", fg="red", bold=True, file = sys.stderr)
            sys.exit(1)

    def write_to_disk(self, output_prefix):
        '''save index file to disk'''

        try:
            self.idx.save(output_prefix)
        except:
            click.secho("[ERROR] saving index to disk failed", fg="red", bold=True, file = sys.stderr)
            sys.exit(1)


@click.command(name = "index")
@click.option('-f', '--fasta', "fasta_file", required=True, type=click.Path(exists=True), help="FASTA file")
@click.option('-n', '--names', "names_file", required=True, type=click.Path(exists=True), help="Names file")
@click.option('-k', '--kmer-size', "kSize", required=True, type=click.IntRange(15, 31, clamp=False), help = "kmer size" )
@click.option('-o', '--output', "output_prefix", required=False, default=None, help = "index output file prefix")
def main(fasta_file, names_file, kSize, output_prefix):
    '''FASTA file indexing'''  
    
    if not output_prefix:
        output_prefix = os.path.basename(fasta_file)
        output_prefix = os.path.splitext(output_prefix)[0]
        output_prefix = "idx" + "_" + output_prefix

    if not kSize % 2:
        click.secho(f"[WARNING] kmer size: {kSize} is even, please enter an odd value", fg="yellow", bold=True, file=sys.stderr)
        sys.exit(1)

    
    idx = Index(fasta_file, names_file)
    idx.validate_names()
    idx.index(kSize)
    idx.write_to_disk(output_prefix)