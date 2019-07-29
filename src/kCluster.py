import sys
import click

from src.click_context import cli
from .index.kp_indexing import kmers as index_kmers       # pylint: disable=relative-beyond-top-level
from .index.kp_indexing import skipmers as index_skipmers       # pylint: disable=relative-beyond-top-level
from .cluster.clustering import main as cluster_main    # pylint: disable=relative-beyond-top-level
from .pairwise.virtualQs import main as pairwise_main   # pylint: disable=relative-beyond-top-level
from .dump.dump import main as dump_main                # pylint: disable=relative-beyond-top-level
from .preprocessing.skipmers_cds import preprocess_cds           # pylint: disable=relative-beyond-top-level



# cli.add_command(index_main, name="index")
cli.add_command(index_kmers, name="index_kmers")
cli.add_command(index_skipmers, name="index_skipmers")
cli.add_command(pairwise_main, name="pairwise")
cli.add_command(cluster_main, name="cluster")
cli.add_command(dump_main, name="dump")
cli.add_command(preprocess_cds, name="preprocess_cds")


if __name__ == '__main__':
    cli()