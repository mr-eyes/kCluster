import sys
import click

from index import main as index_main
from cluster import main as cluster_main
from pairwise import main as pairwise_main
from dump import main as dump_main
from version import __version__


@click.group()
@click.version_option(version=__version__, prog_name="kCluster")
def cli():
    pass

cli.add_command(index_main, name="index")
cli.add_command(pairwise_main, name="pairwise")
cli.add_command(cluster_main, name="cluster")
cli.add_command(dump_main, name="dump")



if __name__ == '__main__':
    KCLUSTER =("""
        ██╗  ██╗ ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ 
        ██║ ██╔╝██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗
        █████╔╝ ██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝
        ██╔═██╗ ██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗
        ██║  ██╗╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║
        ╚═╝  ╚═╝ ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝                                            
    \n\n""")

    click.secho(KCLUSTER, fg='green', bold=True, file=sys.stderr)
    cli()