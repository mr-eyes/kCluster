import sqlite3
import sys
import os
import click
from src.click_context import cli

class Dump:

    def __init__(self, logger_obj, sqlite_file):
        self.Logger = logger_obj
        try:
            self.conn = sqlite3.connect(sqlite_file)
        
        except sqlite3.Error as err:
            print(err)
            _err_msg = f"couldn't connect to {sqlite_file}"
            click.echo(click.style(_err_msg, fg='red', bold=True), file=sys.stderr)
            sys.exit(1)
        
        self.sqlite_file = sqlite_file

    
    def get_table_columns(self, table_name):
        """return sqlite database table columns.
        Args:
            table_name (str): table name.
        Return:
            col_names (list)
        """

        curs = self.conn.execute(f'select * from {table_name}')
        col_names = list(map(lambda x: x[0], curs.description))
        return col_names

    def get_table_data(self, table_name):
        curs = self.conn.execute(f'select * from {table_name}')
        rows = curs.fetchall()
        return rows

    def table_export(self, table_name):
        """export sqlite table to TSV file printed to the stdout"""

        header = self.get_table_columns(table_name)
        rows = self.get_table_data(table_name)
        print("\t".join(header))
        for row in rows:
            row = list(map(str, row))
            print("\t".join(row))

    def simple_export(self):
        header = ["seq_1", "seq_2", "shared_kmers"]
        rows = self.get_table_data("virtualQs")
        print("\t".join(header))
        for row in rows:
            row = map(str, list([row[1], row[2], row[-1]]))
            print("\t".join(row))
        

@cli.command(name = "dump", help_priority=4)
@click.option('-d', '--db', required=True, type=click.Path(exists=True), help="sqlite database file")
@click.option('-t', '--table', required=False, type=click.Choice(['virtualQs', 'meta_info', 'namesmap']), show_default=True, default="virtualQs", help="database table to be exported")
@click.option('--simple', 'simple_output', is_flag=True, required=False, help="export in a tsv output [seq1,seq2,shared] no virtualQs")
@click.pass_context
def main(ctx, db, table, simple_output):
    """Dump sqlite database table to the stdout in TSV format."""
    tsv = Dump(logger_obj = ctx.obj, sqlite_file = db)
    if simple_output:
        tsv.simple_export()
    else:
        tsv.table_export(table)