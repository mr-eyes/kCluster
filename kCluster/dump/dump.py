import sqlite3
import sys
import os
import click

class Dump:

    def __init__(self, sqlite_file):
         
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
        


@click.command("dump")
@click.option('-d', '--db', required=True, type=click.Path(exists=True), help="sqlite database file")
@click.option('-t', '--table', required=False, type=click.Choice(['virtualQs', 'meta_info', 'namesmap']), show_default=True, default="virtualQs", help="database table to be exported")
def main(db, table):
    """Dump sqlite database table to the stdout in TSV format."""

    tsv = Dump(sqlite_file = db)
    tsv.table_export(table)