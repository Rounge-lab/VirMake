import os, sys
import multiprocessing
import subprocess
import click
from get_sample_table import *
from make_config import *

def intilialize_virmake(
    sample_path, outfile="samples.tsv"
):
    """
    Creates a sample file for all samples.
    """

    if os.path.exists(outfile):
            print(f"Output file {outfile} already exists, will not overwrite")
            exit(1)

    write_sample_table(sample_path)

@click.command(
    "init",
    short_help="prepare configuration file and sample table for VirMake run",
)
@click.argument("samples_path", type=click.Path(readable=True))
@click.option(
    "-d",
    "--db-dir",
    default=os.path.join(os.path.realpath("."), "Databases"),
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    show_default=True,
    help="location to store databases",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location for running the application",
    default=".",
)
@click.option(
    "--threads",
    default=8,
    type=int,
    help="number of threads to use per multi-threaded job",
)
def run_init(
    samples_path,
    db_dir,
    working_dir,
    threads,
):
    """
    Make a config file and create a sample table
    """

    #Make sure the directories exist and creates them if not
    os.makedirs(working_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    sample_table = intilialize_virmake(samples_path)
    make_config(db_dir, threads)
