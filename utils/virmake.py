import logging
import os
import subprocess
import sys

import click
from snakemake.io import load_configfile


def get_snakefile(file="workflow/Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), file)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf


# Inspired by ATLAS Metagenome
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    "virmake workflow"


@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Runs the main workflow",
)
@click.argument(
    "workflow",
    type=click.Choice(
        ["qc", "assembly", "identification", "taxonomy", "all", "None"]
    ),
)
@click.option(
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-w",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas.",
    default=".",
)
@click.option(
    "-c",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config-file generated with 'virmake init'",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Test execution.",
)
@click.option(
    "--threads",
    default=16,
    type=int,
    help="Number of threads used on multithreaded jobs",
)
def run_workflow(workflow, dryrun, working_dir, profile, config_file, threads):
    if config_file is None:
        config_file = os.path.join(working_dir, "config.yaml")

    if not os.path.exists(config_file):
        logging.critical(
            f"config-file not found: {config_file}\n"
            "generate one with 'virmake init'"
        )
        exit(1)
    sample_file = os.path.join(working_dir, "samples.tsv")

    if not os.path.exists(sample_file):
        logging.critical(
            f"sample.tsv not found in the working directory. "
            "Generate one with 'virmake init'"
        )
        exit(1)

    snakefile = get_snakefile()
    config = load_configfile(config_file)
    profile = config["workflow_dirs"]["profile_dir"]
    print("Starting workflow...")
    cmd = (
        "time"
        " snakemake --snakefile {snakefile} --directory {working_dir}"
        " --rerun-incomplete "
        "--configfile '{config_file}' --nolock "
        " {profile} --use-conda --use-singularity {dryrun} "
        "--until {target_rule}"
        " -c{threads}"
        " -T 3"
    ).format(
        snakefile=snakefile,
        config_file=config_file,
        profile="" if (profile is None) else "--profile {}".format(profile),
        dryrun="-n" if dryrun else "",
        target_rule=workflow if workflow != "None" else "",
        threads=threads,
        working_dir=working_dir,
    )

    try:
        print(cmd)
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Workflow failed, see log files")
        exit(1)


# Prep For offline
# @cli.command(
#     "prep-offline",
#     context_settings=dict(ignore_unknown_options=True),
#     short_help="Downloads and creates all enviorments needed to run the workflow offline.",
# )
# @click.option(
#     "--threads",
#     default=8,
#     type=int,
#     help="number of threads to use per multi-threaded job",
# )
# def run_prep_offline(threads):
#     """Downloads and creates all enviorments needed to run the workflow offline."""

#     cmd = (
#         "snakemake --snakefile {snakefile}"
#         "--rerun-incomplete "
#         "--conda-frontend mamba"
#         " --nolock  --use-conda --use-singularity --conda-create-envs-only"
#         " --show-failed-logs"
#         " -c{threads}"
#     ).format(snakefile=get_snakefile(), threads=threads)
#     try:
#         subprocess.check_call(cmd, shell=True)
#     except subprocess.CalledProcessError as e:
#         # removes the traceback
#         logging.critical(e)
#         exit(1)


if __name__ == "__main__":
    cli()
