import logging
import pathlib
import subprocess

import click


# Inspired by ATLAS Metagenome
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    "VirMake is a workflow for the analysis of viral metagenomes"


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
    default=8,
    type=int,
    help="Number of threads used on multithreaded jobs",
)
def run_workflow(workflow, dryrun, working_dir, profile, config_file, threads):
    """Runs the main workflow"""
    # load virmake path
    virmake_path = pathlib.Path(__file__).parent

    # load needed paths and check if they exist
    if not config_file:
        config_file = virmake_path / "config.yaml"
    else:
        config_file = pathlib.Path(config_file).resolve()
    if not config_file.exists():
        logging.critical(
            f"config-file not found: {config_file}\n"
            "generate one by running `python setup.py`"
        )
        exit(1)
    if profile:
        profile = f"--profile {pathlib.Path(profile).resolve()} "
        if not profile.exists():
            logging.critical(f"profile not found: {profile}\n")
            exit(1)
    else:
        profile = ""
    if not working_dir:
        working_dir = virmake_path / "working_dir"
    else:
        working_dir = pathlib.Path(working_dir).resolve()

    cmd = (
        "time "
        "snakemake --directory {working_dir} "
        "--rerun-incomplete "
        "--configfile '{config_file}' --nolock "
        "--use-conda --use-singularity {dryrun} "
        "--until {target_rule} "
        "{profile}"
        "-c{threads} -T 3"
    ).format(
        config_file=config_file,
        dryrun="-n" if dryrun else "",
        target_rule=workflow if workflow != "None" else "",
        threads=threads,
        working_dir=working_dir,
        profile=profile,
    )

    print("Starting workflow...")
    print(cmd)
    try:
        subprocess.run(cmd, shell=True)
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
