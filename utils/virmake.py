import logging
import pathlib
import shutil
import subprocess
import os

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
    "-p",
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-W",
    "--working-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run atlas.",
)
@click.option(
    "-C",
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
    "-c",
    "--threads",
    default=24,
    type=int,
    help="Number of threads used on multithreaded jobs",
)
def run_workflow(workflow, dryrun, working_dir, profile, config_file, threads):
    """Runs the main workflow"""
    # load virmake path
    virmake_path = pathlib.Path(__file__).parent

    # load needed paths and check if they exist
    if not config_file:
        config_file = virmake_path / "workflow" / "config.yaml"
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


# # Prepare VirMake for offline use
# @cli.command(
#     "prep",
#     context_settings=dict(ignore_unknown_options=True),
#     short_help="Downloads and creates all enviorments needed to run the workflow offline.",
# )
# @click.option(
#     "--threads",
#     default=24,
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


# Download sample data
@cli.command(
    "get",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Downloads read data from public databases.",
)
@click.argument(
    "database",
    type=click.Choice(["SRA"]),
)
@click.argument(
    "accession",
    type=str,
)
@click.option(
    "-o",
    "--output-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="Location to download data to. Default is 'working_dir/input'.",
    default=pathlib.Path(__file__).parent / "working_dir" / "input",
)
def run_get(database, accession, output_dir):
    """Downloads read data from public databases."""
    db_commands = {
        "SRA": "fasterq-dump {accession} -O {output_dir}".format(
            accession=accession, output_dir=output_dir
        )
    }
    cmd = db_commands[database]
    try:
        print(f"Downloading {accession} from {database}...")
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


# Clean
@cli.command(
    "clean",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Clean VirMake directory.",
)
@click.argument(
    "target",
    type=click.Choice(["all", "databases", "working_dir", "config"]),
)
@click.option(
    "-y",
    default=False,
)
def clean(target, y):
    """clean virmake directory."""
    virmake_path = pathlib.Path(__file__).parent
    if not y:
        click.confirm(
            f"are you sure you want to delete {target}?",
            abort=True,
            default=False,
        )
    if target == "all":
        if [
            (virmake_path / "working_dir").exists(),
            (virmake_path / "databases").exists(),
            (virmake_path / "config.yaml").exists(),
        ]:
            shutil.rmtree(virmake_path / "working_dir")
            shutil.rmtree(virmake_path / "databases")
            os.remove(virmake_path / "config.yaml")
    elif target == "databases":
        if (virmake_path / "databases").exists():
            shutil.rmtree(virmake_path / "databases")
    elif target == "working_dir":
        if (virmake_path / "working_dir").exists():
            dirs = os.listdir(virmake_path / "working_dir")
            if ".snakemake" in dirs:
                logging.warning(
                    "Not removing .snakemake directory!"
                    "This is the snakemake working directory."
                    "If you want to remove this directory, "
                    "by doing so you will corrupt your workflow setup."
                    "To retain your workflow setup, copy it into your new"
                    "working directory."
                )
                dirs.remove(".snakemake")
            for d in dirs:
                shutil.rmtree(virmake_path / "working_dir" / d)
    elif target == "config":
        if (virmake_path / "workflow" / "config.yaml").exists():
            os.remove(virmake_path / "workflow" / "config.yaml")
    else:
        logging.critical(f"unknown target: {target}")
        exit(1)


if __name__ == "__main__":
    cli()
