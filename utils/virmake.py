import logging
import pathlib
import shutil
import subprocess
import os
import yaml

import click




# Inspired by ATLAS Metagenome
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    "VirMake is a workflow for the analysis of viral metagenomes"


# load config file
try:
    with open("config/params.yaml", "r") as cf:
        try:
            config = yaml.safe_load(cf)
        except yaml.YAMLError as ye:
            logging.critical(ye)
            exit(1)
except FileNotFoundError:
    logging.critical("Config file not found: config/params.yaml")
    exit(1)

# load virmake path
virmake_path = pathlib.Path(__file__).parent.absolute()


@cli.command(
    "db",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Downloads databases for steps defined in config.\n"+
        "Databases will be downloaded as needed if not downloaded using this option.",
)
@click.option(
    "-c",
    "--threads",
    default=24,
    type=int,
    help="maximum number of threads used on multithreaded jobs",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="test execution",
)
def get_dbs(threads, dryrun):
    """Downloads all databases needed to run the workflow."""
    cmd = (
        "snakemake --snakefile {snakefile} "
        "--conda-frontend mamba --configfile {config} "
        "--nolock  --use-conda "
        "--show-failed-logs --directory {workflow_dir} "
        "-c{threads} "
        "SETUP_DB "
        "{dryrun} "
    ).format(
        snakefile=virmake_path / "workflow" / "Snakefile",
        threads=threads,
        config=virmake_path / "config" / "params.yaml",
        workflow_dir=virmake_path / "workflow",
        dryrun = "-n" if dryrun else ""
    )
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)

@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True,
                          allow_extra_args=True),
    short_help="Runs the main workflow",
)
@click.argument(
    "workflow",
    type=click.Choice(
        [
            "all",
            "qc",
            "assembly",
            "identification",
            "mapping",
            "taxonomy",
            "function",
            "stats",
        ]
    ),
    default="all",
)
@click.option(
    "-p",
    "--profile",
    default=None,
    help="snakemake profile e.g. for cluster execution.",
)
@click.option(
    "-d",
    "--workflow-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run virmake.",
)
@click.option(
    "-C",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config file generated during virmake setup",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="test execution",
)
@click.option(
    "-c",
    "--threads",
    default=24,
    type=int,
    help="maximum number of threads used on multithreaded jobs",
)
@click.option(
    "-s",
    "--slurm",
    is_flag=True,
    default=False,
    show_default=True,
    help="use slurm cluster to run parallel jobs",
)
@click.option(
    "-T",
    "--jobs_at_once",
    default=3,
    type=int,
    help="number of jobs to add to queue at once",
)
@click.pass_context

#add cli command for the pipeline argument in make_virmake.py
def run_workflow(
    ctx,
    workflow,
    dryrun,
    workflow_dir,
    profile,
    config_file,
    threads,
    slurm,
    jobs_at_once,
):
    """Runs the main workflow"""

    # load needed paths and check if they exist
    if not config_file:
        config_file = virmake_path / "config" / "params.yaml"
    else:
        config_file = pathlib.Path(config_file).resolve()
    if not config_file.exists():
        logging.critical(
            f"config-file not found: {config_file}\n"
            "generate one by running `python setup.py`"
        )
        exit(1)
    if not profile:
        profile = virmake_path / "config"
    else:
        profile = pathlib.Path(profile).resolve()
    if not workflow_dir:
        workflow_dir = virmake_path / "workflow"
    else:
        workflow_dir = pathlib.Path(workflow_dir).resolve()
    
    extra_args = " ".join(ctx.args)

    cmd = (
        "mkdir -p '{benchmark}' && "
        "/usr/bin/time -p -o '{benchmark}/total_time.txt' "
        "snakemake --directory '{workflow_dir}' "
        "--configfile '{config_file}' "
        "--until {target_rule} "
        "-c{threads} -T {jobs_at_once} "
        "{slurm} -j{threads} "
        "{dryrun} "
        "{extra_args}"
    ).format(
        benchmark=pathlib.Path(config["path"]["benchmark"]),
        workflow_dir=workflow_dir,
        config_file=config_file,
        dryrun="-n" if dryrun else "",
        target_rule=workflow.upper(),
        threads=threads,
        slurm=f"--profile {profile} --slurm --default-resources slurm_account={config['slurm_account']}"
        if slurm
        else "--use-conda --rerun-incomplete --nolock --conda-frontend mamba",
        jobs_at_once=jobs_at_once,
        extra_args=extra_args
    )

    print("Starting workflow...")
    print(cmd)
    try:
        subprocess.run(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        logging.critical(f"Workflow failed, see log files")
        exit(1)


# Prepare VirMake for offline use
@cli.command(
    "prep",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Downloads and creates all environments needed to run the workflow offline.",
)
@click.option(
    "-c",
    "--threads",
    default=24,
    type=int,
    help="number of threads to use per multi-threaded job",
)
def run_prep_offline(threads):
    """Downloads and creates all environments needed to run the workflow offline."""
    cmd = (
        "snakemake --snakefile {snakefile} "
        "--conda-frontend mamba --configfile {config} "
        "--nolock  --use-conda --conda-create-envs-only "
        "--show-failed-logs --directory {workflow_dir} "
        "-c{threads}"
    ).format(
        snakefile=virmake_path / "workflow" / "Snakefile",
        threads=threads,
        config=virmake_path / "config" / "params.yaml",
        workflow_dir=virmake_path / "workflow",
    )
    try:
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
    type=click.Choice(["all", "databases", "results", "config"]),
)
@click.option(
    "-y",
    "--yes",
    default=False,
    is_flag=True,
    help="Skip confirmation",
)
def clean(target, yes):
    """clean virmake directory."""
    if not yes:
        click.confirm(
            f"are you sure you want to delete {target}?",
            abort=True,
            default=False,
        )
    if target == "all":
        if (virmake_path / "results").exists():
            shutil.rmtree(virmake_path / "results")
        if (virmake_path / "resources" / "databases").exists():
            shutil.rmtree(virmake_path / "resources" / "databases")
        if (virmake_path / "config").exists():
            shutil.rmtree(virmake_path / "config")
    elif target == "databases":
        if (virmake_path / "resources" / "databases").exists():
            shutil.rmtree(virmake_path / "resources" / "databases")
    elif target == "results":
        if (virmake_path / "results").exists():
            dirs = os.listdir(virmake_path / "results")
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
                shutil.rmtree(virmake_path / "results" / d)
    elif target == "config":
        if (virmake_path / "config").exists():
            shutil.rmtree(virmake_path / "config")
    else:
        logging.critical(f"unknown target: {target}")
        exit(1)


if __name__ == "__main__":
    cli()
