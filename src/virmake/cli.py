"""locate package assets through :mod:`virmake.paths`. ``params.yaml`` should be loaded
per command rather than at import time"""

import pathlib
import shutil
import subprocess

import click
import yaml

from virmake import paths


def resolve_config_file(config_file: str | None) -> pathlib.Path:
    if config_file:
        path = pathlib.Path(config_file)
    else:
        path = pathlib.Path("config") / "params.yaml"
    path = path.resolve()
    if not path.exists():
        raise click.ClickException(
            f"config file not found: {path}. Pass --config-file or run from a "
            "directory containing config/params.yaml."
        )
    return path


def load_params(config_file: pathlib.Path) -> dict:
    """Load and check params.yaml"""
    try:
        with open(config_file) as fh:
            config = yaml.safe_load(fh)
    except yaml.YAMLError as e:
        raise click.ClickException(f"could not parse {config_file}: {e}") from e
    if not isinstance(config, dict) or "path" not in config:
        raise click.ClickException(
            f"{config_file} does not look like a VirMake params file "
            "(missing 'path' section)."
        )
    return config


# Inspired by ATLAS Metagenome
@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    """VirMake is a workflow for the analysis of viral metagenomes from paired-end data.
    It offers taxonomic and functional annotation, supports offline running and has
    support for HPC cluster execution."""


@cli.command(
    "db",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Downloads databases for steps defined in config.\n"
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
@click.option(
    "-C",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="params.yaml to use [default: ./config/params.yaml]",
)
def get_dbs(threads, dryrun, config_file):
    """Downloads all databases needed to run the workflow."""
    config = resolve_config_file(config_file)
    cmd = (
        "snakemake --snakefile {snakefile} "
        "--conda-frontend mamba --configfile {config} "
        "--nolock  --use-conda "
        "--show-failed-logs --directory {workflow_dir} "
        "-c{threads} "
        "SETUP_DB "
        "{dryrun} "
    ).format(
        snakefile=paths.snakefile(),
        config=config,
        workflow_dir=pathlib.Path.cwd(),
        threads=threads,
        dryrun="-n" if dryrun else "",
    )
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise click.ClickException(f"database setup failed: {e}") from e


@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
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
    "--wf_profile",
    type=click.Path(exists=True, resolve_path=True),
    help="directory with config file for command line arguments\n[default: ./config]",
)
@click.option(
    "-d",
    "--workflow-dir",
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    help="location to run virmake [default: current directory]",
)
@click.option(
    "-C",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="config file with virmake parameters\n[default: ./config/params.yaml]",
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
    "--slurm",
    is_flag=True,
    default=False,
    show_default=True,
    help="use slurm cluster to run parallel jobs\n"
    "Parameters in the slurm profile located at config/slurm/config.yaml "
    "may need to be set.",
)
@click.pass_context
def run_workflow(ctx, workflow, dryrun, workflow_dir, wf_profile, config_file, slurm):
    """Runs the main workflow"""

    config_path = resolve_config_file(config_file)
    params = load_params(config_path)

    benchmark = params.get("path", {}).get("benchmark")
    if not benchmark:
        raise click.ClickException(f"config.path.benchmark missing from {config_path}")

    if wf_profile is None:
        wf_profile = pathlib.Path("config").resolve()
    if not wf_profile.is_dir():
        raise click.ClickException(
            f"workflow profile directory not found: {wf_profile} "
            "(expected config/config.yaml)."
        )
    if workflow_dir is None:
        workflow_dir = pathlib.Path.cwd()
    if slurm:
        slurm_profile = (pathlib.Path("config") / "slurm").resolve()
        if not slurm_profile.is_dir():
            raise click.ClickException(
                f"slurm profile directory not found: {slurm_profile} "
                "(expected config/slurm/config.yaml)."
            )
    else:
        slurm_profile = ""

    extra_args = " ".join(ctx.args)

    cmd = (
        "mkdir -p '{benchmark}' && "
        "/usr/bin/time -p -o '{benchmark}/total_time.txt' "
        "snakemake --directory '{workflow_dir}' "
        "--configfile '{config_file}' "
        "--until {target_rule} "
        "--workflow-profile {wf_profile} "
        "{slurm} "
        "{dryrun} "
        "{extra_args}"
    ).format(
        benchmark=benchmark,
        workflow_dir=workflow_dir,
        wf_profile=wf_profile,
        config_file=config_path,
        dryrun="-n" if dryrun else "",
        target_rule=workflow.upper(),
        slurm=f"--profile {slurm_profile} " if slurm else "",
        extra_args=extra_args,
    )

    print("Starting workflow...")
    print(cmd)
    subprocess.run(cmd, shell=True)


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
@click.option(
    "-C",
    "--config-file",
    type=click.Path(exists=True, resolve_path=True),
    help="params.yaml to use [default: ./config/params.yaml]",
)
def run_prep_offline(threads, config_file):
    """Downloads and creates environments needed to run the workflow offline."""
    config = resolve_config_file(config_file)
    cmd = (
        f"snakemake --snakefile {paths.snakefile()} "
        f"--conda-frontend mamba --configfile {config} "
        "--nolock  --use-conda --conda-create-envs-only "
        f"--show-failed-logs --directory {pathlib.Path.cwd()} "
        f"-c{threads}"
    )
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        raise click.ClickException(f"environment creation failed: {e}") from e


@cli.command(
    "clean",
    context_settings=dict(ignore_unknown_options=True),
    short_help="Clean VirMake directories in the current working directory.",
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
    # TODO: remove after implementing run-directory semantics
    if not yes:
        click.confirm(
            f"are you sure you want to delete {target}?",
            abort=True,
            default=False,
        )
    cwd = pathlib.Path.cwd()
    targets = {
        "all": [cwd / "results", cwd / "resources" / "databases", cwd / "config"],
        "databases": [cwd / "resources" / "databases"],
        "results": [cwd / "results"],
        "config": [cwd / "config"],
    }
    for path in targets[target]:
        if path.exists():
            click.echo(f"removing {path}")
            shutil.rmtree(path)


if __name__ == "__main__":
    cli()
