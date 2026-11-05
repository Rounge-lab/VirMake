"""locate package assets through :mod:`virmake.paths`. ``params.yaml`` should be loaded
per command rather than at import time"""

import pathlib
import shutil
import subprocess

import click
import yaml

from virmake import paths
from virmake.config import (
    default_db_dir,
    default_envs_dir,
    make_params,
    make_results_dirs,
    resolve_resource_dir,
    write_params,
    write_profiles,
)
from virmake.samples import write_samples_table
from virmake.testdata import download_test_data


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
    short_help="Prepare a run directory (samples table + configuration).",
)
@click.argument(
    "rundir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
)
@click.option(
    "--reads",
    type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path),
    help="Directory with paired-end reads ({sample}_{frac}.fastq[.gz]).",
)
@click.option(
    "--qc-reads",
    type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path),
    help="Directory with QC'd reads ({sample}_{frac}.fastq[.gz]); skips QC for these samples.",
)
@click.option(
    "--contigs",
    type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path),
    help="Directory with pre-assembled contigs ({sample}.fasta).",
)
@click.option(
    "--test-data",
    is_flag=True,
    default=False,
    help="Download the synthetic test dataset (Zenodo) into <RUNDIR>/test_data and use it as QC reads.",
)
@click.option(
    "--db-dir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
    help="Databases directory [default: per-user data dir, with confirmation].",
)
@click.option(
    "--envs-dir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
    help="Conda environment prefix [default: per-user data dir, with confirmation].",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Regenerate config and samples.tsv in an existing run directory "
    "(results/ is never touched).",
)
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    default=False,
    help="Accept default resource locations without confirmation.",
)
def prep(
    rundir,
    reads,
    qc_reads,
    contigs,
    test_data,
    db_dir,
    envs_dir,
    force,
    yes,
):
    """Prepare a run directory: samples.tsv, config, and results/."""
    rundir = rundir.expanduser()
    if test_data and (reads or qc_reads or contigs):
        raise click.ClickException(
            "--test-data cannot be combined with --reads/--qc-reads/--contigs."
        )
    if rundir.exists() and any(rundir.iterdir()) and not force:
        raise click.ClickException(
            f"{rundir} already exists and is not empty. Use --force to regenerate "
            "config and the samples table (results/ is never touched)."
        )

    db_dir = resolve_resource_dir(
        db_dir, default=default_db_dir(), description="database", assume_yes=yes
    )
    envs_dir = resolve_resource_dir(
        envs_dir,
        default=default_envs_dir(),
        description="conda environment",
        assume_yes=yes,
    )

    # Create roots early so permission problems surface now, not mid-run.
    rundir.mkdir(parents=True, exist_ok=True)
    db_dir.mkdir(parents=True, exist_ok=True)
    envs_dir.mkdir(parents=True, exist_ok=True)

    if test_data:
        click.echo("Downloading test data from Zenodo...")
        download_test_data(rundir / "test_data")
        qc_reads = rundir / "test_data"

    params = make_params(
        rundir,
        db_dir,
        reads=str(reads) if reads else "",
        contigs=str(contigs) if contigs else "",
    )
    write_params(rundir, params)
    write_profiles(rundir, envs_dir)
    make_results_dirs(rundir)

    samples_table = write_samples_table(
        rundir / "samples.tsv",
        reads=str(reads) if reads else "",
        qc_reads=str(qc_reads) if qc_reads else "",
        contigs=str(contigs) if contigs else "",
    )
    if samples_table.empty:
        click.echo(
            "WARNING: no samples found. samples.tsv is header-only. Provide "
            "--reads/--qc-reads/--contigs, use --test-data, or edit "
            f"{rundir / 'samples.tsv'} by hand."
        )

    click.echo(f"Run directory ready: {rundir}")
    click.echo(f"  databases: {db_dir}")
    click.echo(f"  conda envs: {envs_dir}")
    click.echo("Next steps:")
    click.echo(f"  virmake envs {rundir}  # optional: prefetch rule environments")
    click.echo(f"  virmake db {rundir}    # download databases")
    click.echo(f"  virmake run {rundir}   # run the workflow")


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
