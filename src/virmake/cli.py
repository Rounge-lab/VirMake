"""Virmake CLI

Commands operate on the run directory created by ``virmake prep``:

- ``prep <RUNDIR>``  create the run directory (samples table + config)
- ``envs <RUNDIR>``  prefetch the workflow's conda environments
- ``db   <RUNDIR>``  download the workflow's databases
- ``run  <RUNDIR>``  execute the workflow (local or SLURM)

Snakemake is invoked as an argument list (never ``shell=True``) with the
package-resolved Snakefile and the run directory as working directory.
"""

import pathlib
import shlex
import shutil
import subprocess
import sys
import time

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

_REQUIRED_DB_FLAGS = ("virsorter2", "genomad", "checkv")


def resolve_rundir(rundir: str | pathlib.Path) -> pathlib.Path:
    """Validate that ``rundir`` was prepared by ``virmake prep``."""
    rundir = pathlib.Path(rundir).expanduser().resolve()
    if not (rundir / "config" / "params.yaml").is_file():
        raise click.ClickException(
            f"{rundir} is not a prepared run directory (missing config/params.yaml).\n"
            f"Create one with: virmake prep {rundir}"
        )
    return rundir


def load_params(rundir: pathlib.Path) -> dict:
    """Load the run directory's params.yaml (minimal check; full validation is P5)."""
    path = rundir / "config" / "params.yaml"
    try:
        config = yaml.safe_load(path.read_text())
    except yaml.YAMLError as e:
        raise click.ClickException(f"could not parse {path}: {e}") from e
    if not isinstance(config, dict) or "path" not in config:
        raise click.ClickException(
            f"{path} does not look like a VirMake params file (missing 'path' section)."
        )
    return config


def _load_workflow_profile(rundir: pathlib.Path) -> dict:
    path = rundir / "config" / "config.yaml"
    if not path.is_file():
        raise click.ClickException(
            f"{path} is missing. Regenerate it with: virmake prep {rundir} --force"
        )
    profile = yaml.safe_load(path.read_text())
    if not isinstance(profile, dict):
        raise click.ClickException(f"{path} is not a valid profile file.")
    return profile


def _conda_prefix(rundir: pathlib.Path, assume_yes: bool) -> pathlib.Path:
    """The machine-level conda env prefix recorded for this run (plan D5).

    If the profile carries no prefix (hand-edited config), resolve one via
    the resource-dir flow and write it back to the profile.
    """
    profile = _load_workflow_profile(rundir)
    prefix = str(profile.get("conda-prefix") or "")
    if prefix:
        return pathlib.Path(prefix).expanduser()
    prefix = resolve_resource_dir(
        None,
        default=default_envs_dir(),
        description="conda environment",
        assume_yes=assume_yes,
    )
    profile["conda-prefix"] = str(prefix.expanduser().resolve())
    with open(rundir / "config" / "config.yaml", "w") as fh:
        yaml.safe_dump(profile, fh, sort_keys=False)
    return prefix


def _snakemake_base(rundir: pathlib.Path) -> list[str]:
    """Common snakemake invocation prefix (package Snakefile, run dir, params)."""
    return [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(paths.snakefile()),
        "--directory",
        str(rundir),
        "--configfile",
        str(rundir / "config" / "params.yaml"),
    ]


def _run_snakemake(cmd: list[str], description: str) -> float:
    """Run snakemake (list args); returns elapsed wall seconds."""
    click.echo(f"{description.capitalize()}: {shlex.join(cmd)}")
    start = time.perf_counter()
    proc = subprocess.run(cmd)
    elapsed = time.perf_counter() - start
    if proc.returncode != 0:
        raise click.ClickException(
            f"{description} failed with exit code {proc.returncode}"
        )
    return elapsed


def _confirm_remove(path: pathlib.Path, yes: bool) -> None:
    """Remove ``path`` after confirmation; recreate it empty (keeps perms visible)."""
    if not path.exists():
        click.echo(f"Nothing to remove: {path} does not exist.")
        return
    if not yes and not click.confirm(
        f"Delete {path} and ALL of its contents?\n"
        "This location is shared across runs on this machine.",
        default=False,
    ):
        raise click.Abort()
    click.echo(f"Removing {path}")
    shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def _preflight(rundir: pathlib.Path, params: dict, prefix: pathlib.Path) -> None:
    """Non-blocking sanity checks before handing off to snakemake (plan §3.2)."""
    samples = pathlib.Path(params["path"]["samples"])
    if not samples.is_file() or samples.stat().st_size == 0:
        raise click.ClickException(f"samples table missing or empty: {samples}")
    if len([ln for ln in samples.read_text().splitlines() if ln.strip()]) <= 1:
        click.echo(
            f"WARNING: {samples} appears to contain no samples (header only).",
            err=True,
        )
    missing = [
        name
        for name in _REQUIRED_DB_FLAGS
        if not (pathlib.Path(params["path"]["database"][name]) / "flag").is_file()
    ]
    if missing:
        click.echo(
            f"WARNING: databases not set up: {', '.join(missing)}.\n"
            f"         Run `virmake db {rundir}` first; required steps will fail "
            "if they stay missing.",
            err=True,
        )
    if not prefix.is_dir() or not any(prefix.iterdir()):
        click.echo(
            f"WARNING: conda env prefix {prefix} is empty.\n"
            f"         Run `virmake envs {rundir}` to prefetch; environments are "
            "otherwise created lazily on first use.",
            err=True,
        )
    if sys.platform == "darwin":
        click.echo(
            "WARNING: macOS support is best-effort; some workflow environments "
            "may not be available.",
            err=True,
        )


# CLI


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
def cli():
    """VirMake is a workflow for the analysis of viral metagenomes."""


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
    "envs",
    short_help="Prefetch the workflow's conda environments.",
)
@click.argument(
    "rundir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
)
@click.option(
    "-c",
    "--threads",
    default=24,
    type=int,
    show_default=True,
    help="cores used while creating environments",
)
@click.option(
    "--clean",
    is_flag=True,
    default=False,
    help="remove ALL VirMake conda environments on this machine",
)
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    default=False,
    help="skip confirmation for --clean and default locations",
)
def make_envs(rundir, threads, clean, yes):
    """Create all conda environments needed by the workflow (prefetch).

    Environments live in the prefix recorded by `virmake prep` and are
    shared across runs on this machine (content-hashed by snakemake).
    """
    rundir = resolve_rundir(rundir)
    prefix = _conda_prefix(rundir, assume_yes=yes)
    if clean:
        _confirm_remove(prefix, yes)
        return
    cmd = _snakemake_base(rundir) + [
        "--software-deployment-method",
        "conda",
        "--conda-prefix",
        str(prefix),
        "--conda-create-envs-only",
        "--cores",
        str(threads),
    ]
    _run_snakemake(cmd, description="environment creation")


@cli.command(
    "db",
    short_help="Download the workflow's databases.",
)
@click.argument(
    "rundir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
)
@click.option(
    "-c",
    "--threads",
    default=24,
    type=int,
    show_default=True,
    help="cores used on multithreaded jobs",
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
    "--clean",
    is_flag=True,
    default=False,
    help="remove the recorded databases directory",
)
@click.option(
    "-y",
    "--yes",
    is_flag=True,
    default=False,
    help="skip confirmation for --clean and default locations",
)
def setup_dbs(rundir, threads, dryrun, clean, yes):
    """Download all databases needed by the workflow into the recorded db dir."""
    rundir = resolve_rundir(rundir)
    prefix = _conda_prefix(rundir, assume_yes=yes)
    params = load_params(rundir)
    if clean:
        # db root = parent of the checkv entry (the writer guarantees the layout)
        db_root = pathlib.Path(params["path"]["database"]["checkv"]).parent
        _confirm_remove(db_root, yes)
        return
    cmd = _snakemake_base(rundir) + [
        "--software-deployment-method",
        "conda",
        "--conda-prefix",
        str(prefix),
        "--cores",
        str(threads),
        "SETUP_DB",
    ]
    if dryrun:
        cmd.append("-n")
    _run_snakemake(cmd, description="database setup")


@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True, allow_extra_args=True),
    short_help="Run the workflow on a prepared run directory.",
)
@click.argument(
    "rundir",
    type=click.Path(file_okay=False, path_type=pathlib.Path),
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
    "--cores",
    default=1,
    type=int,
    show_default=True,
    help="cores for local execution",
)
@click.option(
    "--slurm",
    is_flag=True,
    default=False,
    show_default=True,
    help="submit jobs via SLURM (profile: config/slurm)",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="test execution",
)
@click.pass_context
def run_workflow(ctx, rundir, workflow, cores, slurm, dryrun):
    """Run the workflow on a prepared run directory."""
    rundir = resolve_rundir(rundir)
    params = load_params(rundir)
    prefix = _conda_prefix(rundir, assume_yes=False)
    _preflight(rundir, params, prefix)

    cmd = _snakemake_base(rundir) + [
        "--workflow-profile",
        str(rundir / "config"),
        workflow.upper(),
        "--cores",
        str(cores),
    ]
    if slurm:
        cmd += ["--profile", str(rundir / "config" / "slurm")]
    if dryrun:
        cmd.append("-n")
    cmd += list(ctx.args)

    click.echo("Starting workflow...")
    elapsed = _run_snakemake(cmd, description="workflow")
    if not dryrun:
        bench = rundir / "results" / "benchmark" / "total_time.txt"
        bench.parent.mkdir(parents=True, exist_ok=True)
        bench.write_text(f"{elapsed:.1f} seconds\n")


if __name__ == "__main__":
    cli()
