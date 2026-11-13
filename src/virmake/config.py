"""Run-directory configuration for VirMake.

Writers for each config file:

- ``config/params.yaml``
- ``config/config.yaml``       — Snakemake workflow profile
- ``config/slurm/config.yaml`` — Snakemake SLURM executor profile

Resource locations for databases and environments are now user-wide with defaults
through platformdirs."""

from pathlib import Path

import click
import platformdirs
import yaml

from virmake import paths


def default_data_root() -> Path:
    """Per-user VirMake data root (XDG data dir or platform equivalent)."""
    return Path(platformdirs.user_data_dir("virmake"))


def default_db_dir() -> Path:
    """Default databases location: ``<data root>/databases``."""
    return default_data_root() / "databases"


def default_envs_dir() -> Path:
    """Default conda env prefix: ``<data root>/envs``."""
    return default_data_root() / "envs"


def resolve_resource_dir(
    requested: str | Path | None,
    *,
    default: Path,
    description: str,
    assume_yes: bool = False,
) -> Path:
    """Fallback to ``default`` requires interactive confirmation unless ``assume_yes``.
    On decline, prompts for a path."""
    if requested is not None:
        return Path(requested).expanduser()
    if assume_yes:
        click.echo(f"Using default {description} location: {default}")
        return default
    click.echo(
        f"Default {description} location is {default}\n"
        "WARNING: this is inside your HOME directory. On many clusters HOME\n"
        "has a small quota, while these resources can be large (tens of GB).\n"
        "Consider a project-shared location."
    )
    if click.confirm(f"Use {default} for {description}?", default=False):
        return default
    return click.prompt(
        f"Enter the {description} location",
        type=click.Path(exists=False, file_okay=False, path_type=Path),
    ).expanduser()


# params.yaml


def write_params(rundir: str | Path, params: dict) -> Path:
    """``<rundir>/config/params.yaml``"""
    path = Path(rundir) / "config" / "params.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        yaml.safe_dump(params, fh, sort_keys=False)
    return path


def make_results_dirs(rundir: str | Path) -> None:
    """``<rundir>/results/{output,log,benchmark,temp}``"""
    rundir = Path(rundir)
    for sub in ("output", "log", "benchmark", "temp"):
        (rundir / "results" / sub).mkdir(parents=True, exist_ok=True)


# Profile files


def _load_template(name: str) -> dict:
    return yaml.safe_load((paths.templates_dir() / name).read_text())


def make_path_block(
    rundir: str | Path,
    db_dir: str | Path,
    reads: str = "",
    contigs: str = "",
) -> dict:
    """Resolves all the required paths in the config"""
    rundir = Path(rundir).expanduser().resolve()
    db_dir = Path(db_dir).expanduser().resolve()

    def _input(value: str) -> str:
        return str(Path(value).expanduser().resolve()) if value else ""

    return {
        "virmake": str(paths.package_root()),
        "envs": str(paths.envs_dir()),
        "output": str(rundir / "results" / "output"),
        "log": str(rundir / "results" / "log"),
        "benchmark": str(rundir / "results" / "benchmark"),
        "temp": str(rundir / "results" / "temp"),
        "scripts": str(paths.scripts_dir()),
        "samples": str(rundir / "samples.tsv"),
        "input_reads": _input(reads),
        "input_contigs": _input(contigs),
        "database": {
            "RefSeq": str(db_dir / "RefSeq" / "viral.1.1.genomic.fna"),
            "virsorter2": str(db_dir / "virsorter2"),
            "genomad": str(db_dir / "genomad"),
            "vibrant": str(db_dir / "vibrant" / "vibrant-1.2.1"),
            "checkv": str(db_dir / "checkv"),
            "INPHARED": str(db_dir / "INPHARED"),
            "vcontact2": str(db_dir / "vcontact2"),
            "vcontact3": str(db_dir / "vcontact3"),
            "DRAM": str(db_dir / "DRAM" / "DRAM_data"),
        },
    }


def make_params(
    rundir: str | Path,
    db_dir: str | Path,
    reads: str = "",
    contigs: str = "",
) -> dict:
    """template + parsed paths"""
    params = _load_template("params.yaml")
    params["path"] = make_path_block(rundir, db_dir, reads, contigs)
    return params


def make_workflow_profile(envs_dir: str | Path) -> dict:
    """``<rundir>/config/config.yaml``"""
    profile = _load_template("workflow-profile.yaml")
    profile["conda-prefix"] = str(Path(envs_dir).expanduser().resolve())
    return profile


def make_slurm_profile() -> dict:
    """``config/slurm/config.yaml``"""
    return _load_template("slurm-profile.yaml")


def write_profiles(rundir: str | Path, envs_dir: str | Path) -> None:
    """Write ``<rundir>/config/config.yaml`` and ``config/slurm/config.yaml``."""
    rundir = Path(rundir)
    slurm_dir = rundir / "config" / "slurm"
    slurm_dir.mkdir(parents=True, exist_ok=True)
    with open(rundir / "config" / "config.yaml", "w") as fh:
        yaml.safe_dump(make_workflow_profile(envs_dir), fh, sort_keys=False)
    with open(slurm_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(make_slurm_profile(), fh, sort_keys=False)
