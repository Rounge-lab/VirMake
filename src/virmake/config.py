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


def make_params(
    rundir: str | Path,
    db_dir: str | Path,
    reads: str = "",
    contigs: str = "",
) -> dict:
    rundir = Path(rundir).expanduser().resolve()
    db_dir = Path(db_dir).expanduser().resolve()

    def _input(value: str) -> str:
        return str(Path(value).expanduser().resolve()) if value else ""

    return {
        "slurm_account": "default",
        "keep_qc_reads": True,
        "assembler": "metaSpades",  # may be changed to megahit
        "identifier": "virsorter2",  # may be changed to genomad
        "trim_percentage": 0.05,
        "min_coverage": 75,
        "min_contig_size": 1000,
        "vibrant": {"is_virome": "no"},
        "virsorter2": {
            "id": {
                "min_length": 3000,
                "min_score": 0.5,
                "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
            },
            "for_dramv": {
                "min_length": 100,
                "min_score": 0.01,
                "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
            },
        },
        "vcontact3": {"db_version": "latest", "exports": "cytoscape"},
        "path": {
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
        },
        "dereplication": {
            "ani": 97,
            "precluster_ani": 95,
            "min_aligned_fraction": 70,
            "vOTU_num_start": 1,
            "vOTU_prefix": "vOTU",
            "vOTU_suffix": "",
            "vOTU_num_len": 5,
        },
        "quality_threshold": "medium",
        "threads": 24,
        "job_type": {
            "small": "normal",
            "normal": "normal",
            "big": "normal",
        },
        "memory": {
            "tiny": 3750,
            "small": 7500,
            "normal": 15000,
            "big": 30000,
            "huge": 60000,
            "gigantic": 120000,
            "enormous": 240000,
        },
        "time": {
            "tiny": "30 m",
            "small": "1 h",
            "normal": "6 h",
            "big": "24 h",
            "long": "120 h",
        },
        "rule_inclusion": {
            "all": {
                "qc": True,
                "assembly": True,
                "metaquast": True,
                "identification": True,
                "taxonomy": True,
                "mapping": True,
                "instrain": True,
                "function": True,
                "stats": True,
            },
            "stats": {
                "metaquast": True,
                "mapping": True,
                "instrain": False,
                "dramv": True,
            },
        },
    }


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


def make_workflow_profile(envs_dir: str | Path) -> dict:
    """``config/config.yaml``"""
    return {
        "software-deployment-method": "conda",
        "conda-prefix": str(Path(envs_dir).expanduser().resolve()),
        "show-failed-logs": True,
        "rerun-incomplete": True,
        "keep-going": True,
        "printshellcmds": True,
        "jobname": "{rule}.{jobid}",
        "max-jobs-per-second": 1,
    }


def make_slurm_profile() -> dict:
    """``config/slurm/config.yaml``"""
    return {
        "executor": "slurm",
        "default-resources": {
            "slurm_account": "",
            "runtime": "120 * attempt",
            "mem_mb": "7500 * attempt",
        },
        "max-status-checks-per-second": 0.1,
        "latency-wait": 60,
        "restart-times": 2,
        "jobs": 200,
        "slurm-efficiency-report": True,
    }


def write_profiles(rundir: str | Path, envs_dir: str | Path) -> None:
    """Write ``<rundir>/config/config.yaml`` and ``config/slurm/config.yaml``."""
    rundir = Path(rundir)
    slurm_dir = rundir / "config" / "slurm"
    slurm_dir.mkdir(parents=True, exist_ok=True)
    with open(rundir / "config" / "config.yaml", "w") as fh:
        yaml.safe_dump(make_workflow_profile(envs_dir), fh, sort_keys=False)
    with open(slurm_dir / "config.yaml", "w") as fh:
        yaml.safe_dump(make_slurm_profile(), fh, sort_keys=False)
