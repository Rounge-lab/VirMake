import importlib.resources
import subprocess
import sys
from pathlib import Path
from uuid import uuid4

import pytest

from virmake.config import make_params, make_results_dirs, write_params
from virmake.samples import write_sample_sheet


def _package_root() -> Path:
    return Path(str(importlib.resources.files("virmake")))


@pytest.fixture(scope="session")
def snakemake_available():
    """Skip when this interpreter cannot run snakemake (e.g. a dev-only box)."""
    proc = subprocess.run(
        [sys.executable, "-m", "snakemake", "--version"],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.skip(
            "snakemake not usable in this environment "
            f"(`python -m snakemake --version` failed: {proc.stderr.strip()})"
        )
    return proc.stdout.strip()


@pytest.fixture
def run_dir(tmp_path):
    """Minimal run directory using the real config and sample-sheet writers."""
    params = make_params(
        rundir=tmp_path,
        db_dir=tmp_path / "databases",
    )
    params["run_id"] = str(uuid4())

    write_params(tmp_path, params)
    make_results_dirs(tmp_path)
    write_sample_sheet(tmp_path / "samples.tsv", ())
    return tmp_path


def test_workflow_parses_from_foreign_cwd(snakemake_available, run_dir):
    snakefile = _package_root() / "workflow" / "Snakefile"
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(snakefile),
        "--directory",
        str(run_dir),
        "--configfile",
        str(run_dir / "config" / "params.yaml"),
        "--list-rules",
    ]
    proc = subprocess.run(cmd, cwd=run_dir, capture_output=True, text=True, timeout=600)
    assert proc.returncode == 0, f"snakemake stderr:\n{proc.stderr}"
    # One rule per module: proves the `use rule *` wiring loaded everything
    # (this also parse-exercises the unreviewed qc/taxonomy/setup_db modules).
    for rule in (
        "metaSpades",
        "combine_coverage",
        "dramv_annotate",
        "get_sample_stats",
    ):
        assert rule in proc.stdout, f"rule {rule} missing from --list-rules output"
