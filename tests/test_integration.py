"""Integration tests: Snakemake dryruns over a synthetic run dir.

- the full `all` DAG closes without any downloaded database (db flags are
  rule outputs — review finding from setup_db.smk);
- `virmake db -n` builds the SETUP_DB sub-DAG under explicit deployment flags;
- the workflow profile is accepted by the pinned snakemake (strict key check);
- legacy `--until TARGET` is DAG-equivalent to explicit `TARGET`.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from virmake import paths
from virmake.cli import cli

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        sys.platform == "win32", reason="snakemake unsupported on Windows"
    ),
]


@pytest.fixture(scope="module")
def rundir(tmp_path_factory):
    base = tmp_path_factory.mktemp("integration")
    reads = base / "reads"
    reads.mkdir()
    for sample in ("A", "B"):
        for frac in ("1", "2"):
            (reads / f"{sample}_{frac}.fastq").write_text("")
    run = base / "run"
    result = CliRunner().invoke(
        cli,
        [
            "prep",
            str(run),
            "--reads",
            str(reads),
            "--db-dir",
            str(base / "db"),
            "--envs-dir",
            str(base / "envs"),
        ],
    )
    assert result.exit_code == 0, result.output
    return run


def _snakemake_dryrun(rundir: Path, *args: str) -> subprocess.CompletedProcess:
    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(paths.snakefile()),
        "--directory",
        str(rundir),
        "--configfile",
        str(rundir / "config" / "params.yaml"),
        "--cores",
        "1",
        *args,
    ]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=600)


def _job_rules(proc: subprocess.CompletedProcess) -> set[str]:
    """Rule names from the dryrun's Job stats table."""
    return set(re.findall(r"^(\S+)\s+\d+\s+\d+", proc.stdout, re.MULTILINE))


class TestCliDryruns:
    def test_run_all_dryrun_closes_dag(self, rundir):
        result = CliRunner().invoke(cli, ["run", str(rundir), "all", "-n"])
        assert result.exit_code == 0, result.output

    def test_run_qc_dryrun_accepts_workflow_profile(self, rundir):
        result = CliRunner().invoke(cli, ["run", str(rundir), "qc", "-n"])
        assert result.exit_code == 0, result.output

    def test_db_dryrun_targets_setup_db(self, rundir):
        result = CliRunner().invoke(cli, ["db", str(rundir), "-n"])
        assert result.exit_code == 0, result.output


class TestUntilEquivalence:
    @pytest.mark.parametrize("target", ["QC", "MAPPING"])
    def test_explicit_target_matches_legacy_until(self, rundir, target):
        direct = _snakemake_dryrun(rundir, target, "-n")
        until = _snakemake_dryrun(rundir, "--until", target, "-n")
        assert direct.returncode == 0, direct.stderr
        assert until.returncode == 0, until.stderr
        assert _job_rules(direct) == _job_rules(until), (
            f"--until {target} and target {target} DAGs differ:\n"
            f"only-target: {sorted(_job_rules(direct))}\n"
            f"only-until:  {sorted(_job_rules(until))}"
        )
