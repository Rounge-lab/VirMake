"""Tests for envs/db/run commands."""

import subprocess
import sys

import pytest
import yaml
from click.testing import CliRunner

from virmake.cli import cli
from virmake.config import make_params, make_results_dirs, write_params, write_profiles


def _combined(result):
    """stdout + stderr across click versions (8.1 mixes; 8.2+ separates)."""
    out = result.output
    try:
        out += result.stderr
    except (ValueError, AttributeError):  # pragma: no cover
        pass
    return out


@pytest.fixture
def rundir(tmp_path):
    rundir = tmp_path / "run"
    write_params(rundir, make_params(rundir, tmp_path / "db"))
    write_profiles(rundir, tmp_path / "envs")
    make_results_dirs(rundir)
    header = "sample_id\tr1\tr2\tqc_r1\tqc_r2\tcontigs"
    row = "A\t/reads/A_1.fastq\t/reads/A_2.fastq\t\t\t"
    (rundir / "samples.tsv").write_text(header + "\n" + row + "\n")
    return rundir


@pytest.fixture
def snakemake_calls(monkeypatch):
    calls = []

    def fake(cmd, **kwargs):
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, returncode=0)

    monkeypatch.setattr("virmake.cli.subprocess.run", fake)
    return calls


class TestEnvs:
    def test_builds_envs_with_create_only(self, rundir, tmp_path, snakemake_calls):
        result = CliRunner().invoke(cli, ["envs", str(rundir)])
        assert result.exit_code == 0, _combined(result)
        cmd = snakemake_calls[0]
        assert cmd[0].endswith("python") or cmd[0] == sys.executable
        assert cmd[1:3] == ["-m", "snakemake"]
        assert "--conda-create-envs-only" in cmd
        assert "--software-deployment-method" in cmd
        i = cmd.index("--conda-prefix")
        assert cmd[i + 1] == str(tmp_path / "envs")
        i = cmd.index("--directory")
        assert cmd[i + 1] == str(rundir)
        i = cmd.index("--cores")
        assert cmd[i + 1] == "24"

    def test_clean_removes_prefix(self, rundir, tmp_path):
        envs_dir = tmp_path / "envs"
        (envs_dir / "abc123").mkdir(parents=True)
        result = CliRunner().invoke(cli, ["envs", str(rundir), "--clean", "--yes"])
        assert result.exit_code == 0, _combined(result)
        assert envs_dir.is_dir() and not any(envs_dir.iterdir())

    def test_clean_aborts_without_yes(self, rundir, tmp_path):
        envs_dir = tmp_path / "envs"
        (envs_dir / "abc").mkdir(parents=True)
        result = CliRunner().invoke(cli, ["envs", str(rundir), "--clean"], input="n\n")
        assert result.exit_code != 0
        assert any(envs_dir.iterdir())


class TestDb:
    def test_targets_setup_db(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["db", str(rundir)])
        assert result.exit_code == 0, _combined(result)
        cmd = snakemake_calls[0]
        assert "SETUP_DB" in cmd
        assert "--conda-create-envs-only" not in cmd
        assert "--conda-prefix" in cmd

    def test_dryrun_flag(self, rundir, snakemake_calls):
        CliRunner().invoke(cli, ["db", str(rundir), "-n"])
        assert "-n" in snakemake_calls[0]

    def test_clean_removes_db_root(self, rundir, tmp_path):
        db_dir = tmp_path / "db"
        (db_dir / "checkv").mkdir(parents=True)
        (db_dir / "checkv" / "flag").write_text("")
        result = CliRunner().invoke(cli, ["db", str(rundir), "--clean", "--yes"])
        assert result.exit_code == 0, _combined(result)
        assert db_dir.is_dir() and not any(db_dir.iterdir())


class TestRun:
    def test_command_shape(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert result.exit_code == 0, _combined(result)
        cmd = snakemake_calls[0]
        assert cmd[1:3] == ["-m", "snakemake"]
        for flag in (
            "--snakefile",
            "--directory",
            "--configfile",
            "--workflow-profile",
        ):
            assert flag in cmd
        i = cmd.index("--workflow-profile")
        assert cmd[i + 1] == str(rundir / "config")
        assert "ALL" in cmd
        i = cmd.index("--cores")
        assert cmd[i + 1] == "1"
        assert (rundir / "results" / "benchmark" / "total_time.txt").is_file()

    def test_target_maps_to_rule(self, rundir, snakemake_calls):
        CliRunner().invoke(cli, ["run", str(rundir), "assembly"])
        assert "ASSEMBLY" in snakemake_calls[0]

    def test_dryrun_skips_benchmark(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["run", str(rundir), "-n"])
        assert result.exit_code == 0, _combined(result)
        assert "-n" in snakemake_calls[0]
        assert not (rundir / "results" / "benchmark" / "total_time.txt").exists()

    def test_slurm_profile(self, rundir, snakemake_calls):
        CliRunner().invoke(cli, ["run", str(rundir), "--slurm"])
        cmd = snakemake_calls[0]
        i = cmd.index("--profile")
        assert cmd[i + 1] == str(rundir / "config" / "slurm")

    def test_extra_args_passed_through(self, rundir, snakemake_calls):
        CliRunner().invoke(cli, ["run", str(rundir), "all", "--", "--touch"])
        assert "--touch" in snakemake_calls[0]

    def test_failure_raises(self, rundir, monkeypatch):
        def failing(cmd, **kwargs):
            return subprocess.CompletedProcess(cmd, returncode=1)

        monkeypatch.setattr("virmake.cli.subprocess.run", failing)
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert result.exit_code != 0
        assert "failed" in _combined(result)

    def test_unprepared_dir_gives_prep_hint(self, tmp_path):
        result = CliRunner().invoke(cli, ["run", str(tmp_path / "nope")])
        assert result.exit_code != 0
        assert "virmake prep" in _combined(result)


class TestPreflight:
    def test_db_warning_when_flags_missing(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert "virmake db" in _combined(result)

    def test_envs_warning_when_prefix_empty(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert "virmake envs" in _combined(result)

    @pytest.mark.skipif(sys.platform == "darwin", reason="darwin warning is expected")
    def test_no_warnings_when_ready(self, rundir, tmp_path, snakemake_calls):
        db_dir, envs_dir = tmp_path / "db", tmp_path / "envs"
        for name in ("virsorter2", "checkv"):  # identifier + quality step
            d = db_dir / name
            d.mkdir(parents=True, exist_ok=True)
            (d / "flag").write_text("")
        (envs_dir / "somehash").mkdir(parents=True, exist_ok=True)
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert result.exit_code == 0
        assert "WARNING" not in _combined(result)

    def test_header_only_samples_warns(self, rundir, snakemake_calls):
        (rundir / "samples.tsv").write_text(
            "sample_id\tr1\tr2\tqc_r1\tqc_r2\tcontigs\n"
        )
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert result.exit_code == 0
        assert "no samples" in _combined(result)

    def test_missing_samples_aborts(self, rundir, snakemake_calls):
        (rundir / "samples.tsv").unlink()
        result = CliRunner().invoke(cli, ["run", str(rundir)])
        assert result.exit_code != 0
        assert "samples table" in _combined(result)


class TestSlurmAccount:
    def test_account_written_to_profile(self, rundir, snakemake_calls):
        result = CliRunner().invoke(
            cli, ["run", str(rundir), "--slurm", "--account", "ec34"]
        )
        assert result.exit_code == 0, _combined(result)
        profile = yaml.safe_load(
            (rundir / "config" / "slurm" / "config.yaml").read_text()
        )
        assert profile["default-resources"]["slurm_account"] == "ec34"
        assert "SLURM account recorded" in _combined(result)

    def test_empty_account_warns_but_submits(self, rundir, snakemake_calls):
        result = CliRunner().invoke(cli, ["run", str(rundir), "--slurm"])
        assert result.exit_code == 0, _combined(result)
        assert "SLURM account" in _combined(result)
        assert "--profile" in snakemake_calls[0]

    def test_recorded_account_silences_warning(self, rundir, snakemake_calls):
        path = rundir / "config" / "slurm" / "config.yaml"
        profile = yaml.safe_load(path.read_text())
        profile["default-resources"]["slurm_account"] = "ec34"
        path.write_text(yaml.safe_dump(profile, sort_keys=False))
        result = CliRunner().invoke(cli, ["run", str(rundir), "--slurm"])
        assert result.exit_code == 0, _combined(result)
        assert "SLURM account" not in _combined(result)

    def test_account_requires_slurm(self, rundir):
        result = CliRunner().invoke(cli, ["run", str(rundir), "--account", "ec34"])
        assert result.exit_code != 0
        assert "requires --slurm" in _combined(result)
