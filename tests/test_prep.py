from pathlib import Path

import pandas as pd
import platformdirs
import pytest
import yaml
from click.testing import CliRunner

from virmake.cli import cli


@pytest.fixture
def patched_defaults(monkeypatch, tmp_path):
    """Point per-user resource defaults at tmp_path/'appdata'."""
    appdata = tmp_path / "appdata"
    monkeypatch.setattr(
        platformdirs, "user_data_dir", lambda appname, **kw: str(appdata)
    )
    return appdata


@pytest.fixture
def reads_dir(tmp_path):
    d = tmp_path / "reads"
    d.mkdir()
    for name in ("A_1.fastq", "A_2.fastq"):
        (d / name).write_text("")
    return d


@pytest.fixture
def contigs_dir(tmp_path):
    d = tmp_path / "contigs"
    d.mkdir()
    (d / "A.fasta").write_text("")
    return d


class TestPrepTree:
    def test_creates_full_tree(self, tmp_path, patched_defaults):
        rundir = tmp_path / "run"
        result = CliRunner().invoke(cli, ["prep", str(rundir), "--yes"])
        assert result.exit_code == 0, result.output
        assert (rundir / "samples.tsv").is_file()
        assert (rundir / "config" / "params.yaml").is_file()
        assert (rundir / "config" / "config.yaml").is_file()
        assert (rundir / "config" / "slurm" / "config.yaml").is_file()
        for sub in ("output", "log", "benchmark", "temp"):
            assert (rundir / "results" / sub).is_dir()
        # defaults recorded (D5)
        params = yaml.safe_load((rundir / "config" / "params.yaml").read_text())
        assert params["path"]["database"]["checkv"] == str(
            patched_defaults / "databases" / "checkv"
        )
        profile = yaml.safe_load((rundir / "config" / "config.yaml").read_text())
        assert profile["conda-prefix"] == str(patched_defaults / "envs")
        # resource roots created eagerly
        assert (patched_defaults / "databases").is_dir()
        assert (patched_defaults / "envs").is_dir()

    def test_hints_printed(self, tmp_path, patched_defaults):
        result = CliRunner().invoke(cli, ["prep", str(tmp_path / "run"), "--yes"])
        assert "virmake envs" in result.output
        assert "virmake db" in result.output
        assert "virmake run" in result.output


class TestPrepForce:
    def test_refuses_existing_nonempty(self, tmp_path, patched_defaults):
        rundir = tmp_path / "run"
        rundir.mkdir()
        (rundir / "config").mkdir()
        result = CliRunner().invoke(cli, ["prep", str(rundir), "--yes"])
        assert result.exit_code != 0
        assert "already exists" in result.output

    def test_force_regenerates_and_never_touches_results(
        self, tmp_path, patched_defaults
    ):
        rundir = tmp_path / "run"
        CliRunner().invoke(cli, ["prep", str(rundir), "--yes"])
        # user edits a parameter; produces a result file
        params_path = rundir / "config" / "params.yaml"
        params = yaml.safe_load(params_path.read_text())
        params["threads"] = 99
        params_path.write_text(yaml.safe_dump(params))
        marker = rundir / "results" / "output" / "precious.txt"
        marker.write_text("do not delete")
        result = CliRunner().invoke(cli, ["prep", str(rundir), "--yes", "--force"])
        assert result.exit_code == 0, result.output
        assert yaml.safe_load(params_path.read_text())["threads"] == 24
        assert marker.read_text() == "do not delete"


class TestPrepResourceDirs:
    def test_explicit_dirs_no_interaction(self, tmp_path, patched_defaults, reads_dir):
        db, envs = tmp_path / "db", tmp_path / "envs"
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--reads",
                str(reads_dir),
                "--db-dir",
                str(db),
                "--envs-dir",
                str(envs),
            ],
        )
        assert result.exit_code == 0, result.output
        params = yaml.safe_load(
            (tmp_path / "run" / "config" / "params.yaml").read_text()
        )
        assert params["path"]["database"]["checkv"] == str(db / "checkv")

    def test_interactive_confirm_defaults(self, tmp_path, patched_defaults):
        result = CliRunner().invoke(
            cli, ["prep", str(tmp_path / "run")], input="y\ny\n"
        )
        assert result.exit_code == 0, result.output
        params = yaml.safe_load(
            (tmp_path / "run" / "config" / "params.yaml").read_text()
        )
        assert params["path"]["database"]["checkv"].startswith(
            str(patched_defaults / "databases")
        )

    def test_interactive_decline_prompts_for_paths(self, tmp_path, patched_defaults):
        db, envs = tmp_path / "db", tmp_path / "envs"
        result = CliRunner().invoke(
            cli,
            ["prep", str(tmp_path / "run")],
            input=f"n\n{db}\nn\n{envs}\n",
        )
        assert result.exit_code == 0, result.output
        assert (db).is_dir()
        assert (envs).is_dir()
        profile = yaml.safe_load(
            (tmp_path / "run" / "config" / "config.yaml").read_text()
        )
        assert profile["conda-prefix"] == str(envs)


class TestPrepInputs:
    def test_reads_and_contigs_populate_samples(
        self, tmp_path, patched_defaults, reads_dir, contigs_dir
    ):
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--reads",
                str(reads_dir),
                "--contigs",
                str(contigs_dir),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert result.exit_code == 0, result.output
        table = pd.read_csv(tmp_path / "run" / "samples.tsv", sep="\t")
        assert table["sample_id"].tolist() == ["A"]
        assert Path(table["r1"].iloc[0]).name == "A_1.fastq"
        assert Path(table["r2"].iloc[0]).name == "A_2.fastq"
        assert table["qc_r1"].isna().all()  # no qc reads given
        assert Path(table["contigs"].iloc[0]).name == "A.fasta"

    def test_missing_reads_dir_rejected(self, tmp_path, patched_defaults):
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--reads",
                str(tmp_path / "nope"),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert result.exit_code != 0

    def test_empty_table_warns(self, tmp_path, patched_defaults):
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert result.exit_code == 0
        assert "header-only" in result.output

    def test_plain_fastq_reads_warn(self, tmp_path, patched_defaults):
        plain = tmp_path / "plain"
        plain.mkdir()
        (plain / "A_1.fastq").write_text("")
        (plain / "A_2.fastq").write_text("")
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--reads",
                str(plain),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert result.exit_code == 0
        assert "not gzipped" in result.output

    def test_gzipped_reads_do_not_warn(self, tmp_path, patched_defaults, reads_dir):
        # fixture writes plain .fastq; gz variant here
        gz = tmp_path / "gz"
        gz.mkdir()
        (gz / "A_1.fastq.gz").write_text("")
        (gz / "A_2.fastq.gz").write_text("")
        result = CliRunner().invoke(
            cli,
            [
                "prep",
                str(tmp_path / "run"),
                "--reads",
                str(gz),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert result.exit_code == 0
        assert "not gzipped" not in result.output


class TestPrepTestData:
    @pytest.fixture
    def fake_download(self, monkeypatch):
        def _fake(destination):
            destination = destination  # Path from cli
            destination.mkdir(parents=True, exist_ok=True)
            for frac in ("1", "2"):
                (destination / f"synthetic_{frac}.fastq.gz").write_text("")
            return []

        monkeypatch.setattr("virmake.cli.download_test_data", _fake)
        return _fake

    def test_test_data_populates_qc_reads(
        self, tmp_path, patched_defaults, fake_download
    ):
        rundir = tmp_path / "run"
        result = CliRunner().invoke(cli, ["prep", str(rundir), "--test-data", "--yes"])
        assert result.exit_code == 0, result.output
        assert (rundir / "test_data").is_dir()
        # qc reads recorded; header check via pandas
        table = pd.read_csv(rundir / "samples.tsv", sep="\t")
        assert table["sample_id"].tolist() == ["synthetic"]
        assert table["qc_r1"].iloc[0].endswith("synthetic_1.fastq.gz")

    def test_test_data_conflicts_with_reads(
        self, tmp_path, patched_defaults, reads_dir, fake_download
    ):
        result = CliRunner().invoke(
            cli,
            ["prep", str(tmp_path / "run"), "--test-data", "--reads", str(reads_dir)],
        )
        assert result.exit_code != 0
        assert "cannot be combined" in result.output
