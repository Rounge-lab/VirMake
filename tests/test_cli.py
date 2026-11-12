"""CLI import-safety, run-dir resolution, and asset-path tests."""

import subprocess
import sys

import click
import pytest
import yaml
from click.testing import CliRunner

from virmake import __version__, paths
from virmake.cli import cli, load_params, resolve_rundir


def test_version():
    result = CliRunner().invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


class TestPaths:
    def test_all_assets_resolve_into_the_package(self):
        assert paths.snakefile().is_file()
        assert (paths.envs_dir() / "DRAMv.post-deploy.sh").is_file()
        assert paths.scripts_dir().is_dir()
        for name in (
            "contamination_dummy.tsv",
            "quality_summary_dummy.tsv",
            "viral-boundary-dummy.tsv",
        ):
            assert (paths.dummy_resources_dir() / name).is_file()


class TestImportSafety:
    """Regression guards for the load-config-at-import bug."""

    def test_group_help_without_config(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)  # deliberately no config/params.yaml
        result = CliRunner().invoke(cli, ["--help"])
        assert result.exit_code == 0
        assert "VirMake" in result.output

    @pytest.mark.parametrize("command", ["prep", "envs", "db", "run"])
    def test_command_help_without_config(self, command, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(cli, [command, "--help"])
        assert result.exit_code == 0

    @pytest.mark.parametrize("module", ["virmake.cli", "virmake.samples"])
    def test_module_import_never_pulls_in_snakemake(self, module, tmp_path):
        # Subprocess guard: the differential-test oracle in test_samples.py
        # legitimately imports snakemake into the pytest process, so this
        # D12 check must run in a clean interpreter.
        code = (
            "import sys\n"
            f"import {module}\n"
            "mods = [m for m in sys.modules if m.startswith('snakemake')]\n"
            "assert not mods, f'snakemake modules imported: {mods}'\n"
        )
        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            cwd=tmp_path,
        )
        assert proc.returncode == 0, proc.stderr


class TestResolveRundir:
    def _minimal_rundir(self, tmp_path):
        rundir = tmp_path / "run"
        (rundir / "config").mkdir(parents=True)
        (rundir / "config" / "params.yaml").write_text("path: {}\n")
        return rundir

    def test_valid_rundir_resolves(self, tmp_path):
        rundir = self._minimal_rundir(tmp_path)
        assert resolve_rundir(rundir) == rundir.resolve()

    def test_missing_rundir_raises_with_prep_hint(self, tmp_path):
        with pytest.raises(click.ClickException, match="virmake prep"):
            resolve_rundir(tmp_path / "nope")

    def test_load_params_roundtrip(self, tmp_path):
        rundir = self._minimal_rundir(tmp_path)
        assert load_params(rundir) == {"path": {}}

    def test_load_params_rejects_non_virmake_yaml(self, tmp_path):
        rundir = tmp_path / "run"
        (rundir / "config").mkdir(parents=True)
        (rundir / "config" / "params.yaml").write_text(yaml.safe_dump({"nope": 1}))
        with pytest.raises(click.ClickException, match="missing 'path'"):
            load_params(rundir)


class TestRunErrorPaths:
    def test_run_fails_cleanly_on_unprepared_dir(self, tmp_path):
        result = CliRunner().invoke(cli, ["run", str(tmp_path / "nope")])
        assert result.exit_code != 0
        assert "virmake prep" in result.output
