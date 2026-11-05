"""CLI import-safety, config resolution, and asset-path tests."""

import subprocess
import sys

import click
import pytest
import yaml
from click.testing import CliRunner

from virmake import paths
from virmake.cli import cli, load_params, resolve_config_file


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

    @pytest.mark.parametrize("command", ["db", "run", "prep", "clean"])
    def test_command_help_without_config(self, command, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(cli, [command, "--help"])
        assert result.exit_code == 0

    @pytest.mark.parametrize("module", ["virmake.cli", "virmake.samples"])
    def test_module_import_never_pulls_in_snakemake(self, module, tmp_path):
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


class TestConfigResolution:
    def test_default_location_missing_raises(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        with pytest.raises(click.ClickException, match="config file not found"):
            resolve_config_file(None)

    def test_explicit_argument_resolves(self, tmp_path):
        config = tmp_path / "elsewhere" / "params.yaml"
        config.parent.mkdir()
        config.write_text(yaml.safe_dump({"path": {}}))
        assert resolve_config_file(str(config)) == config.resolve()

    def test_load_params_roundtrip(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        (tmp_path / "config").mkdir()
        params = {"path": {"benchmark": str(tmp_path / "results" / "benchmark")}}
        (tmp_path / "config" / "params.yaml").write_text(yaml.safe_dump(params))
        assert load_params(resolve_config_file(None)) == params

    def test_load_params_rejects_non_virmake_yaml(self, tmp_path):
        bad = tmp_path / "params.yaml"
        bad.write_text(yaml.safe_dump({"not": "virmake"}))
        with pytest.raises(click.ClickException, match="missing 'path'"):
            load_params(bad)

    def test_load_params_rejects_invalid_yaml(self, tmp_path):
        bad = tmp_path / "params.yaml"
        bad.write_text("foo: [unclosed")
        with pytest.raises(click.ClickException, match="could not parse"):
            load_params(bad)


class TestRunErrorPaths:
    def test_run_fails_cleanly_without_config(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        result = CliRunner().invoke(cli, ["run", "all"])
        assert result.exit_code != 0
        assert "config file not found" in result.output
