from pathlib import Path

import click
import platformdirs
import pytest
import yaml
from click.testing import CliRunner

from virmake import paths
from virmake.config import (
    _load_template,
    default_data_root,
    default_db_dir,
    default_envs_dir,
    make_params,
    make_path_block,
    make_results_dirs,
    make_slurm_profile,
    make_workflow_profile,
    resolve_resource_dir,
    write_params,
)


class TestDefaultDirs:
    def test_defaults_under_user_data_root(self):
        assert default_data_root() == Path(platformdirs.user_data_dir("virmake"))
        assert default_db_dir() == default_data_root() / "databases"
        assert default_envs_dir() == default_data_root() / "envs"


class TestResolveResourceDir:
    def test_explicit_wins_without_interaction(self, monkeypatch):
        monkeypatch.setattr(
            click, "confirm", lambda *a, **k: pytest.fail("must not confirm")
        )
        result = resolve_resource_dir(
            "/opt/dbs", default=Path("/default"), description="databases"
        )
        assert result == Path("/opt/dbs")

    def test_explicit_expands_user(self, monkeypatch, tmp_path):
        monkeypatch.setenv("USERPROFILE", str(tmp_path))
        monkeypatch.setenv("HOME", str(tmp_path))
        result = resolve_resource_dir(
            "~/dbs", default=Path("/default"), description="databases"
        )
        assert result == tmp_path / "dbs"

    def test_assume_yes_uses_default(self, monkeypatch):
        monkeypatch.setattr(
            click, "confirm", lambda *a, **k: pytest.fail("must not confirm")
        )
        result = resolve_resource_dir(
            None, default=Path("/default"), description="databases", assume_yes=True
        )
        assert result == Path("/default")

    def test_interactive_confirm_accepts_default(self, monkeypatch):
        monkeypatch.setattr(click, "confirm", lambda *a, **k: True)
        result = resolve_resource_dir(
            None, default=Path("/default"), description="databases"
        )
        assert result == Path("/default")

    def test_interactive_decline_prompts_for_path(self, monkeypatch):
        monkeypatch.setattr(click, "confirm", lambda *a, **k: False)
        monkeypatch.setattr(click, "prompt", lambda *a, **k: Path("/opt/dbs"))
        result = resolve_resource_dir(
            None, default=Path("/default"), description="databases"
        )
        assert result == Path("/opt/dbs")

    def test_end_to_end_via_stdin(self, tmp_path):
        @click.command()
        def _wrapper():
            resolved = resolve_resource_dir(
                None, default=Path("/default"), description="databases"
            )
            click.echo(f"RESOLVED:{resolved}")

        target = tmp_path / "dbs"
        runner = CliRunner()
        accept = runner.invoke(_wrapper, input="y\n")
        assert accept.exit_code == 0
        assert "RESOLVED:/default" in accept.output.replace("\\", "/")

        decline = runner.invoke(_wrapper, input=f"n\n{target}\n")
        assert decline.exit_code == 0
        assert f"RESOLVED:{target}" in decline.output


class TestTemplates:
    """The shipped templates are the source of truth — pin their schema."""

    def test_params_template_parses(self):
        params = _load_template("params.yaml")
        assert set(params) == {
            "slurm_account",
            "keep_qc_reads",
            "assembler",
            "identifier",
            "trim_percentage",
            "min_coverage",
            "min_contig_size",
            "vibrant",
            "virsorter2",
            "vcontact3",
            "path",
            "dereplication",
            "quality_threshold",
            "threads",
            "job_type",
            "memory",
            "time",
            "rule_inclusion",
        }
        assert params["path"] == {}  # injected by make_params

    def test_rule_inclusion_structure(self):
        # The Snakefile's ALL rule reads rule_inclusion.all.* — pin exactly.
        params = _load_template("params.yaml")
        assert set(params["rule_inclusion"]["all"]) == {
            "qc",
            "assembly",
            "metaquast",
            "identification",
            "taxonomy",
            "mapping",
            "instrain",
            "function",
            "stats",
        }
        assert set(params["rule_inclusion"]["stats"]) == {
            "metaquast",
            "mapping",
            "instrain",
            "dramv",
        }

    def test_profile_templates_parse(self):
        wf = _load_template("workflow-profile.yaml")
        assert set(wf) == {
            "software-deployment-method",
            "conda-prefix",
            "show-failed-logs",
            "rerun-incomplete",
            "keep-going",
            "printshellcmds",
            "jobname",
            "max-jobs-per-second",
        }
        assert wf["software-deployment-method"] == "conda"
        assert wf["conda-prefix"] == ""  # injected by make_workflow_profile
        slurm = _load_template("slurm-profile.yaml")
        assert slurm["executor"] == "slurm"


class TestMakeParams:
    def test_sentinel_defaults_flow_through(self, tmp_path):
        # The template is the source of truth; these pins guard against
        # accidental edits to behaviorally important defaults.
        params = make_params(tmp_path, tmp_path / "db")
        assert params["assembler"] == "metaSpades"
        assert params["identifier"] == "virsorter2"
        assert params["quality_threshold"] == "medium"
        assert params["threads"] == 24
        assert params["keep_qc_reads"] is True

    def test_path_block_golden(self, tmp_path):
        rundir = tmp_path / "run"
        db_dir = tmp_path / "db"
        reads = tmp_path / "reads"
        contigs = tmp_path / "contigs"
        assert make_path_block(rundir, db_dir) == {
            "virmake": str(paths.package_root()),
            "envs": str(paths.envs_dir()),
            "output": str(rundir / "results" / "output"),
            "log": str(rundir / "results" / "log"),
            "benchmark": str(rundir / "results" / "benchmark"),
            "temp": str(rundir / "results" / "temp"),
            "scripts": str(paths.scripts_dir()),
            "samples": str(rundir / "samples.tsv"),
            "input_reads": "",
            "input_contigs": "",
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
        block = make_path_block(rundir, db_dir, reads=str(reads), contigs=str(contigs))
        assert block["input_reads"] == str(reads)
        assert block["input_contigs"] == str(contigs)

    def test_make_params_merges_template_and_path(self, tmp_path):
        template = _load_template("params.yaml")
        params = make_params(tmp_path, tmp_path / "db")
        assert params["path"] == make_path_block(tmp_path, tmp_path / "db")
        without_path = {k: v for k, v in params.items() if k != "path"}
        assert without_path == {k: v for k, v in template.items() if k != "path"}


class TestProfiles:
    def test_conda_prefix_injected(self, tmp_path):
        envs = tmp_path / "envs"
        profile = make_workflow_profile(envs)
        assert profile["conda-prefix"] == str(envs)

    def test_dropped_keys_stay_dropped(self, tmp_path):
        profile = make_workflow_profile(tmp_path / "envs")
        assert "nolock" not in profile
        assert "conda-frontend" not in profile
        assert "use-conda" not in profile

    def test_slurm_profile_is_template(self):
        assert make_slurm_profile() == _load_template("slurm-profile.yaml")


class TestWriters:
    def test_write_params_creates_parents_and_roundtrips(self, tmp_path):
        rundir = tmp_path / "run"
        params = make_params(rundir, tmp_path / "db")
        path = write_params(rundir, params)
        assert path == rundir / "config" / "params.yaml"
        assert yaml.safe_load(path.read_text()) == params

    def test_make_results_dirs(self, tmp_path):
        rundir = tmp_path / "run"
        make_results_dirs(rundir)
        for sub in ("output", "log", "benchmark", "temp"):
            assert (rundir / "results" / sub).is_dir()
