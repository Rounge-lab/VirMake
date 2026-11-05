from pathlib import Path

import click
import platformdirs
import pytest
import yaml
from click.testing import CliRunner

from virmake import paths
from virmake.config import (
    default_data_root,
    default_db_dir,
    default_envs_dir,
    make_params,
    make_results_dirs,
    make_slurm_profile,
    make_workflow_profile,
    resolve_resource_dir,
    write_params,
    write_profiles,
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


class TestMakeParams:
    def test_params_golden(self, tmp_path):
        rundir = tmp_path / "run"
        db_dir = tmp_path / "db"
        params = make_params(rundir, db_dir)
        expected = {
            "slurm_account": "default",
            "keep_qc_reads": True,
            "assembler": "metaSpades",
            "identifier": "virsorter2",
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
            "job_type": {"small": "normal", "normal": "normal", "big": "normal"},
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
        assert params == expected

    def test_input_paths_recorded_absolute(self, tmp_path):
        reads = tmp_path / "reads"
        contigs = tmp_path / "contigs"
        params = make_params(
            tmp_path, tmp_path / "db", reads=str(reads), contigs=str(contigs)
        )
        assert params["path"]["input_reads"] == str(reads)
        assert params["path"]["input_contigs"] == str(contigs)

    def test_relative_input_recorded_absolute(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        params = make_params(tmp_path / "run", tmp_path / "db", reads="reads")
        assert params["path"]["input_reads"] == str(tmp_path / "reads")


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


class TestProfiles:
    def test_workflow_profile_golden(self, tmp_path):
        envs = tmp_path / "envs"
        assert make_workflow_profile(envs) == {
            "software-deployment-method": "conda",
            "conda-prefix": str(envs),
            "show-failed-logs": True,
            "rerun-incomplete": True,
            "keep-going": True,
            "printshellcmds": True,
            "jobname": "{rule}.{jobid}",
            "max-jobs-per-second": 1,
        }

    def test_dropped_keys_stay_dropped(self, tmp_path):
        profile = make_workflow_profile(tmp_path / "envs")
        assert "nolock" not in profile
        assert "conda-frontend" not in profile

    def test_slurm_profile_golden(self):
        assert make_slurm_profile() == {
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

    def test_write_profiles_roundtrip(self, tmp_path):
        rundir = tmp_path / "run"
        envs = tmp_path / "envs"
        write_profiles(rundir, envs)
        wf = yaml.safe_load((rundir / "config" / "config.yaml").read_text())
        slurm = yaml.safe_load(
            (rundir / "config" / "slurm" / "config.yaml").read_text()
        )
        assert wf == make_workflow_profile(envs)
        assert slurm == make_slurm_profile()
