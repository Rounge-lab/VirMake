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
import yaml
from click.testing import CliRunner

from virmake import paths
from virmake.cli import cli

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        sys.platform != "linux",
        reason="snakemake workflow execution is Linux-only",
    ),
]


@pytest.fixture(scope="module")
def rundir(tmp_path_factory):
    base = tmp_path_factory.mktemp("integration")
    reads = base / "reads"
    reads.mkdir()
    for sample in ("A", "B"):
        for frac in ("1", "2"):
            (reads / f"{sample}_{frac}.fastq.gz").write_text("")
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
    """Rule names from Snakemake's dry-run Job stats table."""
    combined = proc.stdout + "\n" + proc.stderr
    lines = combined.splitlines()

    position = next(
        (index for index, line in enumerate(lines) if line.strip() == "Job stats:"),
        None,
    )
    if position is None:
        raise AssertionError(f"No Job stats table in dry-run output:\n{combined}")

    rules = set()
    started = False

    for line in lines[position + 1 :]:
        match = re.fullmatch(
            r"\s*([A-Za-z_][A-Za-z0-9_]*)\s+\d+(?:\s+\d+)*\s*",
            line,
        )
        if match:
            started = True
            name = match.group(1)
            if name == "total":
                break
            rules.add(name)
        elif started:
            break

    assert rules, f"Empty Job stats table in dry-run output:\n{combined}"
    return rules


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
        assert direct.returncode == 0, direct.stdout + direct.stderr
        assert until.returncode == 0, until.stdout + until.stderr
        assert _job_rules(direct) == _job_rules(until), (
            f"--until {target} and target {target} DAGs differ:\n"
            f"only-target: {sorted(_job_rules(direct))}\n"
            f"only-until:  {sorted(_job_rules(until))}"
        )


class TestDeclaredTrimmingInputs:
    @pytest.mark.parametrize(
        ("filenames", "directory_name"),
        [
            (("A_1.fastq", "A_2.fastq"), "plain"),
            (("A_R1.fastq.gz", "A_R2.fastq.gz"), "r-suffix"),
            (("A_1.fastq", "A_2.fastq.gz"), "mixed"),
            (("A_R1.fastq.gz", "A_R2.fastq"), "reads with spaces"),
        ],
    )
    def test_fastp_uses_declared_paths(self, tmp_path, filenames, directory_name):
        reads = tmp_path / directory_name
        reads.mkdir()
        for filename in filenames:
            (reads / filename).write_text("")

        run = tmp_path / "run"
        prepared = CliRunner().invoke(
            cli,
            [
                "prep",
                str(run),
                "--reads",
                str(reads),
                "--db-dir",
                str(tmp_path / "db"),
                "--envs-dir",
                str(tmp_path / "envs"),
            ],
        )
        assert prepared.exit_code == 0, prepared.output

        outputs = [
            str(run / "results/output/reads/reads_qc/A_1.fastq.gz"),
            str(run / "results/output/reads/reads_qc/A_2.fastq.gz"),
            str(run / "results/output/qc/fastp/A.json"),
        ]
        proc = _snakemake_dryrun(run, *outputs, "-n", "--printshellcmds")
        combined = proc.stdout + proc.stderr
        assert proc.returncode == 0, combined
        assert _job_rules(proc) == {"fastp_pe"}

        for filename in filenames:
            assert str(reads / filename) in combined

        assert "--disable_html" in combined
        assert "-h " not in combined
        assert "A.html" not in combined

        if " " in directory_name:
            for filename in filenames:
                assert f"'{reads / filename}'" in combined


class TestSuppliedAssemblyRouting:
    @pytest.mark.parametrize("with_reads", [False, True])
    def test_supplied_contigs_bypass_assembly(self, tmp_path, with_reads):
        contigs = tmp_path / "supplied contigs"
        contigs.mkdir()
        supplied = contigs / "A.fasta"
        original = b">A\nACGT\n"
        supplied.write_bytes(original)

        run = tmp_path / "run"
        args = [
            "prep",
            str(run),
            "--contigs",
            str(contigs),
            "--db-dir",
            str(tmp_path / "db"),
            "--envs-dir",
            str(tmp_path / "envs"),
        ]
        if with_reads:
            reads = tmp_path / "reads"
            reads.mkdir()
            for mate in ("1", "2"):
                (reads / f"A_{mate}.fastq").write_text("")
            args += ["--reads", str(reads)]

        prepared = CliRunner().invoke(cli, args)
        assert prepared.exit_code == 0, prepared.output

        # Isolate assembly completion from optional metaQUAST/database work.
        params_path = run / "config/params.yaml"
        params = yaml.safe_load(params_path.read_text())
        params["rule_inclusion"]["all"]["metaquast"] = False
        params_path.write_text(yaml.safe_dump(params, sort_keys=False))

        proc = _snakemake_dryrun(run, "ASSEMBLY", "-n", "--printshellcmds")
        combined = proc.stdout + proc.stderr
        assert proc.returncode == 0, combined

        rules = _job_rules(proc)
        assert "ASSEMBLY" in rules
        assert rules.isdisjoint({"metaSpades", "run_megahit", "fastp_se", "fastp_pe"})
        assert supplied.read_bytes() == original

        # Prove that ASSEMBLY requires the supplied file, rather than silently
        # falling back to an assembler or ignoring the declaration.
        supplied.unlink()
        try:
            missing = _snakemake_dryrun(run, "ASSEMBLY", "-n", "--printshellcmds")
            missing_output = missing.stdout + missing.stderr

            assert missing.returncode != 0, missing_output
            assert "MissingInputException" in missing_output, missing_output
            assert str(supplied) in missing_output, missing_output
        finally:
            supplied.write_bytes(original)

        assert supplied.read_bytes() == original
