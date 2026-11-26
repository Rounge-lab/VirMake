"""FastQC artifact extraction, independent of external tools."""

import os
import subprocess
import sys
from zipfile import BadZipFile, ZipFile

import pytest

from virmake import paths
from virmake.workflow.scripts.fastqc import extract_fastqc_data


def test_extracts_only_text_artifact(tmp_path):
    archive = tmp_path / "A_fastqc.zip"
    data = b"##FastQC\t0.12.1\n"
    with ZipFile(archive, "w") as handle:
        handle.writestr("A_fastqc/fastqc_data.txt", data)
        handle.writestr("A_fastqc/fastqc_report.html", "<html></html>")
        handle.writestr("A_fastqc/summary.txt", "summary")

    destination = tmp_path / "results/qc/fastqc/A_1_fastqc_data.txt"
    extract_fastqc_data(archive, destination)

    assert destination.read_bytes() == data
    assert list(destination.parent.iterdir()) == [destination]


@pytest.mark.parametrize(
    "members",
    [
        [],
        ["A_fastqc/summary.txt"],
        ["A_fastqc/fastqc_data.txt", "B_fastqc/fastqc_data.txt"],
    ],
)
def test_missing_or_ambiguous_member_preserves_destination(tmp_path, members):
    archive = tmp_path / "fastqc.zip"
    with ZipFile(archive, "w") as handle:
        for member in members:
            handle.writestr(member, "data")

    destination = tmp_path / "fastqc_data.txt"
    destination.write_bytes(b"existing artifact\n")

    with pytest.raises(ValueError, match="expected exactly one"):
        extract_fastqc_data(archive, destination)

    assert destination.read_bytes() == b"existing artifact\n"


def test_corrupt_archive_preserves_destination(tmp_path):
    archive = tmp_path / "fastqc.zip"
    archive.write_bytes(b"not a zip")
    destination = tmp_path / "fastqc_data.txt"
    destination.write_bytes(b"existing artifact\n")

    with pytest.raises(BadZipFile):
        extract_fastqc_data(archive, destination)

    assert destination.read_bytes() == b"existing artifact\n"


def test_archive_paths_are_not_extracted(tmp_path):
    archive = tmp_path / "fastqc.zip"
    with ZipFile(archive, "w") as handle:
        handle.writestr("../../fastqc_data.txt", b"artifact")
        handle.writestr("../../unexpected.txt", b"must not be extracted")

    destination = tmp_path / "published/data.txt"
    extract_fastqc_data(archive, destination)

    assert destination.read_bytes() == b"artifact"
    assert not (tmp_path / "fastqc_data.txt").exists()
    assert not (tmp_path / "unexpected.txt").exists()
    assert list(destination.parent.iterdir()) == [destination]


def test_extraction_cli_runs_without_package_imports(tmp_path):
    archive = tmp_path / "FastQC results.zip"
    data = b"##FastQC\t0.11.9\n"
    with ZipFile(archive, "w") as handle:
        handle.writestr("sample_fastqc/fastqc_data.txt", data)
        handle.writestr("sample_fastqc/fastqc_report.html", "<html></html>")

    destination = tmp_path / "results with spaces/A_1_fastqc_data.txt"
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)

    proc = subprocess.run(
        [
            sys.executable,
            "-I",
            str(paths.scripts_dir() / "fastqc.py"),
            str(archive),
            str(destination),
        ],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert destination.read_bytes() == data
    assert list(destination.parent.iterdir()) == [destination]


def test_extraction_cli_failure_preserves_existing_artifact(tmp_path):
    archive = tmp_path / "invalid.zip"
    with ZipFile(archive, "w") as handle:
        handle.writestr("sample_fastqc/summary.txt", "not the artifact")

    destination = tmp_path / "A_1_fastqc_data.txt"
    original = b"existing artifact\n"
    destination.write_bytes(original)

    proc = subprocess.run(
        [
            sys.executable,
            "-I",
            str(paths.scripts_dir() / "fastqc.py"),
            str(archive),
            str(destination),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert proc.returncode != 0
    assert "expected exactly one fastqc_data.txt member" in proc.stderr
    assert destination.read_bytes() == original
    assert not list(tmp_path.glob(".*.tmp"))
