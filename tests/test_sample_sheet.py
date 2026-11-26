"""Validation contract for the immutable sample-sheet schema."""

import csv
from pathlib import Path

import pytest
from pydantic import ValidationError

from virmake.sample_sheet import SAMPLE_COLUMNS, Sample, load_sample_sheet


def _write_sheet(path: Path, rows: list[dict]) -> Path:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _reads_row(**overrides) -> dict:
    return {
        "sample_id": "A",
        "read_type": "short",
        "read_layout": "single",
        "read_stage": "raw",
        "r1": "reads.fastq",
        **overrides,
    }


def test_relative_paths_resolve_against_sheet_not_cwd(tmp_path, monkeypatch):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    reads = inputs / "explicit name.fastq"
    reads.write_text("")
    sheet = _write_sheet(
        inputs / "samples.tsv",
        [_reads_row(r1=reads.name)],
    )
    monkeypatch.chdir(tmp_path)

    (sample,) = load_sample_sheet(sheet)

    assert sample.r1 == reads.resolve()
    assert sample.has_reads
    assert not sample.is_paired
    assert sample.requires_trimming
    assert not sample.has_supplied_contigs


def test_mixed_layouts_and_compression(tmp_path):
    for name in ("single.fastq", "mate-one.fastq", "mate-two.fastq.gz"):
        (tmp_path / name).write_text("")
    sheet = _write_sheet(
        tmp_path / "samples.tsv",
        [
            _reads_row(r1="single.fastq"),
            _reads_row(
                sample_id="B",
                read_layout="paired",
                read_stage="qc",
                r1="mate-one.fastq",
                r2="mate-two.fastq.gz",
            ),
        ],
    )

    single, paired = load_sample_sheet(sheet)

    assert not single.is_paired
    assert paired.is_paired
    assert not paired.requires_trimming


def test_contigs_only_and_reads_with_contigs(tmp_path):
    (tmp_path / "reads.fastq").write_text("")
    contigs = tmp_path / "assembly.fasta"
    contigs.write_text("")
    sheet = _write_sheet(
        tmp_path / "samples.tsv",
        [
            {"sample_id": "contigs-only", "contigs": contigs.name},
            _reads_row(contigs=contigs.name),
        ],
    )

    contigs_only, combined = load_sample_sheet(sheet)

    assert not contigs_only.has_reads
    assert contigs_only.read_type is None
    assert not contigs_only.requires_trimming
    assert contigs_only.has_supplied_contigs
    assert combined.has_reads
    assert combined.has_supplied_contigs


def test_long_reads_are_recognized(tmp_path):
    (tmp_path / "reads.fastq").write_text("")
    sheet = _write_sheet(
        tmp_path / "samples.tsv",
        [_reads_row(read_type="long")],
    )

    (sample,) = load_sample_sheet(sheet)

    assert sample.read_type == "long"
    assert sample.requires_trimming


@pytest.mark.parametrize(
    "changes",
    [
        {"sample_id": "../A"},
        {"sample_id": "A B"},
        {"sample_id": "_A"},
        {"read_type": "unknown"},
        {"read_stage": "processed"},
        {"read_layout": "unknown"},
        {"read_type": None},
        {"read_layout": None},
        {"read_stage": None},
        {"read_layout": "paired"},
        {"r2": Path("mate.fastq")},
        {"r1": None, "r2": Path("mate.fastq")},
        {"r1": None},
        {"r1": None, "contigs": Path("assembly.fasta")},
    ],
)
def test_invalid_declarations(changes):
    values = _reads_row(r1=Path("reads.fastq"))
    values.update(changes)

    with pytest.raises(ValidationError):
        Sample.model_validate(values)


def test_duplicate_ids_rejected(tmp_path):
    (tmp_path / "reads.fastq").write_text("")
    sheet = _write_sheet(
        tmp_path / "samples.tsv",
        [_reads_row(), _reads_row()],
    )

    with pytest.raises(ValueError, match="duplicate sample_id"):
        load_sample_sheet(sheet)


@pytest.mark.parametrize("directory", [False, True])
def test_missing_file_or_directory_rejected(tmp_path, directory):
    if directory:
        (tmp_path / "reads.fastq").mkdir()
    sheet = _write_sheet(tmp_path / "samples.tsv", [_reads_row()])

    with pytest.raises(ValueError, match="not an existing file"):
        load_sample_sheet(sheet)


@pytest.mark.parametrize(
    "header",
    [
        [],
        ["sample_id", "r1"],
        [*SAMPLE_COLUMNS, "extra"],
        [*SAMPLE_COLUMNS[:-1], "r1"],
    ],
)
def test_invalid_headers_rejected(tmp_path, header):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text("\t".join(header) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="expected exactly these columns"):
        load_sample_sheet(sheet)


@pytest.mark.parametrize("fields", [["A"], ["A"] + [""] * 7])
def test_malformed_rows_rejected(tmp_path, fields):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(
        "\t".join(SAMPLE_COLUMNS) + "\n" + "\t".join(fields) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="wrong number of TSV fields"):
        load_sample_sheet(sheet)


def test_header_only_sheet_allowed(tmp_path):
    sheet = _write_sheet(tmp_path / "samples.tsv", [])

    assert load_sample_sheet(sheet) == ()


def test_records_are_frozen():
    sample = Sample(sample_id="A", contigs=Path("assembly.fasta"))

    with pytest.raises(ValidationError, match="frozen"):
        sample.sample_id = "B"


def test_loading_does_not_modify_supplied_files(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_bytes(b"supplied reads\n")
    sheet = _write_sheet(tmp_path / "samples.tsv", [_reads_row()])
    original_sheet = sheet.read_bytes()
    original_reads = reads.read_bytes()

    load_sample_sheet(sheet)

    assert sheet.read_bytes() == original_sheet
    assert reads.read_bytes() == original_reads
