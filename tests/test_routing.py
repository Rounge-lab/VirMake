"""Pure routing contracts, independent of Snakemake."""

from pathlib import Path

import pytest

from virmake.routing import assembly_input, qc_reads
from virmake.sample_sheet import Sample


def _sample(**changes) -> Sample:
    values = {
        "sample_id": "A",
        "read_type": "short",
        "read_layout": "single",
        "read_stage": "raw",
        "r1": Path("supplied reads.fastq"),
    }
    values.update(changes)
    return Sample.model_validate(values)


@pytest.mark.parametrize("paired", [False, True])
def test_raw_reads_route_to_generated_outputs(tmp_path, paired):
    sample = _sample(
        read_layout="paired" if paired else "single",
        r2=Path("mate.fastq.gz") if paired else None,
    )
    results = tmp_path / "results"

    paths = qc_reads(sample, results)

    expected = (results / "reads/reads_qc/A_1.fastq.gz",)
    if paired:
        expected += (results / "reads/reads_qc/A_2.fastq.gz",)
    assert paths == expected
    assert sample.requires_trimming
    assert not results.exists()


@pytest.mark.parametrize("paired", [False, True])
def test_supplied_qc_reads_keep_original_paths(tmp_path, paired):
    r1 = tmp_path / "arbitrary first mate.fastq"
    r2 = tmp_path / "arbitrary second mate.fastq.gz"
    sample = _sample(
        read_stage="qc",
        read_layout="paired" if paired else "single",
        r1=r1,
        r2=r2 if paired else None,
    )

    paths = qc_reads(sample, tmp_path / "results")

    assert paths == ((r1, r2) if paired else (r1,))
    assert not sample.requires_trimming
    assert list(tmp_path.iterdir()) == []


def test_contigs_only_sample_has_no_qc_reads(tmp_path):
    sample = Sample(sample_id="A", contigs=Path("supplied.fasta"))

    with pytest.raises(ValueError, match="has no reads"):
        qc_reads(sample, tmp_path / "results")


@pytest.mark.parametrize("stage", ["raw", "qc"])
def test_long_read_qc_rejected(tmp_path, stage):
    sample = _sample(read_type="long", read_stage=stage)

    with pytest.raises(ValueError, match="long-read QC is not supported"):
        qc_reads(sample, tmp_path / "results")


@pytest.mark.parametrize("with_reads", [False, True])
def test_supplied_contigs_take_precedence(tmp_path, with_reads):
    contigs = tmp_path / "supplied assembly.fasta"
    sample = (
        _sample(contigs=contigs)
        if with_reads
        else Sample(sample_id="A", contigs=contigs)
    )

    assert assembly_input(sample, "metaSpades", tmp_path / "output") == contigs
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("assembler", ["metaSpades", "megahit"])
def test_generated_assembly_path(tmp_path, assembler):
    output = tmp_path / "results/output"

    assert assembly_input(_sample(), assembler, output) == (
        output / assembler / "A" / "contigs.fasta"
    )
    assert not output.exists()


def test_unknown_generated_assembler_rejected(tmp_path):
    with pytest.raises(ValueError, match="unsupported assembler"):
        assembly_input(_sample(), "unknown", tmp_path)


def test_routing_does_not_modify_supplied_files(tmp_path):
    reads = tmp_path / "reads.fastq"
    contigs = tmp_path / "contigs.fasta"
    reads.write_bytes(b"original reads\n")
    contigs.write_bytes(b"original contigs\n")
    sample = _sample(read_stage="qc", r1=reads, contigs=contigs)

    qc_reads(sample, tmp_path / "results")
    assembly_input(sample, "metaSpades", tmp_path / "output")

    assert reads.read_bytes() == b"original reads\n"
    assert contigs.read_bytes() == b"original contigs\n"
    assert {path.name for path in tmp_path.iterdir()} == {
        "reads.fastq",
        "contigs.fasta",
    }
