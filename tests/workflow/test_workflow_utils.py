"""Record-based workflow input interfaces."""

from pathlib import Path

import pytest

from virmake.sample_sheet import Sample
from virmake.samples import write_sample_sheet
from virmake.workflow.scripts.workflow_utils import (
    load_samples,
    paired_assembly_reads,
    read_sample_ids,
    resolved_assembly,
    resolved_qc_reads,
)


def _reads(sample_id="A", *, paired=False, stage="raw", read_type="short"):
    return Sample(
        sample_id=sample_id,
        read_type=read_type,
        read_layout="paired" if paired else "single",
        read_stage=stage,
        r1=Path("supplied first.fastq"),
        r2=Path("supplied second.fastq") if paired else None,
    )


def test_loader_preserves_validated_records_and_order(tmp_path):
    contigs = tmp_path / "supplied.fasta"
    contigs.write_text("")
    records = (
        Sample(sample_id="B", contigs=contigs.resolve()),
        Sample(sample_id="A", contigs=contigs.resolve()),
    )
    sheet = tmp_path / "samples.tsv"
    write_sample_sheet(sheet, records)

    samples, ids = load_samples(sheet)

    assert ids == ["B", "A"]
    assert list(samples) == ids
    assert tuple(samples.values()) == records


def test_read_selection_excludes_only_contigs_only_samples():
    samples = {
        "contigs": Sample(sample_id="contigs", contigs=Path("assembly.fasta")),
        "short": _reads("short"),
        "long": _reads("long", read_type="long"),
    }

    assert read_sample_ids(samples) == ["short", "long"]


@pytest.mark.parametrize("paired", [False, True])
def test_generated_qc_paths_use_results_root(tmp_path, paired):
    samples = {"A": _reads(paired=paired)}
    results = tmp_path / "results"

    paths = resolved_qc_reads("A", samples, results)

    expected = [str(results / "reads/reads_qc/A_1.fastq.gz")]
    if paired:
        expected.append(str(results / "reads/reads_qc/A_2.fastq.gz"))
    assert paths == expected
    assert not results.exists()


def test_supplied_qc_paths_are_returned_without_copying(tmp_path):
    sample = _reads(paired=True, stage="qc")

    assert resolved_qc_reads("A", {"A": sample}, tmp_path) == [
        str(sample.r1),
        str(sample.r2),
    ]
    assert list(tmp_path.iterdir()) == []


def test_assembly_returns_string_and_prefers_supplied_contigs(tmp_path):
    contigs = tmp_path / "supplied assembly.fasta"
    sample = Sample(sample_id="A", contigs=contigs)

    path = resolved_assembly("A", {"A": sample}, "unknown", tmp_path)

    assert path == str(contigs)
    assert path + ".fai" == str(contigs) + ".fai"
    assert list(tmp_path.iterdir()) == []


def test_paired_assembly_resolves_two_inputs(tmp_path):
    samples = {"A": _reads(paired=True)}

    assert paired_assembly_reads("A", samples, tmp_path) == [
        str(tmp_path / "reads/reads_qc/A_1.fastq.gz"),
        str(tmp_path / "reads/reads_qc/A_2.fastq.gz"),
    ]


def test_single_end_assembly_is_rejected_when_inputs_are_requested(tmp_path):
    samples = {"A": _reads()}

    # Merely loading/selecting samples does not reject unrelated QC work.
    assert read_sample_ids(samples) == ["A"]

    with pytest.raises(ValueError, match="requires paired-end short reads"):
        paired_assembly_reads("A", samples, tmp_path)


@pytest.mark.parametrize("stage", ["raw", "qc"])
def test_selected_long_read_processing_is_rejected(tmp_path, stage):
    samples = {"A": _reads(stage=stage, read_type="long")}

    with pytest.raises(ValueError, match="long-read QC is not supported"):
        resolved_qc_reads("A", samples, tmp_path)
