"""QC verification declarations, independent of external tools."""

from pathlib import Path

import pytest

from virmake.sample_sheet import Sample
from virmake.workflow.scripts.qc_spec import build_qc_spec


def _reads(sample_id="A", *, stage="raw", paired=False, read_type="short"):
    return Sample(
        sample_id=sample_id,
        read_type=read_type,
        read_layout="paired" if paired else "single",
        read_stage=stage,
        r1=Path(f"{sample_id}_1.fastq"),
        r2=Path(f"{sample_id}_2.fastq") if paired else None,
    )


@pytest.mark.parametrize("paired", [False, True])
def test_raw_reads_declare_fastp_and_each_fastqc_mate(paired):
    spec = build_qc_spec((_reads(paired=paired),), "1.2.3")

    assert spec["producer"] == {"name": "virmake", "version": "1.2.3"}
    artifacts = spec["artifacts"]
    expected_files = [
        "qc/fastp/A.json",
        "qc/fastqc/A_1_fastqc_data.txt",
    ]
    if paired:
        expected_files.append("qc/fastqc/A_2_fastqc_data.txt")

    assert [artifact["file"] for artifact in artifacts] == expected_files
    assert [artifact["format"] for artifact in artifacts] == (
        ["bio.fastp"] + ["bio.fastqc"] * (2 if paired else 1)
    )
    assert all(artifact["schema_version"] == "1.0" for artifact in artifacts)


def test_supplied_qc_reads_do_not_claim_fastp_report():
    spec = build_qc_spec((_reads(stage="qc", paired=True),), "1.2.3")

    assert [artifact["file"] for artifact in spec["artifacts"]] == [
        "qc/fastqc/A_1_fastqc_data.txt",
        "qc/fastqc/A_2_fastqc_data.txt",
    ]
    assert all(artifact["format"] == "bio.fastqc" for artifact in spec["artifacts"])


def test_contigs_only_samples_have_no_qc_artifacts():
    samples = (
        Sample(sample_id="contigs", contigs=Path("assembly.fasta")),
        _reads("B", stage="qc"),
    )

    spec = build_qc_spec(samples, "1.2.3")

    assert len(spec["artifacts"]) == 1
    assert spec["artifacts"][0]["file"] == "qc/fastqc/B_1_fastqc_data.txt"


def test_case_distinct_sample_names_have_distinct_addresses():
    spec = build_qc_spec(
        (_reads("A", stage="qc"), _reads("a", stage="qc")),
        "1.2.3",
    )

    artifacts = spec["artifacts"]
    assert len({artifact["id"] for artifact in artifacts}) == 2
    assert len({artifact["type"] for artifact in artifacts}) == 2
    assert [artifact["file"] for artifact in artifacts] == [
        "qc/fastqc/A_1_fastqc_data.txt",
        "qc/fastqc/a_1_fastqc_data.txt",
    ]


@pytest.mark.parametrize("stage", ["raw", "qc"])
def test_long_reads_are_rejected(stage):
    with pytest.raises(ValueError, match="long-read QC is not supported"):
        build_qc_spec((_reads(stage=stage, read_type="long"),), "1.2.3")


@pytest.mark.parametrize(
    "samples",
    [
        (),
        (Sample(sample_id="A", contigs=Path("assembly.fasta")),),
    ],
)
def test_no_read_bearing_samples_are_rejected(samples):
    with pytest.raises(ValueError, match="at least one read-bearing sample"):
        build_qc_spec(samples, "1.2.3")


def test_duplicate_ids_are_rejected():
    sample = _reads()

    with pytest.raises(ValueError, match="duplicate sample_id"):
        build_qc_spec((sample, sample), "1.2.3")


def test_blank_producer_version_is_rejected():
    with pytest.raises(ValueError, match="producer_version"):
        build_qc_spec((_reads(),), " ")
