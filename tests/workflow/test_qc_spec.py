"""QC family declarations, independent of external tools."""

from pathlib import Path

import pytest

from virmake.sample_sheet import Sample
from virmake.workflow.scripts.qc_spec import (
    build_qc_spec,
    qc_artifact_paths,
)


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
    samples = (_reads(paired=paired),)
    spec = build_qc_spec(samples, "1.2.3")

    assert spec["producer"] == {
        "name": "virmake",
        "version": "1.2.3",
    }

    expected_templates = [
        "qc/fastp/{sample}.json",
        "qc/fastqc/{sample}_1_fastqc_data.txt",
    ]
    expected_paths = [
        "qc/fastp/A.json",
        "qc/fastqc/A_1_fastqc_data.txt",
    ]
    if paired:
        expected_templates.append("qc/fastqc/{sample}_2_fastqc_data.txt")
        expected_paths.append("qc/fastqc/A_2_fastqc_data.txt")

    artifacts = spec["artifacts"]
    assert [artifact["file"] for artifact in artifacts] == expected_templates
    assert qc_artifact_paths(samples) == expected_paths
    assert [artifact["format"] for artifact in artifacts] == (
        ["bio.fastp"] + ["bio.fastqc"] * (2 if paired else 1)
    )
    assert all(artifact["schema_version"] == "1.0" for artifact in artifacts)
    assert all(artifact["sanitize"] is True for artifact in artifacts)
    assert all(
        artifact["for_each"]
        == {
            "dimension": "sample",
            "values": ["A"],
        }
        for artifact in artifacts
    )


def test_supplied_qc_reads_do_not_claim_fastp_report():
    samples = (_reads(stage="qc", paired=True),)
    spec = build_qc_spec(samples, "1.2.3")

    assert [artifact["file"] for artifact in spec["artifacts"]] == [
        "qc/fastqc/{sample}_1_fastqc_data.txt",
        "qc/fastqc/{sample}_2_fastqc_data.txt",
    ]
    assert qc_artifact_paths(samples) == [
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
    assert spec["artifacts"][0]["for_each"] == {
        "dimension": "sample",
        "values": ["B"],
    }
    assert qc_artifact_paths(samples) == [
        "qc/fastqc/B_1_fastqc_data.txt",
    ]


def test_family_membership_preserves_sheet_order():
    samples = (
        _reads("B", stage="qc", paired=True),
        _reads("A", stage="raw"),
        _reads("C", stage="raw", paired=True),
    )
    artifacts = build_qc_spec(samples, "1.2.3")["artifacts"]

    assert [artifact["for_each"]["values"] for artifact in artifacts] == [
        ["A", "C"],  # fastp: raw only
        ["B", "A", "C"],  # FastQC mate 1: all reads
        ["B", "C"],  # FastQC mate 2: paired only
    ]


def test_family_identity_templates():
    artifacts = build_qc_spec(
        (_reads(paired=True),),
        "1.2.3",
    )["artifacts"]

    assert [artifact["id"] for artifact in artifacts] == [
        "fastp-{sample}",
        "fastqc-{sample}-r1",
        "fastqc-{sample}-r2",
    ]
    assert [artifact["type"] for artifact in artifacts] == [
        "virmake.qc.fastp-{sample}",
        "virmake.qc.fastqc-{sample}-r1",
        "virmake.qc.fastqc-{sample}-r2",
    ]


def test_case_distinct_values_and_paths_are_preserved_for_coda():
    samples = (
        _reads("A", stage="qc"),
        _reads("a", stage="qc"),
    )
    artifacts = build_qc_spec(samples, "1.2.3")["artifacts"]

    assert len(artifacts) == 1
    assert artifacts[0]["for_each"]["values"] == ["A", "a"]
    assert artifacts[0]["sanitize"] is True
    assert qc_artifact_paths(samples) == [
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


def test_concrete_paths_reject_invalid_qc_selection():
    with pytest.raises(ValueError, match="at least one read-bearing sample"):
        qc_artifact_paths(())


@pytest.mark.parametrize("stage", ["raw", "qc"])
@pytest.mark.parametrize("paired", [False, True])
def test_all_qc_families_have_module_attribution(stage, paired):
    spec = build_qc_spec(
        (_reads(stage=stage, paired=paired),),
        "1.2.3",
    )

    assert spec["producer"] == {
        "name": "virmake",
        "version": "1.2.3",
    }
    assert all(artifact["module_id"] == "qc" for artifact in spec["artifacts"])
