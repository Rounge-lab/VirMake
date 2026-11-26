"""Build workflow-coda QC artifact families and workflow dependencies."""

from collections.abc import Sequence

from virmake.sample_sheet import Sample


def _qc_groups(
    samples: Sequence[Sample],
) -> tuple[list[str], list[str], list[str]]:
    """Return raw, read-bearing, and paired short-read sample IDs."""
    raw = []
    reads = []
    paired = []
    seen: set[str] = set()

    for sample in samples:
        if sample.sample_id in seen:
            raise ValueError(f"duplicate sample_id: {sample.sample_id!r}")
        seen.add(sample.sample_id)

        if not sample.has_reads:
            continue
        if sample.read_type != "short":
            raise ValueError(
                f"sample {sample.sample_id!r}: long-read QC is not supported"
            )

        reads.append(sample.sample_id)

        if sample.requires_trimming:
            raw.append(sample.sample_id)
        if sample.is_paired:
            paired.append(sample.sample_id)

    if not reads:
        raise ValueError("QC requires at least one read-bearing sample")

    return raw, reads, paired


def build_qc_spec(
    samples: Sequence[Sample],
    producer_version: str,
) -> dict:
    """Declare QC families with file paths relative to the spec directory."""
    if not producer_version.strip():
        raise ValueError("producer_version must be non-blank")

    raw, reads, paired = _qc_groups(samples)
    artifacts = []

    if raw:
        artifacts.append(
            {
                "for_each": {
                    "dimension": "sample",
                    "values": raw,
                },
                "sanitize": True,
                "file": "qc/fastp/{sample}.json",
                "id": "fastp-{sample}",
                "type": "virmake.qc.fastp-{sample}",
                "format": "bio.fastp",
                "schema_version": "1.0",
            }
        )

    for mate, members in (("1", reads), ("2", paired)):
        if not members:
            continue

        artifacts.append(
            {
                "for_each": {
                    "dimension": "sample",
                    "values": members,
                },
                "sanitize": True,
                "file": (
                    "qc/fastqc/{sample}_"
                    f"{mate}_fastqc_data.txt"
                ),
                "id": f"fastqc-{{sample}}-r{mate}",
                "type": f"virmake.qc.fastqc-{{sample}}-r{mate}",
                "format": "bio.fastqc",
                "schema_version": "1.0",
            }
        )

    return {
        "producer": {
            "name": "virmake",
            "version": producer_version,
            "module_id": "virmake.qc",
        },
        "artifacts": artifacts,
    }


def qc_artifact_paths(samples: Sequence[Sample]) -> list[str]:
    """Return concrete relative paths required by the Snakemake DAG."""
    raw, reads, paired = _qc_groups(samples)

    return [
        *(f"qc/fastp/{sample}.json" for sample in raw),
        *(
            f"qc/fastqc/{sample}_1_fastqc_data.txt"
            for sample in reads
        ),
        *(
            f"qc/fastqc/{sample}_2_fastqc_data.txt"
            for sample in paired
        ),
    ]
