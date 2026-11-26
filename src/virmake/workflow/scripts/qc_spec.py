"""Build explicit workflow-coda QC verification declarations."""

from collections.abc import Sequence

from virmake.sample_sheet import Sample


def build_qc_spec(
    samples: Sequence[Sample],
    producer_version: str,
) -> dict:
    """Return a spec whose file paths are relative to the results directory.

    Persistent artifact naming contract:
      qc/fastp/<sample>.json
      qc/fastqc/<sample>_<mate>_fastqc_data.txt

    Only raw short-read samples produce fastp reports. Every short-read
    sample produces one FastQC artifact per declared mate.

    Numeric addresses avoid collisions from case-folding or sanitizing sample
    identifiers. Sample identifiers remain unchanged in filesystem paths.
    """
    if not producer_version.strip():
        raise ValueError("producer_version must be non-blank")

    artifacts = []
    seen: set[str] = set()

    for position, sample in enumerate(samples, start=1):
        if sample.sample_id in seen:
            raise ValueError(f"duplicate sample_id: {sample.sample_id!r}")
        seen.add(sample.sample_id)

        if not sample.has_reads:
            continue
        if sample.read_type != "short":
            raise ValueError(
                f"sample {sample.sample_id!r}: long-read QC is not supported"
            )

        address = f"sample-{position}"
        if sample.requires_trimming:
            artifacts.append(
                {
                    "file": f"qc/fastp/{sample.sample_id}.json",
                    "id": f"fastp-{address}",
                    "type": f"virmake.qc.fastp.{address}",
                    "format": "bio.fastp",
                    "schema_version": "1.0",
                }
            )

        mates = ("1", "2") if sample.is_paired else ("1",)
        for mate in mates:
            artifacts.append(
                {
                    "file": (
                        f"qc/fastqc/{sample.sample_id}_{mate}_fastqc_data.txt"
                    ),
                    "id": f"fastqc-{address}-mate-{mate}",
                    "type": f"virmake.qc.fastqc.{address}.mate-{mate}",
                    "format": "bio.fastqc",
                    "schema_version": "1.0",
                }
            )

    if not artifacts:
        raise ValueError("QC requires at least one read-bearing sample")

    return {
        "producer": {
            "name": "virmake",
            "version": producer_version,
        },
        "artifacts": artifacts,
    }
