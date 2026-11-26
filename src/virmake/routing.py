"""Pure input routing; no filesystem mutation or availability checks."""

from pathlib import Path

from virmake.sample_sheet import Sample


def qc_reads(sample: Sample, results_dir: str | Path) -> tuple[Path, ...]:
    """Return supplied QC reads or declared trimming-output paths.

    Generated paths are a naming contract, not evidence that a producer exists
    or that trimming has run. Stage preflight must check those prerequisites.
    """
    if not sample.has_reads:
        raise ValueError(f"sample {sample.sample_id!r} has no reads")
    if sample.read_type != "short":
        raise ValueError(f"sample {sample.sample_id!r}: long-read QC is not supported")

    if not sample.requires_trimming:
        # Sample validation guarantees r1 and, for paired reads, r2.
        assert sample.r1 is not None
        if sample.is_paired:
            assert sample.r2 is not None
            return sample.r1, sample.r2
        return (sample.r1,)

    directory = Path(results_dir) / "reads" / "reads_qc"
    mates = ("1", "2") if sample.is_paired else ("1",)
    return tuple(directory / f"{sample.sample_id}_{mate}.fastq.gz" for mate in mates)


def assembly_input(
    sample: Sample,
    assembler: str,
    output_dir: str | Path,
) -> Path:
    """Return supplied contigs preferentially, or a generated assembly path."""
    if sample.has_supplied_contigs:
        assert sample.contigs is not None
        return sample.contigs

    if not sample.has_reads:
        raise ValueError(f"sample {sample.sample_id!r} has neither contigs nor reads")
    if assembler not in {"metaSpades", "megahit"}:
        raise ValueError(f"unsupported assembler: {assembler!r}")

    return Path(output_dir) / assembler / sample.sample_id / "contigs.fasta"
