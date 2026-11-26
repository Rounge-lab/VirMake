"""Sample-table construction. Construction through stdlib removed coupling between cli
and snakemake"""

import csv
import re
from pathlib import Path

import pandas as pd

from virmake.sample_sheet import SAMPLE_COLUMNS, Sample

_READS_RE = re.compile(r"(?P<sample>[^/]+)_(?P<frac>R?[12])\.fastq(?P<gz>\.gz)?$")
_CONTIGS_RE = re.compile(r"(?P<sample>[^/]+)\.fasta$")


def _match_dir(path: str | Path, pattern: re.Pattern) -> list[tuple[Path, dict]]:
    """Replicate ``glob_wildcards(str(path) + "/{...}")`` over one level."""
    directory = Path(path)
    matches = []
    if not directory.is_dir():
        return matches
    for entry in directory.iterdir():
        # skip hidden files
        if entry.name.startswith("."):
            continue
        # glob matches files only, never directories (os.walk filenames)
        if not entry.is_file():
            continue
        m = pattern.fullmatch(entry.name)
        if m:
            matches.append((entry, m.groupdict()))
    return matches


def get_reads_samples(path: str | Path, qc: bool = False) -> pd.DataFrame:
    """Discover paired reads; reject incomplete or ambiguous pairs."""
    samples: dict[str, dict[str, Path]] = {}
    conventions: dict[str, set[str]] = {}

    for entry, groups in sorted(_match_dir(path, _READS_RE)):
        sample_id = groups["sample"]
        suffix = groups["frac"]
        mate = suffix[-1]
        convention = "R" if suffix.startswith("R") else "numeric"

        mates = samples.setdefault(sample_id, {})
        conventions.setdefault(sample_id, set()).add(convention)

        if mate in mates:
            raise ValueError(
                f"Duplicate mate {mate} for sample {sample_id!r}: "
                f"{mates[mate]} and {entry}"
            )
        mates[mate] = entry.resolve()

    rows = []
    for sample_id, mates in sorted(samples.items()):
        if len(conventions[sample_id]) != 1:
            raise ValueError(
                f"Ambiguous paired-read naming for sample {sample_id!r}: "
                "do not mix _1/_2 with _R1/_R2."
            )

        missing = {"1", "2"} - mates.keys()
        if missing:
            raise ValueError(
                f"Missing mate {', '.join(sorted(missing))} "
                f"for sample {sample_id!r} in {path}."
            )

        rows.append(
            {
                "sample_id": sample_id,
                "r1": str(mates["1"]),
                "r2": str(mates["2"]),
            }
        )

    reads_table = pd.DataFrame(rows, columns=["sample_id", "r1", "r2"])
    if qc:
        reads_table = reads_table.rename(columns={"r1": "qc_r1", "r2": "qc_r2"})
    return reads_table


def get_contigs_samples(path: str | Path) -> pd.DataFrame:
    rows = [
        {"sample_id": groups["sample"], "contigs": str(entry)}
        for entry, groups in _match_dir(path, _CONTIGS_RE)
    ]
    contigs_table = pd.DataFrame(rows, columns=["sample_id", "contigs"])
    contigs_table = contigs_table.drop_duplicates(subset=["sample_id"])
    return contigs_table.sort_values(by=["sample_id"])


def build_samples_table(
    reads: str = "", qc_reads: str = "", contigs: str = ""
) -> pd.DataFrame:
    """Merge reads, QC reads, and contigs tables on sample_id (outer joins)."""
    reads_table = (
        get_reads_samples(reads)
        if reads
        else pd.DataFrame(columns=["sample_id", "r1", "r2"])
    )
    qc_reads_table = (
        get_reads_samples(qc_reads, qc=True)
        if qc_reads
        else pd.DataFrame(columns=["sample_id", "qc_r1", "qc_r2"])
    )
    contigs_table = (
        get_contigs_samples(contigs)
        if contigs
        else pd.DataFrame(columns=["sample_id", "contigs"])
    )
    samples_table = pd.merge(reads_table, qc_reads_table, how="outer", on=["sample_id"])
    return pd.merge(samples_table, contigs_table, how="outer", on=["sample_id"])


def write_samples_table(
    samples_path: str | Path, reads: str = "", qc_reads: str = "", contigs: str = ""
) -> pd.DataFrame:
    """Build the samples table and write it as TSV; returns the table."""
    samples_table = build_samples_table(reads, qc_reads, contigs)
    samples_table.to_csv(samples_path, sep="\t", index=False)
    return samples_table


def discover_samples(
    reads: str = "",
    qc_reads: str = "",
    contigs: str = "",
) -> tuple[Sample, ...]:
    """Build immutable declarations from paired-read and contig discovery.

    Raw and supplied-QC reads cannot both declare the same sample.
    Supplied contigs are independent and may accompany either read stage.
    """
    declarations: dict[str, dict] = {}

    for directory, stage in ((reads, "raw"), (qc_reads, "qc")):
        if not directory:
            continue
        for row in get_reads_samples(directory).to_dict(orient="records"):
            sample_id = row["sample_id"]
            if sample_id in declarations:
                raise ValueError(
                    f"sample {sample_id!r} is declared in both raw and QC reads"
                )
            declarations[sample_id] = {
                "sample_id": sample_id,
                "read_type": "short",
                "read_layout": "paired",
                "read_stage": stage,
                "r1": Path(row["r1"]),
                "r2": Path(row["r2"]),
            }

    if contigs:
        for row in get_contigs_samples(contigs).to_dict(orient="records"):
            sample_id = row["sample_id"]
            declaration = declarations.setdefault(sample_id, {"sample_id": sample_id})
            declaration["contigs"] = Path(row["contigs"]).resolve()

    return tuple(
        Sample.model_validate(declarations[sample_id])
        for sample_id in sorted(declarations)
    )


def write_discovered_sample_sheet(
    samples_path: str | Path,
    reads: str = "",
    qc_reads: str = "",
    contigs: str = "",
) -> tuple[Sample, ...]:
    """Write new-schema input declarations, returning validated records.

    Complete discovery and declaration validation before opening the sheet.
    This is a preparation API, not a workflow-generated-path updater.
    """
    samples = discover_samples(reads, qc_reads, contigs)

    with Path(samples_path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=SAMPLE_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for sample in samples:
            writer.writerow(sample.model_dump(mode="json"))

    return samples
