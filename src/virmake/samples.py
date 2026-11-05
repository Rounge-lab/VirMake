"""Sample-table construction. Construction through stdlib removed coupling between cli
and snakemake"""

import re
from pathlib import Path

import pandas as pd

_READS_RE = re.compile(r"(?P<sample>[^/]+)_(?P<frac>[^/]+)\.fastq(?P<gz>\.gz)?$")
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
    rows = [
        {
            "sample_id": groups["sample"],
            "p": str(entry),
            "frac": groups["frac"],
            "gz": groups["gz"] is not None,
        }
        for entry, groups in _match_dir(path, _READS_RE)
    ]
    plain = {r["sample_id"] for r in rows if not r["gz"]}
    zipped = {r["sample_id"] for r in rows if r["gz"]}
    if plain & zipped:
        raise ValueError(
            "Error: compressed and uncompressed read files from the same sample."
        )
    tmp_df = pd.DataFrame(rows, columns=["sample_id", "p", "frac"])
    r1 = tmp_df.loc[tmp_df["frac"] == "1", ["sample_id", "p"]].rename(
        columns={"p": "r1"}
    )
    r2 = tmp_df.loc[tmp_df["frac"] == "2", ["sample_id", "p"]].rename(
        columns={"p": "r2"}
    )
    reads_table = pd.merge(r1, r2, how="inner", on=["sample_id"])
    reads_table = reads_table.sort_values(by=["sample_id"])
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
    contigs_table = contigs_table.sort_values(by=["sample_id"])
    return contigs_table


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
    return pd.merge(
        samples_table, contigs_table, how="outer", on=["sample_id"]
    )


def write_samples_table(
    samples_path: str | Path, reads: str = "", qc_reads: str = "", contigs: str = ""
) -> pd.DataFrame:
    """Build the samples table and write it as TSV; returns the table."""
    samples_table = build_samples_table(reads, qc_reads, contigs)
    samples_table.to_csv(samples_path, sep="\t", index=False)
    return samples_table
