import pathlib
import click
import pandas as pd
from snakemake.io import glob_wildcards


def get_reads_samples(path, qc=False):

    path = path.rstrip("/")

    ( SAMPLES_fastq, Ffq ) = glob_wildcards(path + "/{sample,[^/]+}_{frac,[^/]+}.fastq")
    ( SAMPLES_gz, Fgz ) = glob_wildcards(path + "/{sample,[^/]+}_{frac,[^/]+}.fastq.gz")

    if set(SAMPLES_fastq).intersection(set(SAMPLES_gz)):
        raise ValueError("Error: compressed and uncompressed read files from the same sample.")

    tmp_df = pd.DataFrame({
        "sample_id": SAMPLES_fastq + SAMPLES_gz,
        "p": [ path+f"/{sample}_{frac}.fastq" for sample, frac in zip(SAMPLES_fastq, Ffq)] +
            [ path+f"/{sample}_{frac}.fastq.gz" for sample, frac in zip(SAMPLES_gz, Fgz)],
        "frac": Ffq + Fgz,
    })

    reads_table = pd.merge(tmp_df.loc[tmp_df["frac"] == "1",["sample_id","p"]].rename(columns={"p": "r1"}),
                           tmp_df.loc[tmp_df["frac"] == "2",["sample_id","p"]].rename(columns={"p": "r2"}),
                           how="inner", on=["sample_id"])

    reads_table.sort_values(by = ["sample_id"],inplace=True)

    if qc:
        reads_table = reads_table.rename(columns={"r1": "qc_r1", "r2": "qc_r2"})

    return reads_table

def get_contigs_samples(path):

    path = path.rstrip("/")

    (
        SAMPLES,
    ) = glob_wildcards(path + "/{sample}.fasta")

    # remove duplicates
    SAMPLES = sorted(list(set(SAMPLES)))

    contigs_table = pd.DataFrame({
        "sample_id": SAMPLES ,
        "contigs": [ path+f"/{sample}.fasta" for sample in SAMPLES ],
    })

    return contigs_table

@click.command()
@click.argument("virmake_path")
@click.option(
    "-w",
    "--work-dir",
    default="results",
    help="Path to working directory.",
)
@click.option(
    "-r",
    "--reads",
    default="",
    help="Location of read files",
)
@click.option(
    "-q",
    "--qc-reads",
    default="",
    help="Location of QC read files",
)
@click.option(
    "-c",
    "--contigs",
    default="",
    help="Location of assembled contigs.",
)
def write_samples_file(virmake_path, work_dir, reads, qc_reads, contigs):

    if not reads == "":
        reads_table = get_reads_samples(reads)
    else:
        reads_table = pd.DataFrame(columns = ["sample_id", "r1", "r2"])

    if not qc_reads == "":
        qc_reads_table = get_reads_samples(qc_reads, qc=True)
    else:
        qc_reads_table = pd.DataFrame(columns = ["sample_id", "qc_r1", "qc_r2"])

    if not contigs == "":
        contigs_table = get_contigs_samples(contigs)
    else:
        contigs_table = pd.DataFrame(columns = ["sample_id", "contigs"])

    virmake_path = pathlib.Path(virmake_path)

    samples_table = pd.merge(left=reads_table, right=qc_reads_table, how="outer", on=["sample_id"])
    samples_table = pd.merge(left=samples_table, right=contigs_table, how="outer", on=["sample_id"])

    samples_table.to_csv(virmake_path / work_dir / "samples.tsv", sep="\t", index=False)

if __name__ == "__main__":
    write_samples_file()
