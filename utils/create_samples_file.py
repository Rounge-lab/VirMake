import pathlib
import click
import pandas as pd
from snakemake.io import glob_wildcards


def get_reads_samples(path):

    path = path.rstrip("/")

    (
        SAMPLES,
        FRAC,
    ) = glob_wildcards(path + "/{sample}_{frac}.fastq.gz")

    # remove duplicates
    SAMPLES = sorted(list(set(SAMPLES)))
    
    reads_table = pd.DataFrame({
        "sample_id": SAMPLES ,
        "r1": [ path+f"/{sample}_1.fastq.gz" for sample in SAMPLES ],
        "r2": [ path+f"/{sample}_2.fastq.gz" for sample in SAMPLES ]
    })

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
    "-c",
    "--contigs",
    default="",
    help="Location of assembled contigs.",
)
def write_samples_file(virmake_path, work_dir, reads, contigs):
    # if reads == "" and contigs == "":
    #     samples_table = pd.DataFrame(columns = ["sample_id", "r1", "r2", "contigs"])
    # else:
    if not reads == "":
        reads_table = get_reads_samples(reads)
    else:
        reads_table = pd.DataFrame(columns = ["sample_id", "r1", "r2"])
    if not contigs == "":
        contigs_table = get_contigs_samples(contigs)
    else:
        contigs_table = pd.DataFrame(columns = ["sample_id", "contigs"])
    
    virmake_path = pathlib.Path(virmake_path)

    samples_table = pd.merge(left=reads_table, right=contigs_table, how="outer", on=["sample_id"])
    samples_table.to_csv(virmake_path / work_dir / "samples.tsv", sep="\t", index=False)

if __name__ == "__main__":
    write_samples_file()
