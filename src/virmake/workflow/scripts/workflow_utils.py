import pandas as pd
from pathlib import Path

from virmake.routing import assembly_input, qc_reads
from virmake.sample_sheet import Sample, load_sample_sheet


def load_samples(path: str | Path) -> tuple[dict[str, Sample], list[str]]:
    """Load validated declarations indexed by ID, preserving sheet order."""
    records = load_sample_sheet(path)
    return (
        {sample.sample_id: sample for sample in records},
        [sample.sample_id for sample in records],
    )


def read_sample_ids(samples: dict[str, Sample]) -> list[str]:
    """Select read-bearing samples without silently excluding long reads.

    Selected-stage input routing must reject unsupported read types.
    """
    return [
        sample_id
        for sample_id, sample in samples.items()
        if sample.has_reads
    ]


def resolved_qc_reads(
    sample_id: str,
    samples: dict[str, Sample],
    results_dir: str | Path,
) -> list[str]:
    """Resolve all declared QC mates for a selected sample."""
    return [str(path) for path in qc_reads(samples[sample_id], results_dir)]


def resolved_assembly(
    sample_id: str,
    samples: dict[str, Sample],
    assembler: str,
    output_dir: str | Path,
) -> str:
    """Resolve supplied or generated contigs, retaining a string interface."""
    return str(assembly_input(samples[sample_id], assembler, output_dir))


def paired_assembly_reads(
    sample_id: str,
    samples: dict[str, Sample],
    results_dir: str | Path,
) -> list[str]:
    """Resolve inputs for the currently paired-only assembly producers.

    Called from a selected job's input function, not during module import.
    Supplied contigs bypass these producers through resolved_assembly.
    """
    paths = resolved_qc_reads(sample_id, samples, results_dir)
    if len(paths) != 2:
        raise ValueError(
            f"sample {sample_id!r}: generated assembly currently requires "
            "paired-end short reads"
        )
    return paths


def get_samples(path):
    sample_table = pd.read_csv(path, delimiter='\t')
    samples = sample_table["sample_id"].to_list()
    return sample_table, samples

def get_qc_reads_loc(wildcards,
                     sample_table,
                     standard_output,
                     r,
                     keep_qc_reads):
    qc_r1_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'qc_r1'].iloc[0]
    qc_r2_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'qc_r2'].iloc[0]
    if pd.isna(qc_r1_loc) or pd.isna(qc_r2_loc) or qc_r1_loc == "" or qc_r2_loc == "":
        if keep_qc_reads:
            return standard_output+"/fastp_pe/qc_reads/"+wildcards.sample+f"_{r}.fastq.gz"
        else:
            return standard_output+"/fastp_pe/"+wildcards.sample+f"_{r}.fastq.gz"
    else:
        if r == "1":
            return qc_r1_loc
        elif r == "2":
            return qc_r2_loc

def get_assembly_loc(wildcards,
                     sample_table,
                     assembler,
                     standard_output):
    assembly_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'contigs'].iloc[0]
    if pd.isna(assembly_loc) or assembly_loc == "":
        return standard_output + "/" + assembler + "/"+wildcards.sample+"/contigs.fasta"
    else:
        return assembly_loc

# GLOBAL FUNCTIONS #
def get_min_quality(threshold):
    """
    Gets the threshold for checkv quality control.
    """
    if threshold.lower() == "complete":
        return ["Complete"]
    elif threshold.lower() == "high":
        return ["Complete", "High-quality"]
    elif threshold.lower() == "medium":
        return ["Complete", "High-quality", "Medium-quality"]
    elif threshold.lower() == "low":
        return ["Complete", "High-quality", "Medium-quality", "Low-quality"]
    elif threshold.lower() == "not-determined":
        return ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]
    else:
        print("Quality threshold not defined properly. Using medium.")
        return ["Complete", "High-quality", "Medium-quality"]

def vibrant_virome(is_virome):
    """
    Gathers if the input is virome.
    """
    if is_virome.lower == "yes":
        return "-virome"
    else:
        return ""
