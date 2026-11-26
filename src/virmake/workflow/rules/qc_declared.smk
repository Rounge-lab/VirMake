"""Declaration-based QC producers and strict artifact verification."""

import re
from pathlib import Path

import yaml

from virmake import __version__
from scripts.qc_spec import build_qc_spec
from scripts.workflow_utils import load_samples, resolved_qc_reads

samples, SAMPLE = load_samples(config["path"]["samples"])

OUTPUT = Path(config["path"]["output"])
READS = str(OUTPUT / "reads" / "reads_qc")
FASTP = str(OUTPUT / "qc" / "fastp")
FASTQC = str(OUTPUT / "qc" / "fastqc")

RAW_SINGLE = [
    sample_id
    for sample_id, sample in samples.items()
    if sample.requires_trimming
    and sample.read_type == "short"
    and not sample.is_paired
]
RAW_PAIRED = [
    sample_id
    for sample_id, sample in samples.items()
    if sample.requires_trimming
    and sample.read_type == "short"
    and sample.is_paired
]


def sample_pattern(ids):
    # A never-matching expression keeps rules with empty selections inactive.
    return "(?:" + "|".join(re.escape(value) for value in ids) + ")" if ids else "(?!)"


def generated_reads(path):
    return path if config["keep_qc_reads"] else temp(path)


def qc_spec():
    # Called only for selected QC jobs: unsupported reads do not block
    # unrelated targets simply because this module is imported.
    return build_qc_spec(tuple(samples.values()), __version__)


def qc_artifacts(wildcards):
    return [
        str(OUTPUT / artifact["file"])
        for artifact in qc_spec()["artifacts"]
    ]


def fastqc_read(wildcards):
    paths = resolved_qc_reads(wildcards.sample, samples, OUTPUT)
    mate = int(wildcards.mate)
    if mate > len(paths):
        raise ValueError(
            f"sample {wildcards.sample!r} has no mate {wildcards.mate}"
        )
    return paths[mate - 1]


def fastqc_staged_read(wildcards):
    source = Path(fastqc_read(wildcards))
    suffix = ".fastq.gz" if source.name.endswith(".gz") else ".fastq"
    return str(
        OUTPUT / "qc" / "fastqc_tmp"
        / f"{wildcards.sample}_{wildcards.mate}{suffix}"
    )


rule QC:
    input:
        str(OUTPUT / "qc_manifest.json"),


rule fastp_se:
    input:
        R1=lambda w: str(samples[w.sample].r1),
    output:
        R1=generated_reads(READS + "/{sample}_1.fastq.gz"),
        json=FASTP + "/{sample}.json",
    wildcard_constraints:
        sample=sample_pattern(RAW_SINGLE),
    log:
        config["path"]["log"] + "/fastp/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/fastp/{sample}.txt",
    conda:
        config["path"]["envs"] + "/fastp.yaml"
    threads: 4
    resources:
        mem_mb=config["memory"]["tiny"],
        runtime=config["time"]["small"],
    shell:
        """
        fastp \
            -i {input.R1:q} \
            -o {output.R1:q} \
            -j {output.json:q} \
            --disable_html \
            -w {threads} \
            > {log:q} 2>&1
        """


rule fastp_pe:
    input:
        R1=lambda w: str(samples[w.sample].r1),
        R2=lambda w: str(samples[w.sample].r2),
    output:
        R1=generated_reads(READS + "/{sample}_1.fastq.gz"),
        R2=generated_reads(READS + "/{sample}_2.fastq.gz"),
        json=FASTP + "/{sample}.json",
    wildcard_constraints:
        sample=sample_pattern(RAW_PAIRED),
    log:
        config["path"]["log"] + "/fastp/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/fastp/{sample}.txt",
    conda:
        config["path"]["envs"] + "/fastp.yaml"
    threads: 4
    resources:
        mem_mb=config["memory"]["tiny"],
        runtime=config["time"]["small"],
    shell:
        """
        fastp \
            -i {input.R1:q} \
            -I {input.R2:q} \
            -o {output.R1:q} \
            -O {output.R2:q} \
            -j {output.json:q} \
            --disable_html \
            -w {threads} \
            > {log:q} 2>&1
        """


rule stage_fastqc_read:
    input:
        reads=fastqc_read,
    output:
        reads=temp(
            str(OUTPUT / "qc" / "fastqc_tmp"
                / "{sample}_{mate}.{extension}")
        ),
    wildcard_constraints:
        sample=sample_pattern([
            sample_id
            for sample_id, sample in samples.items()
            if sample.has_reads and sample.read_type == "short"
        ]),
        mate="[12]",
        extension="fastq(?:\\.gz)?",
    shell:
        """
        ln -s -- {input.reads:q} {output.reads:q}
        """


rule fastqc:
    input:
        reads=fastqc_staged_read,
    output:
        archive=temp(
            str(OUTPUT / "qc" / "fastqc_tmp"
                / "{sample}_{mate}_fastqc.zip")
        ),
        html=temp(
            str(OUTPUT / "qc" / "fastqc_tmp"
                / "{sample}_{mate}_fastqc.html")
        ),
    wildcard_constraints:
        sample=sample_pattern([
            sample_id
            for sample_id, sample in samples.items()
            if sample.has_reads and sample.read_type == "short"
        ]),
        mate="[12]",
    params:
        outdir=str(OUTPUT / "qc" / "fastqc_tmp"),
    log:
        config["path"]["log"] + "/fastqc/{sample}_{mate}.log",
    benchmark:
        config["path"]["benchmark"] + "/fastqc/{sample}_{mate}.txt",
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["normal"],
    shell:
        """
        mkdir -p {params.outdir:q}
        fastqc {input.reads:q} \
            --outdir {params.outdir:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """


rule extract_fastqc:
    input:
        archive=rules.fastqc.output.archive,
    output:
        data=str(
            OUTPUT / "qc" / "fastqc"
            / "{sample}_{mate}_fastqc_data.txt"
        ),
    params:
        extract=config["path"]["scripts"] + "/fastqc.py",
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    threads: 1
    resources:
        mem_mb=config["memory"]["tiny"],
        runtime=config["time"]["tiny"],
    shell:
        """
        python {params.extract:q} {input.archive:q} {output.data:q}
        """


rule write_qc_spec:
    input:
        sheet=config["path"]["samples"],
    output:
        spec=str(OUTPUT / "qc_spec.yaml"),
    params:
        specification=lambda w: qc_spec(),
    run:
        destination = Path(output.spec)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            yaml.safe_dump(params.specification, sort_keys=False),
            encoding="utf-8",
        )


rule verify_qc:
    input:
        spec=rules.write_qc_spec.output.spec,
        artifacts=qc_artifacts,
    output:
        manifest=str(OUTPUT / "qc_manifest.json"),
    log:
        config["path"]["log"] + "/qc_verification.log",
    conda:
        config["path"]["envs"] + "/coda.yaml"
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        workflow-coda verify {output.manifest:q} \
            --spec {input.spec:q} \
            --producer-name virmake \
            --producer-version {__version__:q} \
            --on-failure fail \
            > {log:q} 2>&1
        """
