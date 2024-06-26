from scripts.workflow_utils import get_samples, get_qc_reads_loc

sample_table, SAMPLE = get_samples(config["path"]["samples"])
FRAC = ["1", "2"]

# QUALITY CONTROL #
rule QC:
    input:
        expand(
            config["path"]["output"] + "/fastp_pe/{sample}_1.fastq",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/fastp_pe/{sample}_2.fastq",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/fastp_pe/{sample}.html",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/fastp_pe/{sample}.json",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/fastqc/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
        config["path"]["output"] + "/multiqc/multiqc.html",
    output:
        config["path"]["temp"] + "/finished_QC",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    message:
        "[QC] Finished QC."
    shell:
        """
        touch {output}
        """

rule fastp_pe:
    """
    performs quality control/pre-processing of the raw reads
    """
    input:
        R1=config["path"]["input_reads"] + "/{sample}_1.fastq.gz",
        R2=config["path"]["input_reads"] + "/{sample}_2.fastq.gz",
    output:
        R1=config["path"]["output"] + "/fastp_pe/{sample}_1.fastq",
        R2=config["path"]["output"] + "/fastp_pe/{sample}_2.fastq",
        html=config["path"]["output"] + "/fastp_pe/{sample}.html",
        json=config["path"]["output"] + "/fastp_pe/{sample}.json",
    log:
        config["path"]["log"] + "/fastp_pe/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/fastp_pe/{sample}.txt"
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["normal"],
    conda:
        config["path"]["envs"] + "/fastp.yaml"
    message:
        "[fastp_pe] Executing FASTP quality control/pre-processing (trimming) on raw reads..."
    threads: 1
    shell:
        """
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}\
        -h {output.html} -j {output.json} &> {log}
        """

rule move_qc_for_fastqc:
    input:
        R1=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="1"),
        R2=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="2"),
    output:
        R1 = temp(config["path"]["output"] + "/qc/tmp/{sample}_1.fastq"),
        R2 = temp(config["path"]["output"] + "/qc/tmp/{sample}_2.fastq"),
    shell:
        """
            ln -s {input.R1} {output.R1}
            ln -s {input.R2} {output.R2}
        """

rule fastqc:
    """
    performs quality control of processed QC reads
    """
    input:
        expand(rules.move_qc_for_fastqc.output.R1, sample=SAMPLE),
        expand(rules.move_qc_for_fastqc.output.R2, sample=SAMPLE),
    output:
        report=expand(
            config["path"]["output"] + "/fastqc/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
        dir=directory(config["path"]["output"] + "/fastqc/"),
    threads: config["threads"]
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    message:
        "[fastqc] Executing FASTQC quality control on trimmed reads..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    log:
        config["path"]["log"] + "/fastqc.log",
    benchmark:
        config["path"]["benchmark"] + "/fastqc.txt"
    shell:
        """
        fastqc {input} -o {output.dir} -t {threads} &> {log}
        """

rule multiqc_read_qc:
    """
    Summarize QC after read processing
    """
    input:
        fastqc_output=expand(
            config["path"]["output"] + "/fastqc/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
    output:
        report=config["path"]["output"] + "/multiqc/multiqc.html",
    params:
        read_qc_path=config["path"]["output"] + "/fastqc/",
        out_dir=config["path"]["output"] + "/multiqc/",
        out_name="multiqc"
    threads: config["threads"]
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    message:
        "[multiqc_qc] Executing MULTIQC summary of read quality control..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    log:
        config["path"]["log"] + "/multiqc.log",
    benchmark:
        config["path"]["benchmark"] + "/multiqc.txt"
    shell:
        """
        multiqc --force -o {params.out_dir} -n {params.out_name} {params.read_qc_path} 2> {log}
        """