from scripts.workflow_utils import get_samples

sample_table, SAMPLE = get_samples(config["path"]["samples"])

# ASSEMBLY #

rule ASSEMBLY:
    input:
        expand(
            config["path"]["output"] + "/metaSpades/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/metaSpades/{sample}/scaffolds.fasta",
            sample=SAMPLE,
        ),
        config["path"]["output"] + "/metaQUAST/report.html",
        config["path"]["output"]
        + "/metaQUAST/combined_reference/transposed_report.tsv",
    output:
        config["path"]["temp"] + "/finished_ASSEMBLY",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    message:
        "[ASSEMBLY] Assembly of all samples finished"
    shell:
        """
        touch {output}
        """


rule metaSpades:
    """
    Assembles all sequences with metaSpades
    """
    input:
        R1=config["path"]["output"] + "/fastp_pe/{sample}_1.fastq",
        R2=config["path"]["output"] + "/fastp_pe/{sample}_2.fastq",
    output:
        dir=directory(
            config["path"]["output"] + "/metaSpades/{sample}/",
        ),
        contigs=config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",
        scaffolds=config["path"]["output"] + "/metaSpades/{sample}/scaffolds.fasta",
    params:
        temp_dir=config["path"]["temp"] + "/metaSpades/",
    message:
        "[metaSpades] Performing assembly of paired end reads..."
    conda:
        config["path"]["envs"] + "/metaSpades.yaml"
    log:
        config["path"]["log"] + "/metaSpades/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/metaSpades/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["metaspades"],
        runtime=config["time"]["metaspades"],
    shell:
        """
        metaspades.py -1 {input.R1} -2 {input.R2} -o {output.dir}\
        --tmp-dir {params.temp_dir} -t {threads} -m {resources.mem_mb} &> {log}
        """

rule metaQUAST:
    """
    performs quality control on all assembled contigs
    with the reference database RefSeq Viral
    """
    params:
        reference=config["path"]["database"]["RefSeq"] + "/viral.1.1.genomic.fna",
    input:
        expand(
            config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",
            sample=SAMPLE,
        ),
    output:
        dir=directory(config["path"]["output"] + "/metaQUAST/"),
        report=config["path"]["output"] + "/metaQUAST/report.html",
        transposed_report=config["path"]["output"]
        + "/metaQUAST/combined_reference/transposed_report.tsv",
    message:
        "[metaQUAST] Running metaQUAST quality control on assembled contigs..."
    conda:
        config["path"]["envs"] + "/metaQUAST.yaml"
    log:
        config["path"]["log"] + "/metaQUAST.log",
    benchmark:
        config["path"]["benchmark"] + "/metaQUAST.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["metaquast"],
        runtime=config["time"]["metaquast"],
    shell:
        """
        mkdir -p {output.dir}
        metaquast.py {input} -o {output.dir}\
        -r {params.reference} --threads {threads} --max-ref-number 0 &> {log}
        """