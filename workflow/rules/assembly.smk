from scripts.workflow_utils import get_samples, get_qc_reads_loc, get_assembly_loc

sample_table, SAMPLE = get_samples(config["path"]["samples"])

# ASSEMBLY #

rule ASSEMBLY:
    input:
        expand(config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",sample=SAMPLE),
        config["path"]["output"] + "/metaQUAST/combined_reference/transposed_report.tsv" if config["rule_inclusion"]["all"]["metaquast"] else [],
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
        R1=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="1"),
        R2=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="2"),
    output:
        dir=directory(
            config["path"]["output"] + "/metaSpades/{sample}/",
        ),
        contigs=config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",
        scaffolds=config["path"]["output"] + "/metaSpades/{sample}/scaffolds.fasta",
    params:
        temp_dir=config["path"]["temp"] + "/metaSpades/",
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


rule run_megahit:
    input:
        reads=lambda w: [ get_qc_reads_loc(w, sample_table, config["path"]["output"],r=rd) for rd in ["1", "2"]],
    output:
        contigs=config["path"]["output"] + "/megahit/{sample}/contigs.fasta",
    benchmark:
        config["path"]["benchmark"] + "/megahit/{sample}.txt"
    params:
        # all parameters are optional
        extra="--min-count 10 --k-list 21,29,39,59,79,99,119,141",
    log:
        config["path"]["log"] + "/megahit/{sample}.log",
    threads: 8
    resources:
        mem_mb=250000,
        runtime=1440
    wrapper:
        "v7.4.0/bio/megahit"


rule get_contigs_for_metaquast:
    input:
        contigs=lambda w: get_assembly_loc(w, sample_table, config["assembler"], config["path"]["output"]),
    output:
        contigs=temp(config["path"]["output"]+"/assembly/tmp/{sample}.fasta")
    shell:
        "cp {input.contigs} {output.contigs}"

rule metaQUAST:
    """
    Performs quality control on all assembled contigs
    with the reference database RefSeq Viral
    """
    input:
        contigs=expand(
            config["path"]["output"]+"/assembly/tmp/{sample}.fasta",
            sample=SAMPLE,
        ),
        reference=config["path"]["database"]["RefSeq"],
    output:
        dir=directory(config["path"]["output"] + "/metaQUAST/"),
        report=config["path"]["output"] + "/metaQUAST/report.html",
        transposed_report=config["path"]["output"]+ "/metaQUAST/combined_reference/transposed_report.tsv",
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
        metaquast.py {input.contigs} -o {output.dir}\
            -r {input.reference} --threads {threads} --max-ref-number 0 &> {log}
        """
