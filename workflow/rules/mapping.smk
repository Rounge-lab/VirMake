from scripts.workflow_utils import get_samples

SAMPLE, FRAC = get_samples(config["path"]["input"])

# MAPPING
rule MAPPING:
    input:
        config["path"]["output"] + "/mapping/index",
        expand(
            config["path"]["output"] + "/mapping/BAM/{sample}.map.bam",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/postfilter_base_coverage.txt.gz",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/postfilter_coverage_histogram.txt",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/postfilter_coverage_stats.txt",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/postfilter_coverage_binned.txt",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/trimmed_mean_coverage.tsv",
            sample=SAMPLE,
        ),
        config["path"]["output"] + "/contig_stats/raw_coverage_table.tsv",
        expand(
            config["path"]["output"] + "/instrain/{sample}",
            sample=SAMPLE,
        ),
        config["path"]["output"] + "/instrain/compared_samples",
    output:
        config["path"]["temp"] + "/finished_MAPPING",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    message:
        "[MAPPING] Mapping finished..."
    shell:
        """
        touch {output}
        """


rule build_index:
    """
    Builds an index to prepare for mapping
    """
    input:
        config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
    output:
        index_dir=directory(config["path"]["output"] + "/mapping/index"),
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    message:
        "[build_index] Building index for mapping..."
    benchmark:
        config["path"]["benchmark"] + "/build_index.txt"
    log:
        config["path"]["log"] + "/bowtie2_build_index.log",
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        mkdir -p {output.index_dir}
        bowtie2-build {input} {output.index_dir}/mapping_index &> {log}
        """


rule read_mapping:
    """
    performs read mapping between original sample and vOTUs
    """
    input:
        R1=config["path"]["output"] + "/fastp_pe/{sample}_1.fastq",
        R2=config["path"]["output"] + "/fastp_pe/{sample}_2.fastq",
        index_dir=rules.build_index.output.index_dir,
    output:
        config["path"]["output"] + "/mapping/BAM/{sample}.map.bam",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        bowtie2=config["path"]["log"] + "/bowtie2_mapping/{sample}.log",
        samtools=config["path"]["log"] + "/sam_to_bam/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/bowtie2_mapping/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        bowtie2 -p {threads} -x {input.index_dir}/mapping_index \
        -1 {input.R1} -2 {input.R2} 2> {log.bowtie2} | \
        samtools view -b -o {output} - &> {log.samtools}
        """

rule contig_stats:
    """
    Creates simple coverage statisitcs for each read mapping
    """
    input:
        genomes=config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
        bam=rules.read_mapping.output,
    output:
        basecov=config["path"]["output"]
        + "/contig_stats/{sample}/postfilter_base_coverage.txt.gz",
        covhist=config["path"]["output"]
        + "/contig_stats/{sample}/postfilter_coverage_histogram.txt",
        covstats=config["path"]["output"]
        + "/contig_stats/{sample}/postfilter_coverage_stats.txt",
        bincov=config["path"]["output"]
        + "/contig_stats/{sample}/postfilter_coverage_binned.txt",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        config["path"]["log"] + "/contig_stats/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/contig_stats/{sample}.txt"
    message:
        "[contig_stats] Creating coverage statistics for each sample..."
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        pileup.sh ref={input.genomes} in={input.bam} \
            threads={threads} \
            -Xmx{resources.mem_mb}m \
            covstats={output.covstats} \
            hist={output.covhist} \
            basecov={output.basecov} \
            concise=t \
            secondary=t \
            bincov={output.bincov} &> {log}
        """


rule get_trimmed_coverage:
    """
    Gets the trimmed mean of the coverage
    """
    input:
        basecov=rules.contig_stats.output.basecov,
        covstats=rules.contig_stats.output.covstats,
    output:
        trimmed_mean=config["path"]["output"]
        + "/contig_stats/{sample}/trimmed_mean_coverage.tsv",
    params:
        trim_perc=config["trim_percentage"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    log:
        config["path"]["log"] + "/trimmed_mean/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/trimmed_mean/{sample}.txt"
    message:
        "[get_trimmed_coverage] Getting trimmed mean of coverage..."
    script:
        config["path"]["scripts"] + "/trimmed_mean.py"


rule combine_coverage:
    """
    Combines all the coverages into one file,
    prepeares for making the relative abundance
    """
    input:
        covstats=expand(
            config["path"]["output"]
            + "/contig_stats/{sample}/postfilter_coverage_stats.txt",
            sample=SAMPLE,
        ),
    output:
        abundance_table=config["path"]["output"]
        + "/contig_stats/raw_coverage_table.tsv",
    params:
        min_coverage=config["min_coverage"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    conda:
        config["path"]["envs"] + "/tidyverse.yaml"
    message:
        "[combine_coverage] Combining coverage statistics..."
    script:
        config["path"]["scripts"] + "/combine_coverage.R"


rule instrain_profile:
    """
    Create inStrain profiles
    """
    input:
        genome=config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
        mapping=rules.read_mapping.output,
    output:
        directory(config["path"]["output"] + "/instrain/{sample}/"),
    conda:
        config["path"]["envs"] + "/instrain.yaml"
    log:
        config["path"]["log"] + "/instrain/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/instrain/{sample}.txt"
    message:
        "[instrain_profile] Creating inStrain profiles..."
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["normal"],
    shell:
        """
        inStrain profile {input.mapping} {input.genome} -o {output} &> {log}
        """


rule instrain_compare:
    """
    Compare inStrain profiles
    """
    input:
        expand(
            config["path"]["output"] + "/instrain/{sample}/",
            sample=SAMPLE,
        ),
    output:
        directory(config["path"]["output"] + "/instrain/compared_samples"),
    conda:
        config["path"]["envs"] + "/instrain.yaml"
    log:
        config["path"]["log"] + "/instrain/compare.log",
    benchmark:
        config["path"]["benchmark"] + "/instrain/compare.txt"
    message:
        "[instrain_compare] Comparing inStrain profiles..."
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["normal"],
    shell:
        """
        inStrain compare -i {input} -o {output} &> {log}
        """