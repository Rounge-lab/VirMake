from scripts.workflow_utils import get_samples, get_qc_reads_loc
from scripts.get_sample_stats import process_abundance_data

sample_table, SAMPLE = get_samples(config["path"]["samples"])

# MAPPING
rule MAPPING:
    input:
        config["path"]["output"] + "/mapping/rel_abundance_table.tsv",
        config["path"]["output"] + "/instrain/comparison/output/comparison_comparisonsTable.tsv" if config["rule_inclusion"]["all"]["instrain"] else [],
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
        config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
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
        R1=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="1"),
        R2=lambda w: get_qc_reads_loc(w, sample_table, config["path"]["output"],r="2"),
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

rule flagstat:
    """
    Creates a read mapping summary for each sample
    """
    input:
        bam=rules.read_mapping.output,
    output:
        flagstat=config["path"]["output"]+"/mapping/flagstat/{sample}_flagstat.txt",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        config["path"]["log"] + "/mapping/flagstat/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/mapping/flagstat/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        samtools flagstats \
            {input.bam} \
            -O tsv \
            --threads {threads} \
            > {output.flagstat} 2> {log}
        """

rule pileup:
    """
    Creates simple coverage statisitcs for each sample
    """
    input:
        genomes=config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
        bam=rules.read_mapping.output,
    output:
        basecov=config["path"]["output"]
        + "/mapping/pileup/{sample}/postfilter_base_coverage.txt.gz",
        covhist=config["path"]["output"]
        + "/mapping/pileup/{sample}/postfilter_coverage_histogram.txt",
        covstats=config["path"]["output"]
        + "/mapping/pileup/{sample}/postfilter_coverage_stats.txt",
        bincov=config["path"]["output"]
        + "/mapping/pileup/{sample}/postfilter_coverage_binned.txt",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        config["path"]["log"] + "/mapping/pileup/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/mapping/pileup/{sample}.txt"
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

checkpoint combine_coverage:
    """
    Combines all the coverages into one file,
    prepeares for making the relative abundance
    """
    input:
        covstats=expand(
            config["path"]["output"]
            + "/mapping/pileup/{sample}/postfilter_coverage_stats.txt",
            sample=SAMPLE,
        ),
        binned_coverage=expand(
            config["path"]["output"]
            + "/mapping/pileup/{sample}/postfilter_coverage_binned.txt",
            sample=SAMPLE,
        ),
    output:
        covstats=config["path"]["output"] + "/mapping/covstats.tsv",
        mapped_reads=config["path"]["output"] + "/mapping/mapped_reads_table.tsv",
        rel_abundance=config["path"]["output"] + "/mapping/rel_abundance_table.tsv",
    params:
        min_coverage=config["min_coverage"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    conda:
        config["path"]["envs"] + "/tidyverse.yaml"
    script:
        config["path"]["scripts"] + "/combine_coverage.R"


rule instrain_profile:
    """
    Create inStrain profiles
    """
    input:
        genome=config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
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

def samples_for_instrain(wildcards):
    rel_abundance_file = checkpoints.combine_coverage.get(**wildcards).output.rel_abundance
    abundance_summary = process_abundance_data(rel_abundance_file)

    inst_samples = abundance_summary.loc[abundance_summary["n_present"] > 0, "sample_id"]
    return [config["path"]["output"]+f"/instrain/{sa}/" for sa in inst_samples]


rule instrain_compare:
    """
    Compare inStrain profiles
    """
    input:
        lambda w: samples_for_instrain(w)
    output:
        instrain_genome_summary=config["path"]["output"] + "/instrain/comparison/output/comparison_comparisonsTable.tsv",
    params:
        dir=config["path"]["output"] + "/instrain/comparison",
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
        mkdir -p {params.dir}
        inStrain compare --force_compress -i {input} -o {params.dir} &> {log}
        gzip -d {output.instrain_genome_summary}.gz
        """
