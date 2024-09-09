from scripts.workflow_utils import get_samples

sample_table, SAMPLE = get_samples(config["path"]["samples"])
FRAC = ["1", "2"]

rule STATS:
    input:
        config["path"]["output"] + "/statistics/vOTU_stats.tsv",
        config["path"]["output"] + "/statistics/sample_stats.tsv"
    output:
        config["path"]["temp"] + "/finished_STATS",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        touch {output}
        """

rule get_vOTU_stats:
    """
    Gathers info on vOTUs.
    """
    input:
        taxonomy=config["path"]["output"] + "/graphanalyzer/results_vcontact2_vOTU_results.csv",
        gathered_specs=config["path"]["output"]+"/dereplication/checkV_summary.tsv",
        derep_file=config["path"]["output"]+"/dereplication/galah_clusters.tsv",
        contig_id_file=config["path"]["output"]+"/dereplication/old_to_new_ids.tsv",
        rel_abund=config["path"]["output"] + "/mapping/rel_abundance_table.tsv" if config["rule_inclusion"]["stats"]["mapping"] else [],
        # instrain_by_genome=config["path"]["output"]+"/"
        # DRAM_annotations=config["path"]["output"]+"/DRAMv/annotations/annotations.tsv" if config["run_DRAMv"] else [],
        DRAM_distilled_stats=config["path"]["output"]+"/DRAMv/distilled/vMAG_stats.tsv" if config["rule_inclusion"]["stats"]["dramv"] else [],
    output:
        vOTU_stats=config["path"]["output"] + "/statistics/vOTU_stats.tsv"
    script:
        config["path"]["scripts"] + "/get_vOTU_stats.py"

rule get_sample_stats:
    """
    Gathers info on samples.
    """
    input:
        virus_id_tables=expand(config["path"]["output"]+"/virus_identification/{sample}/gathered_quality_tables.tsv", sample=SAMPLE),
        mq_report=config["path"]["output"]+ "/metaQUAST/combined_reference/transposed_report.tsv" if config["rule_inclusion"]["stats"]["metaquast"] else [],
        rel_abund=config["path"]["output"] + "/mapping/rel_abundance_table.tsv" if config["rule_inclusion"]["stats"]["mapping"] else [],
        flagstat=expand(config["path"]["output"]+"/mapping/flagstat/{sample}_flagstat.txt", sample=SAMPLE) if config["rule_inclusion"]["stats"]["mapping"] else [],
    output:
        sample_stats=config["path"]["output"] + "/statistics/sample_stats.tsv"
    script:
        config["path"]["scripts"] + "/get_sample_stats.py"

