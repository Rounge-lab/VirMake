

# FUNCTIONAL ANALYSIS
rule FUNCTION:
    input:
        config["path"]["output"] + "/DRAMv/distilled/vMAG_stats.tsv",
    output:
        config["path"]["temp"] + "/finished_FUNCTION",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        touch {output}
        """


rule dramv_annotate:
    """
    Performs Functional annotation with DRAMv
    """
    input:
        dir=config["path"]["output"] + "/virsorter_for_dram/",
        DRAM_config=config["path"]["database"]["DRAM"] + "/DRAM.config",
    output:
        dir=directory(config["path"]["output"] + "/DRAMv/annotations/"),
        annotations=config["path"]["output"] + "/DRAMv/annotations/annotations.tsv",
    params:
        min_contig_size=config["min_contig_size"],
    conda:
        config["path"]["envs"] + "/DRAMv.yaml"
    log:
        config["path"]["log"] + "/DRAMv.log",
    benchmark:
        config["path"]["benchmark"] + "/DRAMv.txt"
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["big"],
    threads: config["threads"]
    shell:
        """
        DRAM-setup.py import_config --config_loc {input.DRAM_config}
        rm -rdf {output.dir}
        DRAM-v.py annotate -i {input.dir}/for-dramv/final-viral-combined-for-dramv.fa\
            -v {input.dir}/for-dramv/viral-affi-contigs-for-dramv.tab\
            --output_dir {output.dir}\
            --threads {threads}\
            --min_contig_size {params.min_contig_size}\
            &> {log}
        """


rule dramv_distill:
    """
    Performs the distillation of functional annotation
    """
    input:
        rules.dramv_annotate.output.annotations,
    output:
        dir=directory(config["path"]["output"] + "/DRAMv/distilled/"),
        amg_summary=config["path"]["output"] + "/DRAMv/distilled/amg_summary.tsv",
        DRAM_distilled_stats=config["path"]["output"]+"/DRAMv/distilled/vMAG_stats.tsv"
    log:
        config["path"]["log"] + "/DRAMv_distill.log",
    benchmark:
        config["path"]["benchmark"] + "/DRAMv_distill.txt"
    conda:
        config["path"]["envs"] + "/DRAMv.yaml"
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["normal"],
    shell:
        """
        rm -rd {output.dir}
        DRAM-v.py distill -i {input} \
            --output_dir {output.dir} \
            &> {log}
        """
