from scripts.workflow_utils import get_samples

sample_table, SAMPLE = get_samples(config["path"]["samples"])
FRAC = ["1", "2"]

# from scripts.workflow_utils import get_samples

# STATISTICS

rule STATS:
    input:
        config["path"]["output"] + "/statistics/vOTU_stats.tsv",
        # config["path"]["output"] + "/statistics/Sample_stats_vibrant.tsv",
        # config["path"]["output"] + "/statistics/Sample_stats_virsorter2.tsv",
        # config["path"]["output"] + "/statistics/vOTU_AMGs.tsv",
        # config["path"]["output"] + "/statistics/",
        # config["path"]["temp"] + "/html_files.txt",
        # config["path"]["temp"] + "/tables.txt",
        # config["path"]["output"] + "/complete_output.zip",
        # config["path"]["output"] + "/statistics/comparison_comparisonsTable.tsv",
        # config["path"]["benchmark"] + "/_summary/",
        # config["path"]["benchmark"] + "/_summary/merged_benchmarks.csv",
        # config["path"]["benchmark"] + "/_summary/REPORT.txt",
        # config["path"]["benchmark"] + "/_summary/plots",
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
        # instrain_by_genome=config["path"]["output"]+"/"
        # DRAM_annotations=config["path"]["output"]+"/DRAMv/annotations/annotations.tsv" if config["run_DRAMv"] else [],
        # DRAM_distilled=config["path"]["output"]+"/DRAMv/distilled/amg_summary.tsv" if config["run_DRAMv"] else [],
    output:
        vOTU_stats=config["path"]["output"] + "/statistics/vOTU_stats.tsv"
    script:
        config["path"]["scripts"] + "/get_vOTU_stats.py"


rule get_stats:
    """
    Performs the aggregation of all relevant files into summaries tables and plots
    """
    input:
        vOTU_results=config["path"]["output"] + "/graphanalyzer/csv_edit_vOTU_results.xlsx",
        amg_summary=config["path"]["output"] + "/DRAMv/distilled/amg_summary.tsv",
        virsorter2_summary=config["path"]["output"] + "/checkv/virsorter_for_dram/quality_summary.tsv",
        transposed_report=config["path"]["output"] + "/metaQUAST/combined_reference/transposed_report.tsv",
        abundance_table=config["path"]["output"] + "/contig_stats/raw_coverage_table.tsv",
    output:
        config["path"]["output"] + "/statistics/comparison_comparisonsTable.tsv",
#        config["path"]["output"] + "/statistics/Sample_stats_vibrant.tsv",
        config["path"]["output"] + "/statistics/Sample_stats_virsorter2.tsv",
        config["path"]["output"] + "/statistics/vOTU_AMGs.tsv",
        dir=directory(config["path"]["output"] + "/statistics/"),
    params:
        samples=SAMPLE,
        output_path=config["path"]["output"],
    log:
        config["path"]["log"] + "/get_stats.log",
    benchmark:
        config["path"]["benchmark"] + "/get_stats.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["normal"],
    script:
        config["path"]["scripts"] + "/table_stats.py"

rule find_reports:
    """
    Finds all html files in the output directory
    """
    input:
        config["path"]["temp"] + "/finished_QC",
        config["path"]["temp"] + "/finished_ASSEMBLY",
        config["path"]["temp"] + "/finished_IDENTIFICATION",
        config["path"]["temp"] + "/finished_MAPPING",
        config["path"]["temp"] + "/finished_TAXONOMY",
        config["path"]["temp"] + "/finished_FUNCTION",
    output:
        config["path"]["temp"] + "/html_files.txt",
    params:
        all_output_path=config["path"]["output"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        find {params.all_output_path} -name *.html > {output}
        find {params.all_output_path} -name *.pdf >> {output}
        """


rule find_tables:
    input:
        config["path"]["temp"] + "/finished_QC",
        config["path"]["temp"] + "/finished_ASSEMBLY",
        config["path"]["temp"] + "/finished_IDENTIFICATION",
        config["path"]["temp"] + "/finished_MAPPING",
        config["path"]["temp"] + "/finished_TAXONOMY",
        config["path"]["temp"] + "/finished_FUNCTION",
    output:
        config["path"]["temp"] + "/tables.txt",
    params:
        table_list=config["include_tables"],
        all_output_path=config["path"]["output"],
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        rm -f {output}
        for table in {params.table_list}; do
            find {params.all_output_path} -name $table >> {output}
        done
        """


rule zip_output:
    """
    Copies all html files to the statistics directory
    """
    input:
        reports=rules.find_reports.output,
        tables=rules.find_tables.output,
    output:
        config["path"]["output"] + "/complete_output.zip",
    params:
        temp=config["path"]["temp"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    script:
        config["path"]["scripts"] + "/zip_output.py"


# rule gather_benchmarks:
#     input:
#         config["path"]["benchmark"],
#     output:
#         dir=directory(config["path"]["benchmark"] + "/_summary/"),
#         merged=config["path"]["benchmark"] + "/_summary/merged_benchmarks.csv",
#         report=config["path"]["benchmark"] + "/_summary/REPORT.txt",
#     threads: 1
#     resources:
#         mem_mb=config["memory"]["small"],
#         runtime=config["time"]["tiny"],
#     script:
#         config["path"]["scripts"] + "/gather_benchmarks.py"
# rule create_benchmark_plots:
#     input:
#         rules.gather_benchmarks.output.merged,
#     output:
#         directory(config["path"]["benchmark"] + "/_summary/plots"),
#     threads: 1
#     resources:
#         mem_mb=config["memory"]["small"],
#         runtime=config["time"]["tiny"],
#     script:
#         config["path"]["scripts"] + "/create_benchmark_plots.R"

#    """
#    Copies all html files to the statistics directory
#   """
#    input:
#        reports=rules.find_reports.output,
#        tables=rules.find_tables.output,
#    output:
#        config["path"]["output"] + "/complete_output.zip",
#    params:
#        temp=config["path"]["temp"],
#    threads: 1
#    resources:
#        mem_mb=config["memory"]["small"],
#        runtime=config["time"]["tiny"],
#    script:
#        config["path"]["scripts"] + "/zip_output.py"
