from scripts.workflow_utils import get_samples, get_min_quality, get_assembly_loc

sample_table, SAMPLE = get_samples(config["path"]["samples"])

# VIRAL IDENTIFICATION #

rule IDENTIFICATION:
    input:
        config["path"]["output"]+"/dereplication/old_to_new_ids.tsv"
    output:
        config["path"]["temp"] + "/finished_IDENTIFICATION",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    message:
        "[IDENTIFICATION] Finished identification of viral contigs"
    shell:
        """
        touch {output}
        """

rule virsorter:
    """
    performs viral identification with virsorter2
    """
    input:
        assembly_output = lambda w: get_assembly_loc(w, sample_table, config["path"]["output"]),
        flag=config["path"]["database"]["virsorter2"] + "/flag"
    output:
        dir=directory(config["path"]["output"]+"/virsorter/{sample}/"),
        viral_combined=config["path"]["output"]+"/virsorter/{sample}/final-viral-combined.fa",
        viral_boundary=config["path"]["output"]+"/virsorter/{sample}/final-viral-boundary.tsv",
        virus_predictions=config["path"]["output"]+"/virsorter/{sample}/viruses.fasta",
        virus_table=config["path"]["output"]+"/virsorter/{sample}/virus_table.tsv",
    params:
        cutoff_length=config["virsorter2"]["id"]["min_length"],
        cutoff_score=config["virsorter2"]["id"]["min_score"],
        groups=config["virsorter2"]["id"]["viral_groups"],
        db_dir=config["path"]["database"]["virsorter2"],
    message:
        "[virsorter2] Executing viral identification..."
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    log:
        config["path"]["log"] + "/virsorter/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/virsorter/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["normal"],
    shell:
        """
        virsorter run -w {output.dir} \
            -j {threads} --include-groups "{params.groups}"\
            --min-length {params.cutoff_length} \
            --min-score {params.cutoff_score} \
            --db-dir {params.db_dir}\
            -i {input.assembly_output} &> {log}
        cp {output.viral_combined} {output.virus_predictions}
        cp {output.viral_boundary} {output.virus_table}
        """

rule genomad:
    """
    performs viral identification using geNomad
    """
    input:
        assembly_output = lambda w: get_assembly_loc(w, sample_table, config["path"]["output"]),
        flag=config["path"]["database"]["genomad"] + "/flag",
    output:
        dir=directory(config["path"]["output"] + "/genomad/{sample}/"),
        viruses = config["path"]["output"] + "/genomad/{sample}/{sample}_summary/{sample}_virus.fna",
        virus_tab = config["path"]["output"] + "/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv",
        virus_predictions = config["path"]["output"] + "/genomad/{sample}/viruses.fasta",
        virus_table = config["path"]["output"] + "/genomad/{sample}/virus_table.tsv"
    params:
        db_dir=config["path"]["database"]["genomad"] + "/genomad_db",
    conda:
        config["path"]["envs"] + "/genomad.yaml"
    message:
        "[genomad] Executing viral identification..."
    log:
        config["path"]["log"] + "/genomad/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/genomad/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["normal"],
    shell:
        """
        genomad end-to-end {input.assembly_output} {output.dir} {params.db_dir} &> {log}
        cp {output.viruses} {output.virus_predictions}
        cp {output.virus_tab} {output.virus_table}
        """

# rule vibrant:
#     """
#     performs the first pass of viral identification with VIBRANT
#     """
#     input:
#         assembly_output = lambda w: get_assembly_loc(w, sample_table, config["path"]["output"]),
#         flag=config["path"]["database"]["vibrant"] + "/flag",
#     output:
#     params:
#         db_dir=config["path"]["database"]["vibrant"] + "/databases",
#         files_dir=config["path"]["database"]["vibrant"] + "/files",
#         virome=vibrant_virome(config["vibrant"]["is_virome"]),
#     conda:
#         config["path"]["envs"] + "/vibrant.yaml"
#     message:
#     log:
#         config["path"]["log"] + "/vibrant/{sample}.log",
#     benchmark:
#         config["path"]["benchmark"] + "/vibrant/{sample}.txt"
#     threads: config["threads"]
#     resources:
#         mem_mb=config["memory"]["big"],
#         runtime=config["time"]["normal"],
#     shell:
#         """
#         VIBRANT_run.py -i {input}\
#             -t {threads}\
#             -folder {output.dir}\
#             -d {params.db_dir}\
#             -m {params.files_dir}\
#             {params.virome}\
#             &> {log}
#     """

rule checkv:
    """
    performs Quality control on identified viral sequences
    """
    input:
        virus_predictions = config["path"]["output"] + "/{id_tool}/{sample}/viruses.fasta",
        flag=config["path"]["database"]["checkv"] + "/flag",
    output:
        dir=directory(config["path"]["output"] + "/checkv/{id_tool}/{sample}/"),
        summary=config["path"]["output"] + "/checkv/{id_tool}/{sample}/quality_summary.tsv",
        contamination=config["path"]["output"] + "/checkv/{id_tool}/{sample}/contamination.tsv",
        combined=config["path"]["output"] + "/checkv/{id_tool}/{sample}/combined.fna",
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    message:
        "[checkv_virsorter2] Executing quality control on identified sequences..."
    log:
        config["path"]["log"] + "/checkv/{id_tool}/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/checkv/{id_tool}/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        checkv end_to_end {input.virus_predictions}\
            {output.dir} -t {threads} -d {params.db_dir} &> {log}
        touch {output.dir}/proviruses.fna
        touch {output.dir}/viruses.fna
        touch {output.summary}
        cat {output.dir}/proviruses.fna {output.dir}/viruses.fna > {output.combined}
        """

rule reformat_virus_prediction:
    input:
        virus_id_table = config["path"]["output"] + "/{id_tool}/{sample}/virus_table.tsv",
        checkv_contamination=config["path"]["output"] + "/checkv/{id_tool}/{sample}/contamination.tsv",
        checkv_res = config["path"]["output"] + "/checkv/{id_tool}/{sample}/quality_summary.tsv",
    output:
        reformatted_table = config["path"]["output"]+"/virus_identification/{sample}/{id_tool}_checkv_reformat.tsv"
    conda:
        config["path"]["envs"] + "/tidyverse.yaml"
    threads:
        1
    script:
        "../scripts/process_viral_identification.R"

rule filter_predicted_viruses:
    input:
        virus_pred_tables = config["path"]["output"]+"/virus_identification/{sample}/"+config["identifier"].strip("2")+"_checkv_reformat.tsv"
    output:
        gathered_qual = config["path"]["output"]+"/virus_identification/{sample}/gathered_quality_tables.tsv",
        representative_selection = config["path"]["output"]+"/virus_identification/{sample}/representative_virus_predictions.tsv",
        regions = config["path"]["output"]+"/virus_identification/{sample}/predicted_viruses.bed"
    params:
        checkv_quality = ["Complete", "High-quality", "Medium-quality"],
        overlap_threshold = 0.1,
        length_selection = "min_length", ## min_length or max_length
        tool_combination = "all" ## any or all
    conda:
        config["path"]["envs"] + "/tidyverse.yaml"
    threads:
        1
    script:
        "../scripts/select_viral_predictions.R"

rule extract_predicted_viruses:
    input:
        contigs = lambda w: get_assembly_loc(w, sample_table, config["path"]["output"]),
        regions = config["path"]["output"]+"/virus_identification/{sample}/predicted_viruses.bed"
    output:
        pred_vir_to_rename = temp(config["path"]["output"]+"/virus_identification/{sample}/tmp_pred_vir.fasta"),
        predicted_viruses = config["path"]["output"]+"/virus_identification/{sample}/predicted_viruses.fasta"
    params:
        index_file = lambda w: get_assembly_loc(w, sample_table, config["path"]["output"])+".fai",
    conda:
        config["path"]["envs"] + "/DRAMv.yaml"
    threads:
        1
    shell:
        """
            bedtools getfasta -fi {input.contigs} -bed {input.regions} -name -fo {output.pred_vir_to_rename}
            [ -f {params.index_file} ] && rm {params.index_file}
            awk -F ':' '/^>/ {{ print $1; next }} {{ print }}' {output.pred_vir_to_rename} > {output.predicted_viruses}
        """

rule gather_checkV_summaries:
    input:
        gathered_qual = expand(config["path"]["output"]+"/virus_identification/{sample}/gathered_quality_tables.tsv", sample=SAMPLE),
    output:
        checkV_gathered = config["path"]["output"]+"/dereplication/checkV_summary.tsv",
        checkM_format = config["path"]["output"]+"/dereplication/checkM_summary.tsv"
    conda:
        config["path"]["envs"] + "/tidyverse.yaml"
    script:
        "../scripts/gather_checkV.R"


rule clear_galah_input_folder:
    input:
        checkV_gathered=config["path"]["output"]+"/dereplication/checkV_summary.tsv"
    output:
        clear_flag=temp(config["path"]["output"]+"/dereplication/clear_input.flag")
    params:
        dir_comb=config["path"]["output"]+"/dereplication/split/"
    shell:
        """
            if [ -d {params.dir_comb} ]
            then
                if [ "$(ls -A {params.dir_comb})" ]; then
                for file in {params.dir_comb}*; do rm "$file"; done
                fi
            else
                if [ ! -d {params.dir_comb} ]
                then
                    mkdir {params.dir_comb}
                fi
            fi

            touch {output.clear_flag}
        """

rule split_fasta_for_galah:
    input:
        predicted_viruses = config["path"]["output"]+"/virus_identification/{sample}/predicted_viruses.fasta",
        clear_flag=config["path"]["output"]+"/dereplication/clear_input.flag"
    output:
        split_flag=temp(config["path"]["output"]+"/dereplication/split_flags/{sample}.flag")
    params:
        split_dir=config["path"]["output"]+"/dereplication/split/"
    shell:
        """
        awk -v FOLDER="{params.split_dir}/" '/^>/ {{ file=FOLDER substr($1,2) ".fna" }} {{ print > file }}' {input.predicted_viruses}
        touch {output.split_flag}
        """

rule dereplication:
    input:
        split_flag=expand(config["path"]["output"]+"/dereplication/split_flags/{sample}.flag", sample = SAMPLE),
        checkM_format=config["path"]["output"]+"/dereplication/checkM_summary.tsv"
    output:
        clusters=config["path"]["output"]+"/dereplication/galah_clusters.tsv",
    params:
        dir_comb=config["path"]["output"]+"/dereplication/split/",
        qual_formula = "completeness-5contamination",
        ani = config["dereplication"]["ani"],
        precluster_ani = config["dereplication"]["precluster_ani"],
        min_aligned_fraction = config["dereplication"]["min_aligned_fraction"],
    log:
        config["path"]["log"] + "/dereplication/dereplication.log"
    conda:
        config["path"]["envs"] + "/galah.yaml"
    threads:
        20
    shell:
        """
            galah cluster \
                    --genome-fasta-directory {params.dir_comb} \
                    --output-cluster-definition {output.clusters} \
                    --ani {params.ani} \
                    --precluster-ani {params.precluster_ani} \
                    --min-aligned-fraction {params.min_aligned_fraction} \
                    --checkm-tab-table {input.checkM_format} \
                    --quality-formula {params.qual_formula} \
                    --threads {threads} &> {log}
        """

rule collect_repr_viral_seqs:
    input:
        galah_clusters=config["path"]["output"]+"/dereplication/galah_clusters.tsv"
    output:
        tmp2=temp(config["path"]["output"]+"/dereplication/rename_these.fasta"),
        genomes=config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
        contig_id_file=config["path"]["output"]+"/dereplication/old_to_new_ids.tsv"
    params:
        votu_num_start = config["dereplication"]["vOTU_num_start"],
        vOTU_prefix = config["dereplication"]["vOTU_prefix"],
        vOTU_suffix = config["dereplication"]["vOTU_suffix"],
        vOTU_num_len = config["dereplication"]["vOTU_num_len"],
        script=config["path"]["scripts"]+"/rename_fasta.py"
    shell:
        """
            [ -f {output.genomes} ] && rm {output.genomes}
            [ -f {output.tmp2} ] && rm {output.tmp2}

            cat $(cut -f -1 {input.galah_clusters} | sort -u) >> {output.tmp2}

            python3 {params.script} \
                -i {output.tmp2} \
                --pre "{params.vOTU_prefix}" \
                --pos "{params.vOTU_suffix}" \
                -o {output.genomes} \
                -s {params.votu_num_start} \
                --int_length {params.vOTU_num_len} \
                --contig_id_file {output.contig_id_file}
        """

rule virsorter_for_dram:
    """
    Runs virsorter2 on the vOTUs
    """
    input:
        genomes=config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
        flag=config["path"]["database"]["virsorter2"] + "/flag"
    output:
        dir=directory(config["path"]["output"] + "/virsorter_for_dram/"),
        finished=config["path"]["output"] + "/virsorter_for_dram/finished",
    params:
        cutoff_length=config["virsorter2"]["for_dramv"]["min_length"],
        cutoff_score=config["virsorter2"]["for_dramv"]["min_score"],
        groups=config["virsorter2"]["for_dramv"]["viral_groups"],
        db_dir=config["path"]["database"]["virsorter2"],
    message:
        "[virsorter_for_dram] Running virsorter2 on the vOTUs..."
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    log:
        config["path"]["log"] + "/virsorter_for_dram.log",
    benchmark:
        config["path"]["benchmark"] + "/virsorter_for_dram.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["big"],
        runtime=config["time"]["normal"],
    shell:
        """
        mkdir -p {output.dir}
        virsorter run -w {output.dir} \
            -j {threads} --include-groups "{params.groups}"\
            --seqname-suffix-off\
            --viral-gene-enrich-off\
            --prep-for-dramv\
            --min-length {params.cutoff_length} \
            --min-score {params.cutoff_score} \
            --keep-original-seq all\
            --db-dir {params.db_dir}\
            -i {input.genomes} &> {log}
        touch {output.finished}
        """


rule checkv_vOTU_virsorter2:
    """
    Runs Quality control on the vOTUs after virsorter2
    """
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    input:
        dir=rules.virsorter_for_dram.output.dir,
        finished=rules.virsorter_for_dram.output.finished,
    output:
        dir=directory(config["path"]["output"] + "/checkv/virsorter_for_dram/"),
        summary=config["path"]["output"]+ "/checkv/virsorter_for_dram/quality_summary.tsv",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    message:
        "[checkv_vOTU_virsorter2] Running checkv on the vOTUs after virsorter2..."
    log:
        config["path"]["log"] + "/checkv_virsorter_for_dram.log",
    benchmark:
        config["path"]["benchmark"] + "/checkv_virsorter_for_dram.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
        --db {params.db_dir}/genome_db/checkv_reps &> {log}
        checkv end_to_end {input.dir}/final-viral-combined.fa {output.dir}\
        -t {threads} -d {params.db_dir} &>> {log}
        """
