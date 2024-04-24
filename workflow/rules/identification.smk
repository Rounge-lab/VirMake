from scripts.workflow_utils import get_samples, get_min_quality

SAMPLE, FRAC = get_samples(config["path"]["input"])

# VIRAL IDENTIFICATION #

rule IDENTIFICATION:
    input:
        expand(
            config["path"]["output"] + "/virsorter2/{sample}/",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"] + "/virsorter2/{sample}/finished",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"] + "/checkv/virsorter2/{sample}/",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"]
            + "/checkv/virsorter2/{sample}/quality_summary.tsv",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"]
            + "/checkv/virsorter2/{sample}/combined.fna",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"] + "/filtered_virsorter2/{sample}/",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"]
            + "/filtered_virsorter2/{sample}/filtered_combined.fna",
            sample=SAMPLE,
        )
        + expand(
            config["path"]["output"]
            + "/filtered_virsorter2/{sample}/filtered_contigs",
            sample=SAMPLE,
        ),
        config["path"]["output"] + "/combined_virsorter2/",
        config["path"]["output"] + "/combined_virsorter2/combined_virsorter2.tsv",
        config["path"]["output"] + "/contig_stats/raw_coverage_table.tsv",
        config["path"]["output"] + "/cdhit/",
        config["path"]["output"] + "/cdhit/derep95_combined.fasta",
        config["path"]["output"] + "/vOTU/",
        config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
        config["path"]["output"] + "/virsorter_for_dram/",
        config["path"]["output"] + "/virsorter_for_dram/finished",
        config["path"]["output"] + "/checkv/virsorter_for_dram/",
        config["path"]["output"] + "/checkv/virsorter_for_dram/quality_summary.tsv",
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

rule virsorter2:
    """
    performs viral identification with virsorter2
    """
    input:
        assembly_output = config["path"]["output"] + "/metaSpades/{sample}/contigs.fasta",
    params:
        cutoff_length=config["virsorter2"]["pass1"]["min_lenght"],
        cutoff_score=config["virsorter2"]["pass1"]["min_score"],
        groups=config["virsorter2"]["pass1"]["viral_groups"],
        db_dir=config["path"]["database"]["virsorter2"],
    output:
        dir=directory(config["path"]["output"] + "/virsorter2/{sample}/"),
        final_viral_combined=config["path"]["output"]
        + "/virsorter2/{sample}/final-viral-combined.fa",
        finished=config["path"]["output"] + "/virsorter2/{sample}/finished",
    message:
        "[virsorter2] Executing viral identification..."
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    log:
        config["path"]["log"] + "/virsorter2/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/virsorter2/{sample}.txt"
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
            --keep-original-seq all\
            --db-dir {params.db_dir}\
            -i {input} &> {log}
        touch {output.finished}
        """

rule checkv_virsorter2:
    """
    performs Quality control on identified viral sequences
    """
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    input:
        rules.virsorter2.output.finished,
        dir=rules.virsorter2.output.dir,
        final_viral_combined=rules.virsorter2.output.final_viral_combined,
    output:
        dir=directory(config["path"]["output"] + "/checkv/virsorter2/{sample}/"),
        summary=config["path"]["output"]
        + "/checkv/virsorter2/{sample}/quality_summary.tsv",
        combined=config["path"]["output"]
        + "/checkv/virsorter2/{sample}/combined.fna",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    message:
        "[checkv_virsorter2] Executing quality control on identified sequences..."
    log:
        config["path"]["log"] + "/checkv_virsorter2/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/checkv_virsorter2/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        checkv end_to_end {input.final_viral_combined}\
            {output.dir} -t {threads} -d {params.db_dir} &> {log}
        touch {output.dir}/proviruses.fna
        touch {output.dir}/viruses.fna
        touch {output.summary}
        cat {output.dir}/proviruses.fna {output.dir}/viruses.fna > {output.combined}
        """


rule filter_contigs_virsorter2:
    """
    Gathers only relevant sequences with above threshold score
    """
    params:
        criteria=get_min_quality(config["quality_threshold"]),
    input:
        summary=rules.checkv_virsorter2.output.summary,
        combined=rules.checkv_virsorter2.output.combined,
    output:
        dir=directory(
            config["path"]["output"] + "/filtered_virsorter2/{sample}/",
        ),
        filtered_contigs_viruses=config["path"]["output"]
        + "/filtered_virsorter2/{sample}/filtered_combined.fna",
        filtered_contigs_ID=config["path"]["output"]
        + "/filtered_virsorter2/{sample}/filtered_contigs",
    conda:
        config["path"]["envs"] + "/seqtk.yaml"
    message:
        "[filter_contigs_virsorter2] Filtering contigs with respect to quality scores (quality threshold: {params.criteria})..."
    benchmark:
        config["path"]["benchmark"] + "/filter_contigs_virsorter2/{sample}.txt"
    log:
        config["path"]["log"] + "/filter_contigs_virsorter2/{sample}.log",
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        awk -F '\t' '{params.criteria}  {{print $1}}' {input.summary} > {output.filtered_contigs_ID}
        seqtk subseq {input.combined} {output.filtered_contigs_ID} > {output.filtered_contigs_viruses}
        """

if config["identifier"] == "vibrant":

        rule vibrant:
            """
            performs the first pass of viral identification with VIBRANT
            """
            input:
                assembly_output = rules.metaSpades.output.contigs,
            output:
            params:
                db_dir=config["path"]["database"]["vibrant"] + "/databases",
                files_dir=config["path"]["database"]["vibrant"] + "/files",
                virome=vibrant_virome(config["vibrant"]["is_virome"]),
            conda:
                config["path"]["envs"] + "/vibrant.yaml"
            message:
            log:
                config["path"]["log"] + "/vibrant/{sample}.log",
            benchmark:
                config["path"]["benchmark"] + "/vibrant/{sample}.txt"
            threads: config["threads"]
            resources:
                mem_mb=config["memory"]["big"],
                runtime=config["time"]["normal"],
            shell:
                """
                VIBRANT_run.py -i {input}\
                    -t {threads}\
                    -folder {output.dir}\
                    -d {params.db_dir}\
                    -m {params.files_dir}\
                    {params.virome}\
                    &> {log}
                touch {output.finished}
            """


        rule checkv_vibrant:
            """
            performs Quality control on identified viral sequences
            """
            input:
                rules.vibrant.output.finished,
                dir=rules.vibrant.output.dir,
            output:
                dir=directory(
                    config["path"]["output"] + "/checkv/vibrant/{sample}/",
                ),
                summary=config["path"]["output"]
                + "/checkv/vibrant/{sample}/quality_summary.tsv",
                combined=config["path"]["output"] + "/checkv/vibrant/{sample}/combined.fna",
            params:
                db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
            conda:
                config["path"]["envs"] + "/checkv.yaml"
            message:
                "[checkv_vibrant] Executing quality control on identified sequences..."
            log:
                config["path"]["log"] + "/checkv_vibrant/{sample}.log",
            benchmark:
                config["path"]["benchmark"] + "/checkv_vibrant/{sample}.txt"
            threads: config["threads"]
            resources:
                mem_mb=config["memory"]["normal"],
                runtime=config["time"]["small"],
            shell:
                """
                mkdir -p {output.dir}
                diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
                    --db {params.db_dir}/genome_db/checkv_reps &> {log}
                checkv end_to_end {input.dir}/VIBRANT_contigs/VIBRANT_phages_contigs/contigs.phages_combined.fna\
                {output.dir} -t {threads} -d {params.db_dir} &>> {log}
                cat {output.dir}/proviruses.fna {output.dir}/viruses.fna > {output.combined}
                """


        rule combine_results_vibrant:
            """
            Combine all results from virus identifiaction of all samples with VIBRANT
            """
            input:
                expand(
                    config["path"]["output"]
                    + "/filtered_vibrant/{sample}/filtered_combined.fna",
                    sample=SAMPLE,
                ),
            output:
                dir=directory(config["path"]["output"] + "/combined_vibrant/"),
                combined=config["path"]["output"] + "/combined_vibrant/combined_vibrant.tsv",
            conda:
                config["path"]["envs"] + "/cdhit.yaml"
            message:
                "[combine_results_vibrant] Combining results from all samples..."
            log:
                config["path"]["log"] + "/combine_results_vibrant.log",
            benchmark:
                config["path"]["benchmark"] + "/combine_results_vibrant.txt"
            threads: config["threads"]
            resources:
                mem_mb=config["memory"]["small"],
                runtime=config["time"]["small"],
            shell:
                """
                mkdir -p {output.dir}
                cat {input} > {output.combined}
                """

elif config["identifier"] == "genomad":

        rule genomad:
            """
            performs the first pass of viral identification with geNomad
            """
            input:
                assembly_output = rules.metaSpades.output.contigs,
            output:
                dir=directory(config["path"]["output"] + "/genomad/{sample}/"),
                finished=config["path"]["output"] + "/genomad/{sample}/finished",
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
                genomad end-to-end {input} {output.dir} {params.db_dir} &> {log}
                touch {output.finished}
                """

        rule checkv_genomad:
            """
            performs Quality control on identified viral sequences
            """
            input:
                rules.genomad.output.finished,
                dir=rules.genomad.output.dir,
            output:
                dir=directory(
                    config["path"]["output"] + "/checkv/genomad/{sample}/",
                ),
                summary=config["path"]["output"]
                + "/checkv/genomad/{sample}/quality_summary.tsv",
                combined=config["path"]["output"] + "/checkv/genomad/{sample}/combined.fna",
            params:
                db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
            conda:
                config["path"]["envs"] + "/checkv.yaml"
            message:
                "[checkv_genomad] Executing quality control on identified sequences..."
            log:
                config["path"]["log"] + "/checkv_genomad/{sample}.log",
            benchmark:
                config["path"]["benchmark"] + "/checkv_genomad/{sample}.txt"
            threads: config["threads"]
            resources:
                mem_mb=config["memory"]["normal"],
                runtime=config["time"]["small"],
            shell:
                """
                mkdir -p {output.dir}
                diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
                    --db {params.db_dir}/genome_db/checkv_reps &> {log}
                checkv end_to_end {input.dir}/contigs_summary/contigs_virus.fna\
                {output.dir} -t {threads} -d {params.db_dir} &>> {log}
                cat {output.dir}/proviruses.fna {output.dir}/viruses.fna > {output.combined}
                """

        rule combine_results_genomad:
            """
            Combine all results from virus identifiaction of all samples with genomad
            """
            input:
                expand(
                    config["path"]["output"]
                    + "/filtered_genomad/{sample}/filtered_combined.fna",
                    sample=SAMPLE,
                ),
            output:
                dir=directory(config["path"]["output"] + "/combined_genomad/"),
                combined=config["path"]["output"] + "/combined_genomad/combined_genomad.tsv",
            conda:
                config["path"]["envs"] + "/cdhit.yaml"
            message:
                "[combine_results_genomad] Combining results from all samples..."
            log:
                config["path"]["log"] + "/combine_results_genomad.log",
            benchmark:
                config["path"]["benchmark"] + "/combine_results_genomad.txt"
            threads: config["threads"]
            resources:
                mem_mb=config["memory"]["small"],
                runtime=config["time"]["small"],
            shell:
                """
                mkdir -p {output.dir}
                cat {input} > {output.combined}
                """

rule combine_results_virsorter2:
    """
    Combine all results from virus identifiaction of all samples with virsorter2
    """
    input:
        expand(
            config["path"]["output"]
            + "/filtered_virsorter2/{sample}/filtered_combined.fna",
            sample=SAMPLE,
        ),
    output:
        dir=directory(config["path"]["output"] + "/combined_virsorter2/"),
        combined=config["path"]["output"]
        + "/combined_virsorter2/combined_virsorter2.tsv",
    conda:
        config["path"]["envs"] + "/cdhit.yaml"
    message:
        "[combine_results_virsorter2] Combining results from all samples..."
    log:
        config["path"]["log"] + "/combine_results_virsorter2.log",
    benchmark:
        config["path"]["benchmark"] + "/combine_results_virsorter2.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        cat {input} > {output.combined}
        """

rule dereplication:
    """
    Performs dereplication with clustering on the viral sequences
    """
    input:
        rules.combine_results_virsorter2.output.combined,
    output:
        dir=directory(config["path"]["output"] + "/cdhit/"),
        derep=config["path"]["output"] + "/cdhit/derep95_combined.fasta",
    conda:
        config["path"]["envs"] + "/cdhit.yaml"
    message:
        "[dereplication] Dereplication of viral sequences..."
    log:
        config["path"]["log"] + "/dereplication.log",
    benchmark:
        config["path"]["benchmark"] + "/dereplication.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        cd-hit-est\
        -T {threads}\
        -M {resources.mem_mb}\
        -i {input}\
        -o {output.derep}\
        -c 0.95\
        -aS 0.85\
        -n 9\
        -d 0\
        -p 1\
        -t 4\
        -g 1\
        &> {log}
        """

rule transform_vOTUs:
    """
    Renames viral sequences to unique vOTU_# names
    """
    input:
        rules.dereplication.output.derep,
    output:
        dir=directory(config["path"]["output"] + "/vOTU/"),
        vOTU=config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    conda:
        config["path"]["envs"] + "/seqtk.yaml"
    message:
        "[transform_vOTUs] Renaming viral sequences to vOTU_#..."
    log:
        config["path"]["log"] + "/transform_vOTUs.log",
    benchmark:
        config["path"]["benchmark"] + "/transform_vOTUs.txt"
    shell:
        """
        mkdir -p {output.dir}
        awk '/^>/{{print ">vOTU_" ++i; next}}{{print}}' {input} > {output.vOTU}
        """


rule virsorter_for_dram:
    """
    Runs virsorter2 on the vOTUs
    """
    params:
        cutoff_length=config["virsorter2"]["pass2"]["min_lenght"],
        cutoff_score=config["virsorter2"]["pass2"]["min_score"],
        groups=config["virsorter2"]["pass2"]["viral_groups"],
        db_dir=config["path"]["database"]["virsorter2"],
    input:
        rules.transform_vOTUs.output.vOTU,
    output:
        dir=directory(config["path"]["output"] + "/virsorter_for_dram/"),
        finished=config["path"]["output"] + "/virsorter_for_dram/finished",
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
            -i {input} &> {log}
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
        summary=config["path"]["output"]
        + "/checkv/virsorter_for_dram/quality_summary.tsv",
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
