# PREAMBLE

# set glob wildcards
(
    SAMPLE,
    FRAC,
) = glob_wildcards(config["path"]["input"] + "/{sample}_{frac}.fastq.gz")

# remove duplicates
SAMPLE = sorted(list(set(SAMPLE)))
FRAC = sorted(list(set(FRAC)))


onstart:
    touch(config["path"]["temp"])
    touch(config["path"]["log"])
    print("Samples: " + ", ".join(SAMPLE))


onsuccess:
    print("Workflow finished successfully!")


onerror:
    print("Error has occured. Please, check log files for more details.")


# GLOBAL FUNCTIONS #


def get_min_quality(threshold):
    """
    Gets the threshold for checkv wuality control.
    """
    if threshold.lower() == "complete":
        return "$8~/(Complete)/"
    elif threshold.lower() == "high":
        return "$8~/(High-quality|Complete)/"
    elif threshold.lower() == "medium":
        return "$8~/(Medium-quality|High-quality|Complete)/"
    elif threshold.lower() == "low":
        return "$8~/(Low-quality|Medium-quality|High-quality|Complete)/"
    elif threshold.lower() == "not-determined":
        return "$8~/(Not-determined|Low-quality|Medium-quality|High-quality|Complete)/"
    else:
        return "$8~/(Medium-quality|High-quality|Complete)/"


def vibrant_virome(is_virome):
    """
    Gathers if the input is virome.
    """
    if is_virome.lower == "yes":
        return "-virome"
    else:
        return ""


rule ALL:
    input:
        config["path"]["temp"] + "/finished_QC",
        config["path"]["temp"] + "/finished_ASSEMBLY",
        config["path"]["temp"] + "/finished_IDENTIFICATION",
        config["path"]["temp"] + "/finished_MAPPING",
        config["path"]["temp"] + "/finished_TAXONOMY",
        config["path"]["temp"] + "/finished_FUNCTION",
        config["path"]["temp"] + "/finished_STATS",
    params:
        temp=config["path"]["temp"],
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        rm -rdf {params.temp}
        """


# QUALITY CONTROL #


rule QC:
    input:
        expand(
            config["path"]["output"] + "/fastqc_raw/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
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
            config["path"]["output"] + "/fastqc_qc/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
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


rule fastqc_raw:
    input:
        expand(
            config["path"]["input"] + "/{sample}_{frac}.fastq.gz",
            sample=SAMPLE,
            frac=FRAC,
        ),
    output:
        expand(
            config["path"]["output"] + "/fastqc_raw/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
        dir=directory(config["path"]["output"] + "/fastqc_raw/"),
    threads: config["threads"]
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    message:
        "[fastqc_raw] Executing FASTQC quality control on raw reads..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    log:
        config["path"]["log"] + "/fastqc_raw.log",
    benchmark:
        config["path"]["benchmark"] + "/fastqc_raw.txt"
    shell:
        """
        fastqc {input} -o {output.dir} -t {threads} &> {log}
        """


rule fastp_pe:
    """
    Performes quality control/pre-processing of the raw reads
    """
    input:
        R1=config["path"]["input"] + "/{sample}_1.fastq.gz",
        R2=config["path"]["input"] + "/{sample}_2.fastq.gz",
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


rule fastqc_qc:
    """
    Performes quality control of processed QC reads
    """
    input:
        expand(
            rules.fastp_pe.output.R1,
            sample=SAMPLE,
        ),
        expand(
            rules.fastp_pe.output.R2,
            sample=SAMPLE,
        ),
    output:
        report=expand(
            config["path"]["output"] + "/fastqc_qc/{sample}_{frac}_fastqc.html",
            sample=SAMPLE,
            frac=FRAC,
        ),
        dir=directory(config["path"]["output"] + "/fastqc_qc/"),
    threads: config["threads"]
    conda:
        config["path"]["envs"] + "/fastqc.yaml"
    message:
        "[fastqc_qc] Executing FASTQC quality control on trimmed reads..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    log:
        config["path"]["log"] + "/fastqc_qc.log",
    benchmark:
        config["path"]["benchmark"] + "/fastqc_qc.txt"
    shell:
        """
        fastqc {input} -o {output.dir} -t {threads} &> {log}
        """


# ASSEMBLY #


rule metaSpades:
    """
    Assembles all sequences with metaSpades
    """
    input:
        R1=rules.fastp_pe.output.R1,
        R2=rules.fastp_pe.output.R2,
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


rule megahit:
    """
    Assembles all sequences with megahit
    """
    input:
        R1=rules.fastp_pe.output.R1,
        R2=rules.fastp_pe.output.R2,
    output:
        dir=directory(
            config["path"]["output"] + "/megahit/{sample}/",
        ),
        contigs=config["path"]["output"] + "/megahit/{sample}/contigs.fasta",
    params:
        temp_dir=config["path"]["temp"] + "/megahit/",
    message:
        "[megahit] Performing assembly of paired end reads..."
    conda:
        config["path"]["envs"] + "/megahit.yaml"
    log:
        config["path"]["log"] + "/megahit/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/megahit/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_b=config["memory"]["megahit"],
        runtime=config["time"]["megahit"],
    shell:
        # mkdir and subsequent rm make sure to create output/megahit dir
        # WITHOUT {sample} directory
        """
        mkdir -p {output.dir}
        rm -r {output.dir}
        megahit -1 {input.R1} -2 {input.R2} -o {output.dir}\
        -t {threads} -m {resources.mem_mb_b} &> {log}
        mv {output.dir}/final.contigs.fa {output.contigs}
        """


if config["assembler"].lower() == "metaspades":

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

    rule metaQUAST:
        """
        Performes quality control on all assembled contigs
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

    assembly_output = rules.metaSpades.output.contigs

elif config["assembler"].lower() == "megahit":

    rule ASSEMBLY:
        input:
            expand(
                config["path"]["output"] + "/megahit/{sample}/",
                sample=SAMPLE,
            ),
            expand(
                config["path"]["output"] + "/megahit/{sample}/contigs.fasta",
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

    rule metaQUAST:
        """
        Performes quality control on all assembled contigs
        with the reference database RefSeq Viral
        """
        params:
            reference=config["path"]["database"]["RefSeq"] + "/viral.1.1.genomic.fna",
        input:
            expand(
                config["path"]["output"] + "/megahit/{sample}/contigs.fasta",
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

    assembly_output = rules.megahit.output.contigs


# VIRAL IDENTIFICATION #


rule IDENTIFICATION:
    input:
        expand(
            config["path"]["output"] + "/virsorter2_pass1/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/virsorter2_pass1/{sample}/finished",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/vibrant/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/vibrant/{sample}/finished",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/virsorter2_pass1/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/genomad/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/genomad/{sample}/finished",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/checkv/virsorter2_pass1/{sample}/quality_summary.tsv",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/checkv/virsorter2_pass1/{sample}/combined.fna",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/vibrant/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/vibrant/{sample}/quality_summary.tsv",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/vibrant/{sample}/combined.fna",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/genomad/{sample}/quality_summary.tsv",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/checkv/genomad/{sample}/combined.fna",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/filtered_virsorter2/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/filtered_virsorter2/{sample}/filtered_combined.fna",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/filtered_virsorter2/{sample}/filtered_contigs",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/filtered_vibrant/{sample}/",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"] + "/filtered_vibrant/{sample}/filtered_contigs",
            sample=SAMPLE,
        ),
        expand(
            config["path"]["output"]
            + "/filtered_vibrant/{sample}/filtered_combined.fna",
            sample=SAMPLE,
        ),
        config["path"]["output"] + "/combined_virsorter2/",
        config["path"]["output"] + "/combined_virsorter2/combined_virsorter2.tsv",
        config["path"]["output"] + "/combined_vibrant/",
        config["path"]["output"] + "/combined_vibrant/combined_vibrant.tsv",
        config["path"]["output"] + "/contig_stats/raw_coverage_table.tsv",
        config["path"]["output"] + "/cdhit/",
        config["path"]["output"] + "/cdhit/derep95_combined.fasta",
        config["path"]["output"] + "/vOTU/",
        config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
        config["path"]["output"] + "/virsorter2_pass2/",
        config["path"]["output"] + "/virsorter2_pass2/finished",
        config["path"]["output"] + "/checkv/virsorter2_pass2/",
        config["path"]["output"] + "/checkv/virsorter2_pass2/quality_summary.tsv",
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


rule virsorter2_pass1:
    """
    Performes the first pass of viral identification with virsorter2
    """
    input:
        assembly_output,
    params:
        cutoff_length=config["virsorter2"]["pass1"]["min_lenght"],
        cutoff_score=config["virsorter2"]["pass1"]["min_score"],
        groups=config["virsorter2"]["pass1"]["viral_groups"],
        db_dir=config["path"]["database"]["virsorter2"],
    output:
        dir=directory(config["path"]["output"] + "/virsorter2_pass1/{sample}/"),
        final_viral_combined=config["path"]["output"]
        + "/virsorter2_pass1/{sample}/final-viral-combined.fa",
        finished=config["path"]["output"] + "/virsorter2_pass1/{sample}/finished",
    message:
        "[virsorter2_pass1] Executing viral identification..."
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    log:
        config["path"]["log"] + "/virsorter2_pass1/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/virsorter2_pass1/{sample}.txt"
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


rule vibrant:
    """
    Performes the first pass of viral identification with VIBRANT
    """
    input:
        assembly_output,
    output:
        dir=directory(config["path"]["output"] + "/vibrant/{sample}/"),
        finished=config["path"]["output"] + "/vibrant/{sample}/finished",
    params:
        db_dir=config["path"]["database"]["vibrant"] + "/databases",
        files_dir=config["path"]["database"]["vibrant"] + "/files",
        virome=vibrant_virome(config["vibrant"]["is_virome"]),
    conda:
        config["path"]["envs"] + "/vibrant.yaml"
    message:
        "[vibrant] Executing viral identification..."
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


rule genomad:
    """
    Performes the first pass of viral identification with geNomad
    """
    input:
        assembly_output,
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
        genomad end-to-end {input} {output.dir} {params.db_dir}
        touch {output.finished}
        """


rule checkv_virsorter2:
    """
    Performes Quality control on identified viral sequences
    """
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    input:
        rules.virsorter2_pass1.output.finished,
        dir=rules.virsorter2_pass1.output.dir,
        final_viral_combined=rules.virsorter2_pass1.output.final_viral_combined,
    output:
        dir=directory(config["path"]["output"] + "/checkv/virsorter2_pass1/{sample}/"),
        summary=config["path"]["output"]
        + "/checkv/virsorter2_pass1/{sample}/quality_summary.tsv",
        combined=config["path"]["output"]
        + "/checkv/virsorter2_pass1/{sample}/combined.fna",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    message:
        "[checkv_virsorter2_pass1] Executing quality control on identified sequences..."
    log:
        config["path"]["log"] + "/checkv_virsorter2_pass1/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/checkv_virsorter2_pass1/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["small"],
    shell:
        """
        mkdir -p {output.dir}
        diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
            --db {params.db_dir}/genome_db/checkv_reps &> {log}
        checkv end_to_end {input.final_viral_combined}\
        {output.dir} -t {threads} -d {params.db_dir} &>> {log}
        cat {output.dir}/proviruses.fna {output.dir}/viruses.fna > {output.combined}
        """


rule checkv_vibrant:
    """
    Performes Quality control on identified viral sequences
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


rule checkv_genomad:
    """
    Performes Quality control on identified viral sequences
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


rule filter_contigs_vibrant:
    """
    Gathers only relevant sequences with above threshold score
    """
    params:
        criteria=get_min_quality(config["quality_threshold"]),
    input:
        summary=rules.checkv_vibrant.output.summary,
        combined=rules.checkv_vibrant.output.combined,
    output:
        dir=directory(
            config["path"]["output"] + "/filtered_vibrant/{sample}/",
        ),
        filtered_contigs_viruses=config["path"]["output"]
        + "/filtered_vibrant/{sample}/filtered_combined.fna",
        filtered_contigs_ID=config["path"]["output"]
        + "/filtered_vibrant/{sample}/filtered_contigs",
    conda:
        config["path"]["envs"] + "/seqtk.yaml"
    message:
        "[filter_contigs_vibrant] Filtering contigs with respect to quality scores (quality threshold: {params.criteria})..."
    benchmark:
        config["path"]["benchmark"] + "/filter_contigs_vibrant/{sample}.txt"
    log:
        config["path"]["log"] + "/filter_contigs_vibrant/{sample}.log",
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


rule filter_contigs_genomad:
    """
    Gathers only relevant sequences with above threshold score
    """
    params:
        criteria=get_min_quality(config["quality_threshold"]),
    input:
        summary=rules.checkv_genomad.output.summary,
        combined=rules.checkv_genomad.output.combined,
    output:
        dir=directory(
            config["path"]["output"] + "/filtered_genomad/{sample}/",
        ),
        filtered_contigs_viruses=config["path"]["output"]
        + "/filtered_genomad/{sample}/filtered_combined.fna",
        filtered_contigs_ID=config["path"]["output"]
        + "/filtered_genomad/{sample}/filtered_contigs",
    conda:
        config["path"]["envs"] + "/seqtk.yaml"
    message:
        "[filter_contigs_genomad] Filtering contigs with respect to quality scores (quality threshold: {params.criteria})..."
    benchmark:
        config["path"]["benchmark"] + "/filter_contigs_genomad/{sample}.txt"
    log:
        config["path"]["log"] + "/filter_contigs_genomad/{sample}.log",
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


rule virsorter2_pass2:
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
        dir=directory(config["path"]["output"] + "/virsorter2_pass2/"),
        finished=config["path"]["output"] + "/virsorter2_pass2/finished",
    message:
        "[virsorter2_pass2] Running virsorter2 on the vOTUs..."
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    log:
        config["path"]["log"] + "/virsorter2_pass2.log",
    benchmark:
        config["path"]["benchmark"] + "/virsorter2_pass2.txt"
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
        dir=rules.virsorter2_pass2.output.dir,
        finished=rules.virsorter2_pass2.output.finished,
    output:
        dir=directory(config["path"]["output"] + "/checkv/virsorter2_pass2/"),
        summary=config["path"]["output"]
        + "/checkv/virsorter2_pass2/quality_summary.tsv",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    message:
        "[checkv_vOTU_virsorter2] Running checkv on the vOTUs after virsorter2..."
    log:
        config["path"]["log"] + "/checkv_virsorter2_pass2.log",
    benchmark:
        config["path"]["benchmark"] + "/checkv_virsorter2_pass2.txt"
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


# MAPPING


rule MAPPING:
    input:
        config["path"]["output"] + "/mapping/index",
        expand(
            config["path"]["output"] + "/mapping/SAM/{sample}.map.sam", sample=SAMPLE
        ),
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
        rules.transform_vOTUs.output.vOTU,
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
    Performes read mapping between original sample and vOTUs
    """
    input:
        R1=rules.fastp_pe.output.R1,
        R2=rules.fastp_pe.output.R2,
        index_dir=rules.build_index.output.index_dir,
    output:
        config["path"]["output"] + "/mapping/SAM/{sample}.map.sam",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        config["path"]["log"] + "/bowtie2_mapping/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/bowtie2_mapping/{sample}.txt"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        bowtie2 -p {threads} -x {input.index_dir}/mapping_index\
        -1 {input.R1} -2 {input.R2} -S {output} &> {log}
        """


rule contig_stats:
    """
    Creates simple coverage statisitcs for each read mapping
    """
    input:
        genomes=config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
        sam=rules.read_mapping.output,
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
        pileup.sh ref={input.genomes} in={input.sam} \
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


rule sam_to_bam:
    """
    Converts SAM to BAM
    """
    input:
        rules.read_mapping.output,
    output:
        config["path"]["output"] + "/mapping/BAM/{sample}.map.bam",
    conda:
        config["path"]["envs"] + "/bowtie2.yaml"
    log:
        config["path"]["log"] + "/sam_to_bam/{sample}.log",
    benchmark:
        config["path"]["benchmark"] + "/sam_to_bam/{sample}.txt"
    message:
        "[sam_to_bam] Converting SAM to BAM..."
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        samtools view -bS {input} > {output}
        """


rule instrain_profile:
    """
    Create inStrain profiles
    """
    input:
        genome=rules.transform_vOTUs.output.vOTU,
        mapping=rules.sam_to_bam.output,
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
        inStrain profile {input.mapping} {input.genome} -o {output}
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
        inStrain compare -i {input} -o {output}
        """


# TAXONOMIC ANALYSIS


rule TAXONOMY:
    input:
        config["path"]["output"] + "/prodigal/proteins.faa",
        config["path"]["output"] + "/prodigal/ORFs.genes",
        config["path"]["output"] + "/vcontact2/genes_2_genomes/g2g.csv",
        config["path"]["output"]
        + "/vcontact2/genes_2_genomes/viral_genomes_combined.csv",
        config["path"]["output"] + "/vcontact2/genes_2_genomes/combined_proteins.faa",
        config["path"]["output"] + "/vcontact2/taxonomic_annotation/",
    output:
        config["path"]["temp"] + "/finished_TAXONOMY",
    threads: 1
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        touch {output}
        """


rule prodigal:
    """
    Performes gene prediction on vOTUs
    """
    input:
        rules.transform_vOTUs.output.vOTU,
    output:
        proteins=config["path"]["output"] + "/prodigal/proteins.faa",
        orf=config["path"]["output"] + "/prodigal/ORFs.genes",
    conda:
        config["path"]["envs"] + "/prodigal.yaml"
    log:
        config["path"]["log"] + "/prodigal.log",
    benchmark:
        config["path"]["benchmark"] + "/prodigal.txt"
    threads: config["threads"]
    message:
        "[prodigal] Predicting genes..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        prodigal -i {input} -o {output.proteins} -a {output.orf} -p meta\
        &> {log}
        """


rule gene2genome:
    """
    Performes gene2genome setup for VCONTACT2
    """
    input:
        rules.prodigal.output.orf,
    output:
        config["path"]["output"] + "/vcontact2/genes_2_genomes/g2g.csv",
    conda:
        config["path"]["envs"] + "/vcontact2.yaml"
    log:
        config["path"]["log"] + "/vcontact2_gene2genome.log",
    benchmark:
        config["path"]["benchmark"] + "/vcontact2_gene2genome.txt"
    message:
        "[vcontact2_gene2genome] Setting up gene2genome..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    threads: config["threads"]
    shell:
        """
        vcontact2_gene2genome -p {input} -o {output} -s 'Prodigal-FAA'\
        &> {log}
        """


rule inphared_setup:
    """
    Adds relevant entires from INPHARED into setup files before VCONTACT2
    """
    params:
        inphared_g2g=config["path"]["database"]["INPHARED"]
        + "/vConTACT2_gene_to_genome.csv",
        inphared_proteins=config["path"]["database"]["INPHARED"]
        + "/vConTACT2_proteins.faa",
        simplify_faa=config["path"]["scripts"] + "/simplify_faa-ffn_derep.py",
    input:
        g2g=rules.gene2genome.output,
        proteins=rules.prodigal.output.proteins,
        orf=rules.prodigal.output.orf,
    output:
        combinedg2g=config["path"]["output"]
        + "/vcontact2/genes_2_genomes/viral_genomes_combined.csv",
        combined_proteins=config["path"]["output"]
        + "/vcontact2/genes_2_genomes/combined_proteins.faa",
    benchmark:
        config["path"]["benchmark"] + "/inphared_setup.txt"
    log:
        config["path"]["log"] + "/inphared_setup.log",
    conda:
        config["path"]["envs"] + "/vibrant.yaml"
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    threads: config["threads"]
    shell:
        """
        cat {input.g2g} {params.inphared_g2g} > {output.combinedg2g}
        sed -i 's/,None_provided/,none/g' {output.combinedg2g}
        python3 {params.simplify_faa} {input.orf}
        cat {input.orf}.simple.faa {params.inphared_proteins} > {output.combined_proteins}
        """


rule vcontact2:
    """
    Performs Taxonomic annotation with VCONTACT2
    """
    input:
        proteins=config["path"]["output"]
        + "/vcontact2/genes_2_genomes/combined_proteins.faa",
        g2g=config["path"]["output"]
        + "/vcontact2/genes_2_genomes/viral_genomes_combined.csv",
    output:
        dir=directory(config["path"]["output"] + "/vcontact2/taxonomic_annotation/"),
    benchmark:
        config["path"]["benchmark"] + "/vcontact2.txt"
    conda:
        config["path"]["envs"] + "/vcontact2.yaml"
    log:
        config["path"]["log"] + "/vcontact2.log",
    resources:
        mem_mb=config["memory"]["vcontact2"],
        runtime=config["time"]["vcontact2"],
    threads: config["threads"]
    shell:
        """
        rm -rdf {output.dir}/combined*
        vcontact2 -t {threads} \
            --raw-proteins {input.proteins} \
            --rel-mode 'Diamond' \
            --proteins-fp {input.g2g} \
            --db 'None' \
            --pcs-mode MCL \
            --vcs-mode ClusterONE \
            --output-dir {output.dir} \
            &> {log}
        """


# FUNCTIONAL ANALYSIS


rule FUNCTION:
    input:
        config["path"]["output"] + "/DRAMv/annotations/",
        config["path"]["output"] + "/DRAMv/distilled/",
        config["path"]["output"] + "/DRAMv/distilled/amg_summary.tsv",
        config["path"]["output"] + "/graphanalyzer/",
        config["path"]["output"] + "/graphanalyzer/csv_edit_vOTU_results.xlsx",
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
        rules.virsorter2_pass2.output.dir,
    output:
        dir=directory(config["path"]["output"] + "/DRAMv/annotations/"),
    params:
        DRAM_config=config["path"]["database"]["DRAM"] + "/DRAM.config",
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
        DRAM-setup.py import_config --config_loc {params.DRAM_config}
        rm -rdf {output.dir}
        DRAM-v.py annotate -i {input}/for-dramv/final-viral-combined-for-dramv.fa\
        -v {input}/for-dramv/viral-affi-contigs-for-dramv.tab\
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
        rules.dramv_annotate.output.dir,
    output:
        dir=directory(config["path"]["output"] + "/DRAMv/distilled/"),
        amg_summary=config["path"]["output"] + "/DRAMv/distilled/amg_summary.tsv",
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
        DRAM-v.py distill -i {input}/annotations.tsv\
        --output_dir {output.dir}\
        &> {log}
        cat {output.amg_summary}
        """


rule graphanalyzer:
    """
    Performs the post-processing of VCONTACT2 results
    automatically
    """
    params:
        graph=config["path"]["scripts"] + "/graphanalyzer.py",
        meta=config["path"]["database"]["INPHARED"] + "/data_excluding_refseq.tsv",
    input:
        rules.vcontact2.output.dir,
    output:
        dir=directory(config["path"]["output"] + "/graphanalyzer/"),
        vOTU_results=config["path"]["output"]
        + "/graphanalyzer/csv_edit_vOTU_results.xlsx",
    conda:
        config["path"]["envs"] + "/graphanalyzer.yaml"
    threads: config["threads"]
    resources:
        mem_mb=config["memory"]["normal"],
        runtime=config["time"]["normal"],
    log:
        config["path"]["log"] + "/graphanalyzer.log",
    benchmark:
        config["path"]["benchmark"] + "/graphanalyzer.txt"
    shell:
        """
        python3 {params.graph}\
        --graph {input}/c1.ntw\
        --csv {input}/genome_by_genome_overview.csv\
        --metas {params.meta}\
        --output {output.dir}/\
        --prefix vOTU\
        --suffix vOTU_results\
        --threads {threads}\
        &> {log}
        cat {output.vOTU_results}
        """


# STATISTICS


rule STATS:
    input:
        config["path"]["output"] + "/statistics/Sample_stats_vibrant.tsv",
        config["path"]["output"] + "/statistics/Sample_stats_virsorter2.tsv",
        config["path"]["output"] + "/statistics/vOTU_AMGs.tsv",
        config["path"]["output"] + "/statistics/",
        config["path"]["temp"] + "/html_files.txt",
        config["path"]["temp"] + "/tables.txt",
        config["path"]["output"] + "/complete_output.zip",
        config["path"]["output"] + "/statistics/compared_samples_comparisonsTable.tsv",
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


rule get_stats:
    """
    Performs the aggregation of all relevant files into summaries tables and plots
    """
    input:
        vOTU_results=rules.graphanalyzer.output.vOTU_results,
        amg_summary=rules.dramv_distill.output.amg_summary,
        virsorter2_summary=rules.checkv_vOTU_virsorter2.output.summary,
        transposed_report=rules.metaQUAST.output.transposed_report,
        abundance_table=rules.combine_coverage.output.abundance_table,
    output:
        config["path"]["output"] + "/statistics/compared_samples_comparisonsTable.tsv",
        config["path"]["output"] + "/statistics/Sample_stats_vibrant.tsv",
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
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["tiny"],
    shell:
        """
        echo "" > {output}
        for table in {params.table_list}; do
            find {input} -name $table >> {output}
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
