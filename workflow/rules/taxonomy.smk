

# TAXONOMIC ANALYSIS
rule TAXONOMY:
    input:
        config["path"]["output"] + "/vcontact3/vOTU_assignments.csv",
        # config["path"]["output"] + "/graphanalyzer/results_vcontact2_vOTU_results.csv",
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
    performs gene prediction on vOTUs
    """
    input:
        config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
    output:
        proteins=config["path"]["output"] + "/prodigal/proteins.faa",
        orf=config["path"]["output"] + "/prodigal/ORFs.genes",
    conda:
        config["path"]["envs"] + "/prodigal.yaml"
    log:
        config["path"]["log"] + "/prodigal.log",
    benchmark:
        config["path"]["benchmark"] + "/prodigal.txt"
    threads: 1
    message:
        "[prodigal] Predicting genes..."
    resources:
        mem_mb=config["memory"]["small"],
        runtime=config["time"]["normal"],
    shell:
        """
        prodigal -i {input} -o {output.proteins} -a {output.orf} -p meta\
        &> {log}
        """


rule gene2genome:
    """
    performs gene2genome setup for VCONTACT2
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
    input:
        g2g=rules.gene2genome.output,
        proteins=rules.prodigal.output.proteins,
        orf=rules.prodigal.output.orf,
        inphared_g2g=config["path"]["database"]["INPHARED"] + "/vConTACT2_gene_to_genome.csv",
        inphared_proteins=config["path"]["database"]["INPHARED"] + "/vConTACT2_proteins.faa",
    output:
        combinedg2g=config["path"]["output"] + "/vcontact2/genes_2_genomes/viral_genomes_combined.csv",
        combined_proteins=config["path"]["output"] + "/vcontact2/genes_2_genomes/combined_proteins.faa",
    params:
        simplify_faa=config["path"]["scripts"] + "/simplify_faa-ffn_derep.py",
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
        cat {input.g2g} {input.inphared_g2g} > {output.combinedg2g}
        sed -i 's/,None_provided/,none/g' {output.combinedg2g}
        python3 {params.simplify_faa} {input.orf}
        cat {input.orf}.simple.faa {input.inphared_proteins} > {output.combined_proteins}
        """


rule vcontact2:
    """
    Performs Taxonomic annotation with VCONTACT2
    """
    input:
        proteins=config["path"]["output"]+ "/vcontact2/genes_2_genomes/combined_proteins.faa",
        g2g=config["path"]["output"]+ "/vcontact2/genes_2_genomes/viral_genomes_combined.csv",
        flag=config["path"]["database"]["vcontact2"] + "/flag",
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
        vcontact_res=config["path"]["output"] + "/graphanalyzer/csv_edit_vOTU_results.xlsx",
        vOTU_results=config["path"]["output"] + "/graphanalyzer/results_vcontact2_vOTU_results.csv",
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
        # cat {output.vOTU_results}
        """


rule vcontact3:
    """
    Performs Taxonomic annotation with vConTACT3
    """
    input:
        genomes=config["path"]["output"]+"/dereplication/repr_viral_seqs.fasta",
        db_flag=config["path"]["database"]["vcontact3"] + "/flag",
    output:
        final_assignments=config["path"]["output"] + "/vcontact3/exports/final_assignments.csv",
        vOTU_assignments=config["path"]["output"] + "/vcontact3/vOTU_assignments.csv",
        performance_metrics=config["path"]["output"] + "/vcontact3/exports/performance_metrics.csv",
    params:
        db_domain="prokaryotes",
        db_path=config["path"]["database"]["vcontact3"],
        output_path=config["path"]["output"] + "/vcontact3",
        vOTU_prefix=config["dereplication"]["vOTU_prefix"],
        vOTU_suffix=config["dereplication"]["vOTU_suffix"],
    benchmark:
        config["path"]["benchmark"] + "/vcontact3.txt"
    conda:
        config["path"]["envs"] + "/vcontact3.yaml"
    log:
        config["path"]["log"] + "/vcontact3.log",
    resources:
        mem_mb=config["memory"]["vcontact3"],
        runtime=config["time"]["vcontact3"],
    threads: config["threads"]
    shadow: "minimal"
    shell:
        """

        db_version=$(basename {params.db_path}/*.json .json)

        vcontact3 run \
            --nucleotide {input.genomes} \
            --output {params.output_path} \
            --threads {threads} \
            --db-path {params.db_path} \
            --db-domain {params.db_domain} \
            --db-version "$db_version" \
            --no-progress \
            --force-overwrite \
            --verbose \
            &> {log}
        
        awk -F',' 'NR==1 || $1 ~ /^{params.vOTU_prefix}.*{params.vOTU_suffix}$/' {output.final_assignments} > {output.vOTU_assignments}
        """
