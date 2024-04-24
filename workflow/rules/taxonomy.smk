

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
    performs gene prediction on vOTUs
    """
    input:
        config["path"]["output"] + "/vOTU/vOTU_derep95_combined.fasta",
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
        # cat {output.vOTU_results}
        """
