identifier = config["identifier"]

rule all:
    input:
        config["path"]["database"]["vcontact2"] + "/vcontact2_setup_done.txt",
        config["path"]["database"]["DRAM"],
        config["path"]["database"][identifier],
        config["path"]["database"]["checkv"],
        config["path"]["database"]["INPHARED"] + "/vConTACT2_proteins.faa",
        config["path"]["database"]["INPHARED"] + "/data_excluding_refseq.tsv",
        config["path"]["database"]["INPHARED"] + "/vConTACT2_gene_to_genome.csv",
        config["path"]["database"]["RefSeq"] + "/viral.1.1.genomic.fna"

rule Vcontact2:
    conda:
        config["path"]["envs"] + "/vcontact2.yaml"
    output:
        config["path"]["database"]["vcontact2"] + "/vcontact2_setup_done.txt",
    shell:
        "touch {output}"


rule metaQUAST:
    output:
        reference=config["path"]["database"]["RefSeq"] + "/viral.1.1.genomic.fna",
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
        gunzip viral.1.1.genomic.fna.gz
        mv viral.1.1.genomic.fna {output.reference}
        """


rule DRAMv:
    """
    Downloads DRAMv.
    If Users want to skip this simply add a # symbol in front of:
    directory(config["path"]["database"]["DRAM"]+"/DRAM_data"),
    Keep in mind that the pipeline needs this to run.
    If downloaded independantly, simply add it to the config with same folder structure
    include the DRAM.config from DRAMv here also.
    """
    conda:
        config["path"]["envs"] + "/DRAMv.yaml"
    output:
        directory(config["path"]["database"]["DRAM"]),
    threads: 24
    shell:
        "DRAM-setup.py prepare_databases --output_dir {output} --verbose --skip_uniref --threads {threads}"
        ";"
        "DRAM-setup.py export_config --output_file {output}/DRAM.config"

if identifier == 'vibrant':
    rule vibrant:
        conda:
            config["path"]["envs"] + "/vibrant.yaml"
        output:
            directory(config["path"]["database"]["vibrant"]),
        shell:
            """
            download-db.sh {output}
            """
elif identifier == "gnomad":
    rule genomad:
        output:
            dir=directory(config["path"]["database"]["genomad"]),
            version=config["path"]["database"]["genomad"] + "/genomad_db/version.txt",
        conda:
            config["path"]["envs"] + "/genomad.yaml"
        threads: 24
        shell:
            """
            genomad download-database {output.dir}/.
            """

rule checkv:
    output:
        directory(config["path"]["database"]["checkv"]),
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    log:
        config["path"]["log"] + "/checkv_db/checkv_db.log",
    shell:
        """
        checkv download_database {output}
        
        diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
            --db {params.db_dir}/genome_db/checkv_reps &> {log}
        """

rule inphared:
    output:
        dir=directory(config["path"]["database"]["INPHARED"]),
        g2g=config["path"]["database"]["INPHARED"] + "/vConTACT2_gene_to_genome.csv",
        prot=config["path"]["database"]["INPHARED"] + "/vConTACT2_proteins.faa",
        exclref=config["path"]["database"]["INPHARED"] + "/data_excluding_refseq.tsv",
        version=config["path"]["database"]["INPHARED"] + "/version.txt",
    params:
        version="1May2023",
    shell:
        """
        wget 'https://millardlab-inphared.s3.climb.ac.uk/{params.version}_vConTACT2_gene_to_genome.csv.gz' \
            -O {output.g2g}.gz -nH
        wget 'https://millardlab-inphared.s3.climb.ac.uk/{params.version}_vConTACT2_proteins.faa.gz' \
            -O {output.prot}.gz -nH
        wget 'https://millardlab-inphared.s3.climb.ac.uk/{params.version}_data_excluding_refseq.tsv.gz' \
            -O {output.exclref}.gz -nH
        gunzip {output.dir}/*.gz
        echo "INPHARED VERSION: {params.version}" > {output.version}
        """

rule virsorter2:
    output:
        directory(config["path"]["database"]["virsorter2"]),
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    threads: 24
    shell:
        """
        virsorter setup -d {output}
        """
