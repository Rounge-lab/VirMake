identifier = config["identifier"]

rule SETUP_DB:
    input:
        config["path"]["database"]["RefSeq"] if config["rule_inclusion"]["all"]["metaquast"] else [],
        config["path"]["database"][identifier] + "/flag" if config["rule_inclusion"]["all"]["identification"] else [],
        config["path"]["database"][identifier] + "/flag" if config["rule_inclusion"]["all"]["identification"] else [],
        config["path"]["database"][identifier] + "/flag" if config["rule_inclusion"]["all"]["identification"] else [],
        config["path"]["database"]["checkv"] if config["rule_inclusion"]["all"]["identification"] else [],
        config["path"]["database"]["vcontact2"] + "/flag" if config["rule_inclusion"]["all"]["taxonomy"] else [],
        config["path"]["database"]["INPHARED"] + "/vConTACT2_proteins.faa" if config["rule_inclusion"]["all"]["taxonomy"] else [],
        config["path"]["database"]["INPHARED"] + "/data_excluding_refseq.tsv" if config["rule_inclusion"]["all"]["taxonomy"] else [],
        config["path"]["database"]["INPHARED"] + "/vConTACT2_gene_to_genome.csv" if config["rule_inclusion"]["all"]["taxonomy"] else [],
        config["path"]["database"]["DRAM"] + "/DRAM.config" if config["rule_inclusion"]["all"]["function"] else [],
    output:
        config["path"]["temp"] + "/finished_DB"
    shell:
        "touch {output}"


rule refseq:
    output:
        reference=config["path"]["database"]["RefSeq"],
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
        gunzip viral.1.1.genomic.fna.gz
        mv viral.1.1.genomic.fna {output.reference}
        """

rule vibrant_db:
    conda:
        config["path"]["envs"] + "/vibrant.yaml"
    output:
        db_dir=directory(config["path"]["database"]["vibrant"]),
        flag=config["path"]["database"]["vibrant"] + "/flag",
    shell:
        """
        download-db.sh {output.db_dir}
        touch {output.flag}
        """

rule genomad_db:
    output:
        dir=directory(config["path"]["database"]["genomad"]),
        version=config["path"]["database"]["genomad"] + "/genomad_db/version.txt",
        flag=config["path"]["database"]["genomad"] + "/flag",
    conda:
        config["path"]["envs"] + "/genomad.yaml"
    threads: 24
    shell:
        """
        genomad download-database {output.dir}/.
        touch {output.flag}
        """

rule vs2_db:
    output:
        dir=directory(config["path"]["database"]["virsorter2"]),
        flag=config["path"]["database"]["virsorter2"] + "/flag",
    conda:
        config["path"]["envs"] + "/virsorter2.yaml"
    threads: 24
    shell:
        """
        virsorter setup -d {output.dir}
        touch {output.flag}
        """

rule checkv_db:
    output:
        dir=directory(config["path"]["database"]["checkv"]),
        flag=config["path"]["database"]["checkv"] + "/flag",
    conda:
        config["path"]["envs"] + "/checkv.yaml"
    params:
        db_dir=config["path"]["database"]["checkv"] + "/checkv-db-v1.5",
    log:
        config["path"]["log"] + "/checkv_db/checkv_db.log",
    shell:
        """
        checkv download_database {output.dir}
        diamond makedb --in {params.db_dir}/genome_db/checkv_reps.faa \
            --db {params.db_dir}/genome_db/checkv_reps &> {log}
        touch {output.flag}
        """

rule DRAMv_db:
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
        dram_config=config["path"]["database"]["DRAM"] + "/DRAM.config",
        dram_dir=directory(config["path"]["database"]["DRAM"]),
    threads: 24
    shell:
        """
        # Temporary fix for vogdb
        # Recreate vog.hmm.tar.gz file without the enclosing hmm folder
        mkdir -p {output.dram_dir}/vogdb
        wget -nv "https://fileshare.lisc.univie.ac.at/vog/latest/vog.hmm.tar.gz"
        tar xzf vog.hmm.tar.gz
        cd hmm
        tar czf vog.hmm.tar.gz *.hmm
        mv vog.hmm.tar.gz {output.dram_dir}/vogdb
        cd ..
        rm -rf vog.hmm.tar.gz hmm

        DRAM-setup.py prepare_databases --output_dir {output.dram_dir} \
            --verbose --skip_uniref --threads {threads} \
            --vogdb_loc {output.dram_dir}/vogdb/vog.hmm.tar.gz
        DRAM-setup.py export_config --output_file {output.dram_config}
        """

rule Vcontact2:
    conda:
        config["path"]["envs"] + "/vcontact2.yaml"
    output:
        config["path"]["database"]["vcontact2"] + "/flag",
    shell:
        "touch {output}"

rule inphared_db:
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
