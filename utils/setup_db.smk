database_dir = config["database_dir"]
envs_dir = config["envs_dir"]


rule all:
    input:
        database_dir + "/Vcontact2_setup_done.txt",
        database_dir + "/DRAM/DRAM_data",
        database_dir + "/VIBRANT/vibrant-1.2.1/",
        database_dir + "/checkv/checkv-db-v1.5",
        database_dir + "/INPHARED/1Dec2022_vConTACT2_proteins.faa",
        database_dir + "/INPHARED/1Dec2022_data_excluding_refseq.tsv",
        database_dir + "/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv",
        database_dir + "/RefSeq/viral.1.1.genomic.fna",


rule Vcontact2:
    conda:
        envs_dir + "/vcontact2.yaml"
    output:
        database_dir + "/Vcontact2_setup_done.txt",
    shell:
        "touch {output}"


rule metaQUAST:
    output:
        reference=database_dir + "/RefSeq/viral.1.1.genomic.fna",
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
directory(database_dir+"/DRAM/DRAM_data"),
Keep in mind that the pipeline needs this to run.
If downloaded independantly, simply add it to the config with same folder structure
include the DRAM.config from DRAMv here also.
"""
    conda:
        envs_dir + "/DRAMv.yaml"
    output:
        directory(database_dir + "/DRAM/DRAM_data"),
    threads: 20
    shell:
        "DRAM-setup.py prepare_databases --output_dir {database_dir}/DRAM/DRAM_data --verbose --skip_uniref --threads {threads}"
        ";"
        "DRAM-setup.py export_config --output_file {database_dir}/DRAM/DRAM.config"


rule vibrant:
    conda:
        envs_dir + "/vibrant.yaml"
    output:
        directory(database_dir + "/VIBRANT/vibrant-1.2.1/"),
    shell:
        """
        download-db.sh
        $CONDA_PREFIX/share/vibrant-1.2.1/db/databases/VIBRANT_setup.py -test
        cp -R $CONDA_PREFIX/share/vibrant-1.2.1/* {output}
    """


rule checkv:
    conda:
        envs_dir + "/checkv.yaml"
    output:
        directory(database_dir + "/checkv/checkv-db-v1.5"),
    shell:
        "checkv download_database {database_dir}/checkv/"


rule inphared:
    output:
        database_dir + "/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv",
        database_dir + "/INPHARED/1Dec2022_vConTACT2_proteins.faa",
        database_dir + "/INPHARED/1Dec2022_data_excluding_refseq.tsv",
    shell:
        """
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_vConTACT2_gene_to_genome.csv' -nH
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_vConTACT2_proteins.faa' -nH
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_data_excluding_refseq.tsv' -nH
    mv 1Dec2022_vConTACT2_gene_to_genome.csv {database_dir}/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv
    mv 1Dec2022_vConTACT2_proteins.faa {database_dir}/INPHARED/1Dec2022_vConTACT2_proteins.faa
    mv 1Dec2022_data_excluding_refseq.tsv {database_dir}/INPHARED/1Dec2022_data_excluding_refseq.tsv
    """
