DB_DIR= config["DB_DIR"]

rule all:
    input:
        DB_DIR+"/Vcontact2_setup_done.txt",
        directory(DB_DIR+"/DRAM/DRAM_data"),
        directory(DB_DIR+"/VIBRANT/vibrant-1.2.1/"),
        DB_DIR+"/checkv/checkv-db-v1.5",
        DB_DIR+"/INPHARED/1Dec2022_vConTACT2_proteins.faa",
        DB_DIR+"/INPHARED/1Dec2022_data_excluding_refseq.tsv",
        DB_DIR+"/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv",
        DB_DIR+"/RefSeq/viral.1.1.genomic.fna",

rule Vcontact2:
    conda:
        "../envs/vcontact2.yaml"
    output:
        "{DB_DIR}/Vcontact2_setup_done.txt",
    shell:
        "touch {output}"
rule metaQUAST:
    output:
        reference="{DB_DIR}/RefSeq/viral.1.1.genomic.fna"
    shell:'''
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
        gunzip viral.1.1.genomic.fna.gz
        mv viral.1.1.genomic.fna {output.reference}
        '''
rule DRAMv:
"""
Downloads DRAMv. 
If Users want to skip this simply add a # symbol in front of:
directory(DB_DIR+"/DRAM/DRAM_data"),

Keep in mind that the pipeline needs this to run.
If downloaded independantly, simply add it to the config with same folder structure
include the DRAM.config from DRAMv here also.
"""
    conda:
        "../envs/DRAMv.yaml"
    output:
        directory("{DB_DIR}/DRAM/DRAM_data")
    threads:
        20
    shell:
        "DRAM-setup.py prepare_databases --output_dir {DB_DIR}/DRAM/DRAM_data --verbose --skip_uniref --threads {threads}"
        ";"
        "DRAM-setup.py export_config --output_file {DB_DIR}/DRAM/DRAM.config"

rule vibrant:
    conda:
        "../envs/vibrant.yaml"
    output:
        directory(DB_DIR+"/VIBRANT/vibrant-1.2.1/")
    shell:'''
        download-db.sh
        $CONDA_PREFIX/share/vibrant-1.2.1/db/databases/VIBRANT_setup.py -test
        cp -R $CONDA_PREFIX/share/vibrant-1.2.1/* {output}
    '''

rule checkv:
    conda:
        "../envs/checkv.yaml"
    output:
        directory("{DB_DIR}/checkv/checkv-db-v1.5")
    shell:
        "checkv download_database {DB_DIR}/checkv/"

rule inphared:
    output:
        "{DB_DIR}/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv",
        "{DB_DIR}/INPHARED/1Dec2022_vConTACT2_proteins.faa",
        "{DB_DIR}/INPHARED/1Dec2022_data_excluding_refseq.tsv"
    shell:'''
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_vConTACT2_gene_to_genome.csv' -nH
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_vConTACT2_proteins.faa' -nH 
    wget 'https://millardlab-inphared.s3.climb.ac.uk/1Dec2022_data_excluding_refseq.tsv' -nH
    mv 1Dec2022_vConTACT2_gene_to_genome.csv {DB_DIR}/INPHARED/1Dec2022_vConTACT2_gene_to_genome.csv
    mv 1Dec2022_vConTACT2_proteins.faa {DB_DIR}/INPHARED/1Dec2022_vConTACT2_proteins.faa
    mv 1Dec2022_data_excluding_refseq.tsv {DB_DIR}/INPHARED/1Dec2022_data_excluding_refseq.tsv
    '''
