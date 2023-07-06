# VirMake: a Snakemake pipeline for viral metagenomic data analysis

# About the software:

VirMake is a Snakemake based pipeline that offers viral metagenic data analysis on paired-end data. It offers taxonomic and functional annotation, supports offline running and support for HPC cluster execution. It is made for Linux based systems and has been tested on SLURM cluster execution.

# Usage guide:

## Quick start

`git clone https://github.com/uio-bmi/VirMake.git`

`mamba env create --name virmake --file=virmake.yaml`

`python virmake.py download -d ./databases/`

`python virmake.py prep-offline --threads 8`

Place samples in samples folder, and then run:

`python virmake.py init -d ./databases/ ./samples/`

Get example files by:

`cd ./samples/`

`python ../workflow/scripts/get_example.py`

`cd ../`

Then run:

`python virmake.py run all -c config.yaml --threads 24`


## Installation

Git clone the project/download the workflow from GitHub.
It requires conda/mamba to install the correct environments. When this has been installed, the necessary packages can be installed with the provided YAML file `virmake.yaml`. This can be done with the command:

`conda env create --name virmake --file=virmake.yaml`

Or with mamba:

`mamba env create --name virmake --file=virmake.yaml`

The pipeline requiers two starting files to be initialized before running. These are the `config.yaml` and a `samples.tsv`. The `config.yaml` contains parameters for running the Snakemake, and is where users can customize their analysis. The `samples.tsv` is a file containing the name of each sample, thiscan be user made or created with initialization command.
Note: VirMake follows a precise naming convention for samples, they follow the convention of `SAMPLENAME_R1.fastq.gz` where `SAMPLENAME` is the variable name. All samples must be within the same folder.

An expected structuring of samples and other files looks like this:

        VirMake
            ./databases/        % Contains all databases needed by VirMake
            ./samples/          % Here is the folder for samples
            ./workflow/         % Workflow files
            ./config.yaml       % The generated yaml config file
            ./samples.tsv       % The generated file of all samples

Make sure you are within the directory of VirMake and run:

    Usage: virmake.py init [OPTIONS] SAMPLES_PATH

    python virmake.py init [OPTIONS]
        Options:
      -d, --db-dir PATH       location to store databases
      -w, --working-dir PATH  location for running the application
      --threads INTEGER       number of threads to use per multi-threaded job
      -h, --help              Show this message and exit.

If you followed the expected structure you would run:

`python virmake.py init -d ./databases/ ./samples/`

This will create a basic config file and a sample table which are both needed to run the workflow.

The last required step is to download all needed databases and external files. This command downloads the DRAMv database, it requiers 125 GB minimum RAM and around 35 GB of disk space.
If this creates problems or you want to download it yourself, follow instruction within the `workflow/rules/download.smk`.
Downloading the databases can be done with the command:

    Usage: virmake.py download [OPTIONS]          % Requiers the database location to be provided.

    python virmake.py download [OPTIONS]
    Options:
      -d, --db-dir PATH  location to store databases  [required]
      --threads INTEGER  number of threads to use per multi-threaded job
      -n, --dryrun       Test execution.
      -h, --help         Show this message and exit.

If you have the standard databases location, simply run:

`python virmake.py download -d ./databases/ --threads 16`

Depending on if you run the pipeline locally or on a cluster node, you may need to pre-download all environments for Snakemake.
This can be done with the `prep-offline` command

For a smooth performance and setup for eventual offline running use:

    Usage: virmake.py prep-offline [OPTIONS]

    python virmake.py prep-offline [OPTIONS]
        Options:
      --threads INTEGER  number of threads to use per multi-threaded job
      -h, --help         Show this message and exit.

Simply run:

`python virmake.py prep-offline --threads 8`

This will generate all the needed environments and images to run the pipeline. It can take a while to generate all environments, but it only needs to be done once. This can be beneficial to do if you are running on job clusters that limit internet connections too.

For offline usage simply copy the whole directory over to the machine/server you want to use the pipeline on.
The only thing to watch out for is to edit the config file to have the correct paths to other files that were transferred.
This can easily be done by opening the config file and use search and replace on all earlier paths such as:

Original path on machine with internet:

    c:/home/something/databases/....

To new machine:

    /cluster/project/group/user/databases/...

##  Running the workflow
To run the workflow after the setup steps, simply use the command:

    usage: python virmake.py run [OPTIONS] {qc|assembly|identification|taxonomy|all|None}

    python virmake.py [OPTIONS]
       Options:
      --profile TEXT          snakemake profile e.g. for cluster execution.
      -w, --working-dir PATH  location to run pending.
      -c, --config-file PATH  config-file generated with 'pending init'
      -n, --dryrun            Test execution.
      --threads INTEGER       Number of threads used on multithreaded jobs
      -h, --help              Show this message and exit.

A standard setup will need to run:


`python virmake.py run all -c config.yaml --threads 24`

Or if runing on a Cluster:

`python virmake.py run all --profile Profile/config.yaml -c config.yaml --threads 24`

When running on a cluster, lookup seting up a cluster execution profile file from Snakemake website: [cluster profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
Or look for inspiration within the  `EXAMPLE_PROFILE_CONFIG.yaml`

# Regarding databases.

All databases used can can be downloaded manually and later added to the config, or if some are already downloaded.
Simply edit the config file paths, but follow the structure as in the example config.

# Results explained

The pipeline provides many files and to help navigate these this section will explain what each section provides.
These are all folders within the `results` folder.

## cdhit

The `prep` folder within contains the identified viral contigs from virsorter2, vibrant and a combined file.
The folder contains the cluster file produced by cdhit and the dereplicated file, as well as the renamed to vOTU dereplicated file

## checkv

This folder contains all the checkV resulst grouped by what they were run on, and which sample it is. This includes vibrant, virsorter and the vOTUs. The interesting files here are the `quality_summary.tsv` which is the summarized result of checkv for that run. And in the `filtered` folder, contains two files with only quality controlled contig names within `filtered_contigs` and their fasta sequence in `filtered_combined.fna`


## contig_stats

This folder contains the pileup.sh results and coverage statistics. An intersting file that is used in the aggregation is the `trimmed_mean_coverage.tsv`, which is used when generating the relative abundance file for statistics folder.

## DRAMv

This folder contains the results from DRAMv, both annotate and distilled. The most relevant files can be found within the `distilled` folder. The `amg_summary.tsv` contains all the AMG and functional annotation information. The `product.html` file Is a heatmap of all AMGs, where they exist, their function and how many there are within the vOTU.

## fastqc

Thios folder contains the FastQC results on both the RAW reads and the quality controlled reads. The provided html files for each sequence gives an overview of the quality statistics of each sample.


## graphanalyzer

This folder contains all graphanalyzer results. Most relevant is the folder `single-views_vOTU_results` which contain an interactive plot of the clusters for each vOTU. Another important file is the `results_vcontact2_vOTU_results.csv` which contains the proccessed Vcontact2 output and contains all relevant taxonomic clasification information.

## mapping

Contains all the index and sam files from bowtie2 building and maping. The sam files can be used for further analysis if the users want to.

## metaQUAST

This folder contains the quality controlled reports from all assembled contigs within each sample. The `summary` folder contains summaries of all quality control processes and the `combined_reference` folder contains results pertaining to comparisons towards the reference database of RefSeq Viral.


## metaSpades_assembly

This folder contains all assembled contigs ordered by sample. The most relevant file here is the assembled contig file `contigs.fasta`

## prodigal

This folder contains the results from running prodigal and provides the predicted genes and proteins. The pipeline uses a simplified format of these with the file `orfs.genes.simple.faa`.

## statistics

This folder contains the aggregated statistics and plots for the pipeline.
The Taxonomic annotation information can be found in the three files:
`vOTU_stats_combined.tsv` `vOTU_stats_vibrant.tsv` `vOTU_stats_virsorter2.tsv`. They give insight of the taxonomic classification at Family, Subfamily and Genus level. it also includes the checkv quality score, accession number and if the identified virus is a provirus.

The functional annotation can be found in the file `vOTU_AMGs.tsv`. It provides the protein/gene, origin scaffold, ID and the fucntional description.

Some interesting files for seeing the state of all samples at different stages can be found in: `Sample_stats_vibrant.tsv`, `Sample_stats_virsorter2.tsv` and `Combined_Sample_stats.tsv`.

The file `vOTU_mapped_to_reads.tsv` contains the vOTUs mapped back to their original sequences and if they are lytic or not.
## trimmed

This folder contains the fastp quality controlled raw reads and the relevant reports on all samples.


## vcontact2

This folder contains the VCONTACT2 output. The folder `genes_2_genomes` contains the files used when introducing the INPHARED database to be included in the taxonomic annotation. For further analysis the `c1.clusters` and `c1.ntw` can be used and viewed within [Cytoscape](https://cytoscape.org/).

## vibrant

This folder contains all VIBRANT results grouped by sample and one for the vOTUs. It can be a bit tricky to navigate these but they contain a lot of interesting files. The folder `VIBRANT_results` contains the several tables produced by VIBRANT. The direct viral sequences used by the pipeline is gathered from the folder `VIBRANT_phages/VIBRANT_contigs/contigs.phages_combined.fna` and `VIBRANT_vOTU_derep95_combined/VIBRANT_phages_vOTU_derep95_combined/vOTU_derep95_combined.phages_combined.fna`. The relevant taxonomic information can be found within `VIBRANT_results_contigs/VIBRANT_genome_quality_contigs.tsv` and the different AMG_ files within `VIBRANT_results`

## virsorter2

This folder contains all results from virsorter2. The most relevant file here is the `final-viral-score.tsv` file, which contains the scorings for each contig and what type of virus it was deemed as.
