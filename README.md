# VirMake: a Snakemake pipeline for viral metagenomic data analysis

## Contents

1. [About VirMake](#about-virmake)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Customization](#customization)


## About VirMake

VirMake is a Snakemake based pipeline that offers viral metagenomic data analysis on paired-end data. It offers taxonomic and functional annotation, supports offline running and has support for HPC cluster execution. It is made for Linux based systems and has been tested on x86_64-based Linux.


## Installation

### Prerequisites:

* Git
* Conda (either Miniconda or Anaconda)
* System Requirements: At least 64 GB of memory and 200 GB of free disk space. Additional memory and disk will be needed depending on the data to be analysed.


### Installation steps:

1. Clone the repository using `git clone https://github.com/Rounge-lab/VirMake.git`
2. Run the VirMake setup script `python setup.py`.
   This script will set up a conda environment (in `./venv/`), config files and a table specifying the location of input files.
    * Use the `--working-dir` option to specify the output directory (default is `./results/`).
    * Use the `--reads`, `--qc-reads`, or `--contigs` options to specify the directories where input files are found.

Reads specified with the `--reads` option or the `--qc_reads` option must be paired-end reads and found in the specified folder in files named either `{sample}_1.fastq` and `{sample}_2.fastq` (for uncompressed reads), or `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz` (for gzip compressed reads), where `{sample}` is the sample name. Multiple samples may be present.

If contigs are specified with the `--contigs` option, they must be found in the specified folder in files named `{sample}.fasta`, where `{sample}` is the sample name. Multiple samples may be present.


## Usage

1. Activate the `virmake` conda environment:

```
conda activate ./venv
```

Ensure all commands are run from the root directory of the repository.

2.Run a dry-run of the pipeline:

```
./virmake run -n
```

Optional Help: For additional help and command options, you can access the built-in help system:

```
./virmake -h
```

This will output the steps required for execution.


### Environments preparation

Running VirMake requires the preparation of software and databases, which may take a substantial amount of time. To prepare a run by downloading and setting up environments, use the `prep` command:

```
./virmake prep
```

To prepare a run by downloading and setting up databases, use the `db` command:

```
./virmake db
```

These commands will set up all requirements for running the virmake pipeline, including the prerequisites for any steps specified in the config file (`config/params.yaml`) under `rule_inclusion`.


### Running the workflow

To run the workflow use:

```
./virmake run
```

To run the workflow with more personalized options please use `./virmake run -h` and read the help page.

Please note that the `-c` (or `--threads`) option controls both the number of jobs and the number of threads in each job when running without using Slurm. Multiple jobs are used for assembly and other steps that can be run separately for each sample. For example, running `./virmake run -c 8` will start 8 jobs using 8 threads each (for a total of 64 threads) during assembly.


### Inspecting results

If you are working on a remote server do not forget to **log in with port tunnelling enabled**.
To do this run:

```
ssh -L 8888:localhost:8888 username@server_address
```

Once the workflow is finished you can inspect the results by running:

```
./virmake inspect
```

and click on the link that will be provided in the terminal.
The link should look like `http://127.0.0.1:8888`.


## Customisation

### Workflow parameters

To adjust the workflow settings edit the `params.yaml` file. The file is located in `VirMake/config/params.yaml`.
Parameters include: Assembler choice, Memory allocation, Job execution time, Minimum contig size.
Adjust these settings as needed to suit your specific analysis requirements.


### HPC configuration

You can adjust the HPC profile file to suit your needs. The file is located in `./config/config.yaml`
