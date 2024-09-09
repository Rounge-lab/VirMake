# VirMake: a Snakemake pipeline for viral metagenomic data analysis

## Contents

1. [About VirMake](#about-virmake)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Customization](#customization)
5. [Output explained](#output-explained)

   
## About VirMake

VirMake is a Snakemake based pipeline that offers viral metagenomic data analysis on paired-end data. It offers taxonomic and functional annotation, supports offline running and support for HPC cluster execution. It is made for Linux based systems and has been tested on x86_64-based Linux and with SLURM cluster execution.


## Installation

### Prerequisites:

* Git
* Conda (either Miniconda or Anaconda)
* System Requirements: At least 125 GB of RAM and 180 GB of free disk space.


### Installation steps:

1. Clone the repository using `git clone https://github.com/Rounge-lab/VirMake.git`
2. Run the VirMake setup script `python setup.py`.
   This script will set up a conda environment (./venv/), config files and a table specifying the location of input files.
    * Use the `--working-dir` option to specify the output directory (default is `./results/`).
    * Use the `--reads`, `--qc-reads`, and `--contigs` options to specify the directories where input files are found.

### Directory Structure:
After installation, your directory should look like this:


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
Parameters include: Assembler choice, Memory allocation,Job execution time, Minimum contig size.
Adjust these settings as needed to suit your specific analysis requirements.

### HPC configuration

You can adjust the HPC profile file to suit your needs. The file is located in `./config/config.yaml`


## Output explained

The pipeline provides many files and to help navigate these this section will explain what each section provides.
These are all folders within the `./results/output` folder.


