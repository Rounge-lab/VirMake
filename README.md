# VirMake: a Snakemake pipeline for viral metagenomic data analysis

## Contents

1. [About VirMake](#about-virmake)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Output explained](#output-explained)

## About VirMake

VirMake is a Snakemake based pipeline that offers viral metagenomic data analysis on paired-end data. It offers taxonomic and functional annotation, supports offline running and support for HPC cluster execution. It is made for Linux based systems and has been tested on SLURM cluster execution. It has been tested with x86_64-based Linux.

<!-- ![Flowchart](img/flowchart.png) -->

## Installation

### Prerequisites:
- Git
- Conda (Miniconda or Anaconda) installation with write permissions. Conda must be activated with its base environment.
- At least 125 GB of RAM and 180 GB of free disk space. Additional disk space will be needed for the output files depending on the size and number of your samples!

### To install VirMake follow these steps:

1. Clone the repository using `git clone https://github.com/Rounge-lab/VirMake.git`
2. Run the VirMake setup script `python setup.py -y`. The `-y` flag will automatically trigger installation of the databases. If you want to install the databases yourself, you can omit the flag. Also, **if you encounter any errors** during the installation, please run the setup script again without the flag. We recommend using `screen` or `tmux`to run the setup script.
3. Sit back and relax... The installation will take a while, especially if setting up databases.
4. Check `setup.log` file for any errors. If there are no errors, you are ready to go!

An expected directory structure should look like this (files are omitted for readability):

```
VirMake
├── databases
│   ├── checkv
│   ├── DRAM
│   ├── genomad
│   ├── INPHARED
│   ├── RefSeq
│   ├── vcontact2
│   ├── vibrant
│   └── virsorter2
├── envs
├── img
├── utils
├── venv
├── workflow
│   ├── .snakemake
│   ├── config
│   ├── scripts
│   └── Snakefile
└── working_dir
    └── input
```

## Usage

Before you do anything, make sure that `virmake` conda environment that was created during the installation is activated. You can do this by running:

```
conda activate ./venv
```

Also make sure to run all commands from the root directory of the repository.

Please note that `virmake` has inbuilt help that can be accessed by running:

```
./virmake -h
```

To run a dry-run of the pipeline run the following command:

```
./virmake run -n
```

This should produce a snakemake output with the steps to run.

### Environments preparation

VirMake uses conda environments to run separate rules. These environments will be set up on the first run and will be stored in `./workflow/.snakemake/` folder for subsequent runs. This takes a substantial amount of time on the first run. If you want to set up the environments before running the workflow use:

```
./virmake prep
```

After that you can run the workflow offline.

### Getting samples

To run the workflow you will need to provide input files. The input files should be placed in `./working_dir/input/` folder. These need to be in `.fastq.gz` format. The input files should be named in the following format: `<sample_name>_1.fastq.gz` and `<sample_name>_2.fastq.gz`. VirMake only works for paired-end reads. You need at least two samples to be able to complete the pipeline with a comparison of the samples.

You can also download samples from SRA database by using:

```
./virmake get SRA <accession_number>
```

This command will download the samples from SRA database, place them in `VirMake/working_dir/input/` folder, `gzip` them and rename them accordingly.

You can get example files from [1] by running:

```
./virmake get SRA PRJNA524703
```

[1] Liang, G., Zhao, C., Zhang, H. et al. The stepwise assembly of the neonatal virome is modulated by breastfeeding. Nature 581, 470–474 (2020). [https://doi.org/10.1038/s41586-020-2192-1](https://doi.org/10.1038/s41586-020-2192-1)

### Running the workflow

**Before running the workflow, add all your samples to the working_dir under the input folder.**

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
This will launch shiny app within your browser.

### Workflow params file

To adjust the workflow settings edit the `params.yaml` file. The file is located in `VirMake/workflow/config/params.yaml`.
The file contains the following adjustable settings (default values are provided in square brackets):

```
assembler: [metaSpades]                     # assembler being used
                                            # currently only metaSpades is supported
cd-hit-est:
    coverage: [0.85]                        # coverage threshold for cd-hit-est
    identity_threshold: [0.95]              # identity threshold for cd-hit-est
job_type:
    big: [bigmem]
    normal: [normal]
    small: [normal]
memory:                                     # memory requirements for each job type in megabytes
    big: [32000]
    metaquast: [63000]
    normal: [16000]
    small: [8000]
    tiny: [1000]
    vcontact2: [63000]
min_contig_size: [1000]                     # minimum contig size for DRAMv annotate
min_coverage: [75]                          # minimum coverage for combine_coverage.R script

###############################################################################################

# absolute paths to various folders being used by the pipeline
# we do not recommend changing these unless you know what you are doing!
path:
    benchmark: [/.../VirMake/working_dir/benchmark]
    database:
        DRAM: [/.../VirMake/databases/DRAM]
        INPHARED: [/.../VirMake/databases/INPHARED]
        RefSeq: [/.../VirMake/databases/RefSeq]
        checkv: [/.../VirMake/databases/checkv]
        vcontact2: [/.../VirMake/databases/vcontact2]
        vibrant: [/.../VirMake/databases/vibrant]
        virsorter2: [/.../VirMake/databases/virsorter2]
    envs: [/.../VirMake/envs]
    input: [/.../VirMake/working_dir/input]
    log: [/.../VirMake/working_dir/log]
    output: [/.../VirMake/working_dir/output]
    profile: []
    scripts: [/.../VirMake/workflow/scripts]
    temp: [/.../VirMake/working_dir/temp]
    virmake: [/.../VirMake]

###############################################################################################

quality_threshold: [medium]                 # quality threshold for fastp
threads: [24]                               # minimum number of threads to use for parallelized jobs
time:                                       # time requirements for each job type
    big: [13 h]
    metaquast: [24 h]
    normal: [6 h]
    small: [1 h]
    tiny: [30 min]
    vcontact2: [24 h]
trim_percentage: [0.95]
vibrant:
    is_virome: ['no']                       # is the sample a virome? ('yes' or 'no')
virsorter2:
    pass1:
        min_length: [3000]                  # minimum contig length for virsorter2
                                            # contigs with length < 3000 will be discarded

        min_score: [0.5]                    # minimum quality score for virsorter2

        # viral groups searched for by virsorter2
        viral_groups: [dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae]
    pass2:
        min_length: [1000]
        min_score: [0.5]
        viral_groups: [dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae]
```

### HPC profile file

You can adjust the HPC profile file to suit your needs. The file is located in `./workflow/config/config.yaml`

The profile file should look like this:

```
---
reason: True
show-failed-logs: True
keep-going: True
printshellcmds: True

# Cluster submission

# Provide a custom name for the jobscript that is submitted to the cluster
jobname: "{rule}.{jobid}"

# Maximal number of cluster/drmaa jobs per second, fractions allowed
max-jobs-per-second: 20

# Maximal number of job status checks per second
max-status-checks-per-second: 10

cluster: "sbatch -A [INSTER_ACCOUNT] --output=slurm_out/slurm-%j.out -J {rule}_{wildcards} --mem={resources.mem_mb} --time={resources.runtime} --cpus-per-task={threads} --partition={resources.partition}"
default-resources:
  - mem_mb=4000
  - runtime="0-00:30:00"
  - partition=normal
```

## Output explained

The pipeline provides many files and to help navigate these this section will explain what each section provides.
These are all folders within the `./working_dir/output` folder.


