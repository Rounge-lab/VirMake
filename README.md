# VirMake: a Snakemake pipeline for viral metagenomic data analysis

## Contents

1. [About VirMake](#about-virmake)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Output explained](#output-explained)

## About VirMake

VirMake is a Snakemake based pipeline that offers viral metagenomic data analysis on paired-end data. It offers taxonomic and functional annotation, supports offline running and support for HPC cluster execution. It is made for Linux based systems and has been tested on x86_64-based Linux and with SLURM cluster execution.

<!-- ![Flowchart](img/flowchart.png) -->

## Installation

### Prerequisites:
•  Git
•  Conda (either Miniconda or Anaconda) installation with write permissions. Ensure that Conda's base environment is activated.
•  System Requirements: At least 125 GB of RAM and 180 GB of free disk space.


### To install VirMake follow these steps:

1. Clone the repository using `git clone https://github.com/Rounge-lab/VirMake.git`
2. Run the VirMake setup script `python setup.py`.
   This script will set up a conda environment, config files and a table specifying the location of input files.
    * Use the `--working-dir` option to specify the output directory for configuration files (default is `./working_dir/`).
    * Use the `--reads`, `--reads-qc`, and `--contigs` options to specify the location of any input files.


## Usage

1. Activate the `virmake` conda environment: you do anything:
```
conda activate ./venv
```
Ensure all commands are run from the root directory of the repository.

2.Access help for the `virmake`:

```
./virmake -h
```

3.Run a dry-run of the pipeline:

```
./virmake run -n
```

This will output the steps required for execution.

### Environments preparation

VirMake uses conda environments to run separate rules. These environments will be set up on the first run and will be stored in `./workflow/.snakemake/` folder for subsequent runs. This takes a substantial amount of time on the first run. If you want to set up the environments before running the workflow use:

```
./virmake prep
```

After that you can run the workflow offline.

### Getting samples

To run the workflow you will need to provide input files. The input files should be placed in `./resources/input/` folder. These need to be in `.fastq.gz` format. The input files should be named in the following format: `<sample_name>_1.fastq.gz` and `<sample_name>_2.fastq.gz`. VirMake only works for paired-end reads. You need at least two samples to be able to complete the pipeline with a comparison of the samples.

You can also download samples from SRA database by using:

```
./virmake get SRA <accession_number>
```

This command will download the samples from SRA database, place them in `VirMake/resources/input/` folder, `gzip` them and rename them accordingly.

You can get example files from [1] by running:

```
./virmake get SRA PRJNA524703
```

[1] Liang, G., Zhao, C., Zhang, H. et al. The stepwise assembly of the neonatal virome is modulated by breastfeeding. Nature 581, 470–474 (2020). [https://doi.org/10.1038/s41586-020-2192-1](https://doi.org/10.1038/s41586-020-2192-1)

### Running the workflow

**Before running the workflow, add all your samples to the input folder in the resources folder.**

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

To adjust the workflow settings edit the `params.yaml` file. The file is located in `VirMake/config/params.yaml`.
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
    benchmark: [/.../VirMake/results/benchmark]
    database:
        DRAM: [/.../VirMake/resources/databases/DRAM]
        INPHARED: [/.../VirMake/resources/databases/INPHARED]
        RefSeq: [/.../VirMake/resources/databases/RefSeq]
        checkv: [/.../VirMake/resources/databases/checkv]
        vcontact2: [/.../VirMake/resources/databases/vcontact2]
        vibrant: [/.../VirMake/resources/databases/vibrant]
        virsorter2: [/.../VirMake/resources/databases/virsorter2]
    envs: [/.../VirMake/workflow/envs]
    input: [/.../VirMake/resources/input]
    log: [/.../VirMake/results/log]
    output: [/.../VirMake/results/output]
    profile: []
    scripts: [/.../VirMake/workflow/scripts]
    temp: [/.../VirMake/results/temp]
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

You can adjust the HPC profile file to suit your needs. The file is located in `./config/config.yaml`

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
These are all folders within the `./results/output` folder.


