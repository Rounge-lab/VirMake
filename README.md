# PENDING-Workflow
Creating a snakemake workflow for metagenomic viral identification

# Usage guide:

  
## Installation

Git clone the project/download the workflow from GitHub. 
Make sure you have conda or miniconda installed, preferable with mamba installed.
Create an environment from the enviorment_pending.yaml.
Activate it.
Then make sure you are within the directory of PENDING-WORKFLOW and run:

    python pending.py init
        Options:
      -d, --db-dir PATH       location to store databases
      -w, --working-dir PATH  location for running the application
      --threads INTEGER       number of threads to use per multi-threaded job
      -h, --help              Show this message and exit.

This will create a basic config file and a sample table which are both needed to run the workflow.
  
The last required step is to download all needed databases and external files. This can be done by:

    python pending.py download [OPTIONS]
    Options:
      -d, --db-dir PATH  location to store databases  [required]
      --threads INTEGER  number of threads to use per multi-threaded job
      -n, --dryrun       Test execution.
      -h, --help         Show this message and exit.

The files generated are the `config.yaml` file which contains all locations and setting parameters for several tools in the workflow. And the `sample.tsv` file which gathers all names of the samples within your designated sample folder to be used by the workflow. The tool currently supports samples within one folder that have the format `samplenName_R1.fastq.gz` where sampleName can be whatever you like but the patterns recognized on the ending of files is `_R1.fastq.gz` or `_R2.fastq.gz`.

For a smooth performance and setup for eventual offline running use:

    python pending.py prep-offline [OPTIONS]
        Options:
      --threads INTEGER  number of threads to use per multi-threaded job
      -h, --help         Show this message and exit.

This will generate all the needed environments and images to run the workflow. It can take a while to generate all environments, but it only needs to be done once. This can be beneficial to do if you are running on job clusters that limit internet connections too.

For offline usage simply copy the whole directory over to the machine/server you want to use the workflow on.

The only thing to watch out for is to edit the config file to have the correct paths to other files that were transferred.

This can easily be done by opening the config file and use search and replace on all earlier paths such as:

Original path on machine with internet:

    c:/home/something/databases/....

To new machine:

    /cluster/project/group/user/databases/...

##  Running the workflow
To run the workflow after the setup steps, simply use the command:

    python pending.py [OPTIONS] {qc|assembly|identification|taxonomy|all|None}
       Options:
      --profile TEXT          snakemake profile e.g. for cluster execution.
      -w, --working-dir PATH  location to run pending.
      -c, --config-file PATH  config-file generated with 'pending init'
      -n, --dryrun            Test execution.
      --threads INTEGER       Number of threads used on multithreaded jobs
      -h, --help              Show this message and exit.

























# Developer stuff 
link to refseq:
https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
## Setup for running the snakemake.
Use this on the cluster before submitting jobs:

snakemake --use-conda --conda-create-envs-only -c3

It creates the local envs that snakemake will use.
To be sure no files are locked use:

snakemake --use-conda --conda-create-envs-only -c3 --unlock

try this for running it:
snakemake -s workflow/Snakefile --use-conda --conda-create-envs-only -c4 --profile profile/

snakemake --profile profile/ --use-conda --use-singularity --conda-create-envs-only -c5
snakemake --profile profile/ --use-conda --use-singularity -c5

snakemake --profile profile/ -n

DRAMV setup:

DRAM-setup.py prepare_databases --output_dir DRAM_data --verbose --skip_uniref --threads 15 

snakemake -s workflow/Snakefile --use-conda --use-singularity -c8

snakemake -s workflow/Snakefile --dag

Weird bug with PATH when using python3 and simplify script in Inphared setup

code tried:
    echo "before"
    export PYTHONPATH="{params.python_version}"
    echo $PYTHONPATH
    echo "After"
    which python
    $PYTHONPATH 
    which python
    echo $PATH


