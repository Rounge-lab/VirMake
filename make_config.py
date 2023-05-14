import multiprocessing
import os
import sys
import tempfile
from snakemake.io import load_configfile
from snakemake.utils import update_config as snakemake_update_config
import yaml 

def make_default_config():
    config = {}
    config["Database_dir"] = {
        "main_dir": "Databases",
        "DRAM": "DRAM",
        "CheckV": "checkv/checkv-db-v1.5",
        "VirSorter2": "VirSorter2/db",
        "Inphared": "INPHARED",
        "RefSeq": "RefSeq/viral.1.1.genomic.fna"
    }
    config["sample_table"] = os.getcwd()+"/samples.tsv"
    config["Workflow_dirs"] = {
        "Profile_dir": "profile/",
        "Result_dir": "results/",
        "Log_dir": "logs/",
        "Samples_dir": "samples/",
        "working_dir": os.getcwd()
    }
    config["Virsorter2_settings"] = {
        "pass1": {
            "min_lenght": 3000,
            "min_score": 0.5,
        },
        "pass2": {
            "min_lenght": 1000,
            "min_score": 0.5,
        }
    }
    config["Vibrant_settings"] = {
        "is_virome": "no",
        "cutoff_length": 1000
    }
    config["Cd-hit-est_criteria"] = {
        "ANI": 0.95,
        "coverage": 0.85
    }
    config["CheckV_threshold"] = "Medium"
    config["Threads"] = 24
    config["job_type"] = {
        "small": "normal",
        "normal": "normal",
        "big": "bigmem"
    }
    config["tiny_mem"] = 1000
    config["small_mem"] = 8000
    config["normal_mem"] = 16000
    config["big_mem"] = 32000
    config["vcontact2_mem"] = 63000
    config["metaquast_mem"] = 63000

    config["tiny_time"] = "0-00:30:00"
    config["small_time"] = "0-01:00:00"
    config["normal_time"] = "0-06:00:00"
    config["big_time"] = "0-13:00:00"
    config["vcontact2_time"] = "0-24:00:00"
    config["metaquast_time"] = "0-24:00:00"

    return config

def update_config(config):
    """
    Creates a default config and updates the values based on input config
    """

    default_config = make_default_config()
    snakemake_update_config(default_config, config)

    return default_config


def make_config(
    database_dir,
    threads=8,
    config="config.yaml"
):
    new_config = make_default_config()

    new_config["Threads"] = threads
    new_config["Database_dir"]["main_dir"] = database_dir
    new_config["Database_dir"]["DRAM"] = database_dir +"/"+ new_config["Database_dir"]["DRAM"]
    new_config["Database_dir"]["CheckV"] = database_dir +"/"+ new_config["Database_dir"]["CheckV"]
    new_config["Database_dir"]["VirSorter2"] = database_dir +"/"+ new_config["Database_dir"]["VirSorter2"]
    new_config["Database_dir"]["Inphared"] = database_dir +"/"+ new_config["Database_dir"]["Inphared"]
    new_config["Database_dir"]["RefSeq"] = database_dir +"/"+ new_config["Database_dir"]["RefSeq"]

    print(new_config)
    if os.path.exists(config):
            print(f"Config file {config} already exists")
    else:
        with open(config, "w") as f:
            yaml.dump(new_config, f)