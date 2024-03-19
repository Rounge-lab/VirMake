import yaml
import pathlib
import sys
import logging
import click


# CLI command tool to choose the viral identifier

def make_config():
    """Creates a default config structure."""
    config = {}
    db_path = virmake_path / "databases"

    config["slurm_account"] = "default"
    config["assembler"] = "metaspades"
    config["identifier"] = "virsorter2"
    config["trim_percentage"] = 0.05
    config["min_coverage"] = 75
    config["min_contig_size"] = 1000

    logging.info("Viral identifier chosen is "
                 + config["identifier"]
                 + "vibrant, applying appropriate config settings")

    if config["identifier"] == 'vibrant':
           config["vibrant"] = {
            "is_virome": "no",
            }  
           
           config["path"]["database"] = { "vibrant": str(db_path / "vibrant" / "vibrant-1.2.1") }

    elif config["identifier"] == 'genomad':
           config["path"]["database"] = { "genomad": str(db_path / "genomad") }
    else:
        config["virsorter2"] = {
            "pass1": {
                "min_lenght": 3000,
                "min_score": 0.5,
                "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
            },
            "pass2": {
                "min_lenght": 1000,
                "min_score": 0.5,
                "viral_groups": "dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae",
            },
        }
 

    config["path"] = {
        "virmake": str(virmake_path),
        "envs": str(virmake_path / "envs"),
        "input": str(virmake_path / "working_dir" / "input"),
        "output": str(virmake_path / "working_dir" / "output"),
        "log": str(virmake_path / "working_dir" / "log"),
        "benchmark": str(virmake_path / "working_dir" / "benchmark"),
        "temp": str(virmake_path / "working_dir" / "temp"),
        "scripts": str(virmake_path / "workflow" / "scripts"),
        "database": {
            "DRAM": str(db_path / "DRAM" / "DRAM_data"),
            "checkv": str(db_path / "checkv"),
            "virsorter2": str(db_path / "virsorter2/db"),
            "INPHARED": str(db_path / "INPHARED"),
            "RefSeq": str(db_path / "RefSeq/viral.1.1.genomic.fna"),
            "vcontact2": str(db_path / "vcontact2"),
        },
    }


    config["cd-hit-est"] = {
        "identity_threshold": 0.95,
        "coverage": 0.85,
    }
    config["quality_threshold"] = "medium"
    config["threads"] = 24
    config["job_type"] = {
        "small": "normal",
        "normal": "normal",
        "big": "bigmem",
    }
    config["memory"] = {
        "tiny": 1000,
        "small": 8000,
        "normal": 16000,
        "big": 32000,
        "vcontact2": 63000,
        "metaquast": 63000,
        "metaspades": 63000,
    }
    config["time"] = {
        "tiny": "30 m",
        "small": "1 h",
        "normal": "6 h",
        "big": "13 h",
        "vcontact2": "72 h",
        "metaquast": "24 h",
        "metaspades": "2 h",
    }
    config["include_tables"] = [
        "raw_coverage_table.tsv",
        "trimmed_mean_coverage.tsv",
        "amg_summary.tsv",
        "quality_summary.tsv",
        "quality_summary.tsv",
        "quality_summary.tsv",
        "contig_plasmid_summary",
        "results_vcontact2",
        "compared_samples_comparisonsTable.tsv",
        "transposed_report.tsv",
        "viral_genomes_combined.csv",
        "derep95_combined.fasta.clstr",
    ]

    return config


def main():
    """Saves the config to a file."""
    virmake_path = pathlib.Path(sys.argv[1])
    config_path = virmake_path / "workflow" / "config" / "params.yaml"
    config = make_config(virmake_path)
    if config_path.exists():
        print(f"Config file {config_path} already exists.")
        print("Skipping...")
    else:
        with open(config_path, "w+") as f:
            yaml.dump(config, f)


if __name__ == "__main__":
    main()
