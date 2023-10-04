import yaml
import pathlib
import sys


def make_config(virmake_path):
    """Creates a default config structure."""
    config = {}
    db_path = virmake_path / "databases"

    config["assembler"] = "metaspades"
    config["trim_percentage"] = 0.05
    config["min_coverage"] = 75
    config["min_contig_size"] = 1000
    config["path"] = {
        "virmake": str(virmake_path),
        "envs": str(virmake_path / "envs"),
        "input": str(virmake_path / "working_dir" / "input"),
        "output": str(virmake_path / "working_dir" / "output"),
        "log": str(virmake_path / "working_dir" / "log"),
        "benchmark": str(virmake_path / "working_dir" / "benchmark"),
        "temp": str(virmake_path / "working_dir" / "temp"),
        "scripts": str(virmake_path / "workflow" / "scripts"),
        "profile": "",
        "database": {
            "DRAM": str(db_path / "DRAM" / "DRAM_data"),
            "checkv": str(db_path / "checkv"),
            "virsorter2": str(db_path / "virsorter2/db"),
            "INPHARED": str(db_path / "INPHARED"),
            "RefSeq": str(db_path / "RefSeq/viral.1.1.genomic.fna"),
            "vcontact2": str(db_path / "vcontact2"),
            "vibrant": str(db_path / "vibrant" / "vibrant-1.2.1"),
        },
    }
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
    config["vibrant"] = {
        "is_virome": "no",
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
        "megahit": 0.8,
    }
    config["time"] = {
        "tiny": "30 m",
        "small": "1 h",
        "normal": "6 h",
        "big": "13 h",
        "vcontact2": "24 h",
        "metaquast": "24 h",
        "megahit": "24 h",
    }

    return config


def main():
    """Saves the config to a file."""
    virmake_path = pathlib.Path(sys.argv[1])
    config_path = virmake_path / "workflow" / "config.yaml"
    config = make_config(virmake_path)
    if config_path.exists():
        print(f"Config file {config_path} already exists.")
        print("Skipping...")
    else:
        with open(config_path, "w+") as f:
            yaml.dump(config, f)


if __name__ == "__main__":
    main()
