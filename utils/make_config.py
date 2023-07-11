import yaml
import pathlib
import sys


def make_config(virmake_path):
    """Creates a default config structure."""
    config = {}
    db_path = virmake_path / "databases"

    config["assembler"] = "metaSpades"
    config["path"] = {
        "virmake": str(virmake_path),
        "envs": str(virmake_path / "envs"),
        "input": str(virmake_path / "working_dir" / "input"),
        "output": str(virmake_path / "working_dir" / "output"),
        "log": str(virmake_path / "working_dir" / "log"),
        "benchmark": str(virmake_path / "working_dir" / "benchmark"),
        "temp": str(virmake_path / "working_dir" / "temp"),
        "database": {
            "parent_dir": str(db_path),
            "DRAM": str(db_path / "DRAM"),
            "checkv": str(db_path / "checkv/checkv-db-v1.5"),
            "virsorter2": str(db_path / "virsorter2/db"),
            "INPHARED": str(db_path / "INPHARED"),
            "RefSeq": str(db_path / "RefSeq/viral.1.1.genomic.fna"),
        },
    }
    config["virsorter2"] = {
        "pass1": {
            "min_lenght": 3000,
            "min_score": 0.5,
        },
        "pass2": {
            "min_lenght": 1000,
            "min_score": 0.5,
        },
    }
    config["vibrant"] = {"is_virome": "no", "cutoff_length": 1000}
    config["cd-hit-est"] = {
        "identity_threshold": 0.95,
        "coverage": 0.85,
    }
    config["checkv_threshold"] = "Medium"
    config["threads"] = 24
    config["job_type"] = {
        "small": "normal",
        "normal": "normal",
        "big": "bigmem",
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
