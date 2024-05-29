import yaml
import pathlib
import sys


def make_config():
    """Creates a default config structure."""
    config = {}
    config["snakefile"] = "Snakefile"
    config["latency-wait"] = 60
    config["reason"] = True
    config["show-failed-logs"] = True
    config["rerun-incomplete"] = True
    config["keep-going"] = True
    config["nolock"] = True
    config["printshellcmds"] = True
    config["use-conda"] = True
    config["jobname"] = "{rule}.{jobid}"
    config["max-jobs-per-second"] = 1
    config["max-status-checks-per-second"] = 10
    config["jobs"] = 200
    config[
        "cluster"
    ] = 'sbatch --output=\\"jobs/{rule}/slurm_%x_%j.out\\" --error=\\"jobs/{rule}/slurm_%x_%j.log\\" --mem={resources.mem_mb} --time={resources.runtime}'
    config["conda-frontend"] = "mamba"
    config["set-resources"] = [
        "out_of_memory:mem_mb=50",
        "out_of_memory:runtime=00:01:00",
    ]
    config["default-resources"] = [
        "mem_mb=20000",
        "runtime=02:00:00",
    ]
    config["set-threads"] = ["out_of_memory=1"]
    return config


def main():
    """Saves the config to a file."""
    virmake_path = pathlib.Path(sys.argv[1])
    config_path = virmake_path / "config" / "config.yaml"
    config = make_config()
    if config_path.exists():
        print(f"Config file {config_path} already exists.")
        print("Skipping...")
    else:
        with open(config_path, "w+") as f:
            yaml.dump(config, f)


if __name__ == "__main__":
    main()
