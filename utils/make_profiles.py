import yaml
import pathlib
import sys


def make_slurm_config():
    """Creates a default slurm profile structure."""
    config = {}
    config["executor"] = "slurm"
    config["default-resources"] = {
        "slurm_account": "",
        "runtime": 120,
        "mem_mb": 20000,
    }
    config["max-status-checks-per-second"] = 0.1
    config["latency-wait"] = 60
    config["set-resources"] = {
        "out_of_memory": {
            "mem_mb": 50,
            "runtime": 1
        },
    }
    return config

def make_workflow_config():
    """Creates a default config structure."""
    config = {}
    # config["snakefile"] = "Snakefile"
    config["show-failed-logs"] = True
    config["rerun-incomplete"] = True
    config["keep-going"] = True
    config["nolock"] = True
    config["printshellcmds"] = True
    config["use-conda"] = True
    config["jobname"] = "{rule}.{jobid}"
    config["max-jobs-per-second"] = 1
    config["jobs"] = 200
    config["conda-frontend"] = "mamba"
    return config


def main():
    """Saves the config to a file."""
    virmake_path = pathlib.Path(sys.argv[1])
    
    workflow_config_path = virmake_path / "config" / "config.yaml"
    workflow_config = make_workflow_config()
    
    slurm_config_path = virmake_path / "config" / "slurm" / "config.yaml"
    slurm_config = make_slurm_config()


    if workflow_config_path.exists():
        print(f"Config file {workflow_config_path} already exists.")
        print("Skipping...")
    else:
        with open(workflow_config_path, "w+") as f:
            yaml.dump(workflow_config, f)
    if slurm_config_path.exists():
        print(f"Config file {slurm_config_path} already exists.")
        print("Skipping...")
    else:
        with open(slurm_config_path, "w+") as f:
            yaml.dump(slurm_config, f)


if __name__ == "__main__":
    main()
