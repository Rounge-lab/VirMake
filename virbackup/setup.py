# This script is used to setup the environment for VirMake.

import logging
import os
import pathlib
import subprocess
import sys


# internal functions
def strip_stout(stout):
    """Strip stout and return as string."""
    return stout.decode("utf-8").strip()


# determine absolute path to VirMake
virmake_path = pathlib.Path(__file__).parent.absolute()

print(f"\nVirMake setup started at {virmake_path}\n")

# check if conda is installed
cmd = "conda --version"
try:
    conda_ver = subprocess.run(cmd.split(), capture_output=True)
except FileNotFoundError:
    exit("Conda not found! Please install conda and try again.\n")
print(f"Using {strip_stout(conda_ver.stdout)}\n")

# check if conda is initialized
cmd = "conda activate"
conda_activate_base = subprocess.run(
    cmd.split(), capture_output=True, shell=True
)
if conda_activate_base.returncode != 0:
    exit(strip_stout(conda_activate_base.stderr))

# check if mamba is installed
cmd = "mamba --version"
try:
    mamba_ver = subprocess.run(cmd.split(), capture_output=True)
except FileNotFoundError:
    # if not install it using conda
    print("\nInstalling mamba...\n")
    cmd = "conda install -y -c conda-forge mamba"
    subprocess.run(cmd.split())
    print("\nMamba installed successfully...\n")
print(f"Using {strip_stout(mamba_ver.stdout.splitlines()[0])}")

# check if mamba is initialized
cmd = "mamba activate"
mamba_activate_base = subprocess.run(
    cmd.split(), capture_output=True, shell=True
)
if mamba_activate_base.returncode != 0:
    exit(strip_stout(mamba_activate_base.stderr))

# create virmake environment
cmd = "mamba env list"
venv_list = subprocess.run(cmd.split(), capture_output=True)
if "virmake" not in strip_stout(venv_list.stdout):
    print("\nPreparing VirMake virtual environment...\n")
    cmd = f"mamba env create -f {virmake_path / 'venv.yaml'}"
    subprocess.run(cmd.split())

# create config.yaml
print(f"\nCreating configuration file...\n")
cmd = f"conda run -n virmake python utils/make_config.py {virmake_path}"
subprocess.run(cmd.split())

# create worktree structure
print("Creating working directory structure...")
os.makedirs(virmake_path / "working_dir", exist_ok=True)
os.makedirs(virmake_path / "working_dir" / "input", exist_ok=True)
os.makedirs(virmake_path / "working_dir" / "temp", exist_ok=True)
os.makedirs(virmake_path / "working_dir" / "output", exist_ok=True)

# setup databases
if len(sys.argv) == 1:
    print("\n[ DATABASE SETUP ]\n")
    print(
        "VirMake can automatically download and setup databases for you.",
        "Since VirMake uses DRAM, minimum 125GB of RAM is required.",
        "Around 35GB of free disk space is required to store databases.",
        "This is recommended for first time users.",
        "If you choose not to setup databases now, ",
        "you will have to do it manually later.",
        sep="\n",
    )
    setup_db = input("\nDo you want to setup databases automatically? [Y/n]\n")
else:
    if sys.argv[1].lower() == "-y":
        setup_db = "y"
while True:
    if setup_db.lower() in ["y", ""]:
        print("\nSetting up databases...\n")
        os.makedirs(virmake_path / "databases", exist_ok=True)
        db_files = os.listdir(virmake_path / "databases")
        if db_files != []:
            print("Some files were found in databases directory.")
            overwrite_db = input("Do you wish to overwrite them? [Y/n]\n")
            if overwrite_db.lower() in ["y", ""]:
                for f in db_files:
                    print("Overwriting...")
                    os.remove(virmake_path / "databases" / f)
            elif overwrite_db.lower() == "n":
                print("Skipping database setup...")
                break
            else:
                pass
        # cmd = f"conda run -n virmake python utils/setup_db.py {virmake_path}"
        # subprocess.run(cmd.split())
        cmd = (
            "conda run -n virmake --no-capture-output "
            "snakemake --snakefile utils/setup_db.smk --cores 1 "
            "--rerun-incomplete --conda-frontend mamba "
            f"--config database_dir={virmake_path / 'databases'} "
            f"envs_dir={virmake_path / 'envs'} "
        )
        subprocess.run(cmd.split())

        break
    elif setup_db.lower() == "n":
        print("\nSkipping database setup...\n")
        break
    else:
        pass
