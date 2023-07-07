# This script is used to setup the environment for VirMake.
# TODO introduce logging module

import os
import pathlib
import shutil
import subprocess
import sys


# internal functions
def strip_stdout(stdout):
    """Strip stdout and return as string."""
    return stdout.decode("utf-8").strip()


# determine absolute path to VirMake
virmake_path = pathlib.Path(__file__).parent.absolute()

print(f"\nVirMake setup started at {virmake_path}\n")

# check if conda is installed
try:
    cmd = "conda --version"
    conda_ver = subprocess.run(cmd.split(), capture_output=True)
except FileNotFoundError:
    exit("Conda not found! Please install conda and try again.\n")
print(f"Using {strip_stdout(conda_ver.stdout)}\n")

# check if conda is initialized
cmd = "conda activate base"
conda_activate_base = subprocess.run(
    cmd.split(), capture_output=True, shell=True
)
if conda_activate_base.returncode != 0:
    exit(strip_stdout(conda_activate_base.stderr))

# check if mamba is installed in base
try:
    cmd = "mamba --version"
    mamba_ver = subprocess.run(cmd.split(), capture_output=True)
except FileNotFoundError:
    # if not install it using conda
    print("Installing mamba...\n")
    cmd = "conda install -y -c conda-forge mamba"
    subprocess.run(cmd.split())
    print("Mamba installed successfully.\n")
    cmd = "mamba --version"
    mamba_ver = subprocess.run(cmd.split(), capture_output=True)
print(f"Using {strip_stdout(mamba_ver.stdout.splitlines()[0])}")

# check if mamba is initialized
cmd = "mamba activate base"
mamba_activate_base = subprocess.run(
    cmd.split(), capture_output=True, shell=True
)
if mamba_activate_base.returncode != 0:
    exit(strip_stdout(mamba_activate_base.stderr))

# create virmake environment
cmd = "mamba env list"
venv_list = subprocess.run(cmd.split(), capture_output=True)
if "virmake" not in strip_stdout(venv_list.stdout):
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
                shutil.rmtree(virmake_path / "databases")
                os.makedirs(virmake_path / "databases", exist_ok=True)
            elif overwrite_db.lower() == "n":
                print("Skipping database setup...")
                break
            else:
                pass
        cmd = (
            "conda run -n virmake --no-capture-output "
            "snakemake --snakefile utils/setup_db.smk --cores 8 "
            f"--config database_dir={virmake_path / 'databases'} "
            f"envs_dir={virmake_path / 'envs'} --use-conda --nolock"
        )
        subprocess.run(cmd.split())

        break
    elif setup_db.lower() == "n":
        print("\nSkipping database setup...\n")
        break
    else:
        pass

# prepare main script
print(f"\nPreparing main script at {virmake_path / 'virmake'}\n")
cmd = "conda run -n virmake which python"
virmake_python_path = subprocess.run(cmd.split(), capture_output=True)
virmake_python_path = strip_stdout(virmake_python_path.stdout)
if (virmake_path / "virmake").exists():
    os.remove(virmake_path / "virmake")
with open(virmake_path / "virmake", "w+") as o:
    o.write(f"#!{virmake_python_path}\n\n")
    with open(virmake_path / "utils" / "virmake.py") as i:
        o.write(i.read())
cmd = "chmod +x virmake"
subprocess.run(cmd.split())
