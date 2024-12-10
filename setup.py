# This script is used to setup the environment for VirMake.

import logging
import os
import pathlib
import subprocess
import argparse


def strip_stdout(stdout):
    """strip stdout and return as string"""
    return stdout.decode("utf-8").strip()


def set_logger():
    """logger config"""
    logger = logging.getLogger("setup_logger")
    logger.setLevel(logging.INFO)
    f_handler = logging.FileHandler("setup.log")
    f_handler.setLevel(logging.DEBUG)
    f_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s\n")
    )
    c_handler = logging.StreamHandler()
    c_handler.setLevel(logging.INFO)
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)
    return logger


def check_conda(logger):
    """check if conda is installed and initialized"""
    # check if conda is installed
    try:
        cmd = "conda --version"
        conda_ver = subprocess.run(cmd.split(), capture_output=True)
    except FileNotFoundError:
        logger.error(
            "Conda not found! Please install conda and try again.\n"
            "See setup.log for full traceback.\n"
        )
        logger.critical(strip_stdout(conda_ver.stderr))
        exit(1)
    logger.info(f"Using {strip_stdout(conda_ver.stdout)}\n")

    # check if conda is initialized
    cmd = "conda activate base"
    conda_activate_base = subprocess.run(
        cmd.split(), capture_output=True, shell=True
    )
    if conda_activate_base.returncode != 0:
        logger.error(
            "Conda not initialized! See setup.log for full traceback.\n"
        )
        logger.critical(strip_stdout(conda_activate_base.stderr))
        exit(1)

# Change os that virmake env is created, activated, and a venv built based off that, using venv
def create_venv(logger, virmake_path):
    """create virmake environment"""
    logger.info("\nPreparing VirMake conda env...\n")
    logger.info("\nAttempting to use mamba...\n")
    try:
        cmd = (
            f"mamba env create -f {virmake_path / 'workflow' / 'envs' / 'virmake.yaml'} -p ./venv"
        )
        subprocess.run(cmd.split())
    except FileNotFoundError:
        logger.error(
            "Mamba was not available in base environment, using conda.\n"
        )
        cmd = (
            f"conda env create -f {virmake_path / 'workflow' / 'envs' / 'virmake.yaml'} -p ./venv"
        )
        subprocess.run(cmd.split())


def create_virmake_config(logger, virmake_path, work_dir, reads, contigs):
    """create params.yaml"""
    logger.info(f"\nCreating configuration file...\n")
    cmd = f"mkdir {virmake_path}/config"
    subprocess.run(cmd.split())
    cmd = f"conda run -p venv/ python utils/make_virmake_config.py {virmake_path} -w {work_dir}"
    if not reads == "":
        cmd += f" -r {reads}"
    if not contigs == "":
        cmd += f" -c {contigs}"
    subprocess.run(cmd.split())


def create_slurm_profile(logger, virmake_path):
    """create config.yaml"""
    logger.info(f"\nCreating SLURM profile...\n")
    cmd = (
        f"conda run -p venv/ python utils/make_slurm_profile.py {virmake_path}"
    )
    subprocess.run(cmd.split())


def create_dirs(logger, virmake_path, results_dir):
    """create results directory"""
    logger.info("Creating results directory...")
    os.makedirs(virmake_path / results_dir, exist_ok=True)
    os.makedirs(virmake_path / "resources" / "input", exist_ok=True)


def prep_sample_table(logger, virmake_path, work_dir, reads, qc_reads, contigs):
    """create samples.tsv"""
    logger.info(f"\nCreating samples table...\n")
    cmd = (
        f"conda run -p venv/ python utils/create_samples_file.py "
        f"{virmake_path} -w {work_dir}"
        )
    if not reads == "":
        cmd += f" -r {reads}"
    if not qc_reads == "":
        cmd += f" -q {qc_reads}"
    if not contigs == "":
        cmd += f" -c {contigs}"
    subprocess.run(cmd.split())


def prep_script(logger, virmake_path):
    """prepare main script"""
    logger.info(f"\nPreparing main script at {virmake_path / 'virmake'}\n")
    cmd = "conda run -p venv/ which python"
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


def main():

    parser = argparse.ArgumentParser(description='Prepare for running VirMake by setting input paths and results dir.')
    parser.add_argument('--work-dir', type=str, default="results", help="Path to results directory.")
    parser.add_argument('--reads', type=str, default="", help="Path where reads can be found.")
    parser.add_argument('--qc-reads', type=str, default="", help="Path where QC reads can be found.")
    parser.add_argument('--contigs', type=str, default="", help="Path where contigs can be found.")

    args = parser.parse_args()

    """Run setup"""
    logger = set_logger()

    # determine absolute path to VirMake
    virmake_path = pathlib.Path(__file__).parent.absolute()
    logger.info(f"\nVirMake setup started at {virmake_path}\n")

    # run all steps
    check_conda(logger)
    create_venv(logger, virmake_path)
    create_virmake_config(logger, virmake_path, args.work_dir, args.reads, args.contigs)
    create_slurm_profile(logger, virmake_path)
    create_dirs(logger, virmake_path, args.work_dir)
    prep_sample_table(logger, virmake_path, args.work_dir, args.reads, args.qc_reads, args.contigs)
    prep_script(logger, virmake_path)

    # remove setup.log on success
    logger.info("Success!\n")


if __name__ == "__main__":
    main()
