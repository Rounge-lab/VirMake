# This script is used to setup the environment for VirMake.

import logging
import os
import pathlib
import shutil
import subprocess
import sys


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


def install_check_mamba(logger):
    """
    check if mamba is installed and if not, install it,
    then check if mamba is initialized
    """
    # check if mamba is installed in base
    try:
        cmd = "mamba --version"
        mamba_ver = subprocess.run(cmd.split(), capture_output=True)
    except FileNotFoundError:
        # if not install it using conda
        logger.info("Installing mamba...\n")
        cmd = "conda install -y -c conda-forge mamba"
        subprocess.run(cmd.split())
        logger.info("Mamba installed successfully.\n")
        cmd = "mamba --version"
        mamba_ver = subprocess.run(cmd.split(), capture_output=True)
    logger.info(f"Using {strip_stdout(mamba_ver.stdout.splitlines()[0])}")

    # check if mamba is initialized
    cmd = "mamba activate base"
    mamba_activate_base = subprocess.run(
        cmd.split(), capture_output=True, shell=True
    )
    if mamba_activate_base.returncode != 0:
        logger.error(
            "Mamba not initialized! See setup.log for full traceback.\n"
        )
        logger.critical(strip_stdout(mamba_activate_base.stderr))
        exit(1)


def create_venv(logger, virmake_path):
    """create virmake environment"""
    cmd = "mamba env list"
    venv_list = subprocess.run(cmd.split(), capture_output=True)
    if "virmake" not in strip_stdout(venv_list.stdout):
        logger.info("\nPreparing VirMake virtual environment...\n")
        cmd = f"mamba env create -f {virmake_path / 'venv.yaml'}"
        subprocess.run(cmd.split())


def create_config(logger, virmake_path):
    """create config.yaml"""
    logger.info(f"\nCreating configuration file...\n")
    cmd = f"conda run -n virmake python utils/make_config.py {virmake_path}"
    subprocess.run(cmd.split())


def create_working_dir(logger, virmake_path):
    """create working_dir structure"""
    logger.info("Creating working directory structure...")
    os.makedirs(virmake_path / "working_dir", exist_ok=True)
    os.makedirs(virmake_path / "working_dir" / "input", exist_ok=True)


def setup_db(logger, virmake_path):
    """setup databases"""
    if len(sys.argv) == 1:
        logger.info("\n[ DATABASE SETUP ]\n")
        logger.info(
            "VirMake can automatically download and setup databases for you.\n"
            "Since VirMake uses DRAM, minimum 125GB of RAM is required.\n"
            "Around 35GB of free disk space is required to store databases.\n"
            "This is recommended for first time users.\n"
            "If you choose not to setup databases now,\n"
            "you will have to do it manually later!"
        )
        setup_db = input(
            "\nDo you want to setup databases automatically? [Y/n]\n"
        )
        logger.info(
            f"> {setup_db}",
        )
    else:
        if sys.argv[1].lower() == "-y":
            setup_db = "y"
    while True:
        if setup_db.lower() in ["y", ""]:
            logger.info("\nSetting up databases...\n")
            os.makedirs(virmake_path / "databases", exist_ok=True)
            db_files = os.listdir(virmake_path / "databases")
            if db_files != []:
                logger.info("Some files were found in databases directory.")
                overwrite_db = input("Do you wish to overwrite them? [Y/n]\n")
                if overwrite_db.lower() in ["y", ""]:
                    shutil.rmtree(virmake_path / "databases")
                    os.makedirs(virmake_path / "databases", exist_ok=True)
                elif overwrite_db.lower() == "n":
                    logger.info("\nSkipping database setup...\n")
                    break
                else:
                    pass
            cmd = (
                "conda run -n virmake --no-capture-output "
                "snakemake --snakefile utils/setup_db.smk --cores 8 "
                f"--config database_dir={virmake_path / 'databases'} "
                f"envs_dir={virmake_path / 'envs'} --use-conda --nolock"
            )
            db_workflow = subprocess.run(cmd.split(), capture_output=True)
            if db_workflow.returncode != 0:
                logger.error(
                    "Snakemake error! See setup.log for full traceback.\n"
                )
                logger.critical(strip_stdout(db_workflow.stderr))
                exit(1)

            break
        elif setup_db.lower() == "n":
            logger.info("\nSkipping database setup...\n")
            break
        else:
            pass


def prep_script(logger, virmake_path):
    """prepare main script"""
    logger.info(f"Preparing main script at {virmake_path / 'virmake'}\n")
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


def main():
    """Run setup"""
    logger = set_logger()

    # determine absolute path to VirMake
    virmake_path = pathlib.Path(__file__).parent.absolute()
    logger.info(f"\nVirMake setup started at {virmake_path}\n")

    # run all steps
    check_conda(logger)
    install_check_mamba(logger)
    create_venv(logger, virmake_path)
    create_config(logger, virmake_path)
    create_working_dir(logger, virmake_path)
    setup_db(logger, virmake_path)
    prep_script(logger, virmake_path)

    # remove setup.log on success
    os.remove(f"{virmake_path / 'setup.log'}")
    logger.info("Success!\n")


if __name__ == "__main__":
    main()
