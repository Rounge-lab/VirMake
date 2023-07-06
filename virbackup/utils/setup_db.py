import subprocess
import logging
import pathlib
import sys

virmake_path = pathlib.Path(sys.argv[1])

cmd = (
    "conda run -n virmake "
    "snakemake --snakefile utils/setup_db.smk --cores 1 "
    "--rerun-incomplete --conda-frontend mamba "
    f"--config database_dir={virmake_path / 'databases'} "
    f"envs_dir={virmake_path / 'envs'} "
)
print("Executing: " + cmd)
try:
    subprocess.run(cmd.split(), shell=True)
except subprocess.CalledProcessError as e:
    # removes the traceback
    logging.critical(e)
    exit(1)
