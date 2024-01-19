import shutil
import pathlib


with open(str(snakemake.input.reports), "r") as f:
    report_paths = f.readlines()

with open(str(snakemake.input.tables), "r") as f:
    table_paths = f.readlines()

for p in report_paths:
    ori = p.replace("\n", "")
    dst = (
        str(snakemake.params.temp)
        + "/complete_output/reports/"
        + ori.split("output/")[1]
    )
    dst_dir = "/".join(dst.split("/")[:-1])
    pathlib.Path(dst_dir).mkdir(parents=True, exist_ok=True)
    shutil.copy(ori, dst)

for p in table_paths:
    ori = p.replace("\n", "")
    dst = (
        str(snakemake.params.temp)
        + "/complete_output/tables/"
        + ori.split("output/")[1]
    )
    dst_dir = "/".join(dst.split("/")[:-1])
    pathlib.Path(dst_dir).mkdir(parents=True, exist_ok=True)
    shutil.copy(ori, dst)

complete_output = str(snakemake.params.temp) + "/complete_output"
archive_name = str(snakemake.output).split(".")[0]
shutil.make_archive(archive_name, "zip", complete_output)
