import shutil
import pathlib


with open(str(snakemake.input.html_files), "r") as f:
    paths = f.readlines()

for p in paths:
    ori = p.replace("\n", "")
    dst = (
        str(snakemake.params.temp)
        + "/complete_output/"
        + ori.split("output/")[1]
    )
    dst_dir = "/".join(dst.split("/")[:-1])
    pathlib.Path(dst_dir).mkdir(parents=True, exist_ok=True)
    shutil.copy(ori, dst)

complete_output = str(snakemake.params.temp) + "/complete_output"
shutil.copytree(
    str(snakemake.input.tables), pathlib.Path(complete_output) / "tables"
)
archive_name = str(snakemake.output).split(".")[0]
shutil.make_archive(archive_name, "zip", complete_output)
