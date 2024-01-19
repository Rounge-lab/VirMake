import pandas as pd
from pathlib import Path
from os import walk
import shutil

INPUT = Path(str(snakemake.input)).absolute()
OUTPUT = INPUT / "_summary"


def get_names_paths():
    names = []
    paths = []
    for prefix, dirname, suffixes in walk(INPUT):
        for suffix in suffixes:
            if suffix == "total_time.txt":
                continue
            path = Path(prefix) / suffix
            paths.append(path)
            names.append(
                str(path).replace(f"{str(INPUT)}/", "").replace(".txt", "")
            )
    return names, paths


def merge_benchmarks():
    merged = df = pd.DataFrame(
        columns=[
            "s",
            "h:m:s",
            "max_rss",
            "max_vms",
            "max_uss",
            "max_pss",
            "io_in",
            "io_out",
            "mean_load",
            "cpu_time",
            "rule_name",
        ]
    )
    for name, path in zip(*get_names_paths()):
        df = pd.read_csv(path, sep="\t")
        df["rule_name"] = name
        merged = pd.merge(merged, df, how="outer")
    return merged


def create_report(merged):
    merged_d = merged.drop(columns=["rule_name", "h:m:s"])
    with open(OUTPUT / "REPORT.txt", "w") as f:
        f.write(f"BENCHMARK REPORT\n\n")
        f.write(f"Summed:\n{str(merged_d.sum())}\n\n")
        f.write(f"Mean:\n{str(merged_d.mean())}\n\n")
        f.write(f"Min:\n{str(merged_d.min())}\n\n")
        f.write(f"Max:\n{str(merged_d.max())}\n\n")
        f.write(
            "Refer to: https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake"
        )


if __name__ == "__main__":
    if not OUTPUT.exists():
        OUTPUT.mkdir(parents=True)
    else:
        shutil.rmtree(OUTPUT)
        OUTPUT.mkdir(parents=True)
    merged = merge_benchmarks()
    merged.to_csv(OUTPUT / "merged_benchmarks.csv", index=False)
    create_report(merged)
