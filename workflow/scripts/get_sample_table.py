import os
import sys
import pandas as pd
import logging
from collections import defaultdict


def get_sample_names(sample_path):
    """
    Function that checks if all samples have paired-end files
    and correct extension. If not it does not add them.
    """
    sample_names = []
    arr = os.listdir(sample_path)
    for files in arr:
        if files.endswith(".fastq.gz"):
            sample_names.append(files)

    sample_names_cleared = []
    for x in sample_names:
        if not (
            (x[:-11] + "R1.fastq.gz") in sample_names
            and (x[:-11] + "R2.fastq.gz") in sample_names
        ):
            print("There is a missing instance of sample: ", x)
            print("This sample will not be added to the workflow list")
        else:
            sample_names_cleared.append(x[:-12])

    sample_names_cleared = sorted(list(set(sample_names_cleared)))
    return sample_names_cleared


def create_sample_table(table):
    """
    Function that creates a dataframe from an input table
    """
    df = pd.DataFrame(table)
    return df


def write_sample_table(sample_path):
    """
    Function that writes the sample table to file.
    """
    names = get_sample_names(sample_path)
    df = create_sample_table(names)
    df.columns = ["sample"]
    df.to_csv("samples.tsv", sep="\t", encoding="utf-8")


##Cite atlas
# Smith, J (2011) Metagenome-Atlas code snipper (Version 2.9.0) [Source code].
# https://github.com/metagenome-atlas/atlas/blob/master/atlas/init/create_sample_table.py
def get_sample_names_test(sample_path, outfile="samples_test.tsv"):
    """
    Function that  gathers the  base names of the sequences within the samples/ folder.
    """
    samples = defaultdict(dict)
    seen = set()
    for dir_name, sub_dirs, files in os.walk(os.path.abspath(sample_path)):
        for fname in files:
            if ".fastq" in fname or ".fq" in fname:
                sample_id = fname.split(".fastq")[0].split(".fq")[0]

                sample_id = (
                    sample_id.replace("_R1", "")
                    .replace("_r1", "")
                    .replace("_R2", "")
                    .replace("_r2", "")
                )
                sample_id = sample_id.replace("_", "-").replace(" ", "-")

                fq_path = os.path.join(dir_name, fname)

                if fq_path in seen:
                    continue

                if "_R2" in fname or "_r2" in fname:
                    if "R2" in samples[sample_id]:
                        logging.error(
                            f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}"
                        )

                    samples[sample_id]["R2"] = fq_path
                else:
                    if "R1" in samples[sample_id]:
                        logging.error(
                            f"Duplicate sample {sample_id} was found after renaming; skipping... \n Samples: \n{samples}"
                        )

                    samples[sample_id]["R1"] = fq_path

    samples = pd.DataFrame(samples).T
    if samples.isna().any().any():
        logging.error(f"Missing files:\n {samples}")

    if os.path.exists(outfile):
        logging.error(
            f"Output file {outfile} already exists I don't dare to overwrite it."
        )
    else:
        samples.to_csv(outfile, sep="\t")
    print(samples)
    return samples


def load_sample_table(sample_table="samples.tsv"):
    sampleTable = pd.read_csv(sample_table, index_col=0, sep="\t")
    return sampleTable
