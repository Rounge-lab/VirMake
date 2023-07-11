#!/usr/bin/env python3
import os, sys

"""
Script that downloads the example files,
same as those used within the Thesis Paper.
"""


def download_example():
    """
    Function Handles all samples and addresses to be downloaded
    """
    urlDict = {
        "SRR8653245": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/005/SRR8653245/SRR8653245_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/005/SRR8653245/SRR8653245_2.fastq.gz",
        },
        "SRR8653218": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/008/SRR8653218/SRR8653218_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/008/SRR8653218/SRR8653218_2.fastq.gz",
        },
        "SRR8653221": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/001/SRR8653221/SRR8653221_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/001/SRR8653221/SRR8653221_2.fastq.gz",
        },
        "SRR8653248": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/008/SRR8653248/SRR8653248_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/008/SRR8653248/SRR8653248_2.fastq.gz",
        },
        "SRR8653247": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/007/SRR8653247/SRR8653247_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/007/SRR8653247/SRR8653247_2.fastq.gz",
        },
        "SRR8653084": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/004/SRR8653084/SRR8653084_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/004/SRR8653084/SRR8653084_2.fastq.gz",
        },
        "SRR8652914": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/004/SRR8652914/SRR8652914_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/004/SRR8652914/SRR8652914_2.fastq.gz",
        },
        "SRR8652969": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/009/SRR8652969/SRR8652969_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/009/SRR8652969/SRR8652969_2.fastq.gz",
        },
        "SRR8652861": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/001/SRR8652861/SRR8652861_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/001/SRR8652861/SRR8652861_2.fastq.gz",
        },
        "SRR8653090": {
            "R1": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/000/SRR8653090/SRR8653090_1.fastq.gz",
            "R2": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR865/000/SRR8653090/SRR8653090_2.fastq.gz",
        },
    }
    for key in urlDict:
        download(key, urlDict[key]["R1"], "R1")
        download(key, urlDict[key]["R2"], "R2")


def download(key, URL, end):
    """
    Function downloads and saves the requested sequence
    in current directory
    """
    destination = "./" + key + "_" + end + ".fatq.gz"
    cmd = "wget " + "--quiet " + "-O " + destination + " " + URL
    try:
        os.system(cmd)
        print("Finished downloading: " + key + "_" + end + ".fatq.gz")
    except Exception as e:
        print("Error downloading file: %s" % e)
        sys.exit(1)


download_example()
