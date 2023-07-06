#!/usr/bin/env python

import numpy as np
import re
import gzip

"""
Script to calculate the trimmed mean of contig coverage
for bbmap output file base_coverage.
"""


def get_contig_lengths(contig_summary_file):
    with open(contig_summary_file) as csf:
        lines = csf.readlines()
        cont_id = [x.replace("\n", "").split("\t")[0] for x in lines][1:]
        cont_len = [x.replace("\n", "").split("\t")[2] for x in lines][1:]

        return cont_id, cont_len


def calc_trimmed_mean(contig_basecov, cont_len, trim_perc):
    if len(contig_basecov) > 0:
        cov = [x[1] for x in contig_basecov]
        zeros = int(cont_len) - len(cov)
        cov = [0] * zeros + cov
        trimmed_cov = np.sum(sorted(cov)[: len(cov) - int(len(cov) * trim_perc)]) / int(
            len(cov) * (1 - 2 * trim_perc)
        )
    else:
        trimmed_cov = 0
    return trimmed_cov


def get_trimmed_mean(filename, contig_summary_file, trim_perc=0.05):
    contig_basecov = []
    contig_id = ""

    contig_cov = []

    count = 0

    cont_ids, cont_lengths = get_contig_lengths(contig_summary_file)

    with gzip.open(filename, "rt") as fi:
        while True:
            count += 1
            line = fi.readline().replace("\n", "")

            if not line:
                cont_len = cont_lengths[cont_ids.index(contig_id)]
                trimmed_cov = calc_trimmed_mean(contig_basecov, cont_len, trim_perc)

                contig_cov.append([contig_id, trimmed_cov])
                break

            if re.search("#vOTU_", line):
                if contig_id == "":
                    contig_id = line.replace("#", "")
                else:
                    cont_len = cont_lengths[cont_ids.index(contig_id)]
                    trimmed_cov = calc_trimmed_mean(contig_basecov, cont_len, trim_perc)

                    contig_cov.append([contig_id, trimmed_cov])

                    contig_id = line.replace("#", "")
                    contig_basecov = []
            else:
                contig_basecov.append(
                    [int(x) for x in line.replace("\n", "").split("\t")]
                )
    return contig_cov


def write_coverage_to_file(contig_coverage, filename):
    file = open(filename, "w")
    file.write("contig\tcoverage\n")
    for cont in contig_coverage:
        file.write(cont[0] + "\t" + str(cont[1]) + "\n")

    file.close()


filename = snakemake.input["basecov"]
outfile = snakemake.output["trimmed_mean"]
contig_summary_file = snakemake.input["coverage_stats"]
trim_perc = snakemake.params["trim_perc"]

trimmed_contig_coverage = get_trimmed_mean(
    filename=filename, contig_summary_file=contig_summary_file, trim_perc=trim_perc
)

write_coverage_to_file(trimmed_contig_coverage, outfile)
