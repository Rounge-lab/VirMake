#!/usr/bin/env python3
"""
script to calculate the aggregated results, statistics and simple plots for the
pipeline
"""
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import re


samples = pd.read_table("./samples.tsv", sep="\t", usecols=["sample"])
out_dir = "results/statistics/plots/"
st_dir = "results/statistics/"


def create_sample_stats_virsorter2():
    """
    Creates a table containig necessary stats from virsorter2 results,
    sample stats and checkv
    """
    column_names = [
        "Sample",
        "Raw_reads",
        "QC_reads",
        "Contigs",
        "QC_contigs",
        "checkv_Complete",
        "checkv_High-quality",
        "checkv_Medium-quality",
        "checkv_Low-quality",
        "checkv_Not-determined",
        "VS_dsDNAphage",
        "VS_ssDNA",
        "VS_NCLDV",
        "VS_RNA",
        "VS_lavidaviridae",
    ]

    standard = ["NaN"] * len(column_names)

    df = pd.DataFrame(columns=column_names)

    metaQuast_reads = pd.read_table(
        "results/metaQUAST/combined_reference/transposed_report.tsv",
        usecols=["Assembly", "# contigs (>= 0 bp)"],
    )
    metaQuast_reads["Assembly"] = metaQuast_reads["Assembly"].str.replace(
        "_contigs", ""
    )
    metaQuast_reads = dict(
        zip(metaQuast_reads["Assembly"], metaQuast_reads["# contigs (>= 0 bp)"])
    )

    for x in samples["sample"]:
        checkV = getcheckv(
            "results/checkv/virsorter2/" + x + "/quality_summary.tsv"
        )
        virSorter2 = getVirSorter(
            "results/virsorter2/" + x + "/final-viral-score.tsv"
        )

        raw_reads, filtered_reads = get_trimmed_report(x)
        row = standard
        # Sampe name
        row[0] = x
        # Amount of reads filtered
        row[1] = raw_reads
        row[2] = filtered_reads
        row[3] = metaQuast_reads[x]
        row[4] = getMetaSpades(
            "results/checkv/virsorter2/" + x + "/filtered/filtered_combined.fna"
        )
        # checkv Score
        row[5] = checkV["Complete"]
        row[6] = checkV["High-quality"]
        row[7] = checkV["Medium-quality"]
        row[8] = checkV["Low-quality"]
        row[9] = checkV["Not-determined"]

        row[10] = virSorter2["dsDNAphage"]
        row[11] = virSorter2["ssDNA"]
        row[12] = virSorter2["NCLDV"]
        row[13] = virSorter2["RNA"]
        row[14] = virSorter2["lavidaviridae"]
        df.loc[len(df.index)] = row
        df.to_csv("./results/statistics/Sample_stats_virsorter2.tsv", sep="\t")


def create_sample_stats_vibrant():
    """
    Creates a table containig necessary stats from vibrant results,
    sample stats and checkv
    """
    column_names = [
        "Sample",
        "Raw_reads",
        "QC_reads",
        "Contigs",
        "QC_contigs",
        "checkv_Complete",
        "checkv_High-quality",
        "checkv_Medium-quality",
        "checkv_Low-quality",
        "checkv_Not-determined",
        "fragment",
        "low quality draft",
        "medium quality draft",
        "high quality draft",
        "complete circular",
        "lytic",
        "lysogenic",
    ]

    standard = ["NaN"] * len(column_names)

    df = pd.DataFrame(columns=column_names)

    metaQuast_reads = pd.read_table(
        "results/metaQUAST/combined_reference/transposed_report.tsv",
        usecols=["Assembly", "# contigs (>= 0 bp)"],
    )
    metaQuast_reads["Assembly"] = metaQuast_reads["Assembly"].str.replace(
        "_contigs", ""
    )
    metaQuast_reads = dict(
        zip(metaQuast_reads["Assembly"], metaQuast_reads["# contigs (>= 0 bp)"])
    )

    for x in samples["sample"]:
        checkV = getcheckv(
            "results/checkv/vibrant/" + x + "/quality_summary.tsv"
        )
        vibrant = get_vibrant_quality(
            "results/vibrant/"
            + x
            + "/VIBRANT_contigs/VIBRANT_results_contigs/VIBRANT_genome_quality_contigs.tsv"
        )
        raw_reads, filtered_reads = get_trimmed_report(x)
        row = standard
        # Sampe name
        row[0] = x
        # Amount of reads filtered
        row[1] = raw_reads
        row[2] = filtered_reads
        row[3] = metaQuast_reads[x]
        row[4] = getMetaSpades(
            "results/checkv/vibrant/" + x + "/filtered/filtered_combined.fna"
        )
        # checkv Score
        row[5] = checkV["Complete"]
        row[6] = checkV["High-quality"]
        row[7] = checkV["Medium-quality"]
        row[8] = checkV["Low-quality"]
        row[9] = checkV["Not-determined"]

        row[10] = vibrant["fragment"]
        row[11] = vibrant["low quality draft"]
        row[12] = vibrant["medium quality draft"]
        row[13] = vibrant["high quality draft"]
        row[14] = vibrant["complete circular"]
        row[15] = vibrant["lytic"]
        row[16] = vibrant["lysogenic"]

        df.loc[len(df.index)] = row
        df.to_csv("./results/statistics/Sample_stats_vibrant.tsv", sep="\t")


def get_trimmed_report(sample):
    """
    Returns the fastp quality controlled number of sequences,
    before and after.
    """
    df = pd.read_json("results/trimmed/report/" + sample + ".json")
    raw = df["summary"]["before_filtering"]["total_reads"]
    qc = df["summary"]["after_filtering"]["total_reads"]
    return raw, qc


def getcheckv(input):
    """
    Helper fucntion for gathering the checkv results,
    counts instances of each quality threshold.
    """
    quality_summary = pd.read_table(
        input, usecols=["contig_id", "checkv_quality"]
    )
    counts = quality_summary["checkv_quality"].value_counts()
    index = [
        "Not-determined",
        "Low-quality",
        "Medium-quality",
        "High-quality",
        "Complete",
        "Total-checked",
    ]
    data = ["NaN", "NaN", "NaN", "NaN", "NaN", "NaN"]
    df = pd.DataFrame(data, columns=["counts"], index=index)

    for x in index:
        try:
            df["counts"][x] = counts[x]
        except:
            df["counts"][x] = 0

    featurescheckv = {
        "Not-determined": df["counts"]["Not-determined"],
        "Low-quality": df["counts"]["Low-quality"],
        "Medium-quality": df["counts"]["Medium-quality"],
        "High-quality": df["counts"]["High-quality"],
        "Complete": df["counts"]["Complete"],
        "Total-checked": counts.sum(),
    }

    return featurescheckv


def get_vibrant_quality(quality_file):
    """
    Helper function that counts each instance of each VIBRANT result from input.
    """
    dct = {
        "fragment": 0,
        "low quality draft": 0,
        "medium quality draft": 0,
        "high quality draft": 0,
        "complete circular": 0,
        "lytic": 0,
        "lysogenic": 0,
    }

    df = pd.read_table(quality_file)
    dct.update(df.groupby("type").count()["scaffold"].to_dict())
    dct.update(df.groupby("Quality").count()["scaffold"].to_dict())

    return dct


def fastq_reads(input):
    """
    Helper fucntion for counting amount of sequences within a FASTQ file.
    """
    with open(input, "r") as fp:
        x = len(fp.readlines())
        return x / 4


def getMetaSpades(input):
    """
    Helper fucntion for counting amount of sequences within a metaSpades assembly.
    """
    contig_count = 0
    with open(input, "r") as f:
        for line in f.readlines():
            if ">" in line:
                contig_count += 1
    return contig_count


def getVirSorter(input):
    """
    Helper fucntion for gathering and counting all features gathered
    from virsorter2
    """
    featuresvirsorter2 = {
        "dsDNAphage": 0,
        "ssDNA": 0,
        "NCLDV": 0,
        "RNA": 0,
        "lavidaviridae": 0,
        "Total-Viruses": 0,
    }
    with open(input, "r") as f:
        for line in f.readlines():
            featuresvirsorter2["Total-Viruses"] += 1
            if "dsDNAphage" in line:
                featuresvirsorter2["dsDNAphage"] += 1
            if "ssDNA" in line:
                featuresvirsorter2["ssDNA"] += 1
            if "NCLDV" in line:
                featuresvirsorter2["NCLDV"] += 1
            if "RNA" in line:
                featuresvirsorter2["RNA"] += 1
            if "lavidaviridae" in line:
                featuresvirsorter2["lavidaviridae"] += 1
    new_dict = {k: v - 1 for k, v in featuresvirsorter2.items()}

    return new_dict


def get_vOTU_coverage_Individual(input):
    """
    Helper function for gathering the vOTU coverage
    """
    postfilter_coverage = pd.read_table(input, usecols=["Covered_percent"])
    return postfilter_coverage


def get_graphAnalyzer():
    """
    Helper function for gathering results from graphanalyzer
    """
    gA_df = pd.read_table(
        "./results/graphanalyzer/results_vcontact2_vOTU_results.csv",
        sep=",",
        usecols=[
            "Scaffold",
            "Accession",
            "Status",
            "Family",
            "Subfamily",
            "Genus",
        ],
    )
    gA_df["indexNumber"] = [int(i[5:]) for i in gA_df["Scaffold"]]
    gA_df.sort_values(
        by=["indexNumber"], ascending=[True], inplace=True, ignore_index=True
    )
    gA_df.drop("indexNumber", 1, inplace=True)
    return gA_df


def get_DRAMv(input):
    """
    Helper function for creating dataframe from DRAMv results
    """
    DRAMv = pd.read_table(
        input, sep="\t", usecols=["contig_id", "provirus", "checkv_quality"]
    )


def stats_vOTUs_virsorter2():
    """
    Function that creates the stats table of the vOTUs based on virsorter2 rerun.
    """
    column_names = [
        "vOTU_ID",
        "provirus",
        "Checkv_quality",
        "Accession",
        "Status",
        "Host",
        "Subfamily",
        "Genus",
    ]
    standard = ["NaN"] * len(column_names)
    vOTU_stats = pd.DataFrame(columns=column_names)
    gA_df = get_graphAnalyzer()
    cV_df = pd.read_table(
        "results/checkv/vOTU/virsorter2/quality_summary.tsv",
        sep="\t",
        usecols=["contig_id", "provirus", "checkv_quality"],
    )
    cV_df = cV_df.rename(columns={"contig_id": "Scaffold"})
    result = pd.merge(gA_df, cV_df, on="Scaffold", how="outer")
    result.to_csv("./results/statistics/vOTU_Stats_virsorter2.tsv", sep="\t")


def stats_vOTUs_vibrant():
    """
    Function that creates the stats table of the vOTUs based on vibrant rerun.
    """
    column_names = [
        "vOTU_ID",
        "provirus",
        "Checkv_quality",
        "Accession",
        "Status" "Host",
        "Subfamily",
        "Genus",
    ]
    standard = ["NaN"] * len(column_names)
    vOTU_stats = pd.DataFrame(columns=column_names)
    gA_df = get_graphAnalyzer()
    cV_df = pd.read_table(
        "results/checkv/vOTU/vibrant/quality_summary.tsv",
        sep="\t",
        usecols=["contig_id", "provirus", "checkv_quality"],
    )
    cV_df = cV_df.rename(columns={"contig_id": "Scaffold"})
    result = pd.merge(gA_df, cV_df, on="Scaffold", how="outer")
    result.to_csv("./results/statistics/vOTU_Stats_vibrant.tsv", sep="\t")


def vOTU_AMG_stats():
    """
    Function that creates the  AMG table and aggregates the functional annotation
    from VIBRANT and DRAMv
    """
    column_names = ["protein/gene", "scaffold", "ID", "Description"]
    df = pd.DataFrame(columns=column_names)
    vibrant_amgs = pd.read_table(
        "results/vibrant/vOTU/VIBRANT_vOTU_derep95_combined/VIBRANT_results_vOTU_derep95_combined/VIBRANT_AMG_individuals_vOTU_derep95_combined.tsv",
        sep="\t",
        usecols=["protein", "scaffold", "AMG KO", "AMG KO name"],
    )
    dramv_amgs = pd.read_table(
        "results/DRAMv/distilled/amg_summary.tsv",
        sep="\t",
        usecols=["gene", "scaffold", "gene_id", "gene_description"],
    )
    for index, x in vibrant_amgs.iterrows():
        row = {"protein/gene": "", "scaffold": "", "ID": "", "Description": ""}

        row["protein/gene"] = x["protein"]
        row["scaffold"] = x["scaffold"]
        row["ID"] = x["AMG KO"]
        row["Description"] = x["AMG KO name"]
        df = df.append(row, ignore_index=True)
    for index, x in dramv_amgs.iterrows():
        row = {"protein/gene": "", "scaffold": "", "ID": "", "Description": ""}

        row["protein/gene"] = x["gene"]
        row["scaffold"] = x["scaffold"]
        row["ID"] = x["gene_id"]
        row["Description"] = x["gene_description"]
        df = df.append(row, ignore_index=True)

    df.to_csv("./results/statistics/vOTU_AMGs.tsv", sep="\t")


def combine_vOTU_stats():
    """
    Function that creates the stats table of the vOTUs based on
    both virsorter2 rerun and VIBRANT.
    """
    df_v = pd.read_table(
        "./results/statistics/vOTU_Stats_virsorter2.tsv", sep="\t"
    )
    df_vs2 = pd.read_table(
        "./results/statistics/vOTU_Stats_vibrant.tsv", sep="\t"
    )
    df = df_v.set_index("Scaffold").combine_first(df_vs2.set_index("Scaffold"))
    df.drop(df.filter(regex="Unnamed"), axis=1, inplace=True)
    df.sort_values(
        by="Scaffold",
        inplace=True,
        key=lambda x: x.str.extract("(\d+)").squeeze().astype(int),
    )
    df.to_csv("./results/statistics/vOTU_stats_combined.tsv", sep="\t")


def create_relative_Abundance():
    """
    Function that creates the relative abundance table.
    It calculates the relative abundance based on the coverage file.
    """
    df = pd.read_table("./results/contig_stats/raw_coverage_table.tsv")
    ra_df = pd.DataFrame(columns=df.columns, index=df.index, data=None)
    ra_df["ID"] = df["ID"]
    for column in df.columns[1:]:
        total_instances = df[column].sum(axis=0)
        for index_r, v in enumerate(df[column]):
            relative_abundance = (v / total_instances) * 100
            ra_df[column][index_r] = round(relative_abundance, 2)
    ra_df = ra_df.set_index("ID")
    ra_df.to_csv("./results/statistics/vOTU_Relative_Abundance.tsv", sep="\t")


def combine_sample_stats():
    """
    Function that creates the Combined stats table. Combines results gathered from the
    stats tables of virsorter2 and VIBRANT
    """
    df1 = pd.read_table("./results/statistics/Sample_stats_virsorter2.tsv")
    df2 = pd.read_table("./results/statistics/Sample_stats_vibrant.tsv")

    my_df1 = []
    for index, row in df1.iterrows():
        qc = int(row["QC_contigs"])
        d = {
            "Sample": row["Sample"],
            "Contigs": (row["Contigs"]),
            "Viral Contigs": (
                int(
                    row["checkv_Complete"]
                    + row["checkv_High-quality"]
                    + row["checkv_Medium-quality"]
                    + row["checkv_Low-quality"]
                    + row["checkv_Not-determined"]
                )
            ),
            "QC Viral Contigs": qc,
        }
        my_df1.append(d)
    my_df2 = []
    for index, row in df2.iterrows():
        qc = int(row["QC_contigs"])
        d = {
            "Sample": "",
            "Contigs": 0,
            "Viral Contigs": (
                int(
                    row["checkv_Complete"]
                    + row["checkv_High-quality"]
                    + row["checkv_Medium-quality"]
                    + row["checkv_Low-quality"]
                    + row["checkv_Not-determined"]
                )
            ),
            "QC Viral Contigs": qc,
        }
        my_df2.append(d)
    new_df1 = pd.DataFrame(
        my_df1,
        columns=["Sample", "Contigs", "Viral Contigs", "QC Viral Contigs"],
    )
    new_df2 = pd.DataFrame(
        my_df2,
        columns=["Sample", "Contigs", "Viral Contigs", "QC Viral Contigs"],
    )
    combined_df = new_df1.add(new_df2, fill_value=0)
    combined_df.set_index("Sample")
    combined_df.to_csv(
        "./results/statistics/Combined_Sample_stats.tsv", sep="\t", index=False
    )


def plot_sampleStats():
    """
    Function that plots the samples stats on the contig progression
    """
    data = []
    df = pd.read_table(
        "./results/statistics/Combined_Sample_stats.tsv", sep="\t"
    )
    for index, row in df.iterrows():
        data.append(row.to_list())
    df2 = pd.DataFrame(
        data, columns=["Sample", "Contigs", "Viral Contigs", "QC Viral Contigs"]
    ).set_index("Sample")
    # create figure and axis
    fig, ax = plt.subplots(figsize=(20, 20))
    # setting the axis' labels
    ax.set_ylabel("Number of Contigs", fontsize=32)
    ax.set_xlabel("Pipeline Intervals", fontsize=32)
    plt.rcParams.update({"font.size": 32})
    plt.xticks(fontsize=25, rotation=90)
    plt.yticks(fontsize=25)
    # transposing (switching rows and columns) of DataFrame df and
    # plot a line for each column on the axis ax, which was created previously
    df2.T.plot(ax=ax)
    plt.savefig("./results/statistics/plots/Combined_Sample_stats_all.png")


def plot_sampleStats_viral():
    """
    Function that plots the samples stats on the contig progression,
    but only from identified viral contigs to QC contigs
    """
    data = []
    df = pd.read_table(
        "./results/statistics/Combined_Sample_stats.tsv",
        sep="\t",
        usecols=["Sample", "Viral Contigs", "QC Viral Contigs"],
    )
    for index, row in df.iterrows():
        data.append(row.to_list())
    df2 = pd.DataFrame(
        data, columns=["Sample", "Viral Contigs", "QC Viral Contigs"]
    ).set_index("Sample")
    # create figure and axis
    fig, ax = plt.subplots(figsize=(20, 20))
    # setting the axis' labels
    ax.set_ylabel("Number of Contigs", fontsize=32)
    ax.set_xlabel("Pipeline Intervals", fontsize=32)
    plt.rcParams.update({"font.size": 32})
    plt.xticks(fontsize=25, rotation=90)
    plt.yticks(fontsize=25)
    # transposing (switching rows and columns) of DataFrame df and
    # plot a line for each column on the axis ax, which was created previously
    df2.T.plot(ax=ax)
    plt.xticks(ticks=[0, 1], labels=["Viral Contigs", "QC Viral Contigs"])
    plt.savefig(
        "./results/statistics/plots/Combined_Sample_stats_viral_contigs.png"
    )


def relative_abundace_plot():
    """
    Function that plots the relative abundance as a boxplot
    """
    df = pd.read_table(
        "./results/statistics/vOTU_Relative_Abundance.tsv", sep="\t"
    )
    df = df.drop("ID", axis=1)
    df.columns = ["S" + str(i + 1) for i in range(len(df.columns))]
    b_plot = df.boxplot(column=list(df.columns.values))
    plt.title("Relative Abundance")
    b_plot.plot()
    fig = b_plot.get_figure()
    plt.xlabel("Samples")
    plt.ylabel("Abundance in %")
    fig.savefig(
        "./results/statistics/plots/relative_abundance_plot.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_families_plot():
    """
    Function that plots a barchart of the taxonomic annotation of the vOTUs
    at the family level
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv", sep="\t")

    plot_df = df[["Scaffold", "Family", "Status"]]

    counts_table = plot_df[["Family"]].value_counts()
    b_plot = counts_table.plot(kind="barh")
    fig = b_plot.get_figure()
    plt.ylabel("Family name")
    plt.xlabel("Count")
    plt.xticks(
        rotation=45,
        horizontalalignment="right",
        fontweight="light",
        fontsize="small",
    )
    fig.savefig(
        "./results/statistics/plots/Family_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_Subfamily_plot():
    """
    Function that plots a barchart of the taxonomic annotation of the vOTUs
    at the subfamily level
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv", sep="\t")

    plot_df = df[["Scaffold", "Subfamily", "Status"]]

    counts_table = plot_df[["Subfamily"]].value_counts()
    b_plot = counts_table.plot(kind="barh")
    fig = b_plot.get_figure()
    plt.ylabel("Subfamily name")
    plt.xlabel("Count")
    plt.xticks(
        rotation=45,
        horizontalalignment="right",
        fontweight="light",
        fontsize="small",
    )
    fig.savefig(
        "./results/statistics/plots/Subfamily_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_Genus_plot():
    """
    Function that plots a barchart of the taxonomic annotation of the vOTUs
    at the genus level
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv", sep="\t")

    plot_df = df[["Scaffold", "Genus", "Status"]]

    counts_table = plot_df[["Genus"]].value_counts()
    b_plot = counts_table.plot(kind="barh")
    fig = b_plot.get_figure()
    plt.ylabel("Genus name")
    plt.xlabel("Count")
    plt.xticks(
        rotation=45,
        horizontalalignment="right",
        fontweight="light",
        fontsize="small",
    )
    fig.savefig(
        "./results/statistics/plots/Genus_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_checkv_plot():
    """
    Function that plots a barchart of the vOTU checkv quality
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv")

    plot_df = df[["Scaffold", "checkv_quality", "Status"]]

    counts_table = plot_df[["checkv_quality"]].value_counts()
    b_plot = counts_table.plot(kind="barh")
    fig = b_plot.get_figure()
    plt.xlabel("checkv Quality")
    plt.ylabel("# vOTU in category")
    plt.xticks(
        rotation=45,
        horizontalalignment="right",
        fontweight="light",
        fontsize="small",
    )
    fig.savefig(
        "./results/statistics/plots/checkv_quality_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_to_reads_mapping():
    """
    Function that creates table which contains the vOTU and its
    corresponding actual contig before renaming.
    """
    contig_ID = []
    vOTU_ID = []
    with open("./results/cdhit/derep95_combined.fasta", "r") as sf:
        for ln in sf.readlines():
            if ln.startswith(">"):
                contig_ID.append(ln[1:].strip("\n"))

    with open("./results/cdhit/vOTU_derep95_combined.fasta", "r") as sf:
        for ln in sf.readlines():
            if ln.startswith(">"):
                vOTU_ID.append(ln[1:].strip("\n"))

    combined_df = pd.DataFrame()
    combined_df["vOTU_ID"] = vOTU_ID
    combined_df["Contig_ID"] = contig_ID

    cV_combined_vibrant = pd.DataFrame(columns=["contig_id", "provirus"])
    for x in samples["sample"]:
        cV_df = pd.read_table(
            "./results/checkv/vibrant/" + x + "/quality_summary.tsv", sep="\t"
        )
        cV_combined_vibrant = pd.concat([cV_combined_vibrant, cV_df])

    cV_combined_virsorter2 = pd.DataFrame(columns=["contig_id", "provirus"])
    for x in samples["sample"]:
        cV_df = pd.read_table(
            "./results/checkv/virsorter2/" + x + "/quality_summary.tsv",
            sep="\t",
        )
        cV_combined_virsorter2 = pd.concat([cV_combined_virsorter2, cV_df])

    cV_provirus = pd.concat(
        [
            cV_combined_vibrant.loc[cV_combined_vibrant["provirus"] == "Yes"],
            cV_combined_virsorter2.loc[
                cV_combined_virsorter2["provirus"] == "Yes"
            ],
        ]
    )

    outer = re.compile("/(\|\|.+|_fragment.+)/gm")

    cV_provirus = cV_provirus["contig_id"].to_list()
    cV_provirus = list(
        map(lambda s: re.sub(r"(\|\|.+|_fragment.+)", "", s), cV_provirus)
    )

    for index, row in combined_df.iterrows():
        check = re.sub(r"(\|\|.+|_fragment.+)", "", str(row["Contig_ID"]))
        if check in cV_provirus:
            combined_df.at[index, "Provirus"] = "Yes"
        else:
            combined_df.at[index, "Provirus"] = "No"
    combined_df.to_csv(
        "./results/statistics/vOTU_mapped_to_reads.tsv", sep="\t"
    )


def gather_lytic():
    """
    Function that updates the vOTU stats table to include lytic/lysogenic classification
    """
    provirus_df = pd.read_table(
        "./results/statistics/vOTU_mapped_to_reads.tsv", sep="\t"
    )
    combined_df = pd.read_table(
        "./results/statistics/vOTU_stats_combined.tsv", sep="\t"
    )
    lysogenic_df = pd.read_table(
        "./results/vibrant/vOTU/VIBRANT_vOTU_derep95_combined/VIBRANT_results_vOTU_derep95_combined/VIBRANT_genome_quality_vOTU_derep95_combined.tsv",
        sep="\t",
    )
    combined_df["provirus"] = provirus_df["Provirus"]
    for index, row in combined_df.iterrows():
        lysogenic = lysogenic_df.loc[
            lysogenic_df["scaffold"] == str(row["Scaffold"])
        ]
        if lysogenic.empty:
            combined_df.at[index, "Type"] = "n.a."
        else:
            if len(lysogenic) > 1:
                combined_df.at[index, "Type"] = str(lysogenic.iloc[0]["type"])
            else:
                combined_df.at[index, "Type"] = str(lysogenic["type"].item())
    combined_df.to_csv("./results/statistics/vOTU_stats_combined.tsv", sep="\t")


def vOTU_provirus_plot():
    """
    Function that plots a barchart of if a vOTU is a provirus or not
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv")

    plot_df = df[["Scaffold", "provirus"]]
    print(plot_df)
    plot_df = plot_df[["provirus"]].replace(
        {"Yes": "Provirus", "No": "Non-Provirus"}
    )
    counts_table = plot_df[["provirus"]].value_counts()
    b_plot = counts_table.plot(kind="bar")
    fig = b_plot.get_figure()
    plt.xlabel("Is Provirus or Non-Provirus")
    plt.ylabel("# vOTU")
    plt.xticks(rotation=45, fontweight="light", fontsize="medium")
    fig.savefig(
        "./results/statistics/plots/provirus_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


def vOTU_lysogenic_plot():
    """
    Function that plots a barchart of if a vOTU is lysogenic or lytic
    """
    df = pd.read_table("./results/statistics/vOTU_stats_combined.tsv")

    plot_df = df[["Scaffold", "Type"]]

    counts_table = plot_df[["Type"]].value_counts()
    b_plot = counts_table.plot(kind="bar")
    fig = b_plot.get_figure()
    plt.xlabel("Lysogenic / Lytic")
    plt.ylabel("# vOTU")
    plt.xticks(rotation=45, fontweight="light", fontsize="medium")
    fig.savefig(
        "./results/statistics/plots/lysogenic_vOTU_BarChart.png",
        bbox_inches="tight",
    )
    plt.clf()


if not os.path.exists(st_dir):
    print("Missing output directory, creating now ....")
    os.makedirs(st_dir)
if not os.path.exists(out_dir):
    print("Missing output directory, creating now ....")
    os.makedirs(out_dir)

# Runs all small functions to create statistics, tables and plots.
create_sample_stats_virsorter2()
create_sample_stats_vibrant()
stats_vOTUs_virsorter2()
stats_vOTUs_vibrant()
combine_vOTU_stats()
vOTU_AMG_stats()
create_relative_Abundance()
combine_sample_stats()
vOTU_to_reads_mapping()
gather_lytic()

relative_abundace_plot()
vOTU_families_plot()
vOTU_Subfamily_plot()
vOTU_Genus_plot()
vOTU_checkv_plot()
vOTU_provirus_plot()
vOTU_lysogenic_plot()
plot_sampleStats()
plot_sampleStats_viral()
