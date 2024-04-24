"""
script to calculate the aggregated results and statistics for the pipeline
"""
import pandas as pd
import os
import re


output_path = snakemake.params.output_path
samples = snakemake.params.samples
st_dir = snakemake.output.dir


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
        output_path + "/metaQUAST/combined_reference/transposed_report.tsv",
        usecols=["Assembly", "# contigs (>= 0 bp)"],
    )
    metaQuast_reads["Assembly"] = metaQuast_reads["Assembly"].str.replace(
        "_contigs", ""
    )
    metaQuast_reads = dict(
        zip(
            metaQuast_reads["Assembly"], metaQuast_reads["# contigs (>= 0 bp)"]
        )
    )

    for x in samples:
        checkV = getcheckv(
            output_path
            + "/checkv/virsorter2/"
            + x
            + "/quality_summary.tsv"
        )
        virSorter2 = getVirSorter(
            output_path + "/virsorter2/" + x + "/final-viral-score.tsv"
        )

        raw_reads, filtered_reads = get_trimmed_report(x)
        row = standard
        # Sampe name
        row[0] = x
        # Amount of reads filtered
        row[1] = raw_reads
        row[2] = filtered_reads
        try:
            row[3] = metaQuast_reads[x]
        except KeyError:
            row[3] = metaQuast_reads["contigs"]
        row[4] = getMetaSpades(
            output_path
            + "/filtered_virsorter2/"
            + x
            + "/filtered_combined.fna"
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
        df.to_csv(
            output_path + "/statistics/Sample_stats_virsorter2.tsv", sep="\t"
        )


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
        output_path + "/metaQUAST/combined_reference/transposed_report.tsv",
        usecols=["Assembly", "# contigs (>= 0 bp)"],
    )
    metaQuast_reads["Assembly"] = metaQuast_reads["Assembly"].str.replace(
        "_contigs", ""
    )
    metaQuast_reads = dict(
        zip(
            metaQuast_reads["Assembly"], metaQuast_reads["# contigs (>= 0 bp)"]
        )
    )

    for x in samples:
        checkV = getcheckv(
            output_path + "/checkv/vibrant/" + x + "/quality_summary.tsv"
        )
        vibrant = get_vibrant_quality(
            output_path
            + "/vibrant/"
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
        try:
            row[3] = metaQuast_reads[x]
        except KeyError:
            row[3] = metaQuast_reads["contigs"]
        row[4] = getMetaSpades(
            output_path + "/filtered_vibrant/" + x + "/filtered_combined.fna"
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
        df.to_csv(
            output_path + "/statistics/Sample_stats_vibrant.tsv", sep="\t"
        )


def get_trimmed_report(sample):
    """
    Returns the fastp quality controlled number of sequences,
    before and after.
    """
    df = pd.read_json(output_path + "/fastp_pe/" + sample + ".json")
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
            df.loc[x, "counts"] = counts[x]
        except:
            df.loc[x, "counts"] = 0

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
        output_path + "/graphanalyzer/results_vcontact2_vOTU_results.csv",
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


def vOTU_AMG_stats():
    """
    Function that creates the AMG table and aggregates the functional annotation
    from VIBRANT and DRAMv
    """
    column_names = ["protein/gene", "scaffold", "ID", "Description"]
    dramv_amgs = pd.read_table(
        output_path + "/DRAMv/distilled/amg_summary.tsv",
        sep="\t",
        usecols=["gene", "scaffold", "gene_id", "gene_description"],
    )
    rows_list = []
    for _, row in dramv_amgs.iterrows():
        row_dict = {
            "protein/gene": row["gene"],
            "scaffold": row["scaffold"],
            "ID": row["gene_id"],
            "Description": row["gene_description"],
        }
        rows_list.append(row_dict)
    
    df = pd.DataFrame(rows_list, columns=column_names)

    df.to_csv(output_path + "/statistics/vOTU_AMGs.tsv", sep="\t", index=False)


def create_relative_Abundance():
    """
    Function that creates the relative abundance table.
    It calculates the relative abundance based on the coverage file.
    """
    df = pd.read_table(output_path + "/contig_stats/raw_coverage_table.tsv")
    ra_df = pd.DataFrame(columns=df.columns, index=df.index, data=None)
    ra_df["ID"] = df["ID"]
    for column in df.columns[1:]:
        total_instances = df[column].sum(axis=0)
        for index_r, v in enumerate(df[column]):
            relative_abundance = (v / total_instances) * 100
            ra_df.loc[index_r, column] = round(relative_abundance, 2)
    ra_df = ra_df.set_index("ID")
    ra_df.to_csv(
        output_path + "/statistics/vOTU_Relative_Abundance.tsv", sep="\t"
    )


def combine_sample_stats():
    """
    Function that creates the Combined stats table. Combines results gathered from the
    stats tables of virsorter2 and VIBRANT
    """
    df1 = pd.read_table(
        output_path + "/statistics/Sample_stats_virsorter2.tsv"
    )
    df2 = pd.read_table(output_path + "/statistics/Sample_stats_vibrant.tsv")

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
        output_path + "/statistics/Combined_Sample_stats.tsv",
        sep="\t",
        index=False,
    )


def vOTU_to_reads_mapping():
    """
    Function that creates table which contains the vOTU and its
    corresponding actual contig before renaming.
    """
    contig_ID = []
    vOTU_ID = []
    with open(output_path + "/cdhit/derep95_combined.fasta", "r") as sf:
        for ln in sf.readlines():
            if ln.startswith(">"):
                contig_ID.append(ln[1:].strip("\n"))

    with open(output_path + "/vOTU/vOTU_derep95_combined.fasta", "r") as sf:
        for ln in sf.readlines():
            if ln.startswith(">"):
                vOTU_ID.append(ln[1:].strip("\n"))

    combined_df = pd.DataFrame()
    combined_df["vOTU_ID"] = vOTU_ID
    combined_df["Contig_ID"] = contig_ID

#    cV_combined_vibrant = pd.DataFrame(columns=["contig_id", "provirus"])
#    for x in samples:
#        cV_df = pd.read_table(
#            output_path + "/checkv/vibrant/" + x + "/quality_summary.tsv",
#            sep="\t",
#        )
#        cV_combined_vibrant = pd.concat([cV_combined_vibrant, cV_df])

    cV_combined_virsorter2 = pd.DataFrame(columns=["contig_id", "provirus"])
    for x in samples:
        cV_df = pd.read_table(
            output_path
            + "/checkv/virsorter2/"
            + x
            + "/quality_summary.tsv",
            sep="\t",
        )
        cV_combined_virsorter2 = pd.concat([cV_combined_virsorter2, cV_df])

    cV_provirus = pd.concat(
        [
#            cV_combined_vibrant.loc[cV_combined_vibrant["provirus"] == "Yes"],
            cV_combined_virsorter2.loc[
                cV_combined_virsorter2["provirus"] == "Yes"
            ],
        ]
    )

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
        output_path + "/statistics/vOTU_mapped_to_reads.tsv", sep="\t"
    )


def copy_instrain_compare_output():
    """
    Function that copies the instrain compare output to the statistics folder
    """
    os.system(
        "cp "
        + output_path
        + "/instrain/compared_samples/output/compared_samples_comparisonsTable.tsv "
        + output_path
        + "/statistics/compared_samples_comparisonsTable.tsv "
    )


if not os.path.exists(st_dir):
    print("Missing output directory, creating now ....")
    os.makedirs(st_dir)

# Runs all small functions to create statistics and tables.
if __name__ == "__main__":
    create_sample_stats_virsorter2()
#    create_sample_stats_vibrant()
    vOTU_AMG_stats()
    create_relative_Abundance()
#    combine_sample_stats()
    vOTU_to_reads_mapping()
    copy_instrain_compare_output()
