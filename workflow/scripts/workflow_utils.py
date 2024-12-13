import pandas as pd

def get_samples(path):
    sample_table = pd.read_csv(path, delimiter='\t')
    samples = sample_table["sample_id"].to_list()
    return sample_table, samples

def get_qc_reads_loc(wildcards,
                     sample_table,
                     standard_output,
                     r):
    qc_r1_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'qc_r1'].iloc[0]
    qc_r2_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'qc_r2'].iloc[0]
    if pd.isna(qc_r1_loc) or pd.isna(qc_r2_loc) or qc_r1_loc == "" or qc_r2_loc == "":
        return standard_output+"/fastp_pe/"+wildcards.sample+f"_{r}.fastq"
    else:
        if r == "1":
            return qc_r1_loc
        elif r == "2":
            return qc_r2_loc

def get_assembly_loc(wildcards,
                     sample_table,
                     standard_output):
    assembly_loc = sample_table.loc[sample_table["sample_id"] == wildcards.sample, 'contigs'].iloc[0]
    if pd.isna(assembly_loc) or assembly_loc == "":
        return standard_output + "/metaSpades/"+wildcards.sample+"/contigs.fasta"
    else:
        return assembly_loc

# GLOBAL FUNCTIONS #
def get_min_quality(threshold):
    """
    Gets the threshold for checkv quality control.
    """
    if threshold.lower() == "complete":
        return ["Complete"]
    elif threshold.lower() == "high":
        return ["Complete", "High-quality"]
    elif threshold.lower() == "medium":
        return ["Complete", "High-quality", "Medium-quality"]
    elif threshold.lower() == "low":
        return ["Complete", "High-quality", "Medium-quality", "Low-quality"]
    elif threshold.lower() == "not-determined":
        return ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]
    else:
        print("Quality threshold not defined properly. Using medium.")
        return ["Complete", "High-quality", "Medium-quality"]

def vibrant_virome(is_virome):
    """
    Gathers if the input is virome.
    """
    if is_virome.lower == "yes":
        return "-virome"
    else:
        return ""
