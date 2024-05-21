import pandas as pd

def get_samples(path):
    sample_table = pd.read_csv(path, delimiter='\t')
    samples = sample_table["sample_id"].to_list()
    return sample_table, samples

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
    Gets the threshold for checkv wuality control.
    """
    if threshold.lower() == "complete":
        return "$8~/(Complete)/"
    elif threshold.lower() == "high":
        return "$8~/(High-quality|Complete)/"
    elif threshold.lower() == "medium":
        return "$8~/(Medium-quality|High-quality|Complete)/"
    elif threshold.lower() == "low":
        return "$8~/(Low-quality|Medium-quality|High-quality|Complete)/"
    elif threshold.lower() == "not-determined":
        return "$8~/(Not-determined|Low-quality|Medium-quality|High-quality|Complete)/"
    else:
        return "$8~/(Medium-quality|High-quality|Complete)/"

def vibrant_virome(is_virome):
    """
    Gathers if the input is virome.
    """
    if is_virome.lower == "yes":
        return "-virome"
    else:
        return ""