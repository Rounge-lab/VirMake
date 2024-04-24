from snakemake.io import glob_wildcards

def get_samples(path):
    (
        SAMPLE,
        FRAC,
    ) = glob_wildcards(path + "/{sample}_{frac}.fastq.gz")

    # remove duplicates
    SAMPLE = sorted(list(set(SAMPLE)))
    FRAC = sorted(list(set(FRAC)))
    return SAMPLE, FRAC

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