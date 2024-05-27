"""
script to retrieve sample statistics
"""
import pandas as pd
import re


def calculate_n50(lengths):
    """Calculate the N50 for a sequence of lengths."""
    lengths_sorted = sorted(lengths, reverse=True)
    csum = pd.Series(lengths_sorted).cumsum()
    half_total = csum.iloc[-1] / 2
    n50_index = csum[csum >= half_total].index[0]
    return lengths_sorted[n50_index]

def get_virus_id_stats(virus_id_file, sample_id):
    virus_id_stats = pd.read_csv(virus_id_file, sep="\t", index_col=False)
    stats = {
        'sample_id': sample_id,
        'n_genomes': len(virus_id_stats),
        'genome_median_length': virus_id_stats["length"].median(),
        'genome_min_length': virus_id_stats["length"].min(),
        'genome_max_length': virus_id_stats["length"].max(),
        'genomes_N50': calculate_n50(virus_id_stats["length"]),
        'n_provirus': ((virus_id_stats["provirus"] == "Yes") | (virus_id_stats["vir_id_provirus_assignment"] == "TRUE")).sum(),
        'checkv_complete': (virus_id_stats["checkv_quality"] == "Complete").sum(),
        'checkv_high_quality': (virus_id_stats["checkv_quality"] == "High-quality").sum(),
        'checkv_medium_quality': (virus_id_stats["checkv_quality"] == "Medium-quality").sum(),
        'checkv_low_quality': (virus_id_stats["checkv_quality"] == "Low-quality").sum(),
    }
    return pd.DataFrame([stats])

def summarize_virus_id_stats(stats_files):
    
    sample_pattern = r'virus_identification/(.*?)/gathered_quality_tables\.tsv'
    sample_ids = [re.search(sample_pattern, s).group(1) for s in stats_files]
    
    virus_id_summaries = []

    for vir_id_file, sample_id in zip(stats_files, sample_ids):
        virus_id_summaries.append(get_virus_id_stats(vir_id_file, sample_id))
    
    virus_id_summary = pd.concat(virus_id_summaries, ignore_index=True)

    return virus_id_summary

def process_abundance_data(abundance_file):
    abundance_data = pd.read_csv(abundance_file, sep="\t", index_col=0)
    
    presence = (abundance_data > 0)

    n_presence = presence.sum(axis=0)
    mean_present = abundance_data.where(presence).mean(axis=0)

    abundance_summary = pd.DataFrame({
        'sample_id': abundance_data.columns,
        'n_present': n_presence.values,
        'mean_present': mean_present.values.round(2),
    })

    return(abundance_summary)

def main():
    
    samplewise_summary = summarize_virus_id_stats(snakemake.input.virus_id_tables)
    
    if snakemake.input.rel_abund:
        abundance_stats = process_abundance_data(abundance_file=snakemake.input.rel_abund)
        samplewise_summary = pd.merge(left=samplewise_summary, right=abundance_stats, how="outer", on=["sample_id"])

    samplewise_summary.to_csv(snakemake.output.sample_stats, sep="\t", index=False)

if __name__ == "__main__":
    main()