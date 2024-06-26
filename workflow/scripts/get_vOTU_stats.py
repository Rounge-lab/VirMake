"""
script to retrieve statistics for the viral genomes
"""
import pandas as pd
from pathlib import Path


def process_taxonomy_info(taxonomy_file):
    taxonomy_cols = ["Scaffold","Closer","Accession","Status","VC","Level",
                     "Host","BaltimoreGroup","Realm","Kingdom","Phylum","Class",
                     "Order","Family","Subfamily","Genus"]

    taxonomy = pd.read_csv(taxonomy_file, index_col=0)
    taxonomy = taxonomy[taxonomy_cols]
    taxonomy = taxonomy.rename(columns={"Scaffold": "vOTU"})

    return taxonomy

def process_derep_virus_specs(specs_file, derep_file, id_file):
    specs_cols = ["virus_id","vir_id_tool","length",
                  "start","end","provirus","vir_id_provirus_assignment",
                  "checkv_quality","miuvig_quality","completeness"]

    specs = pd.read_csv(specs_file, sep='\t')
    specs = specs[specs_cols]

    derep = pd.read_csv(derep_file, header=None, names=["rep", "id"], sep='\t')
    derep["old_id"]=derep["rep"].apply(lambda path: Path(path).stem)
    derep["virus_id"]=derep["id"].apply(lambda path: Path(path).stem)
    derep = derep[["virus_id","old_id"]]

    contig_ids = pd.read_csv(id_file, sep='\t')
    contig_ids = contig_ids.rename(columns={"new_id": "vOTU"})

    derep = pd.merge(left=derep, right=contig_ids, how="outer", on="old_id")
    # derep = derep["virus_id", "vOTU"]

    combined_virus_pred_data = pd.merge(left=specs, right=derep, how="left", on="virus_id")

    derep_proc_specs = combined_virus_pred_data.groupby("vOTU").apply(lambda grp: pd.Series({
        'n_genomes': len(grp),
        'vOTU_length': grp.loc[grp["virus_id"] == grp["old_id"].iloc[0], "length"].iloc[0],
        'vOTU_provirus': ((grp.loc[grp["virus_id"] == grp["old_id"].iloc[0], "provirus"].iloc[0] == "Yes") |
                          (grp.loc[grp["virus_id"] == grp["old_id"].iloc[0], "vir_id_provirus_assignment"].iloc[0])),
        'pangenome_median_length': grp["length"].median(),
        'pangenome_min_length': grp["length"].min(),
        'pangenome_max_length': grp["length"].max(),
        'pangenome_n_provirus': ((grp["provirus"] == "Yes") | (grp["vir_id_provirus_assignment"] == "TRUE")).sum(),

    })).reset_index()

    return derep_proc_specs

def process_abundance_data(abundance_file):
    abundance_data = pd.read_csv(abundance_file, sep="\t", index_col=0)

    presence = (abundance_data > 0)

    n_presence = presence.sum(axis=1)
    mean_present = abundance_data.where(presence).mean(axis=1)

    abundance_summary = pd.DataFrame({
        'vOTU': abundance_data.index,
        'n_present': n_presence.values,
        'mean_present': mean_present.values.round(2),
    })

    return(abundance_summary)

def get_dramv_summary(dramv_summary_file):
    dramv_cols = ["vOTU", "Gene count", "Strand switches", "potential AMG count",
                  "Transposase present", "Possible Non-Viral Contig"]#,
                #   "Viral genes with unknown function", "Viral hypothetical genes",
                #   "Viral genes with viral benefits", "Viral genes with host benefits",
                #   "Viral structure genes"]
    dramv_summary = pd.read_csv(dramv_summary_file, sep="\t", index_col=0)
    dramv_summary.index.name = "vOTU"
    dramv_summary.reset_index(inplace=True)
    dramv_summary["vOTU"] = dramv_summary["vOTU"].str.replace(r'-cat_\d+', '', regex=True)

    dramv_summary = dramv_summary[dramv_cols]

    return dramv_summary

def main():

    taxonomy = process_taxonomy_info(taxonomy_file=snakemake.input.taxonomy)

    derep_viral_specs = process_derep_virus_specs(specs_file=snakemake.input.gathered_specs,
                                                  derep_file=snakemake.input.derep_file,
                                                  id_file=snakemake.input.contig_id_file)


    combined_table = pd.merge(left=derep_viral_specs, right=taxonomy, how="outer", on=["vOTU"])

    if snakemake.input.rel_abund:
        abundance_stats = process_abundance_data(abundance_file=snakemake.input.rel_abund)
        combined_table = pd.merge(left=combined_table, right=abundance_stats, how="outer", on=["vOTU"])

    if snakemake.input.DRAM_distilled_stats:
        abundance_stats = get_dramv_summary(dramv_summary_file=snakemake.input.DRAM_distilled_stats)
        combined_table = pd.merge(left=combined_table, right=abundance_stats, how="outer", on=["vOTU"])

    combined_table.to_csv(snakemake.output.vOTU_stats, sep="\t", index=False)


if __name__ == "__main__":
    main()
