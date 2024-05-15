import pandas as pd

def create_relative_Abundance(raw_cov_file, rel_abund_file):
    """
    Function that creates the relative abundance table.
    It calculates the relative abundance based on the coverage file.
    """
    df = pd.read_table(raw_cov_file)
    ra_df = pd.DataFrame(columns=df.columns, index=df.index, data=None)
    ra_df["ID"] = df["ID"]
    for column in df.columns[1:]:
        total_instances = df[column].sum(axis=0)
        for index_r, v in enumerate(df[column]):
            relative_abundance = (v / total_instances) * 100
            ra_df.loc[index_r, column] = round(relative_abundance, 2)
    ra_df = ra_df.set_index("ID")
    ra_df.to_csv(rel_abund_file, sep="\t")

if __name__ == "__main__":
    create_relative_Abundance(raw_cov_file=snakemake.input.raw,
                              rel_abund_file=snakemake.output.rel)
