#!/usr/bin/env Rscript


# env ---------------------------------------------------------------------

library(tidyverse)

# script ------------------------------------------------------------------

lapply(snakemake@input[["covstats"]], function(covstats) {
  suppressMessages(read_tsv(covstats)) %>% print() %>%
    mutate(filtered_coverage = case_when(
      Covered_percent < as.double(snakemake@params[["min_coverage"]]) ~ 0,
      Covered_percent >= as.double(snakemake@params[["min_coverage"]]) ~ Avg_fold
    )) %>%
    rename(ID = 1) %>%
    select(ID, filtered_coverage) %>%
    mutate(sample_id = str_extract(covstats, "(?<=contig_stats/).*(?=/postfilter_coverage_stats)"))
}) %>%
  bind_rows() %>%
  pivot_wider(names_from = sample_id, values_from = filtered_coverage) %>%
  write_tsv(snakemake@output[["abundance_table"]])