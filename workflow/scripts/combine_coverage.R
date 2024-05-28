#!/usr/bin/env Rscript


# env ---------------------------------------------------------------------

suppressPackageStartupMessages(library(tidyverse))

cov_stats <-
  lapply(snakemake@input[["covstats"]], function(covstats) {
    read_tsv(covstats, show_col_types=FALSE) %>%
      rename(vOTU = 1) %>%
      mutate(sample_id = str_extract(covstats, "(?<=pileup/).*(?=/postfilter_coverage_stats)"))
  }) %>%
  bind_rows() %>%
  select(sample_id, everything())


binned_stats <-
  lapply(snakemake@input[["binned_coverage"]], function(covstats) {
    read_tsv(covstats, show_col_types=FALSE, skip=2) %>%
      rename(vOTU = 1) %>%
      mutate(sample_id = str_extract(covstats, "(?<=pileup/).*(?=/postfilter_coverage_binned)"))
  }) %>%
  bind_rows() %>%
  group_by(sample_id, vOTU) %>%
  summarize(median_binned = median(Cov)) %>%
  ungroup() 

cov_stats <-
  cov_stats %>%
  left_join(binned_stats, by = c("sample_id", "vOTU")) 


cov_stats %>%
  write_tsv(snakemake@output[["covstats"]])

cov_stats %>%
  mutate(filtered_reads = case_when(
        Covered_percent < as.double(snakemake@params[["min_coverage"]]) ~ 0,
        Covered_percent >= as.double(snakemake@params[["min_coverage"]]) ~ Plus_reads + Minus_reads
      )) %>%
  select(sample_id, vOTU, filtered_reads) %>%
  pivot_wider(names_from = sample_id, values_from = filtered_reads) %>%
  write_tsv(snakemake@output[["mapped_reads"]])

cov_stats %>%
  mutate(filtered_median_binned = case_when(
        Covered_percent < as.double(snakemake@params[["min_coverage"]]) ~ 0,
        Covered_percent >= as.double(snakemake@params[["min_coverage"]]) ~ median_binned
      )) %>%
  group_by(sample_id) %>%
  mutate(rel_median = round(filtered_median_binned/sum(filtered_median_binned), 5)) %>%
  ungroup() %>%
  select(sample_id, vOTU, rel_median) %>%
  pivot_wider(names_from = sample_id, values_from = rel_median) %>%
  write_tsv(snakemake@output[["rel_abundance"]])
