#!/usr/bin/env Rscript

library(tidyverse)

outdir <- "/cluster/home/paulai/VirMake/working_dir/human/output"
# outdir <- snakemake@params[["outdir"]]

virsorter_round2_scores <-
  read_tsv(paste(outdir, "virsorter2_pass2/final-viral-score.tsv", sep="/")) %>%
  mutate(vOTU = str_extract(seqname, "^vOTU[:digit:]*")) %>%
  select(vOTU, everything())

print(virsorter_round2_scores)

dramv_summary <-
  read_tsv(paste(outdir, "DRAMv/distilled/vMAG_stats.tsv", sep="/")) %>%
  rename(scaffold = 1) %>%
  mutate(vOTU = str_extract(scaffold, "^vOTU[:digit:]*")) %>%
  select(vOTU, everything())

print(dramv_summary)