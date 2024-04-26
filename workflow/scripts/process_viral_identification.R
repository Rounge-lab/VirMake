
library(tidyverse)

id_tool <- snakemake@wildcards[["id_tool"]]


if (id_tool == "virsorter") {
  vir_id_res <- 
    read_tsv(snakemake@input[["virus_id_table"]]) %>% 
    rename(contig_id = seqname,
           vir_id_name = seqname_new) %>%
    mutate(vir_id_provirus_assignment = str_detect(vir_id_name, "partial")) %>% 
    select(contig_id, 
           start_position = trim_bp_start, 
           end_position = trim_bp_end,
           vir_id_name,
           vir_id_provirus_assignment)
}

if (id_tool == "genomad") {
  vir_id_res <- 
    read_tsv(snakemake@input[["virus_id_table"]]) %>% 
    rename(vir_id_name = seq_name) %>%
    mutate(contig_id = str_remove(vir_id_name, "\\|provirus_[:digit:]*_[:digit:]*$")) %>% 
    mutate(start_position = case_when(topology %in% "Provirus" ~ str_extract(coordinates, "^[:digit:]*") %>% as.integer(),
                                    !topology %in% "Provirus" ~ 1 %>% as.integer()),
           end_position = case_when(topology %in% "Provirus" ~ str_extract(coordinates, "[:digit:]*$") %>% as.integer(),
                                    !topology %in% "Provirus" ~ length %>% as.integer())) %>% 
    mutate(vir_id_provirus_assignment = topology %in% "Provirus") %>% 
    select(contig_id, 
           start_position,
           end_position,
           vir_id_name,
           vir_id_provirus_assignment)
}


checkv_res_cont <- 
  read_tsv(snakemake@input[["checkv_contamination"]]) %>% 
  select(provirus, region_types, region_coords_bp, vir_id_name = contig_id)

checkv_res <- 
  read_tsv(snakemake@input[["checkv_res"]]) %>% 
  rename(vir_id_name = contig_id)

checkv_res_cont <-
  checkv_res_cont %>% 
  filter(provirus %in% "Yes") %>% 
  separate_rows(region_types, region_coords_bp, sep = ",") %>% 
  filter(str_detect(region_types, "viral")) %>% 
  mutate(start = str_extract(region_coords_bp, "^[:digit:]*") %>% as.integer(),
         end = str_extract(region_coords_bp, "[:digit:]*$") %>% as.integer(),
         length = end-start+1) %>% 
  select(vir_id_name, start, end, provirus, length) %>% 
  left_join(checkv_res %>% select(-c(contig_length, proviral_length)), by = c("vir_id_name", "provirus")) %>%
  bind_rows(checkv_res %>% 
              filter(provirus %in% "No") %>% 
              mutate(length = contig_length,
                     start = 1,
                     end = length) %>% 
              select(-c(contig_length, proviral_length))) %>% 
  left_join(vir_id_res, by = "vir_id_name") %>% 
  mutate(start = start_position + start -1,
         end = start_position + end -1) %>% 
  select(vir_id_name, contig_id, everything(), vir_id_start_pos = start_position, vir_id_end_pos = end_position) %>% 
  mutate(vir_id_tool = id_tool) %>% 
  group_by(vir_id_name) %>% 
  mutate(n_v = sum(vir_id_name == vir_id_name)) %>% 
  mutate(vir_id_name = case_when(n_v > 1 ~ paste(vir_id_name, seq(max(n_v)), sep = "_"),
                                 TRUE ~ vir_id_name)) %>% 
  ungroup()

checkv_res_cont %>% 
  write_tsv(snakemake@output[["reformatted_table"]])
  
