
library(tidyverse)

full_overlap_frac_threshold <- snakemake@params[["overlap_threshold"]] ## 0-1
accepted_quals <- snakemake@params[["checkv_quality"]] ## vector of checkv qualities
overlap_selection <- snakemake@params[["length_selection"]] ## min_length or max_length
tool_combination <- snakemake@params[["tool_combination"]] ## all or any

sample_id <- snakemake@wildcards[["sample"]]

viral_predictions <- 
  lapply(snakemake@input[["virus_pred_tables"]], function(x) {
    read_tsv(x, col_types = cols(start = "d", end = "d"))
  }) %>%
  bind_rows() %>%
  filter(checkv_quality %in% accepted_quals) %>%
  mutate(id = paste(vir_id_tool, vir_id_name, sep = "_"))

if (snakemake@params[["remove_lt2gene"]]) {
  viral_predictions <-
    viral_predictions %>%
    filter(!str_detect(vir_id_name, "lt2gene"))
}

if (nrow(viral_predictions) == 0) {
  detected_overlaps <-
    viral_predictions %>%
    mutate(value = logical())
} else {
  detected_overlaps <-
    viral_predictions %>%
    group_by(contig_id) %>%
    group_split() %>%
    lapply(function(tmp) {
      if (nrow(tmp) < 2) {
        ## Only one viral prediction tool used, or no predictions of high enough quality
        tmp
      } else {
        ## More than one prediction tool used: Identify overlapping predictions
        tmp %>%
          group_by(vir_id_tool) %>%
          mutate(overlap_start = sapply(start, function(x) {
            tmp_ <-
              c(tmp$id[tmp$start <= x & tmp$end >= x &
                         !tmp$vir_id_tool %in% vir_id_tool[1]],
                tmp$id[tmp$start == x &
                         tmp$vir_id_tool %in% vir_id_tool[1]]) %>%
              sort() %>%
              unique() %>%
              paste(collapse=":")
            ifelse(tmp_ == "" | !str_detect(tmp_, ":"), NA, tmp_)
          }),
          overlap_end = sapply(end, function(x) {
            tmp_ <-
              c(tmp$id [tmp$start <= x & tmp$end >= x &
                          !tmp$vir_id_tool %in% vir_id_tool[1]],
                tmp$id[tmp$end == x &
                         tmp$vir_id_tool %in% vir_id_tool[1]]) %>%
              sort() %>%
              unique() %>%
              paste(collapse = ":")
            ifelse(tmp_ == "" | !str_detect(tmp_, ":"), NA, tmp_)
          })) %>%
          ungroup() %>%
          pivot_longer(c(overlap_start, overlap_end)) %>%
          group_by(vir_id_name) %>%
          filter((any(is.na(value)) & any(!is.na(value)) & !is.na(value)) | ## Remove those with overlaps in the other end
                   all(is.na(value)) | ## And keep any with either no overlaps
                   all(!is.na(value))) %>% ## Or overlaps in both ends
          ungroup() %>%
          distinct(vir_id_name, value, .keep_all = TRUE) %>%
          select(-name) %>%
          group_by(vir_id_name) %>%
          group_split() %>%
          lapply(function(su){
            if (all(is.na(su$value))) {
              su
            } else {
              ids <- str_split(su$value, ":")[[1]]
              su <-
                su %>%
                mutate(full_overlap_length = (min(tmp$end[ tmp$id %in% ids]) - max(tmp$start[ tmp$id %in% ids])),
                       full_length = (max(tmp$end[ tmp$id %in% ids]) - min(tmp$start[ tmp$id %in% ids])),
                       full_overlap_frac = full_overlap_length/full_length)
              if (length(ids == 2)) {
                tmp_range = c(tmp$start[ tmp$id %in% ids & !tmp$id %in% su$id],
                              tmp$end[ tmp$id %in% ids & !tmp$id %in% su$id])
                su %>%
                  mutate(frac_overlap = full_overlap_length/(tmp_range[2]-tmp_range[1]))
              } else {
                su %>%
                  mutate(frac_overlap = NA)
              }
            }
          }) %>%
          bind_rows()
      }
    }) %>%
    bind_rows()
}



if (length(unique(viral_predictions$vir_id_tool)) > 1) {
  ## Select which of the overlapping viral predictions to use
  set.seed(1)
  selected_predictions <-
    detected_overlaps %>%
    filter(!is.na(value)) %>%
    filter(full_overlap_frac > full_overlap_frac_threshold) %>%
    pull(value) %>%
    unique() %>%
    lapply(function(overlapping_preds) {
      tmp_o_p <- strsplit(overlapping_preds, ":")[[1]]
      tmp_overlapping_preds <-
        detected_overlaps %>%
        filter(id %in% tmp_o_p,
               checkv_quality %in% accepted_quals) %>%
        mutate(n_pred = length(tmp_o_p)) %>%
        mutate(length = end - start + 1)

      tmp_overlapping_preds <-
        tmp_overlapping_preds %>%
        mutate(value = case_when(is.na(value) ~ value[!is.na(value)][1], TRUE ~ value)) %>%
        mutate(full_overlap_length = case_when(is.na(full_overlap_length) ~ full_overlap_length[!is.na(full_overlap_length)][1],
                                               TRUE ~ full_overlap_length)) %>%
        mutate(full_length = case_when(is.na(full_length) ~ length,
                                       TRUE ~ full_length)) %>%
        mutate(full_overlap_frac = case_when(is.na(full_overlap_frac) ~ full_length/length,
                                             TRUE ~ full_overlap_frac)) %>%
        mutate(frac_overlap = case_when(is.na(frac_overlap) ~ 1,
                                        TRUE ~ frac_overlap))

      if (overlap_selection == "min_length" & nrow(tmp_overlapping_preds) > 0) {
        tmp_overlapping_preds %>%
          slice_sample(n = nrow(.)) %>%
          slice_min(length, with_ties = FALSE)
      } else if (overlap_selection == "max_length" & nrow(tmp_overlapping_preds) > 0) {
        tmp_overlapping_preds %>%
          slice_sample(n = nrow(.)) %>%
          slice_max(length, with_ties = FALSE)
      }
    }) %>%
    bind_rows()

  if (tool_combination == "any") {
    selected_predictions <-
      selected_predictions %>%
      bind_rows(detected_overlaps %>%
                  filter(is.na(value) | full_overlap_frac < full_overlap_frac_threshold))
  }
} else {
  if (nrow(detected_overlaps > 0)) {
    selected_predictions <-
      detected_overlaps %>%
      mutate(value = id)
  } else {
    selected_predictions <-
      viral_predictions %>%
      mutate(value = "")
  }
  
}

## Rename viral predictions
set.seed(1)
tmp_renaming <-
  selected_predictions %>%
  mutate(virus_id = paste(sample_id, "virus", sample(10^(7):(10^8-1), nrow(.)), sep = "_"))

## Write mapping file
tmp_renaming %>%
  select(virus_id, representative = id, predicted_virus = value) %>%
  separate_rows(predicted_virus, sep = ":", ) %>%
  mutate(representative_prediction_tool = str_extract(representative, "^(genomad|virsorter)"),
         prediction_tool = str_extract(predicted_virus, "^(genomad|virsorter)"),
         across(.cols = c(representative, predicted_virus), .fns = function(x) str_remove(x, "^(genomad|virsorter)_"))) %>%
  write_tsv(snakemake@output[["representative_selection"]])

## Write quality details
tmp_renaming %>%
  select(virus_id, everything()) %>%
  rename(representative = id, predicted_virus = value) %>%
  mutate(across(.cols = c(representative, predicted_virus), .fns = function(x) str_remove(x, "^(genomad|virsorter)_"))) %>%
  write_tsv(snakemake@output[["gathered_qual"]])

## Write bed
tmp_renaming %>%
  select(contig_id, start, end, virus_id) %>%
  mutate(across(c(start, end), .fns = function(x) x-1)) %>%
  write_tsv(snakemake@output[["regions"]], col_names = FALSE)
