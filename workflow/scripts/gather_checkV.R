#!/usr/bin/Rscript


checkV_format <- do.call("rbind", lapply(snakemake@input[["gathered_qual"]], function(f) {
  if (!file.info(f)$size == 0) {
    read.delim(f, stringsAsFactors = FALSE)
  }
}))

new_headers <- c("Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Strain heterogeneity")

checkM_format <- do.call("rbind", lapply(snakemake@input[["gathered_qual"]], function(f) {
  if (!file.info(f)$size == 0) {
    tmp <- read.delim(f, stringsAsFactors = FALSE)
    if (nrow(tmp) > 0) {
      tmp$Completeness <- tmp$completeness
      ## Virus sequences extracted; quality sorted by completeness
      tmp$Contamination <- rep(0, nrow(tmp))
      tmp$`Bin Id` <- tmp$virus_id
      
      if (nrow(tmp) == 1) {
        tmp <- cbind(tmp, matrix(sapply(new_headers, function(x) rep(0,nrow(tmp))),
                                 nrow = 1, dimnames = list(NULL, new_headers)))
      } else {
        tmp <- cbind(tmp, sapply(new_headers, function(x) rep(0,nrow(tmp))))
      }
      tmp[ , c("Bin Id", "Marker lineage", "# genomes", "# markers", "# marker sets", "0", "1", "2", "3", "4", "5+", "Completeness", "Contamination", "Strain heterogeneity")]
    } else {
      tmp
    }
  }
}))

write.table(checkV_format, file = snakemake@output[["checkV_gathered"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(checkM_format, file = snakemake@output[["checkM_format"]], sep = "\t", quote = FALSE, row.names = FALSE)
