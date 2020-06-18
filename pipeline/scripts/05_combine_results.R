suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

path1 <- args[1]

chr <- 1
res <- fread(paste0(path1, "chr_", chr, "_caf_annotated.csv"))
result <- res

#chr <- 2
for (chr in 2:22) {
  res <- fread(paste0(path1, "chr_", chr, "_caf_annotated.csv"))
  result <- bind_rows(result, res)
}

fwrite(result, paste0("all_chr_caf_annotated.csv"), row.names = FALSE)
