suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

path1 <- args[1]
path2 <- args[2]

chr <- 1
res <- fread(paste0(path1, "chr_", chr, ".csv"))
caf <- fread(paste0(path2, "CAFs.chr", chr, ".csv"))
res <- res %>%
  inner_join(caf, by = c("chr" = "CHROM", "pos" = "POS")) %>%
  as.tbl()
write.csv(res, paste0("chr", chr, ".csv"), row.names = FALSE)
result <- res

#chr <- 2
for (chr in 2:22) {
  res <- fread(paste0(path1, "chr_", chr, ".csv"))
  caf <- fread(paste0(path2, "CAFs.chr", chr, ".csv"))
  res <- res %>%
    inner_join(caf, by = c("chr" = "CHROM", "pos" = "POS")) %>%
    as.tbl()
  write.csv(res, paste0("chr", chr, ".csv"), row.names = FALSE)
  results <- bind_rows(result, res)
}

write.csv(results, paste0("all_chr", ".csv"), row.names = FALSE)
