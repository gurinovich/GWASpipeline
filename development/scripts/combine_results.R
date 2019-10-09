library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#data.dir <- "./data/"
data.dir <- args[1]
#res.file <- "./results/result_file"
res.file <- args[2]

chr <- 1
res <- fread(paste0(res.file, chr, ".csv"))
caf <- fread(paste0(data.dir, "CAFs.chr", chr, ".csv"))
res <- res %>%
  inner_join(caf, by = c("chr" = "CHROM", "pos" = "POS")) %>%
  as.tbl()
write.csv(res, paste0(res.file, chr, ".csv"), row.names = FALSE)
result <- res

#chr <- 2
for (chr in 2:22) {
  res <- fread(paste0(res.file, chr, ".csv"))
  caf <- fread(paste0(data.dir, "CAFs.chr", chr, ".csv"))
  res <- res %>%
    inner_join(caf, by = c("chr" = "CHROM", "pos" = "POS")) %>%
    as.tbl()
  write.csv(res, paste0(res.file, chr, ".csv"), row.names = FALSE)
  results <- bind_rows(result, res)
}

write.csv(results, paste0(res.file, ".csv"), row.names = FALSE)
