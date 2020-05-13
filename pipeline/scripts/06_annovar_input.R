suppressPackageStartupMessages(library(dplyr))

args<-commandArgs(TRUE)
chr <- args[1]

snps <- read.csv(paste0("chr", chr, ".csv"))
snps <- as.tbl(snps)
snps

snps <- snps %>%
  mutate(pos2 = pos, A1 = 0, A2 = 0) %>%
  select(chr, pos, pos2, A1, A2, ID)

write.table(snps, paste0("chr",chr,"_snps_input.txt"), quote = F, sep = "\t", row.names = F, col.names = F)
