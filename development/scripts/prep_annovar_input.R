library(dplyr)

args<-commandArgs(TRUE)
chr <- args[1]

res <-  paste("./results/result_file",chr,sep="",".csv")
snps <- read.csv(res)
snps <- as.tbl(snps)
snps

snps <- snps %>%
  mutate(pos2 = pos, A1 = 0, A2 = 0) %>%
  select(chr, pos, pos2, A1, A2, ID)
av <- paste("./tmp/result_file",chr,sep="","_snps_input.txt")
write.table(snps, av, quote = F, sep = "\t", row.names = F, col.names = F)
