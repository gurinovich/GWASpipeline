library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

#data.dir <- "./data/"
#data.dir <- args[1]
###res.dir = "../results/"
res.dir = args[1]
###res.file <- "../results/result_file"
#new.res <- "./result_file"
res.file <- args[2]

##chr7_EL_GWAS.hg19_multianno.csv
## subset fields by annot.fields <- subset(combined[,c("ill_chrom","ill_pos","ill_pos","B_allele_freq","log_r_ratio")])

chr <- 1
res <- fread(paste0(res.file, chr, ".csv"))
annot <- fread(paste0(res.dir, "chr", chr, "_EL_GWAS.hg19_multianno.csv"))
annot.sub <- subset(annot[,c("Chr","Start","Func.refGene","Gene.refGene")])
res <- res %>%
  inner_join(annot.sub, by = c("chr" = "Chr", "pos" = "Start")) %>%
  as.tbl()
write.csv(res, paste0(res.file, chr, ".csv"), row.names = FALSE)
#write.csv(res, paste0( new.res,chr,".csv"), row.names = FALSE)
result <- res

#chr <- 2
for (chr in 2:22) {
  res <- fread(paste0(res.file, chr, ".csv"))
annot <- fread(paste0(res.dir, "chr", chr, "_EL_GWAS.hg19_multianno.csv"))
annot.sub <- subset(annot[,c("Chr","Start","Func.refGene","Gene.refGene")])
#  caf <- fread(paste0(data.dir, "CAFs.chr", chr, ".csv"))
  res <- res %>%
         inner_join(annot.sub, by = c("chr" = "Chr", "pos" = "Start")) %>%
         as.tbl()
  #write.csv(res, paste0(new.res, chr, ".csv"), row.names = FALSE)
  write.csv(res, paste0(res.file, chr, ".csv"), row.names = FALSE)
  results <- bind_rows(result, res)
}

write.csv(results, paste0(res.file, ".csv"), row.names = FALSE)
