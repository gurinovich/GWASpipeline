suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)

chr = args[1]
ref = args[2]

chr <- 1
res <- fread(paste0("chr", chr, ".csv"))
annot <- fread(paste0("chr", chr, "_annotation.", ref,"_multianno.csv"))
annot.sub <- subset(annot[,c("Chr","Start","Func.refGene","Gene.refGene")])
res <- res %>%
  inner_join(annot.sub, by = c("chr" = "Chr", "pos" = "Start")) %>%
  as.tbl()
write.csv(res, paste0("chr", chr, "_annotation.csv"), row.names = FALSE)
result <- res

#chr <- 2
for (chr in 2:22) {
res <- fread(paste0("chr", chr, ".csv"))
annot <- fread(paste0("chr", chr, "_annotation.", ref,"_multianno.csv"))
annot.sub <- subset(annot[,c("Chr","Start","Func.refGene","Gene.refGene")])
  res <- res %>%
         inner_join(annot.sub, by = c("chr" = "Chr", "pos" = "Start")) %>%
         as.tbl()
  write.csv(res, paste0("chr", chr, "_annotation.csv"), row.names = FALSE)
  results <- bind_rows(result, res)
}

write.csv(results, paste0("All_chr_annotation.csv"), row.names = FALSE)
