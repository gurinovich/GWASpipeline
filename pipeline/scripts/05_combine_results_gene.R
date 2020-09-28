#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path1 <- args[1]

sink('combine_results.log', append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

chr <- 1
if("chr_1.csv"%in%dir(path1)){
	res <- fread(paste0(path1, "chr_", chr, ".csv"))
    result <- res
}else{
	result <- matrix()
}

#chr <- 2
for (chr in 2:22) {
  if(paste0("chr_", chr, ".csv")%in%dir(path1)){
  	res <- fread(paste0(path1, "chr_", chr, ".csv"))
  	result <- bind_rows(result, res)}
}

result <- result %>%
  rename(pval = contains("pval"))

result <- result %>%
	select ("chr", "start", "end", "width", "strand", "gene_id", "n.site", "n.alt", "n.sample.alt", contains("Score"), contains("Wald"), "pval")

fwrite(result[order(result$pval),], paste0("all_chr.csv"), row.names = FALSE)

date()
sink()
