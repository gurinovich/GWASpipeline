#!/usr/bin/env Rscript
sink("merge_gds.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))

cat("\n####seqMerge starts\n")
seqMerge(paste0("chr_",1:22,".gds"), "merged.gds")
cat("####seqMerge ends\n\n")

date()
sink()
