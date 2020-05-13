#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SeqArray))

args <- commandArgs(trailingOnly=TRUE)

seqMerge(paste0("chr_",1:22,".gds"), "merged.gds")
