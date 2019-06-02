#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SeqArray))

args <- commandArgs(trailingOnly=TRUE)

seqMerge(args, "merged.gds")
