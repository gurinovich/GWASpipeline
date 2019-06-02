#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(SeqArray))

args <- commandArgs(trailingOnly=TRUE)

vcf.file <- args[1]
gds.file <- args[2]

seqVCF2GDS(vcf.file, gds.file, verbose=FALSE)
