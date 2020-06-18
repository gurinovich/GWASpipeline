suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
annot <- args[2]
nullmod <- args[3]
test <- args[4]
result.file <- args[5]

####Open GDS
gds <- seqOpen(gds.file)

#save vector with snps rs ids
snps <- data.frame(variant.id=seqGetData(gds, "variant.id"), snpID=seqGetData(gds, "annotation/id"))

####Create a SeqVarData object
annot <- readRDS(annot)
seqData <- SeqVarData(gds, sampleData=annot)

####Null model
nullmod <- readRDS(nullmod)

####GWAS
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, test=test, imputed=T, verbose=FALSE)
assoc <- left_join(assoc, snps, by="variant.id")

fwrite(assoc, file = result.file, quote=FALSE, row.names=FALSE)

