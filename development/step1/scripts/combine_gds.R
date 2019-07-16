library(SeqArray)

args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
  
seqMerge(paste0(gds.file,c(1:22),".gds"), paste0(gds.file,".gds"))
