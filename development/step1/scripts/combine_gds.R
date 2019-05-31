library(SeqArray)

vcf.file <- args[1]
gds.file <- args[2]
  
seqMerge(paste0("vcf.file",c(1:22),".gds"), paste0(gds.file,".gds"))
