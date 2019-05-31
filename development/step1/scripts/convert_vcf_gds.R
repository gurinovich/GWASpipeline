library(SeqArray)

vcf.file <- args[1]
gds.file <- args[2]

seqVCF2GDS(vcf.file, gds.file, verbose=FALSE)
  
