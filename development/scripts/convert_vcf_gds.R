library(SeqArray)

args = commandArgs(trailingOnly=TRUE)
#vcf.file <- "./data/vcf_file1_qc2.vcf.gz"
vcf.file <- args[1]
#gds.file <- "./data/gds_file1.gds"
gds.file <- args[2]

seqVCF2GDS(vcf.file, gds.file, verbose=FALSE)