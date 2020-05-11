suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SNPRelate))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
pheno.file <- args[2]
phenotypes <- args[3]
covariates <- unlist(strsplit(args[4], ","))
pc_df <- args[5]
grm <- args[6]
nullmod <- args[7]
test <- args[8]
result.file <- args[9]

####Open GDS
gds <- seqOpen(gds.file)

#save vector with snps rs ids
snps <- seqGetData(gds, "annotation/id")

####Create a SeqVarData object
pheno.dat <- read.csv(pheno.file,stringsAsFactors=F,header=T)
pheno.dat$sample.id <- as.character(pheno.dat$ID)

for(i in 1:32){
if(sum(colnames(pheno.dat)==paste0("PC",i))==1){
colnames(pheno.dat)[colnames(pheno.dat)==paste0("PC",i)] <- paste0("PC",i,".pheno")
}
}

gds.sample.id <- data.frame(sample.id= seqGetData(gds, "sample.id"),stringsAsFactors=F)

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to include : ",length(analysis.sample.id)))

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####PCA
pc.df <- readRDS(pc_df)

####add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

####covariance matrix from pcrelate output
grm <- readRDS(grm)

####Null model
nullmod <- readRDS("nullmod.rds")

####GWAS
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod, test=test, imputed=T, verbose=FALSE)
assoc$variant.id <- snps

assoc <- assoc %>%
  select(-4) %>%
  rename(rs.id = variant.id) %>%
  arrange(Score.pval)

write.csv(assoc, file = result.file, quote=FALSE, row.names=FALSE)

