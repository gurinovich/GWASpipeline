library(SeqArray)
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(dplyr)
library(SNPRelate)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
#gds.file <- "./data/gds_file"
gds.file <- args[1]
#pheno.file <- "./data/pheno_file.csv"
pheno.file <- args[2]
#phenotypes <- "outcome"
phenotypes <- args[3]
#num_covariates <- 6
num_covariates <- args[4]
#covariates <- c("age", "sex", "PC1", "PC2", "PC3", "PC4")
covariates <- vector(mode = "character", length = num_covariates)
for (i in 1:num_covariates) {
  covariates[i] <- args[4+i]
}
model <- args[5+i]

####Open GDS
gds <- seqOpen(paste0(gds.file,".gds"))

####Create a SeqVarData object
pheno.dat <- read.csv(pheno.file,stringsAsFactors=F,header=T,na.strings=c(NA,""))
pheno.dat$sample.id <- as.character(pheno.dat$ID)
colnames(pheno.dat)[colnames(pheno.dat)=="PC1"] <- "PC1.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC2"] <- "PC2.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC3"] <- "PC3.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC4"] <- "PC4.pheno"

gds.sample.id <- data.frame(sample.id= seqGetData(gds, "sample.id"),stringsAsFactors=F)

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)], "PC1.pheno","PC2.pheno","PC3.pheno","PC4.pheno")]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to include : ",length(analysis.sample.id)))

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####PCA
pc.df <- readRDS("./data/pc.df.rds")

####add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

####covariance matrix from pcrelate output
grm <- readRDS("./data/grm.rds")

####Null model
if(model=="linear"){
	model.switch <- "gaussian"
}
if(model=="logistic"){
	model.switch <- "binomial")
}
nullmod <- fitNullModel(seqData, outcome=phenotypes, 
                        covars=covariates,
                        cov.mat=grm,
                        family=model.switch, verbose=T)

saveRDS(nullmod,"./data/nullmod.rds")

