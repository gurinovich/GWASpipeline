date()

n.args <- length(commandArgs())
opt.file <- commandArgs()[n.args-1]
gds.file <- commandArgs()[n.args]
source(opt.file)

####Open GDS
library(SeqArray)
gds <- seqOpen(gds.file)

####Create a SeqVarData object
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(dplyr)

pheno.dat <- read.csv(pheno.file,stringsAsFactors=F,header=T,na.strings=c(NA,""))
pheno.dat$sample.id <- paste0(pheno.dat$FID,"_",pheno.dat$IID)
colnames(pheno.dat)[colnames(pheno.dat)=="PC1"] <- "PC1.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC2"] <- "PC2.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC3"] <- "PC3.pheno"
colnames(pheno.dat)[colnames(pheno.dat)=="PC4"] <- "PC4.pheno"
pheno.dat$sample.id <- as.character(pheno.dat$sample.id)

gds.sample.id <- data.frame(sample.id= seqGetData(gds, "sample.id"),stringsAsFactors=F)
annot <- left_join(gds.sample.id, pheno.dat, by="sample.id")
annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)], "PC1.pheno","PC2.pheno","PC3.pheno","PC4.pheno")]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to include : ",length(analysis.sample.id)))

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####Read in pruned set
library(SNPRelate)
pruned <- readRDS("pruned.rds")

####KING
king <- readRDS("king.rds")
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)

####PC-AiR
pcs <- readRDS("pcs.rds")

library(dplyr)
library(ggplot2)
pc.df <- readRDS("pc.df.rds")

####PC-Relate
seqSetFilter(seqData, variant.id=pruned)
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)
pcrel <- pcrelate(iterator, pcs=pcs$vectors[,1:2], training.set=pcs$unrels)
seqResetFilter(seqData, verbose=FALSE)
saveRDS(pcrel,"pcrel.rds")

kinship <- pcrel$kinBtwn
png("kinship.png")
ggplot(kinship, aes(k0, kin)) +
    geom_hline(yintercept=2^(-seq(3,9,2)/2), linetype="dashed", color="grey") +
    geom_point(alpha=0.5) +
    ylab("kinship estimate") +
    theme_bw()
dev.off()

####add PCs to sample annotation in SeqVarData object
annot <- AnnotatedDataFrame(pc.df)
sampleData(seqData) <- annot

####covariance matrix from pcrelate output
grm <- pcrelateToMatrix(pcrel, scaleKin=2)
saveRDS(grm,"grm.rds")

date()
