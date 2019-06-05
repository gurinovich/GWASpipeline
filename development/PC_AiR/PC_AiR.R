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

####LD pruning to get variant set
if(is.null(snpset)){
library(SNPRelate)
snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e7, 
                          ld.threshold=sqrt(0.1))
pruned <- unlist(snpset, use.names=FALSE)
saveRDS(pruned, "pruned.rds")
}else{
  pruned <- unlist(read.table(snpset.file))
}

####KING
king <- snpgdsIBDKING(gds, sample.id=analysis.sample.id, snp.id=pruned, verbose=FALSE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(king,"king.rds")

####PC-AiR
pcs <- pcair(seqData, kinobj=kingMat, kin.thresh=2^(-11/2),
                      divobj=kingMat, div.thresh=-2^(-11/2),
             sample.include=analysis.sample.id,
             snp.include=pruned)
saveRDS(pcs,"pcs.rds")

library(ggplot2)
pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pData(annot), pc.df, by="sample.id")
saveRDS(pc.df,"pc.df.rds")

png("PC1vsPC2.png")
ggplot(pc.df, aes(PC1, PC2)) + 
geom_point()
dev.off()

png("PC3vsPC4.png")
ggplot(pc.df, aes(PC3, PC4)) +
geom_point()
dev.off()

date()
