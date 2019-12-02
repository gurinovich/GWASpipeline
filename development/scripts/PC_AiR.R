library(SeqArray)
library(GENESIS)
library(Biobase)
library(SeqVarTools)
library(dplyr)
library(SNPRelate)
library(ggplot2)
library(data.table)

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
#snpset.file <- NA
snpset.file <- args[5+i]

# print(gds.file)
# print(pheno.file)
# print(phenotypes)
# print(num_covariates)
# print(covariates)
# print(snpset.file)


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

####LD pruning to get variant set
snp.dat <- data.frame(variant.id = seqGetData(gds, "variant.id"), chr = seqGetData(gds, "chromosome"), pos = seqGetData(gds, "position"))
snp.dat$chr_pos <- paste0(snp.dat$chr, ":", snp.dat$pos)
if(is.na(snpset.file)){
  snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e7, ld.threshold=sqrt(0.1))
  pruned <- unlist(snpset, use.names=FALSE)
  saveRDS(pruned, "./data/pruned.rds")
  pruned.dat <- snp.dat[snp.dat$variant.id%in%pruned,]
  fwrite(pruned.dat[,c("chr", "pos")], "./data/snpset.txt", quote=FALSE, row.names=FALSE, )
}else{
  snpset.dat <- fread(snpset.file,stringsAsFactors=F,header=T,na.strings=c(NA,""))
  snpset.dat$chr_pos <- paste0(snpset.dat$chr, ":", snpset.dat$pos)
  snp.intersect.dat <- snp.dat[snp.dat$chr_pos%in% snpset.dat$chr_pos,]
  pruned <- unlist(snp.intersect.dat$variant.id)
}

####KING
king <- snpgdsIBDKING(gds, sample.id=analysis.sample.id, snp.id=pruned, verbose=FALSE)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)
saveRDS(king,"./data/king.rds")

####PC-AiR
pcs <- pcair(seqData, kinobj=kingMat, kin.thresh=2^(-11/2),
                      divobj=kingMat, div.thresh=-2^(-11/2),
             sample.include=analysis.sample.id,
             snp.include=pruned)
saveRDS(pcs,"./data/pcs.rds")

pc.df <- as.data.frame(pcs$vectors)
names(pc.df) <- paste0("PC", 1:ncol(pcs$vectors))
pc.df$sample.id <- row.names(pcs$vectors)
pc.df <- left_join(pData(annot), pc.df, by="sample.id")
saveRDS(pc.df,"./data/pc.df.rds")

png("./figures/PC1vsPC2.png")
ggplot(pc.df, aes(PC1, PC2)) + geom_point()
dev.off()

png("./figures/PC3vsPC4.png")
ggplot(pc.df, aes(PC3, PC4)) + geom_point()
dev.off()

