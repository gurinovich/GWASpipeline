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

num_covariates <- args[4]

covariates <- vector(mode = "character", length = num_covariates)
for (i in 1:num_covariates) {
  covariates[i] <- args[4+i]
}
#snpset.file <- NA
snpset.file <- args[5+i]

####Open GDS
gds <- seqOpen(gds.file)

####Create a SeqVarData object
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat$ID)

for(i in 1:32){
if(!is.na(colnames(pheno.dat)==paste0("PC",i)){
colnames(pheno.dat)[colnames(pheno.dat)==paste0("PC",i)] <- paste0("PC",i,".pheno")
}

gds.sample.id <- data.frame(sample.id=seqGetData(gds, "sample.id"),stringsAsFactors=F)

annot <- left_join(gds.sample.id, pheno.dat)

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to include : ",length(analysis.sample.id)))
saveRDS(analysis.sample.id, "./data/analysis.sample.id.rds")

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
  saveRDS(pruned, "./data/pruned.rds")
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
