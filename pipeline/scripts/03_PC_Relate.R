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
snpset.file <- args[5]
pruned <- args[6]
king <- args[7]
pcs <- args[8]
pc_df <- args[9]

####Open GDS
gds <- seqOpen(gds.file)

####Create a SeqVarData object
pheno.dat <- read.csv(pheno.file,stringsAsFactors=F,header=T)
pheno.dat$sample.id <- as.character(pheno.dat$ID)

for(i in 1:32){
if(sum(colnames(pheno.dat)==paste0("PC",i))==1){
colnames(pheno.dat)[colnames(pheno.dat)==paste0("PC",i)] <- paste0("PC",i,".pheno")
}
}

gds.sample.id <- data.frame(sample.id= seqGetData(gds, "sample.id"),stringsAsFactors=F)

annot <- left_join(gds.sample.id, pheno.dat, by="sample.id")

annot <- annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])]

analysis.sample.id <- na.omit(annot[,c("sample.id",colnames(annot)[colnames(annot)%in%c(phenotypes,covariates)])])$sample.id
print(paste0("number of individuals to include : ",length(analysis.sample.id)))
saveRDS(analysis.sample.id, "analysis.sample.id.rds")

metadata <- data.frame(labelDescription=colnames(annot),row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)
all.equal(annot$sample.id,seqGetData(gds, "sample.id"))
seqData <- SeqVarData(gds, sampleData=annot)

####Read in pruned set
pruned <- readRDS(pruned)

####KING
king <- readRDS(king)
kingMat <- king$kinship
dimnames(kingMat) <- list(king$sample.id, king$sample.id)

####PC-AiR
pcs <- readRDS(pcs)

pc.df <- readRDS(pc_df)

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

