#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
gds.file <- args[1]
pheno.file <- args[2]
analysis.sample.id <- args[3]
model <- args[4]
phenotypes <- args[5]  
out.file <- args[6]
log.file <- args[7]

sink(log.file, append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(reshape2))

####Open GDS
gds <- seqOpen(gds.file)

####group
pheno.dat <- read.csv(pheno.file, stringsAsFactors=F, header=T)
pheno.dat$sample.id <- as.character(pheno.dat[,1])
analysis.sample.id <- readRDS(analysis.sample.id)
pheno.dat <- pheno.dat[pheno.dat$sample.id%in%analysis.sample.id,]

table(pheno.dat$group)

if(sum(colnames(pheno.dat)%in%c("group"))==1){
	group_names <- names(table(pheno.dat$group))
	out.caf <- c()
	out.dosage <- c()
	for(i in 1:length(group_names)){
		group.id <- pheno.dat[pheno.dat$group==group_names[i],]$sample.id
		seqSetFilter(gds, sample.id = group.id)
		genotype.dat <- seqGetData(gds, "genotype")
		genotype <- genotype.dat[1,,]+genotype.dat[2,,]
		caf <- apply(genotype,2,sum)/dim(genotype)[1]/2
		dosage.dat <- seqGetData(gds, "annotation/format/DS")$data
		dosage <- apply(dosage.dat,2,sum)/dim(dosage.dat)[1]/2
		seqResetFilter(gds)
		out.caf <- cbind(out.caf, caf)
		out.dosage <- cbind(out.dosage, dosage)
	}
	out.caf <- as.data.frame(out.caf)
	out.dosage <- as.data.frame(out.dosage)
	colnames(out.caf) <- paste0(group_names, ".caf")
	colnames(out.dosage) <- paste0(group_names, ".dosage")
	out.caf.dosage <- cbind(out.caf, out.dosage)
	out.caf.dosage$chr <- seqGetData(gds,"chromosome")
	out.caf.dosage$pos <- seqGetData(gds,"position")
}else{
	group_names <- c()
	print("group variable not found")
	}

if (model == "logistic") {
	case.id <- pheno.dat[pheno.dat[,phenotypes]==1,]$sample.id
	seqSetFilter(gds, sample.id = case.id)
	genotype.dat <- seqGetData(gds,"genotype")
	genotype <- genotype.dat[1,,]+genotype.dat[2,,]
	case.caf <- apply(genotype,2,sum)/dim(genotype)[1]/2
	dosage.dat <- seqGetData(gds, "annotation/format/DS")$data
	case.dosage <- apply(dosage.dat,2,sum)/dim(dosage.dat)[1]/2
	seqResetFilter(gds)
	
	control.id <- pheno.dat[pheno.dat[,phenotypes]==0,]$sample.id
	seqSetFilter(gds, sample.id = control.id)
	genotype.dat <- seqGetData(gds,"genotype")
	genotype <- genotype.dat[1,,]+genotype.dat[2,,]
	control.caf <- apply(genotype,2,sum)/dim(genotype)[1]/2
	dosage.dat <- seqGetData(gds, "annotation/format/DS")$data
	control.dosage <- apply(dosage.dat,2,sum)/dim(dosage.dat)[1]/2
	seqResetFilter(gds)

	out.caf.dosage$n.case <- length(case.id)
	out.caf.dosage$n.control <- length(control.id)
	out.caf.dosage$case.caf <- case.caf
	out.caf.dosage$control.caf <- control.caf
	out.caf.dosage$case.dosage <- case.dosage
	out.caf.dosage$control.dosage <- control.dosage
}


seqSetFilter(gds, sample.id = analysis.sample.id)
genotype.dat <- seqGetData(gds, "genotype")
genotype <- genotype.dat[1,,]+genotype.dat[2,,]
caf <- apply(genotype,2,sum)/dim(genotype)[1]/2
dosage.dat <- seqGetData(gds, "annotation/format/DS")$data
dosage <- apply(dosage.dat,2,sum)/dim(dosage.dat)[1]/2
seqResetFilter(gds)

out.caf.dosage$caf <- caf
out.caf.dosage$dosage <- dosage

allele <- seqGetData(gds, "allele")
all.dat <- colsplit(allele, ",", c("REF", "ALT"))

out <- cbind(out.caf.dosage, all.dat)

write.csv(out, file = out.file, quote=FALSE, row.names=FALSE)

date()
sink()
