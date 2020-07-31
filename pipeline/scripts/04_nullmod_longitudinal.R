#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
phenotypes <- args[1]
covariates <- unlist(strsplit(args[2], ","))
model <- args[3]
analysis.sample.id <- args[4]
pc_df <- args[5]
grm <- args[6]

sink("nullmod_longitudinal.log", append=FALSE, split=TRUE)
date()
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GMMAT))

####PCA
pc.df <- readRDS(pc_df)
analysis.sample.id <- readRDS(analysis.sample.id)
pc.df <- pc.df[pc.df$sample.id%in%analysis.sample.id,]

####covariance matrix from pcrelate output
grm <- readRDS(grm)

####Mixed effect null model

if(model=="linear"){
	model.switch <- gaussian(link = "identity")
}
if(model=="logistic"){
	model.switch <- binomial(link = "logit")
}

fix.eff=paste(phenotypes,"~ 1")
if(!is.null(covariates)){
	for(covi in covariates)
		fix.eff=paste(fix.eff,"+",covi)
}
fix.eff=formula(fix.eff)

cat("\n####glmmkin starts\n")
nullmod <- glmmkin(fix.eff, data=pc.df, kins=as.matrix(grm), id="sample.id", family = model.switch)
cat("####glmmkin ends\n\n")

saveRDS(nullmod,"nullmod_longitudinal.rds")

date()
sink()
