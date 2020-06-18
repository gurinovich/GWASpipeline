suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

args = commandArgs(trailingOnly=TRUE)
result.file <- args[1]
caf.file <- args[2] 
model <- args[3]
combine.file <- args[4]


result.dat <- fread(result.file, stringsAsFactors = F, header = T)
caf.dat <- fread(caf.file, stringsAsFactors = F, header = T)

combine.dat <- left_join(result.dat, caf.dat, by=c("chr","pos"))
maf <- ifelse(combine.dat$caf>0.5, 1-combine.dat$caf, combine.dat$caf)
out.dat <- combine.dat[maf*2*combine.dat$n.obs>0.99,]

if (model == "logistic") {
	case.maf <- ifelse(combine.dat$case.caf>0.5, 1-combine.dat$case.caf, combine.dat$case.caf)
	control.maf <- ifelse(combine.dat$control.caf>0.5, 1-combine.dat$control.caf, combine.dat$control.caf)

	out.dat <- out.dat[(case.maf*2*out.dat$n.case>0.99)&(control.maf*2*out.dat$n.control>0.99),]
}

fwrite(out.dat, combine.file, quote=F, row.names=F)
