suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))

args = commandArgs(trailingOnly=TRUE)
pheno.file <- args[1]
model <- args[2] 
path <- args[3]

pheno <- read.csv(pheno.file, stringsAsFactors = F)
pheno <- as.tbl(pheno)
group_names <- names(table(pheno$group))

if (model == "logistic") {
	group_names <- c(group_names, "cases", "controls")
}
	
caf <- function(x) {
  sum(x) / (2*length(x))
}

groups <- as.vector(group_names)

for (g in groups) {
  for (c in 1:22) {
    dosf <- fread(paste0(path, g, ".txt.", c, ".dosages"))
    dosf <- as.tbl(dosf)
    names(dosf)[1:7] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
    doscaf <- map_dbl(as_data_frame(t(dosf[8:ncol(dosf)])), caf)
 #   head(as_data_frame(doscaf))
    dosf <- dosf %>%
      select(-QUAL, -V8:-ncol(.)) %>%
      mutate(CAF_dos = doscaf)
  write.csv(dosf, file = paste0(g, ".", c, ".CAF_dosages.csv"), row.names = F)
  }
}
