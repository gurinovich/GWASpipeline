suppressPackageStartupMessages(library(dplyr))

args = commandArgs(trailingOnly=TRUE)
pheno.file <- args[1]
model <- args[2]  
analysis <- args[3]

pheno <- read.csv(pheno.file, stringsAsFactors = F)
analysis.sample.id <- readRDS(analysis)
pheno <- pheno[pheno$ID%in%analysis.sample.id,]
pheno <- as.tbl(pheno)

group_names <- names(table(pheno$group))

for (g in group_names) {
  pheno %>%
    filter(group == g) %>%
    select(ID) %>%
    write.table(file = paste0(g, ".txt"), quote = F, row.names = F, col.names = F)
}

# add cases and controls groups if model == "logistic"
if (model == "logistic") {
  pheno %>%
    filter(outcome == 0) %>%
    select(ID) %>%
    write.table(file = "controls.txt", quote = F, row.names = F, col.names = F)
  
  pheno %>%
    filter(outcome == 1) %>%
    select(ID) %>%
    write.table(file = "cases.txt", quote = F, row.names = F, col.names = F)
}
