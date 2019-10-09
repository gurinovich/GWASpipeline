library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#pheno.file <- "./data/pheno_file.csv"
pheno.file <- args[1]
  
pheno <- read.csv(pheno.file, stringsAsFactors = F)
pheno <- as.tbl(pheno)

group_names <- names(table(pheno$group))

# g <- group_names[1]
for (g in group_names) {
  pheno %>%
    filter(group == g) %>%
    select(ID) %>%
    write.table(file = paste0(dirname(pheno.file), "/", g, ".txt"), quote = F, row.names = F, col.names = F)
}

write.table(as.data.frame(group_names), file = paste0(dirname(pheno.file), "/groups.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
