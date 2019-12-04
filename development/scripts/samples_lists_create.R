library(dplyr)

args = commandArgs(trailingOnly=TRUE)
#pheno.file <- "./data/pheno_file_logistic.csv"
pheno.file <- args[1]
#model <- "logistic"
model <- args[2] 
  
pheno <- read.csv(pheno.file, stringsAsFactors = F)
pheno <- as.tbl(pheno)

group_names <- names(table(pheno$group))

# g <- group_names[1]
for (g in group_names) {
  pheno %>%
    filter(group == g) %>%
    select(ID) %>%
    write.table(file = paste0("./tmp/", g, ".txt"), quote = F, row.names = F, col.names = F)
}

# add cases and controls groups if model == "logistic"
if (model == "logistic") {
  pheno %>%
    filter(outcome == 0) %>%
    select(ID) %>%
    write.table(file = "./tmp/controls.txt", quote = F, row.names = F, col.names = F)
  
  pheno %>%
    filter(outcome == 1) %>%
    select(ID) %>%
    write.table(file = "./tmp/cases.txt", quote = F, row.names = F, col.names = F)
  
  group_names <- c("cases", "controls", group_names)
}

write.table(as.data.frame(group_names), file = "./tmp/groups.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
