library(dplyr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
#res.file <- "./results/result_file"
res.file <- args[1]

res.df <- fread(paste0(res.file, "1.csv"))
res.df <- as.tbl(res.df)
res.df

#i <- 2
for (i in 2:22) {
  temp <- fread(paste0(res.file, i, ".csv"))
  temp <- as.tbl(temp)
  res.df <- res.df %>%
    bind_rows(temp)
}

rm(temp)

res.df <- res.df %>%
  arrange(Score.pval)

write.csv(res.df, file = paste0(res.file, ".csv"), quote=FALSE, row.names=FALSE)
