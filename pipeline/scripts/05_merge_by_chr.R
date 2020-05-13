suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))

args = commandArgs(trailingOnly=TRUE)
pheno.file <- args[1]
model <- args[2] 
path1 <- args[3]
path2 <- args[4]

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

#merge files per choromosome:

for(chr in 1:22) {
  gr <- groups[1]
  genAF1 <- fread(paste0(path1, gr, ".txt.", chr, ".afreq"), stringsAsFactors = F, header = T)
  genAF1 <- as.tbl(genAF1)
  names(genAF1)[c(1, 5:6)] <- c("CHROM", paste0("CAF_geno_", gr) , paste0(names(genAF1)[6], "_", gr))
  genAF2 <- fread(paste0(path2, gr, ".", chr, ".CAF_dosages.csv"), stringsAsFactors = F, header = T)
  genAF2 <- as.tbl(genAF2)
    names(genAF2)[7] <- paste0(names(genAF2)[7], "_", gr)
  data1 <- genAF1 %>%
    inner_join(genAF2, by = c("CHROM", "ID", "REF", "ALT")) %>%
    select(CHROM, POS, ID:ALT, FILTER, starts_with("OBS_CT_"), starts_with("CAF_geno_"), starts_with("CAF_dos_"))
  data2 <- data1
#gr <- groups[2]
  for (gr in groups[2:length(groups)]) {
    genAF1 <- fread(paste0(path1, gr, ".txt.", chr, ".afreq"), stringsAsFactors = F, header = T)
    genAF1 <- as.tbl(genAF1)
    names(genAF1)[c(1, 5:6)] <- c("CHROM", paste0("CAF_geno_", gr) , paste0(names(genAF1)[6], "_", gr))
    genAF2 <- fread(paste0(path2, gr, ".", chr, ".CAF_dosages.csv"), stringsAsFactors = F, header = T)
    genAF2 <- as.tbl(genAF2)
    names(genAF2)[7] <- paste0(names(genAF2)[7], "_", gr)
    data1 <- genAF1 %>%
      inner_join(genAF2, by = c("CHROM", "ID", "REF", "ALT")) %>%
      select(CHROM, POS, ID:ALT, FILTER, starts_with("OBS_CT_"), starts_with("CAF_geno_"), starts_with("CAF_dos_"))
    data2 <- bind_cols(data2, data1[7:9])  
  }
  write.csv(data2, file = paste0("CAFs.chr", chr, ".csv"), row.names = F)
}
