library(dplyr)
library(data.table)
library(purrr)

args = commandArgs(trailingOnly=TRUE)

#work.dir <- "./data/"
work.dir <- args[1]

caf <- function(x) {
  sum(x) / (2*length(x))
}

#groups <- c("NECS.old", "NECS.new", "Illumina", "cases", "controls", "NECS.old.controls")
groups <- read.table(paste0(work.dir, "groups.txt"))
groups <- as.vector(groups$V1)

# g <- groups[1]
for (g in groups) {
# c <- 1
  for (c in 1:22) {
    dosf <- fread(paste0(work.dir, g, ".", c, ".dosages"))
    dosf <- as.tbl(dosf)
    names(dosf)[1:7] <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
    doscaf <- map_dbl(as_data_frame(t(dosf[8:ncol(dosf)])), caf)
 #   head(as_data_frame(doscaf))
    dosf <- dosf %>%
      select(-QUAL, -V8:-ncol(.)) %>%
      mutate(CAF_dos = doscaf)
  write.csv(dosf, file = paste0(work.dir, g, ".", c, ".CAF_dosages.csv"), row.names = F)
  }
}

#merge files per choromosome:

for(chr in 1:22) {
  gr <- groups[1]
  genAF1 <- fread(paste0(work.dir, gr, ".", chr, ".afreq"), stringsAsFactors = F, header = T)
  genAF1 <- as.tbl(genAF1)
  names(genAF1)[c(1, 5:6)] <- c("CHROM", paste0("CAF_geno_", gr) , paste0(names(genAF1)[6], "_", gr))
  genAF2 <- fread(paste0(work.dir, gr, ".", chr, ".CAF_dosages.csv"), stringsAsFactors = F, header = T)
  genAF2 <- as.tbl(genAF2)
    names(genAF2)[7] <- paste0(names(genAF2)[7], "_", gr)
  data1 <- genAF1 %>%
    inner_join(genAF2, by = c("CHROM", "ID", "REF", "ALT")) %>%
    select(CHROM, POS, ID:ALT, FILTER, starts_with("OBS_CT_"), starts_with("CAF_geno_"), starts_with("CAF_dos_"))
  data2 <- data1
#gr <- groups[2]
  for (gr in groups[2:length(groups)]) {
    genAF1 <- fread(paste0(work.dir, gr, ".", chr, ".afreq"), stringsAsFactors = F, header = T)
    genAF1 <- as.tbl(genAF1)
    names(genAF1)[c(1, 5:6)] <- c("CHROM", paste0("CAF_geno_", gr) , paste0(names(genAF1)[6], "_", gr))
    genAF2 <- fread(paste0(work.dir, gr, ".", chr, ".CAF_dosages.csv"), stringsAsFactors = F, header = T)
    genAF2 <- as.tbl(genAF2)
    names(genAF2)[7] <- paste0(names(genAF2)[7], "_", gr)
    data1 <- genAF1 %>%
      inner_join(genAF2, by = c("CHROM", "ID", "REF", "ALT")) %>%
      select(CHROM, POS, ID:ALT, FILTER, starts_with("OBS_CT_"), starts_with("CAF_geno_"), starts_with("CAF_dos_"))
    data2 <- bind_cols(data2, data1[7:9])  
  }
  write.csv(data2, file = paste0(work.dir, "CAFs.chr", chr, ".csv"), row.names = F)
}
