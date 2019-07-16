#!/bin/bash

module load R/3.5.1
module load gcc/7.4.0

source $1

#convert VCF files to GDS files
for chr in {1..22}
do
	Rscript ./scripts/convert_vcf_gds.R ${vcf_file}${chr}".vcf.gz" ${gds_file}${chr}".gds"
done

wait

#combine gds files per chromosome into one file
Rscript ./scripts/combine_gds.R $gds_file

#wait

#PC_AiR step
Rscript ./scripts/PC_AiR.R $gds_file $pheno_file $phenotypes $num_covariates ${covariates[@]} $snpset_file
