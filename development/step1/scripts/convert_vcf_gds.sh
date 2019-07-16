#!/bin/bash

module load R/3.5.0
module load gcc/7.2.0

vcf_file=$1
gds_file=$2
pheno_file=$3

#convert VCF files to GDS files
for chr in $(seq 1 22); do
	Rscript convert_vcf_gds.R $vcf_file$chr.vcf.gz $gds_file$chr.gds
done

wait

#combine gds files per chromosome into one file
Rscript combine_gds.R $vcf_file $gds_file

#PCA_AiR step
Rscript PCA_AiR.R $gds_file $pheno_file