#!/bin/bash

vcf_file=$1".vcf.gz"

#only keep variants that have been successfully genotyped in 50% of individuals and remove SNPs not in HWE & remove SNPs with missing rate > 1%
vcftools --gzvcf $vcf_file --max-missing 0.5 --hwe 0.000001 --max-missing 0.01 --recode --stdout | gzip -c > $1"_qc.vcf.gz"


