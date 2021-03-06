#!/bin/bash

vcf_file=$1".vcf.gz"

# remove SNPs with > 1% missing rate
vcftools --gzvcf $vcf_file --max-missing 0.99 --recode --stdout | gzip -c > $1"_qc.vcf.gz"

wait

# remove monomorphic SNPs:
bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $1"_qc.vcf.gz" -Oz -o $1"_qc2.vcf.gz"

# remove $1"_qc.vcf.gz"
rm -f $1"_qc.vcf.gz"
