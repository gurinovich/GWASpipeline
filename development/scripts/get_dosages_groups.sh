#!/bin/bash

data_folder=$1
vcf_file=$2

file="./tmp/groups.txt"

groups=()
while IFS= read -r line; do
  groups+=("$line")
done < $file

for i in "${groups[@]}"
do
  for j in {1..22}
  do
  bcftools view -S "./tmp/"${i}.txt ${vcf_file}${j}"_qc2.vcf.gz" -o "./tmp/"${i}.$j.temp
  wait
  rm -f ${vcf_file}${j}"_qc2.vcf.gz"
  wait
  plink2 --vcf "./tmp/"${i}.$j.temp --double-id --freq --out "./tmp/"${i}.$j
  wait
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%DS]\n" "./tmp/"${i}.$j.temp -o "./tmp/"${i}.$j.dosages
  done
done