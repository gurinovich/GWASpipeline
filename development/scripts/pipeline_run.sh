#!/bin/bash

module load R/3.6.0
module load vcftools
module load bcftools

source $1

#QC
for chr in {1..22}
do
	sh ./scripts/QC1.sh ${vcf_file}${chr} &> "./log/qc_chr"$chr
done

wait

#convert VCF files to GDS files
for chr in {1..22}
do
	Rscript ./scripts/convert_vcf_gds.R ${vcf_file}${chr}"_qc.vcf.gz" ${gds_file}${chr}".gds"
done

wait

#combine gds files per chromosome into one file
Rscript ./scripts/combine_gds.R $gds_file

wait

#PC_AiR step
Rscript ./scripts/PC_AiR.R $gds_file $pheno_file $phenotypes $num_covariates ${covariates[@]} $snpset_file

wait

#GRM step
Rscript ./scripts/PC_Relate.R $gds_file $pheno_file $phenotypes $num_covariates ${covariates[@]} $snpset_file

wait

#generate Null model
Rscript ./scripts/assocTestSingle_nullmod.R $gds_file $pheno_file $phenotypes $num_covariates ${covariates[@]}

wait

#GWAS
for chr in {1..22}
do
	Rscript ./scripts/assocTestSingle.R ${gds_file}${chr}".gds" $pheno_file $phenotypes $num_covariates ${covariates[@]} ${result_file}${chr}".csv"
done

wait

#summary
find ./results -maxdepth 1 -name '*.csv' | xargs -n 1 tail -n +2 | awk -F ',' '{print $2 "," $3 "," $5 "," $6 "," $7 "," $8 "," $9 "," $10}' > ${result_file}".csv"

wait

sed -i '1i chr,pos,n.obs,freq,Score,Score.SE,Score.Stat,Score.pval' ${result_file}".csv"

wait

#summary plots
Rscript ./scripts/qqplot_manhattanplot_MAF.R ${result_file}".csv"

