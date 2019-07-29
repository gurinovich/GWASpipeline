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

