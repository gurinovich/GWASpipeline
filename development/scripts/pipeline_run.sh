#!/bin/bash

module load R/3.6.0
module load vcftools
module load bcftools
module load plink/2.00a1LM
module load annovar/2018apr

source $1
num_covariates="${#covariates[@]}"
mkdir ./tmp


#QC
for chr in {1..22}
do
	sh ./scripts/QC1.sh ${vcf_file}${chr} &> "./log/qc_chr"$chr
done

wait

#convert VCF files to GDS files
for chr in {1..22}
do
	Rscript ./scripts/convert_vcf_gds.R ${vcf_file}${chr}"_qc2.vcf.gz" ${gds_file}${chr}".gds"
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
Rscript ./scripts/assocTestSingle_nullmod.R $gds_file $pheno_file $phenotypes $num_covariates ${covariates[@]} $model

wait

#GWAS
for chr in {1..22}
do
	Rscript ./scripts/assocTestSingle.R ${gds_file}${chr}".gds" $pheno_file $phenotypes $num_covariates ${covariates[@]} ${result_file}${chr}".csv" $test
done

wait

#summary
find ./results -maxdepth 1 -name '*.csv' | xargs -n 1 tail -n +2 | awk -F ',' '{print $1 "," $2 "," $3 "," $5 "," $6 "," $7 "," $8 "," $9 "," $10}' > ${result_file}".csv"

wait

sed -i '1i rs.id,chr,pos,n.obs,freq,Est,Est.SE,Stat,pval' ${result_file}".csv"

wait

#summary plots
Rscript ./scripts/qqplot_manhattanplot_MAF.R ${result_file}".csv"

wait

### calculate and add CAFs to the results

#To create files for each of the subgroups in the phenotype file in the column "group" (FOR NOW ASSUME it is provided - EDIT later if not provided, or indicate that it has to be provided and the same value can be in it if no different groups to compare) AND allele frequencies for cases and controls if model == "logistic":
Rscript ./scripts/samples_lists_create.R $pheno_file $model

wait

# calculate allele frequencies for each group of subjects
file=$data_folder"groups.txt"

groups=()
while IFS= read -r line; do
  groups+=("$line")
done < $file

for i in "${groups[@]}"
do
  for j in {1..22}
  do
  bcftools view -S ${data_folder}${i}.txt ${vcf_file}${j}"_qc2.vcf.gz" -o ${data_folder}${i}.$j.temp
  wait
  plink2 --vcf ${data_folder}${i}.$j.temp --double-id --freq --out ${data_folder}${i}.$j
  wait
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%DS]\n" ${data_folder}${i}.$j.temp -o ${data_folder}${i}.$j.dosages
  done
done

wait

#calculate CAFs from dosages:
Rscript ./scripts/calc_cafs_dosages.R $data_folder

wait

# add info on AFs and others to the results files (updates the files)
Rscript ./scripts/combine_results.R $data_folder $result_file

wait

#create input files for annovar
for chr in {1..22}
do
        Rscript ./scripts/prep_annovar_input.R $chr
done

wait

#Download humandb with specified version

annotate_variation.pl -downdb -buildver $refver -webfrom annovar refGene ./tmp/humandb/
wait

# Run annovar

for chr in {1..22}
do
table_annovar.pl -build $refver ./tmp/"result_file"${chr}"_snps_input.txt" ./tmp/humandb/ -out ./results/"chr"${chr}"_EL_GWAS" -remove -protocol refGene -operation g -nastring . -csvout
done

##Remove tmp dir
\rm -rf ./tmp


