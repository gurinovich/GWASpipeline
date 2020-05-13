#!/usr/bin/env nextflow

log.info """\
-

G W A S  ~  P I P E L I N E

================================
indir     : $params.indir
outdir    : $params.outdir

vcf       : $params.vcf
pheno     : $params.pheno
snpset    : $params.snpset

phenotypes: $params.phenotypes
covars    : $params.covars
model     : $params.model
test      : $params.test
ref       : $params.ref

-
"""

// Read vcf files
Channel
  .fromPath(params.vcf).splitCsv(header: true)
  .map {row -> [row.id, row.file_path]}
  .ifEmpty {error "File ${params.vcf} not parsed properly"}
  .set {vcf_files}

// --------------------------------------

/*
** STEP 1 QC
*/
process qc_miss {
  publishDir "${params.outdir}/qc1", mode: 'copy'
  tag "$id"

  input:
  set id, vcf from vcf_files

  output:
  file '*vcf.gz' into qc1

  script:
  """
  vcftools --gzvcf $PWD/${vcf} --max-missing 0.99 --recode --stdout | gzip -c > ${id}_qc1.vcf.gz
  """
}

process qc_mono {
  publishDir "${params.outdir}/qc2", mode: 'copy'
  
  input:
  each x from 1..22
  file(vcf) from qc1.collect()

  output:
  file '*qc2.vcf.gz' into qc2_1
  file '*qc2.vcf.gz' into qc2_2

  script:
  """
  bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' chr_${x}_qc1.vcf.gz -Oz -o chr_${x}_qc2.vcf.gz
  """
}


/*
** STEP 2 vcf_to_gds
*/
process vcf_to_gds {
  publishDir "${params.outdir}/gds_files", mode: 'copy'
  
  input:
  each x from 1..22
  file(vcf) from qc2_1.collect()

  output:
  file '*' into gds_files1
  file '*' into gds_files2

  script:
  """
  Rscript $PWD/scripts/02_vcf_to_gds.R chr_${x}_qc2.vcf.gz chr_${x}.gds
  """
}
process merge_gds {
  publishDir "${params.outdir}/gds_merged", mode: 'copy'

  input:
  file '*' from gds_files1.collect()

  output:
  file 'merged.gds' into gds_merged

  script:
  """
  Rscript $PWD/scripts/02_merge_gds.R
  """
}

/*
** STEP 3 PCA and GRM
*/
process pcair {
  publishDir "${params.outdir}/pcair", mode: 'copy'
  
  input:
  file "merged.gds" from gds_merged.collect()

  output:
  file '*' into pcair

  script:
  """
  Rscript $PWD/scripts/03_PC_AiR.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.snpset}
  """
}

process pcrelate {
  publishDir "${params.outdir}/pcrelate", mode: 'copy'
  
  input:
  file "merged.gds" from gds_merged.collect()
  file '*' from pcair.collect()

  output:
  file '*' into pcrelate

  script:
  """
  Rscript $PWD/scripts/03_PC_Relate.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.snpset} pruned.rds king.rds pcs.rds pc.df.rds
  """
}

/*
** STEP 4 nullmod and GWAS
*/
process nullmod {
  publishDir "${params.outdir}/nullmod", mode: 'copy'
  
  input:
  file "merged.gds" from gds_merged.collect()
  file '*' from pcair.collect()
  file '*' from pcrelate.collect()

  output:
  file '*' into nullmod

  script:
  """
  Rscript $PWD/scripts/04_nullmod.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.model} pc.df.rds grm.rds
  """
}

process gwas {
  publishDir "${params.outdir}/gwas", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from gds_files2.collect()
  file '*' from pcair.collect()
  file '*' from pcrelate.collect()
  file '*' from nullmod.collect()

  output:
  file '*' into gwas

  script:
  """
  Rscript $PWD/scripts/04_GWAS.R chr_${x}.gds ${params.pheno} ${params.phenotypes} ${params.covars} pc.df.rds grm.rds nullmod.rds ${params.test} chr_${x}.csv 
  """
}

/*
** STEP 5 summary and plot
*/

process samples_by_group {
  publishDir "${params.outdir}/samples_by_group", mode: 'copy'
  
  input:
  file '*' from pcair.collect()
  
  output:
  file '*' into samples_by_group

  script:
  """
  Rscript $PWD/scripts/05_samples_lists_create.R ${params.pheno} ${params.model} analysis.sample.id.rds
  """
}

group_datasets = Channel.fromPath(params.group)

group_datasets.into {group_datasets1; group_datasets2; group_datasets3}

process maf_by_group1 {
  publishDir "${params.outdir}/maf_by_group1", mode: 'copy'
  
  input:
  each i from 1..22
  file group from group_datasets1
  file '*' from qc2_2.collect()
  
  output:
  file "${group}.${i}.temp" into maf_by_group1

  script:
  """
  bcftools view -S ${group} chr_${i}_qc2.vcf.gz -o ${group}.${i}.temp
  """
}

process maf_by_group2 {
  publishDir "${params.outdir}/maf_by_group2", mode: 'copy'
  
  input:
  each i from 1..22
  file group from group_datasets2
  file '*' from maf_by_group1.collect()
  
  output:
  file '*' into maf_by_group2

  script:
  """
  plink2 --vcf ${group}.${i}.temp --double-id --freq --out ${group}.${i}
  bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER[\t%DS]\n" ${group}.${i}.temp -o ${group}.${i}.dosages
  """
}

process caf_dosage {
  publishDir "${params.outdir}/caf_dosage", mode: 'copy'
  
  input:
  file '*' from maf_by_group2.collect()
  
  output:
  file '*' into caf_dosage

  script:
  """
  Rscript $PWD/scripts/05_calc_cafs_dosages.R ${params.pheno} ${params.model} $PWD/results/maf_by_group2/
  """
}

process merge_by_chr {
  publishDir "${params.outdir}/merge_by_chr", mode: 'copy'
  
  input:
  file '*' from caf_dosage.collect()
  
  output:
  file '*' into merge_by_chr

  script:
  """
  Rscript $PWD/scripts/05_merge_by_chr.R ${params.pheno} ${params.model} $PWD/results/maf_by_group2/ $PWD/results/caf_dosage/
  """
}

process add_to_results {
  publishDir "${params.outdir}/add_to_results", mode: 'copy'
  
  input:
  file '*' from merge_by_chr.collect()
  
  output:
  file '*' into add_to_results1
  file '*' into add_to_results2
  file '*' into add_to_results3

  script:
  """
  Rscript $PWD/scripts/05_combine_results.R $PWD/results/gwas/ $PWD/results/merge_by_chr/
  """
}

process plot {
  publishDir "${params.outdir}/qq_manhattan", mode: 'copy'
  
  input:
  file '*' from add_to_results1.collect()
  
  output:
  file '*' into qq_manhattan

  script:
  """
  Rscript $PWD/scripts/05_qqplot_manhattanplot.R all_chr.csv
  """
}

/*
** STEP 6 annotation
*/

process annovar_input {
  publishDir "${params.outdir}/annovar_input", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from add_to_results2.collect()
  
  output:
  file '*' into annovar_input

  script:
  """
  Rscript $PWD/scripts/06_annovar_input.R ${x}
  """
}

process annovar_ref {  
  input:
  
  output:
  
  script:
  """
  annotate_variation.pl --downdb --buildver ${params.ref} --webfrom annovar refGene $PWD/results/humandb
  """
}

process annovar {
  publishDir "${params.outdir}/annovar", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from annovar_input.collect()

  output:
  file '*' into annovar

  script:
  """
  table_annovar.pl -build ${params.ref} chr${x}_snps_input.txt $PWD/results/humandb/ -out chr${x}_annotation -remove -protocol refGene -operation g -nastring . -csvout
  """
}

process add_annovar {
  publishDir "${params.outdir}/annotated_results", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from add_to_results3.collect()
  file '*' from annovar.collect()

  output:
  file '*' into annotated_results

  script:
  """
  Rscript $PWD/scripts/06_add_anno_results.R ${x} ${params.ref}
  """
}

