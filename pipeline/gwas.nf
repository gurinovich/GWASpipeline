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
  file '*' into gds_files3

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
  file 'merged.gds' into gds_merged1
  file 'merged.gds' into gds_merged2
  file 'merged.gds' into gds_merged3

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
  file "merged.gds" from gds_merged1.collect()

  output:
  file '*' into pcair1
  file '*' into pcair2
  file '*' into pcair3

  script:
  """
  Rscript $PWD/scripts/03_PC_AiR.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.snpset}
  """
}

process pcrelate {
  publishDir "${params.outdir}/pcrelate", mode: 'copy'
  
  input:
  file "merged.gds" from gds_merged2.collect()
  file '*' from pcair1.collect()

  output:
  file '*' into pcrelate

  script:
  """
  Rscript $PWD/scripts/03_PC_Relate.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.snpset} analysis.sample.id.rds annot.rds pruned.rds king.rds pcs.rds pc.df.rds
  """
}

/*
** STEP 4 nullmod and GWAS
*/
process nullmod {
  publishDir "${params.outdir}/nullmod", mode: 'copy'
  
  input:
  file "merged.gds" from gds_merged3.collect()
  file '*' from pcair2.collect()
  file '*' from pcrelate.collect()

  output:
  file '*' into nullmod

  script:
  """
  Rscript $PWD/scripts/04_nullmod.R merged.gds ${params.phenotypes} ${params.covars} ${params.model} analysis.sample.id.rds annot.rds pc.df.rds grm.rds
  """
}

process gwas {
  publishDir "${params.outdir}/gwas", mode: 'copy'
    
  input:
  each x from 1..22
  file '*' from gds_files2.collect()
  file '*' from nullmod.collect()

  output:
  file '*' into gwas1
  file '*' into gwas2

  script:
  """
  Rscript $PWD/scripts/04_GWAS.R chr_${x}.gds annot_pc.rds nullmod.rds ${params.test} chr_${x}.csv 
  """
}

/*
** STEP 5 summary and plot
*/

process caf_by_group {
  publishDir "${params.outdir}/caf_by_group", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from gds_files3.collect()
  file '*' from pcair3.collect()
  
  output:
  file '*' into caf_by_group

  script:
  """
  Rscript $PWD/scripts/05_caf_by_group.R chr_${x}.gds ${params.pheno} analysis.sample.id.rds ${params.model} ${params.phenotypes} chr_${x}_caf_by_group.csv
  """
}

process merge_by_chr {
  publishDir "${params.outdir}/merge_by_chr", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from gwas1.collect()
  file '*' from caf_by_group.collect()
  
  output:
  file '*' into merge_by_chr

  script:
  """
  Rscript $PWD/scripts/05_merge_by_chr.R chr_${x}.csv chr_${x}_caf_by_group.csv ${params.model} chr_${x}_caf_annotated.csv
  """
}

process combine_results {
  publishDir "${params.outdir}/combined_results", mode: 'copy'
  
  input:
  file '*' from merge_by_chr.collect()
  
  output:
  file '*' into combined_results1
  file '*' into combined_results2

  script:
  """
  Rscript $PWD/scripts/05_combine_results.R $PWD/results/merge_by_chr/
  """
}

process plot {
  publishDir "${params.outdir}/qq_manhattan", mode: 'copy'
  
  input:
  file '*' from combined_results1.collect()
  
  output:
  file '*' into qq_manhattan

  script:
  """
  Rscript $PWD/scripts/05_qqplot_manhattanplot.R all_chr_caf_annotated.csv
  """
}

/*
** STEP 6 annotation
*/

process annovar_input {
  publishDir "${params.outdir}/annovar_input", mode: 'copy'
    
  input:
  file '*' from combined_results2.collect()
  
  output:
  file '*' into annovar_input1
  file '*' into annovar_input2

  script:
  """
  Rscript $PWD/scripts/06_annovar_input.R
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
  file '*' from annovar_input1.collect()

  output:
  file '*' into annovar

  script:
  """
  table_annovar.pl -build ${params.ref} top_snps_input.txt $PWD/results/humandb/ -out top_annotation -remove -protocol refGene -operation g -nastring . -csvout
  """
}

process add_annovar {
  publishDir "${params.outdir}/annotated_results", mode: 'copy'
  
  input:
  file '*' from annovar_input2.collect()
  file '*' from annovar.collect()

  output:
  file '*' into annotated_results

  script:
  """
  Rscript $PWD/scripts/06_add_anno_results.R
  """
}

