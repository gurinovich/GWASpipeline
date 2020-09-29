#!/usr/bin/env nextflow


/*
USER INPUT PARAMETERS
*/
date = new Date().format( 'yyyyMMdd' )

params.vcf_list    = null
params.pheno       = null

params.phenotypes  = null
params.covars      = null

params.pca_grm     = false
params.snpset      = null
params.grm         = null

params.gwas        = false
params.model       = "linear"
params.test        = "Score"
params.imputed     = false

params.gene_based  = false
params.max_maf     = 0.01
params.method      = "Burden"

params.longitudinal = false
params.random_slope = null

params.group               = "null"
params.min_maf             = 0.1
params.max_pval_manhattan  = 1

params.max_pval    = 1
params.ref_genome  = "hg19"

params.help        = null

println()

/*
OUTPUT DIRECTORY
*/

params.outdir = "$PWD/Analysis_Results-${date}"

/*
HELP MESSAGE
*/

if(params.help){
log.info """
USAGE:

The command for runing the pipeline is:

 for gwas:
 nextflow gwas.nf --vcf_list $PWD/data/toy_vcf.csv --pheno $PWD/data/pheno_file_logistic.csv --snpset $PWD/data/snpset.txt --phenotypes outcome --covars age,sex,PC1,PC2,PC3,PC4 --ref_genome hg19 --gwas --model linear --test Wald

 for Gene-based Test:

 for longitudinal analysis:

Mandatory arguments:
--vcf_list                 String        Name of the two-column csv mapping file: id , file_path 
--pheno                    String        Name of the phenotype file
--phenotype                String        Name of the phenotype column


Optional arguments:
--snpset                   String       Name of the two column txt file separated by comma: chr,pos
--covars                   String       Name of the covariates to include in analysis model separated by comma
--ref_genome               String       Name of the reference genome: hg19 or hg38
--gwas                     Logical      If true, run gwas
--model                    String       Name of regression model for gwas: linear or logistic
--test                     String       Name of statistical test for significance: Wald or Score

--outdir


"""
  exit 1
}

log.info """\
-

G W A S  ~  P I P E L I N E

================================
outdir    : $params.outdir

vcf       : $params.vcf_list
pheno     : $params.pheno
snpset    : $params.snpset

phenotypes: $params.phenotypes
covars    : $params.covars
model     : $params.model
test      : $params.test
ref       : $params.ref_genome

-
"""

// Read vcf files
Channel
  .fromPath(params.vcf_list).splitCsv(header: true)
  .map {row -> [row.id, row.file_path]}
  .ifEmpty {error "File ${params.vcf} not parsed properly"}
  .set {vcf_files}

// --------------------------------------

/*
** STEP 0 Configuration Check
*/
if(!(params.pca_grm|params.gwas|params.gene_based|params.longitudinal)){
  log.info """
  EXIT: NO ANALYSIS SELECTED
  """
  exit 1
}
if(params.gwas&params.gene_based|params.gwas&params.longitudinal|params.gene_based&params.longitudinal|params.gwas&params.gene_based&params.longitudinal){
  log.info """
  EXIT: MORE THAN ONE ANLYSES SELECTED
  """
  exit 1
}


/*
** STEP 1 QC
*/
process qc_miss {
  publishDir "${params.outdir}/QC/qc_miss", mode: 'copy'
  tag "$id"
  
  input:
  set id, vcf from vcf_files

  output:
  file '*vcf.gz' into qc1

  script:
  """
  vcftools --gzvcf ${vcf} --max-missing 0.99 --recode --stdout | gzip -c > ${id}_qc1.vcf.gz
  """
}

process qc_mono {
  publishDir "${params.outdir}/QC/qc_mono", mode: 'copy'
  
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
  publishDir "${params.outdir}/GDS/gds_files", mode: 'copy'
    
  input:
  each x from 1..22
  file(vcf) from qc2_1.collect()

  output:
  file '*' into gds_files1
  file '*.gds' into gds_files2
  file '*.gds' into gds_files3

  script:
  """
  Rscript $PWD/scripts/02_vcf_to_gds.R chr_${x}_qc2.vcf.gz chr_${x}.gds vcf_to_gds_chr${x}.log
  """
}

process merge_gds {
  publishDir "${params.outdir}/GDS/gds_merged", mode: 'copy'

  input:
  file '*' from gds_files1.collect()

  output:
  file '*' into gds_merged1
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
if(params.pca_grm){
  process pcair {
  publishDir "${params.outdir}/PCA_GRM/pcair", mode: 'copy'
      
  input:
  file "*" from gds_merged1.collect()

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
  publishDir "${params.outdir}/PCA_GRM/pcrelate", mode: 'copy'
  
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
}

/*
** STEP 4 nullmod and gwas/gene-based/longitudinal analysis
*/
if((params.gwas|params.gene_based)&params.pca_grm){
  process nullmod {
    publishDir "${params.outdir}/Association_Test/nullmod", mode: 'copy'
  
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
}

if((params.gwas|params.gene_based)&!params.pca_grm){
  process nullmod_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/nullmod", mode: 'copy'
  
    input:
    file "merged.gds" from gds_merged3.collect()

    output:
    file '*' into nullmod
    file '*' into nullmod1

    script:
    """
    Rscript $PWD/scripts/04_nullmod_skip_pca_grm.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.model} ${params.grm}
    """
  }
}

if(params.gwas&params.pca_grm){
  process gwas {
    publishDir "${params.outdir}/Association_Test/gwas", mode: 'copy'
    
    input:
    each x from 1..22
    file '*' from gds_files2.collect()
    file '*' from nullmod.collect()
 
    output:
    file '*' into gwas1
    file '*.csv' into gwas2

    script:
    """
    Rscript $PWD/scripts/04_gwas.R chr_${x}.gds annot_pc.rds nullmod.rds ${params.test} ${params.imputed} chr_${x}.csv chr${x}_gwas.log
    """
  }
}

if(params.gwas&!params.pca_grm){ 
  process gwas_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/gwas", mode: 'copy'
    
    input:
    each x from 1..22
    file '*' from gds_files2.collect()
    file '*' from nullmod.collect()
 
    output:
    file '*' into gwas1
    file '*.csv' into gwas2

    script:
    """
    Rscript $PWD/scripts/04_gwas.R chr_${x}.gds annot.rds nullmod.rds ${params.test} ${params.imputed} chr_${x}.csv chr${x}_gwas.log
    """
  }
}

if(params.gene_based&params.pca_grm){
  process gene_based {
    publishDir "${params.outdir}/Association_Test/gene_based", mode: 'copy'
    
    input:
    each x from 1..22
    file '*' from gds_files2.collect()
    file '*' from nullmod.collect()
 
    output:
    file '*' into gene_based
  
    script:
    """
    Rscript $PWD/scripts/04_gene_based.R chr_${x}.gds annot_pc.rds nullmod.rds ${params.max_maf} ${params.method} chr_${x}.csv chr_${x}.rds chr${x}_gene_based.log
    """
  }
}

if(params.gene_based&!params.pca_grm){
  process gene_based_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/gene_based", mode: 'copy'
    
    input:
    each x from 1..22
    file '*' from gds_files2.collect()
    file '*' from nullmod.collect()
 
    output:
    file '*' into gene_based
 

    script:
    """
    Rscript $PWD/scripts/04_gene_based.R chr_${x}.gds annot.rds nullmod.rds ${params.max_maf} ${params.method} chr_${x}.csv chr_${x}.rds chr${x}_gene_based.log
    """
  }
}

if(params.longitudinal&params.pca_grm){
  process nullmod_longitudinal {
    publishDir "${params.outdir}/Association_Test/nullmod_longitudinal", mode: 'copy'
  
    input:
    file '*' from pcair2.collect()
    file '*' from pcrelate.collect()

    output:
    file '*' into nullmod_longitudinal

    script:
    """
    Rscript $PWD/scripts/04_nullmod_longitudinal.R ${params.pheno} ${params.phenotypes} ${params.covars} ${params.model} analysis.sample.id.rds pc.df.rds grm.rds ${params.random_slope}
    """
  }
}

if(params.longitudinal&!params.pca_grm){
  process nullmod_longitudinal_skip_pca_grm {
    publishDir "${params.outdir}/Association_Test/nullmod_longitudinal", mode: 'copy'
  
    input:
    file "merged.gds" from gds_merged3.collect()

    output:
    file '*' into nullmod_longitudinal
    file '*' into nullmod1

    script:
    """
    Rscript $PWD/scripts/04_nullmod_longitudinal_skip_pca_grm.R merged.gds ${params.pheno} ${params.phenotypes} ${params.covars} ${params.model} ${params.grm} ${params.random_slope}
    """
  }
}

if(params.longitudinal){
  process gwas_longitudinal {
    publishDir "${params.outdir}/Association_Test/gwla", mode: 'copy'
    
    input:
    each x from 1..22
    file '*' from gds_files2.collect()
    file '*' from nullmod_longitudinal.collect()
 
    output:
    file '*' into gwas_longitudinal
   
    script:
    """
    Rscript $PWD/scripts/04_gwas_longitudinal.R chr_${x}.gds nullmod_longitudinal.rds chr_${x}.txt chr${x}_gwas_longitudinal.log
    """
  }
}

/*
** STEP 5 summary and plot
*/

if(params.gene_based){
  process combine_results_gene {
    publishDir "${params.outdir}/Summary_Plot/combined_results", mode: 'copy'
  
    input:
    file '*' from gene_based.collect()
  
    output:
    file '*' into combined_results

    script:
    """
    Rscript $PWD/scripts/05_combine_results_gene.R ${params.outdir}/Association_Test/gene_based/
    """
  }

  process plot_gene {
    publishDir "${params.outdir}/Summary_Plot/qq_plot", mode: 'copy'
  
    input:
    file '*' from combined_results.collect()
  
    output:
    file '*' into qq_plot

    script:
    """
    Rscript $PWD/scripts/05_qqplot_gene.R all_chr.csv
    """
  }
}


if((params.gwas|params.longitudinal)&params.pca_grm){
 process caf_by_group {
  publishDir "${params.outdir}/Summary_Plot/caf_by_group", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from gds_files3.collect()
  file '*' from pcair3.collect()
  
  output:
  file '*' into caf_by_group

  script:
  """
  Rscript $PWD/scripts/05_caf_by_group.R chr_${x}.gds ${params.pheno} analysis.sample.id.rds ${params.model} ${params.phenotypes} ${params.group} chr_${x}_caf_by_group.csv chr${x}_caf_by_group.log
  """
 }
}

if((params.gwas|params.longitudinal)&!params.pca_grm){
 process caf_by_group_skip_pca_grm {
  publishDir "${params.outdir}/Summary_Plot/caf_by_group", mode: 'copy'
  
  input:
  each x from 1..22
  file '*' from gds_files3.collect()
  file '*' from nullmod1.collect()
  
  output:
  file '*' into caf_by_group

  script:
  """
  Rscript $PWD/scripts/05_caf_by_group.R chr_${x}.gds ${params.pheno} analysis.sample.id.rds ${params.model} ${params.phenotypes} ${params.group} chr_${x}_caf_by_group.csv chr${x}_caf_by_group.log
  """
 }
}

if(params.longitudinal){
  process coincide_gwas {
    publishDir "${params.outdir}/Summary_Plot/gwla", mode: 'copy'
  
    input:
    each x from 1..22
    file '*' from gwas_longitudinal.collect()
  
    output:
    file '*' into gwas1

    script:
    """
    Rscript $PWD/scripts/05_coincide_by_chr.R chr_${x}.txt chr_${x}.csv chr${x}_coincide.log
    """
  }
}

if(params.gwas){
  process merge_by_chr {
    publishDir "${params.outdir}/Summary_Plot/merge_by_chr", mode: 'copy'
  
    input:
    each x from 1..22
    file '*' from gwas1.collect()
    file '*' from caf_by_group.collect()
  
    output:
    file '*' into merge_by_chr

    script:
    """
    Rscript $PWD/scripts/05_merge_by_chr.R chr_${x}.csv chr_${x}_caf_by_group.csv ${params.model} chr_${x}_caf_annotated.csv chr${x}_merge.log
    """
  }
}

if(params.longitudinal){
  process merge_by_chr_longitudinal {
    publishDir "${params.outdir}/Summary_Plot/merge_by_chr", mode: 'copy'
  
    input:
    each x from 1..22
    file '*' from gwas1.collect()
    file '*' from caf_by_group.collect()
  
    output:
    file '*' into merge_by_chr

    script:
    """
    Rscript $PWD/scripts/05_merge_by_chr.R chr_${x}.csv chr_${x}_caf_by_group.csv ${params.model} chr_${x}_caf_annotated.csv chr${x}_merge.log
    """
  }
}

if(params.gwas|params.longitudinal){
process combine_results {
  publishDir "${params.outdir}/Summary_Plot/combined_results", mode: 'copy'
  
  input:
  file '*' from merge_by_chr.collect()
  
  output:
  file '*' into combined_results1
  file '*' into combined_results2

  script:
  """
  Rscript $PWD/scripts/05_combine_results.R ${params.outdir}/Summary_Plot/merge_by_chr/
  """
}

process plot {
  publishDir "${params.outdir}/Summary_Plot/qq_manhattan", mode: 'copy'
  
  input:
  file '*' from combined_results1.collect()
  
  output:
  file '*' into qq_manhattan

  script:
  """
  Rscript $PWD/scripts/05_qqplot_manhattanplot.R all_chr_caf_annotated.csv ${params.min_maf} ${params.max_pval_manhattan}
  """
}

/*
** STEP 6 annotation
*/

process annovar_input {
  publishDir "${params.outdir}/Annotation/annovar_input", mode: 'copy'
    
  input:
  file '*' from combined_results2.collect()
  
  output:
  file '*' into annovar_input1
  file '*' into annovar_input2

  script:
  """
  Rscript $PWD/scripts/06_annovar_input.R ${params.max_pval}
  """
}

process annovar_ref {  
  publishDir "${params.outdir}/Annotation/", mode: 'copy'
  
  script:
  """
  mkdir -p ${params.outdir}/Annotation
  annotate_variation.pl --downdb --buildver ${params.ref_genome} --webfrom annovar refGene ${params.outdir}/Annotation/humandb
  """
}

process annovar {
  publishDir "${params.outdir}/Annotation/annovar", mode: 'copy'
    
  input:
  file '*' from annovar_input1.collect()

  output:
  file '*' into annovar

  script:
  """
  table_annovar.pl -build ${params.ref_genome} top_snps_input.txt ${params.outdir}/Annotation/humandb/ -out top_annotation -remove -protocol refGene -operation g -nastring . -csvout
  """
}

process add_annovar {
  publishDir "${params.outdir}/Annotation/annotated_results", mode: 'copy'
  
  input:
  file '*' from annovar_input2.collect()
  file '*' from annovar.collect()

  output:
  file '*' into report

  script:
  """
  Rscript $PWD/scripts/06_add_anno_results.R ${params.ref_genome}
  """
}
}

/*
** STEP 7 Report
*/

if(params.gwas|params.longitudinal){
  process report {
    publishDir "${params.outdir}/Report", mode: 'copy'
  
    input:
    file '*' from report.collect()
  
    output:
  
    script:
    """
    mkdir -p ${params.outdir}/Report
    Rscript -e 'rmarkdown::render("$PWD/scripts/07_report.Rmd")' ${params.outdir}
    mv $PWD/scripts/07_report.html ${params.outdir}/Report/Report.html
    """
  }
}

if(params.gene_based){
  process report_gene {
    publishDir "${params.outdir}/Report", mode: 'copy'
  
    input:
    file '*' from qq_plot.collect()
  
    output:
  
    script:
    """
    mkdir -p ${params.outdir}/Report
    Rscript -e 'rmarkdown::render("$PWD/scripts/07_report_gene.Rmd")' ${params.outdir} ${params.max_pval}
    mv $PWD/scripts/07_report_gene.html ${params.outdir}/Report/07_report_gene.html
    """
  }
}
