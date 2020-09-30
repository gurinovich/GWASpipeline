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
  .map {row -> tuple(row.Chromosome, file(row.VCF))}
  .ifEmpty {error "File ${params.vcf_list} not parsed properly"}
  .set {vcf_files}

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
  tag "$chr"
  publishDir "${params.outdir}/QC/qc_miss", mode: 'copy'
  
  input:
  set val(chr), file(vcf) from vcf_files

  output:
  set val(chr), file("*vcf.gz") into qc1

  script:
  """
  vcftools --gzvcf $vcf --max-missing 0.99 --recode --stdout | gzip -c > qc1.vcf.gz
  """
}

process qc_mono {
  tag "$chr"
  publishDir "${params.outdir}/QC/qc_mono", mode: 'copy'
  
  input:
  set val(chr), file(vcf) from qc1

  output:
  set val(chr), file("*vcf.gz") into qc2_1, qc2_2

  script:
  """
  bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $vcf -Oz -o qc2.vcf.gz
  """
}

