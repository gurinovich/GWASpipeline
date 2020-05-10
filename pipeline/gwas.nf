#!/usr/bin/env nextflow

log.info """\
-

G W A S  ~  P I P E L I N E

================================
indir     : $params.indir
outdir    : $params.outdir
vcf       : $params.vcf

-
"""

// Read vcf files
Channel
  .fromPath(params.vcf).splitCsv(header: true)
  .map {row -> [row.id, row.file_path]}
  .ifEmpty {error "File ${params.vcf} not parsed properly"}
  .set {vcf_files}
fix = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
fix.into {fix1; fix2; fix3}

// --------------------------------------

/*
** STEP 1 QC
*/
process qc_miss {
  publishDir "${params.outdir}/temp1", mode: 'copy'
  tag "$id"

  input:
  set id, vcf from vcf_files

  output:
  file '*vcf.gz' into temp1

  script:
  """
  vcftools --gzvcf $PWD/${vcf} --max-missing 0.99 --recode --stdout | gzip -c > ${id}_qc1.vcf.gz
  """
}

process qc_mono {
  publishDir "${params.outdir}/temp2", mode: 'copy'
  
  input:
  val x from fix1
  file(vcf) from temp1.collect()

  output:
  file '*qc2.vcf.gz' into temp2

  script:
  """
  bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="AR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' chr_${x}_qc1.vcf.gz -Oz -o chr_${x}_qc2.vcf.gz
  """
}


/*
** STEP 2
*/
process vcf_to_gds {
  publishDir "${params.outdir}/gds_files", mode: 'copy'
  
  input:
  val x from fix2
  file(vcf) from temp2.collect()

  output:
  file '*gds' into gds_files

  script:
  """
  Rscript $PWD/scripts/02_vcf_to_gds.R chr_${x}_qc2.vcf.gz chr_${x}.gds
  """
}
process merge_gds {
  publishDir "${params.outdir}/gds_merged", mode: 'copy'

  input:
  file(gds) from gds_files.collect()

  output:
  file 'merged.gds' into merged_gds

  script:
  """
  Rscript $PWD/scripts/02_merge_gds.R ${gds}
  """
}
