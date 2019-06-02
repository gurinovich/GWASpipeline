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

// --------------------------------------

/*
** STEP 1
*/
process vcf_to_gds {
  publishDir "${params.outdir}/gds_files", mode: 'copy'
  tag "$id"

  input:
  set id, vcf from vcf_files

  output:
  file '*gds' into gds_files

  script:
  """
  Rscript $PWD/scripts/01_vcf_to_gds.R ${vcf} ${id}.gds
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
  Rscript $PWD/scripts/01_merge_gds.R ${gds}
  """
}
