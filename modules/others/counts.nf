#!/usr/bin/env nextflow

process REDUCE_COUNTS {
  
  publishDir "${out}/counts", mode: 'copy'
  
  input:
  path samples
  val out

  output:
  path "counts.txt", emit: matrix

  script:
  """
  ../../../bin/make_counts.py --samples ${samples}
  """
}
