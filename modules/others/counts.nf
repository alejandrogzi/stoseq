#!/usr/bin/env nextflow

process REDUCE_COUNTS {

  debug true

  input:
  path samples

  output:
  path "counts.txt", emit: matrix

  script:
  """
  ../../../bin/make_counts.py --samples ${samples}
  """
}
