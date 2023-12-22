#!/usr/bin/env nextflow

process MAKE_TRANSCRIPTOME {

  input:
  file fasta
  file gtf

  output:
  path("*.fa"), emit: transcriptome

  script:
  """
  cargo install to-trans

  to-trans --fasta $fasta --gtf $gtf
  """
}

