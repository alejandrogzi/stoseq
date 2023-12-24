#!/usr/env/bin nextflow

process RIBODETECTOR {

  publishDir "${out}/ribodetector/${sample}", mode: 'copy'
  cpus 16

  input:
  tuple val(sample), path(fastqs)
  val out

  output:
  tuple val(sample), path("*.fastq.gz"), emit: reads

  script:
  """
  ribodetector_cpu \\
  -t ${task.cpus} \\
  -l 50 \\
  -i ${fastqs} \\
  -e rrna \\
  -o ${sample}.non_rRNA_1.fastq ${sample}.non_rRNA_2.fastq

  gzip ${sample}.non_rRNA_*.fastq
  """
}
