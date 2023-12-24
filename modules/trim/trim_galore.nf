


process TRIM_GALORE {

  publishDir "${out}/trim/${sample}", mode: 'copy', overwrite: 'false'
  cpus 8

  input:
  tuple val(sample), path(fastqs)
  val out

  output:
  tuple val(sample), path("*.f*q.gz"), emit: trim_fq

  script:
  """
  trim_galore \\
  --cores ${task.cpus} \\
  --quality 30 \\
  --paired \\
  --gzip \\
  --illumina \\
  --length 30 \\
  $fastqs
  """
}
