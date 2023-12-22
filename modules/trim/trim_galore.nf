


process TRIM_GALORE {

  publishDir "${out}/trim/${sample}", mode: 'copy', overwrite: 'false'
  cpus 8

  input:
  tuple val(sample), path(fastqs)

  output:
  tuple val(sample), path("*.fastq.gz"), emit: trim_fq

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
