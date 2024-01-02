#!/usr/env/bin nextflow


process RRNA {

  publishDir "${out}/RRNA_reads/${sample}", mode: 'copy'
  cpus 10 

  input:
  tuple val(sample), path(fastqs)
  path db
  val out

  output:
    tuple val(sample), path("*non_rRNA.fastq.gz"), emit: reads
    tuple val(sample), path("*.log")     , emit: log

  script:
  """
  sortmerna \\
  ${'--ref '+db.join(' --ref ')} \\
  --reads ${fastqs[0]} \\
  --reads ${fastqs[1]} \\
  --threads $task.cpus \\
  --workdir . \\
  --aligned rRNA_reads \\
  --fastx \\
  --other non_rRNA_reads \\
  --paired_in \\
  --out2

  mv non_rRNA_reads_fwd.f*q.gz ${sample}_1.non_rRNA.fastq.gz
  mv non_rRNA_reads_rev.f*q.gz ${sample}_2.non_rRNA.fastq.gz
  mv rRNA_reads.log ${sample}.sortmerna.log
  """
}

