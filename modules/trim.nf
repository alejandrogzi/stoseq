#!/usr/env/bin nextflow


process TRIMMOMATIC {

    publishDir "${out}/trimmomatic/${sample}", mode: 'copy', overwrite: 'false'
    executor 'local'
    cpus 6
    debug true

    input:
        tuple val(sample), path(fastq)
        val out

    output:
        tuple val(sample), path("*.paired.trim*.fastq.gz"), emit: trim_paired

    script:
    """
    mkdir -p ${out}/trimmomatic/${sample}

    trimmomatic \\
    PE \\
    -threads ${task.cpus} \\
    ${fastq[0]} ${fastq[1]} \\
    ${sample}.paired.trim_1.fastq.gz ${sample}.unpaired.trim_1.fastq.gz \\
    ${sample}.paired.trim_2.fastq.gz ${sample}.unpaired.trim_2.fastq.gz \\
    -phred33 \\
    SLIDINGWINDOW:4:20 MINLEN:36

    rm ${sample}.unpaired.trim_1.fastq.gz ${sample}.unpaired.trim_2.fastq.gz
    rm ${params.dir}/${fastq[0]} ${params.dir}/${fastq[1]}
    """
}