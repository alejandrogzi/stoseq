#!/usr/env/bin nextflow


process TRIMMOMATIC {

    publishDir "${out}/trim/${sample}", mode: 'copy', overwrite: 'false'
    cpus 6

    input:
        tuple val(sample), path(fastq)
        val out
        val dir

    output:
        tuple val(sample), path("*.paired.trim*.fastq.gz"), emit: trim_paired

    script:
    """
    mkdir -p ${out}/trim/${sample}

    trimmomatic \\
    PE \\
    -threads ${task.cpus} \\
    ${fastq[0]} ${fastq[1]} \\
    ${sample}.paired.trim_1.fastq.gz ${sample}.unpaired.trim_1.fastq.gz \\
    ${sample}.paired.trim_2.fastq.gz ${sample}.unpaired.trim_2.fastq.gz \\
    -phred33 \\
    SLIDINGWINDOW:4:20 MINLEN:36

    rm ${sample}.unpaired.trim_1.fastq.gz ${sample}.unpaired.trim_2.fastq.gz
    #rm ${dir}/${fastq[0]} ${dir}/${fastq[1]}

    #find ../../ -name ${fastq[0]} -type f -delete
    #find ../../ -name ${fastq[1]} -type f -delete
    """
}
