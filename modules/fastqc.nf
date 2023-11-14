#!/usr/env/bin nextflow

process FASTQC {

    publishDir "${out}/fastqc/${sample}", mode: 'copy', overwrite: 'false'
    cpus 4

    input:
        tuple val(sample), file(fastq)
        val(out)
    
    output:
        path("*.paired.trim*_fastqc.zip"), emit: fastqc_zip
        path("*.paired.trim*_fastqc.html")

    script:
    """
    mkdir -p ${outdir}/fastqc/${sample}

    fastqc ${fastq[0]} ${fastq[1]} \\
    --threads ${task.cpus} \\
    --kmers 7 --quiet --extract -f fastq
    """
}